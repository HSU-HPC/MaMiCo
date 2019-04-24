#!/usr/bin/env python3
import argparse
import hashlib
import logging
import os
import subprocess
import json
import copy
import shutil

template = json.loads(r'''
{
  "project_name": "NavierStokes",
  "paths": {
    "peano_kernel_path": "./Peano/",
    "exahype_path": "./ExaHyPE",
    "output_directory": "./ApplicationExamples/CompressibleNavierStokes"
  },
  "computational_domain": {
    "dimension": 2,
    "end_time": 0.55,
    "offset": [
      0.0,
      0.0
    ],
    "width": [
      10.0,
      10.0
    ]
  },
  "shared_memory": {
    "cores": 28,
    "background_job_consumers": 14,
    "properties_file": "sharedmemory.properties",
    "autotuning_strategy": "dummy",
    "manual_pinning": false
  },
  "distributed_memory": {
      "ranks_per_node": 1,
      "load_balancing_strategy": "greedy_naive",
      "load_balancing_type": "static",
      "node_pool_strategy": "fair",
      "timeout": 6000,
      "buffer_size": 1600
  },
  "optimisation": {
    "fuse_algorithmic_steps": true,
    "fuse_algorithmic_steps_factor": 0.99,
    "spawn_predictor_as_background_thread": true,
    "spawn_amr_background_threads": true,
    "disable_vertex_exchange_in_time_steps": true,
    "time_step_batch_factor": 0.0,
    "disable_metadata_exchange_in_batched_time_steps": true,
    "double_compression": 0.0,
    "spawn_double_compression_as_background_thread": true
  },
  "solvers": [
    {
      "type": "ADER-DG",
      "name": "NavierStokesSolver_ADERDG",
      "order": 1,
      "maximum_mesh_size": 0.37,
      "maximum_mesh_depth": 10,
      "time_stepping": "global",
      "aderdg_kernel": {
        "language": "C",
        "nonlinear": true,
        "terms": [
          "flux",
          "viscous_flux",
          "source"
        ],
        "space_time_predictor": {},
        "optimised_terms": [],
        "optimised_kernel_debugging": [],
        "implementation": "generic"
      },
      "point_sources": 0,
      "variables": [
        {
          "name": "rho",
          "multiplicity": 1
        },
        {
          "name": "j",
          "multiplicity": 2
        },
        {
          "name": "E",
          "multiplicity": 1
        }
      ],
      "parameters": {
        "viscosity": 0.1,
        "scenario": "convergence"
      },
      "plotters": [
        {
          "type": "user::defined ErrorWriter",
          "name": "ErrorWriter",
          "time": 0.5,
          "repeat": 0.5,
          "output": "./results/solution",
          "variables": 4
        }
      ]
    }
  ]
}

''')

job_template = \
'''#@ job_type     = parallel
#@ class        = test
#@ node         = {nodes}
#@ tasks_per_node = {tasks_per_node}
#@ island_count = 1
#@ network.MPI = sn_all,not_shared,us 
##@ energy_policy_tag = ExaHyPE_Euler_energy_tag
##@ minimize_time_to_solution = yes
#@ wall_clock_limit = 0:30:00
#@ job_name = {job_name}
#@ initialdir = {working_dir}
#@ error  =  {log_dir}/job.$(schedd_host).$(jobid).err
#@ output =  {log_dir}/job.$(schedd_host).$(jobid).out
#@ notification=error
#@ notify_user=lukas@krenz.land
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh
module switch intel/18.0
module switch tbb/2018
module switch gcc/5

export OMP_NUM_THREADS={threads_per_task}
export MP_TASK_AFFINITY=core:{threads_per_task}

poe {exahype_bin} {exahype_config_file}
'''


def get_meshsize(factor):
    # We need to make sure we are slightly larger than actual mesh-size
    # otherwise wrong mesh size might be chosen.
    eps = 10e-3
    return 10/3**factor + eps

def get_exahype_root():
    return os.path.realpath(os.path.dirname(__file__) + "/../..") 

def get_application_path():
    return os.path.realpath(os.path.dirname(__file__))

def render_template(template, config, file_name):
    template['solvers'][0]['order'] = config['order']
    template['solvers'][0]['maximum_mesh_size'] = config['mesh_size']
    template['shared_memory']['cores'] = config['threads_per_task']
    template['shared_memory']['background_job_consumers'] = config['background_threads_per_task']
    template['distributed_memory']['ranks_per_node'] = config['tasks_per_node']

    #results_dir = '{results_dir}/order_{order}_{factor}'.format(
    #    results_dir=config['results_dir'],
   #     order=config['order'],
    #    factor=config['factor'])

    template['solvers'][0]['plotters'][0]['output'] = config['results_dir'] + "/solution"
    rendered_template = json.dumps(template, indent=4)

    # Quick hack to ensure semi-unique name.
    template_hash = hashlib.sha512(rendered_template.encode('utf-8')).hexdigest()[:8]
    config['template_hash'] = template_hash
    
    file_name = os.path.realpath(file_name.format(template_hash=template_hash))
    print(file_name)
    logging.info("Created file {file_name}".format(file_name=file_name))

    with open(file_name, 'w') as f:
        f.write(rendered_template)

    config['exahype_config_file'] = file_name

def render_jobscript(template, config, file_name):
    rendered_jobscript = template.format(**config)

    # Quick hack to ensure semi-unique name.
    #jobscript_hash = hashlib.sha512(rendered_jobscript.encode('utf-8')).hexdigest()[:8]
    file_name = os.path.abspath(file_name.format(template_hash=config['template_hash']))
    logging.info("Created file {file_name}".format(file_name=file_name))

    with open(file_name, 'w') as f:
        f.write(rendered_jobscript)

    config['jobscript_path'] = file_name

def run(template, my_env, config):
    logging.info("Starting running {config}".format(config=config))                                                  
    config['job_name'] = "exahype_o{order}_f{factor}".format(order=config['order'],                                  
                                                             factor=config['factor'])                                

    config['exahype_bin'] = '{bin_dir}/exahype_order_{order}'.format(                                                
        bin_dir=config['bin_dir'],
        order=config['order'])

    config['results_dir'] = '{dir}/order_{order}_factor_{factor}'.format(                                            
        dir=config['results_dir'],
        order=config['order'],
        factor=config['factor'])

    # Render exahype file.
    os.makedirs(config['script_dir'], exist_ok=True)
    template_file_name = config['script_dir'] + '/rendered_template_{template_hash}.exahype2'                        
    render_template(template, config, template_file_name)

    # Create dirs for output and logs.
    config['log_dir'] = config['tmp_dir'] + '/logs/'
    config['working_dir'] = config['tmp_dir'] + '/workingdirs/' + config['template_hash'] + '/'                      

    os.makedirs(config['log_dir'], exist_ok=True)
    os.makedirs(config['working_dir'], exist_ok=True)
    os.makedirs(config['results_dir'], exist_ok=True)
    os.makedirs(config['script_dir'], exist_ok=True)

    # Copy log filter to working dir
    shutil.copy2(get_application_path() + '/exahype.log-filter',                                                     
                 config['working_dir'])

    # Render jobscript.
    jobscript_file_name = config['script_dir'] + '/rendered_jobscript_{template_hash}.cmd'                           
    render_jobscript(job_template, config, jobscript_file_name)                                                      

    # Submit jobscript.
    subprocess.run(['llsubmit', config['jobscript_path']], env=my_env)                                               

    #subprocess.run([exahype_bin, template_path], env=my_env)                                                        
    logging.info("Submitted jobscript {jobscript_file_name} for config {template_file_name}.".format(                
        jobscript_file_name="",
        template_file_name=config['exahype_config_file']))

def build(template, my_env, config):
    os.makedirs(config['bin_dir'], exist_ok=True)
    os.makedirs(config['script_dir'], exist_ok=True)

    os.chdir(get_application_path()) # move to correct directory                                                     

    toolkit_bin = get_exahype_root() + '/Toolkit/toolkit.sh'                                                         

    template_file_name = config['script_dir'] + '/rendered_template_{template_hash}.exahype2'                        
    render_template(template, config, template_file_name)

    logging.info("Started building with {}".format(config))

    # Clean files from prior build.
    subprocess.run(['make', 'clean' ,'-j{}'.format(os.cpu_count())], env=my_env)                                     
    logging.info("Cleaned.")

    # Run toolkit for currenct settings.
    subprocess.run([toolkit_bin, config['exahype_config_file']], env=my_env)                                         
    logging.info("Ran toolkit.")

    # Compile program.
    subprocess.run(['make', '-j{}'.format(os.cpu_count())], env=my_env)                                              
    logging.info("Compiled.")

    # Move binary for later use.
    compiled_name = 'ExaHyPE-{}'.format(template['project_name'])
    new_name =  '/exahype_order_{}'.format(config['order'])
    new_path = os.path.realpath(config['bin_dir'] + new_name)

    shutil.copy2(compiled_name, new_path)
    logging.info("Copied file to {}".format(new_path))

    return new_path
   
def main():
    parser = argparse.ArgumentParser(description='Run simulations for various orders and mesh-sizes.')               
    parser.add_argument('--build', '-b', action='store_true')                                                        
    parser.add_argument('--run', '-r', action='store_true')

    args = parser.parse_args()

    print(__file__)
    print(get_exahype_root())

    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)                                        

    # Set up environment for compiling.
    # TODO: Maybe extract to shell script?
    my_env = os.environ.copy()
    my_env['SHAREDMEM'] = 'TBB'
    my_env['DISTRIBUTEDMEM'] = 'None'
    my_env['MODE'] = 'Release'

    user_config = {'nodes': 1,
                   'tasks_per_node': 1,
                   'threads_per_task': 28,
                   'background_threads_per_task': 14,

                   'tmp_dir': os.path.expanduser('/gss/scratch/pr83no/ga24dib2/exahype/'),                           
                   'results_dir': os.path.expanduser('/gpfs/work/pr83no/ga24dib2/exahype/results/'),                 
                   'bin_dir': os.path.expanduser('/gpfs/work/pr83no/ga24dib2/exahype/bin/'),                         
                   'script_dir': os.path.expanduser('/gss/scratch/pr83no/ga24dib2/exahype/templates/')               
    }
    user_config['total_tasks'] = user_config['nodes'] * user_config['tasks_per_node']                                

    order_grid = [1,2,3,4,5,6]
    max_factor = 4
    factor_grid = range(1, max_factor+1)

    # Build for various orders.
    if args.build:
        logging.info("Start compiling.")
        for order in order_grid:
            run_config = {'order': order,
                    'factor': max_factor,
                    'mesh_size': get_meshsize(max_factor)}
            config = {**user_config, **run_config}
            build(template=template, my_env=my_env, config=config)                                                   

    if args.run:
        logging.info("Start running.")
        for order in order_grid:
            for factor in factor_grid:
                run_config = {'order': order,
                        'factor': factor,
                        'mesh_size': get_meshsize(factor)}
                config = {**user_config, **run_config}
                run(template=template, my_env=my_env, config=config)                                                 

if __name__ == '__main__':
    main()
