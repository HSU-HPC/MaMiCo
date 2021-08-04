# -*- coding: utf-8 -*-
'''
This file can be used to train and save networks in the form of SavedModel files to later be imported by the NeuralNetJunctor
'''



#this suppresses info messages in console (2 suppresses warnings, 3 suppresses errors)
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'  # or any {'0', '1', '2', '3'}

import time

#suppresses deprecation warnings
import tensorflow.python.util.deprecation as deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False

import tensorflow as tf
import tensorflow_probability as tfp

import numpy as np
import pandas as pd
tf.keras.backend.set_floatx('float32')
from matplotlib import pyplot as plt
import math

tfd = tfp.distributions

class timecallback(tf.keras.callbacks.Callback):
    def __init__(self):
        self.times = []
        # use this value as reference to calculate cummulative time taken
        self.current = time.process_time()
    def on_epoch_end(self,epoch,logs = {}):
        self.times.append(time.process_time() - self.current)
        self.current = time.process_time()
    def on_train_end(self,logs = {}):
        print("Average Time: "+str(np.sum(self.times)/len(self.times)))
        print("Total Time: "+str(np.sum(self.times)))
        print(len(self.times))
        
class LearningRateLoggingCallback(tf.keras.callbacks.Callback):

    def on_epoch_end(self, epoch,logs={}):
        print("lr: ")
        print(self.model.optimizer._decayed_lr(tf.float32).numpy())
        if(self.model.optimizer._decayed_lr(tf.float32).numpy()<2e-6):
            self.model.optimizer.lr.decay_rate=1
            self.model.optimizer.lr.initial_learning_rate=2e-6

def prior_trainable(kernel_size, bias_size=0, dtype=None):
  n = kernel_size + bias_size
  c = np.log(np.expm1(1))
  return tf.keras.Sequential([
      tfp.layers.VariableLayer(2*n, dtype=dtype),
      tfp.layers.DistributionLambda(lambda t: tfd.Independent(
          tfd.Normal(loc=t[...,:n], scale=tf.nn.softplus(c+t[..., n:])),
          reinterpreted_batch_ndims=1)),
  ])

def posterior_mean_field(kernel_size, bias_size=0, dtype=None):
  n = kernel_size + bias_size
  c = np.log(np.expm1(1))#e-7
  return tf.keras.Sequential([
      tfp.layers.VariableLayer(2*n, dtype=dtype),
      tfp.layers.DistributionLambda(lambda t: tfd.Independent(
          tfd.Normal(loc=t[..., :n]+0.2, scale=0.00001*tf.nn.softplus(c + t[..., n:])),
          reinterpreted_batch_ndims=1)),
  ])

#loads all csv files specified in the lists and arranges them by timestep
def get_real_trainig_data(datasets_in, datasets_label):
    dflist_in = []
    dflist_label = []
    
    if(len(datasets_in)!=len(datasets_label)):
        print("size of inputs does not match size of labels")
        exit(1)
    tmp=len(datasets_label[0])
    for i in range(1,len(datasets_label)):
        if(len(datasets_label[i])!=tmp):
            print("not all scenarios are the same length")
            exit(2)
    
    for scenario_in in datasets_in:
        for dataset_in in scenario_in:
            dflist_in.append(pd.read_csv(dataset_in, sep=";", header=None))
        
    for scenario_label in datasets_label:
        for dataset_label in scenario_label:
            dflist_label.append(pd.read_csv(dataset_label, sep=";", header=None))
        
    for i in range(len(dflist_in)):
        dflist_in[i]=dflist_in[i][(dflist_in[i][1]!=0) & (dflist_in[i][1]!=13) & (dflist_in[i][2]!=0) & (dflist_in[i][2]!=13) & (dflist_in[i][3]!=0) & (dflist_in[i][3]!=13)]
        dflist_in[i]=dflist_in[i][[0,8]]
        dflist_in[i].columns = range(dflist_in[i].shape[1])
        dflist_in[i] = dflist_in[i].reset_index(drop=True)
    
    for i in range(len(dflist_label)):
        dflist_label[i]=dflist_label[i][[0,8,14]]
        dflist_label[i].columns = range(dflist_label[i].shape[1])
        #this was needed for some weird issue, where some masses are set to zero (might have been a problem with using the wrong savepoint for md)
        dflist_label[i][1].replace({0 : 1}, inplace=True)

    label_array = np.empty((0, 216), float)
    input_array = np.empty((0, 1512), float)
    
    n_scenarios = len(datasets_label)
    n_datasets = len(datasets_label[0])
    
    for i in range(1, int(dflist_in[0].iloc[[-1]][0])+1):
        for j in range(n_scenarios):
            tmp=j*n_datasets
            for k in range(n_datasets):
                input_array = np.append(input_array, np.array([dflist_in[j][dflist_in[j][0]==i].transpose().iloc[1].values.tolist()]), axis=0)
                label_array = np.append(label_array, np.true_divide(np.array([dflist_label[tmp+k][dflist_label[tmp+k][0]==i].transpose().iloc[2].values.tolist()]), np.array([dflist_label[tmp+k][dflist_label[tmp+k][0]==i].transpose().iloc[1].values.tolist()])), axis=0)

    return input_array, label_array, n_scenarios, n_datasets;

def act(x):
    return x*x*x+100*x

act_vectorize=np.vectorize(act)
act_convert= lambda x: act_vectorize(x).astype(np.float32)

def act_cubic(arg):
    return tf.convert_to_tensor(arg, dtype=tf.float32)

def compile_prob_model(size_in, size_out, n_datasets, kl_mod):
    model = tf.keras.models.Sequential([
      tf.keras.layers.InputLayer(input_shape=(size_in,), name="input"),
      tfp.layers.DenseVariational(size_out*2, posterior_mean_field, prior_trainable, kl_weight=kl_mod,
                                   kl_use_exact=True, activation="relu"),
      #tf.keras.layers.Dense(2*size_out, activation=act_cubic, kernel_initializer=tf.keras.initializers.RandomUniform(minval=1/1500, maxval=2/1500),),
      tfp.layers.DistributionLambda(lambda t: tfd.Normal(loc=t[..., :size_out], scale=tf.nn.softplus(t[...,size_out:]))),
    ])
    
    lr = 0.000005
    
    lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
        lr,
        decay_steps=100,
        decay_rate=0.95,
        staircase=False)
    
    opt=tf.keras.optimizers.Adam(
        learning_rate=lr,
        name='adam_opt'
    )

    
    negloglik = lambda y, y_pred: -y_pred.log_prob(y)
    
    model.compile(optimizer=opt,
                  loss=negloglik,
                  metrics=[tf.metrics.MeanAbsoluteError()])
    
    return model

def compile_std_model(size_in, size_out):
    
    model = tf.keras.models.Sequential([
      tf.keras.layers.InputLayer(input_shape=size_in),
      tf.keras.layers.Dense(units=size_out/4, activation="relu"),
      tf.keras.layers.Dense(units=size_out/4, activation="relu"),
      tf.keras.layers.Dense(units=size_out)
    ])

    
    loss_fn = tf.keras.losses.MeanSquaredError(reduction="auto")
    lr = 0.0000005
    
    lrs= tf.keras.optimizers.schedules.ExponentialDecay(
        lr,
        decay_steps=1000,
        decay_rate=0.96,
        staircase=False)
    
    opt=tf.keras.optimizers.Adam(
        learning_rate=lr,
        name='adam_opt'
    )
    
    model.compile(optimizer=opt,
                  loss=loss_fn,
                  metrics=[tf.metrics.MeanAbsoluteError()])
    
    return model

def prob_run(datasets_in, datasets_label, epochs, batch_sizes, kl_mod):
    if(len(ep)!=len(batch_sizes)):
        print("ep and batch_size do not match")
        return
    
    input_array, label_array, n_scenarios, n_datasets = get_real_trainig_data(datasets_in, datasets_label)

    model = compile_prob_model(np.shape(input_array)[1], np.shape(label_array)[1], n_datasets, kl_mod)
    
    timetaken = timecallback()
    lrcall=LearningRateLoggingCallback()
    
    loss=[]
    mae=[]

    for i in range(len(epochs)):
        temp=model.fit(input_array, label_array, batch_size=batch_sizes[i], epochs=epochs[i], shuffle=True,callbacks = [lrcall,timetaken])
        run_comparison_prob(model,batch_sizes[i],kl_mod)
        
        loss.extend(temp.history["loss"])
        mae.extend(temp.history["mean_absolute_error"])
    
    tf.keras.models.save_model(model, 'model_prob_s'+str(n_scenarios)+'_b'+str(n_datasets)+"_final_1")

def std_run(datasets_in, datasets_label, epochs, batch_sizes):
    input_array, label_array, n_scenarios, n_datasets = get_real_trainig_data(datasets_in, datasets_label)
    
    timetaken = timecallback()
    lrcall=LearningRateLoggingCallback()

    model = compile_std_model(np.shape(input_array)[1], np.shape(label_array)[1])
    loss=[]
    mae=[]
    kl_mod=0
    for i in range(len(epochs)):
        temp=model.fit(input_array, label_array, batch_size=batch_sizes[i], epochs=epochs[i], shuffle=True,callbacks = [lrcall,timetaken])
        run_comparison_std(model,batch_sizes[i],kl_mod)
        loss.extend(temp.history["loss"])
        mae.extend(temp.history["mean_absolute_error"])
    
    tf.keras.models.save_model(model, 'model_std_s'+str(n_scenarios)+'_b'+str(n_datasets)+"_sp")
    




#the following are 2d lists containing paths to input and label datasets
#each csv file is a writetofile output
#an input dataset can be provided multiple label datasets as indicated
#all of these need to contain the same amount of timesteps
#also each set of datasets(i.e. sets of labels assigned to the same input data) should be equally sized
#otherwise, the batch size will not function as intended
#please make sure that the indices used in get_real_trainig_data match the columns in the csv containing the momentum and mass data

datasets_in = [["../clean_data/writer2_clean.csv"],
                        #["input2"]
                        ]

datasets_label = [["../clean_data/writer1_clean.csv",
                               #"label2 for input1"
                            ],
                            #["label1 for input2",
                            #"label2 for input2"]
                        ]


#the epoch and batch size lists can contain multiple entries but need to be of equal length
epochs=[2500]
batch_sizes=[100]

#this scales the error of the prob network, shold be 1/dataset_size, i.e. total number of timesteps over all datasets
#to clarify, when providing three csv files as label data each containing 250 time steps, this should be 1/750
kl_mod=1/(3000)

#std_run(datasets_in, datasets_label, epochs, batch_sizes)
#prob_run(datasets_in, datasets_label, epochs, batch_sizes, kl_mod)



