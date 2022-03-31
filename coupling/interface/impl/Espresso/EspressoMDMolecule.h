// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDMOLECULE_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDMOLECULE_H_

#include "coupling/interface/Molecule.h"
#include "domain_decomposition.hpp"
#include "energy.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "tarch/la/Vector.h"
#include "utils.hpp"

namespace coupling {
namespace interface {
class EspressoMDMolecule;
}
} // namespace coupling

/** interface for espresso molecule access.
 *  @author Rahul Arora, Philipp Neumann
 */
class coupling::interface::EspressoMDMolecule : public coupling::interface::Molecule<3> {
public:
  EspressoMDMolecule(Particle *myMolecule) : _myMolecule(myMolecule) {}
  EspressoMDMolecule() : _myMolecule(NULL) {}
  virtual ~EspressoMDMolecule() {}

  void setMolecule(Particle *newMolecule) { _myMolecule = newMolecule; }

  /** returns/ sets the velocity of the molecule
        The velocity of the particles have been scaled in Espresso with
     time_step */
  virtual tarch::la::Vector<3, double> getVelocity() const {
    tarch::la::Vector<3, double> velocity(0.0);
    for (int i = 0; i < 3; i++) {
      velocity[i] = _myMolecule->m.v[i] / time_step;
    }
    return velocity;
  }

  virtual void setVelocity(const tarch::la::Vector<3, double> &velocity) {
    for (int i = 0; i < 3; i++) {
      _myMolecule->m.v[i] = time_step * velocity[i];
    }
  }

  /** returns/ sets the position of the molecule */
  virtual tarch::la::Vector<3, double> getPosition() const {
    tarch::la::Vector<3, double> position(0.0);
    for (int i = 0; i < 3; i++) {
      position[i] = _myMolecule->r.p[i];
    }
    return position;
  }

  virtual void setPosition(const tarch::la::Vector<3, double> &position) {
    for (int i = 0; i < 3; i++) {
      _myMolecule->r.p[i] = position[i];
    }
  }

  /** sets the force acting on this molecule. This function is only called in
   * the USHER
   *  scheme so far If you want to set the force of a newly created molecule,
   *  you need to implement this function.
   */

  virtual void setForce(const tarch::la::Vector<3, double> &force) {
    for (int i = 0; i < 3; i++) {
      _myMolecule->f.f[i] = force[i];
    }
  }

  virtual tarch::la::Vector<3, double> getForce() const {
    tarch::la::Vector<3, double> force(0.0);
    for (int i = 0; i < 3; i++) {
      force[i] = _myMolecule->f.f[i];
    }
    return force;
  }

  /** returns/ sets the potential energy of the molecule */

  virtual double getPotentialEnergy() const {

    // Determine the cell in which the particle is located
    ParticleList *pl = NULL, *tmp;
    int ind, c;

    for (c = 0; c < local_cells.n; c++) {
      tmp = local_cells.cell[c];
      ind = _myMolecule - tmp->part;
      if (ind >= 0 && ind < tmp->n) {
        pl = tmp;
        break;
      }
    }

    if (!pl) {
      std::cout << "ERROR: The particle is loacted outside the simulation domian" << std::endl;
      exit(EXIT_FAILURE);
      return 0.0;
    }

    // pl is the cell which contains the particle _myMolecule
    IA_Neighbor *neighbor;
    Particle *p2;
    unsigned int np2;
    double potentialEnergy = 0.0;
    double dist2, vec21[3];

    // loop over all the neighbor cells and calculate potential energy
    // contribution from each molecule in the cell
    for (unsigned int n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      neighbor = &dd.cell_inter[c].nList[n];
      p2 = neighbor->pList->part;
      np2 = neighbor->pList->n;
      for (unsigned int j = 0; j < np2; j++) {
        dist2 = distance2vec(_myMolecule->r.p, p2[j].r.p, vec21);
        IA_parameters *ia_params = get_ia_param(_myMolecule->p.type, p2->p.type);
        Particle *temp = &(p2[j]);
        potentialEnergy += calc_non_bonded_pair_energy(_myMolecule, temp, ia_params, vec21, sqrt(dist2), dist2);
      }
    }
    return potentialEnergy;
  }
  virtual void setPotentialEnergy(const double &potentialEnergy) { _potentialEnergy = potentialEnergy; }

private:
  Particle *_myMolecule;
  double _potentialEnergy;
};

#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDMOLECULE_H_
