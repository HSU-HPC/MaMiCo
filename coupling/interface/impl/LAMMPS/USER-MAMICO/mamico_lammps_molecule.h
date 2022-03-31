// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef MAMICO_LMP_LAMMPS_MOLECULE_H
#define MAMICO_LMP_LAMMPS_MOLECULE_H

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/interface/Molecule.h"
#include "coupling/interface/MamicoInterfaceProvider.h"
#include "mamico_lammps_md_solver_interface.h"
#include "mamico_cell.h"

namespace LAMMPS_NS {

/** gives access to one particular LAMMPS molecule/atom.
 *  @author Philipp Neumann
 */
template <unsigned int dim>
class MamicoLammpsMolecule : public coupling::interface::Molecule<dim> {
public:
  MamicoLammpsMolecule(double **x, double **v, double **f, int n, double cutoff)
      : coupling::interface::Molecule<dim>(), _x(x), _v(v), _f(f), _n(n),
        _cutOffRadiusSquared(cutoff * cutoff), _tolerance(1.0e-8) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Allocate MamicoLammpsMolecule for molecule " << _n
              << std::endl;
#endif
  }

  // default constructor; no access possible, only null pointers
  MamicoLammpsMolecule()
      : coupling::interface::Molecule<dim>(), _x(NULL), _v(NULL), _f(NULL),
        _n(-1) {}

  virtual tarch::la::Vector<dim, double> getVelocity() const {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    if (_v == NULL) {
      std::cout << "ERROR getVelocity(): _v==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    tarch::la::Vector<dim, double> vel(0.0);
    for (int d = 0; d < dim; d++)
      vel[d] = _v[_n][d];
    return vel;
  }
  virtual void setVelocity(const tarch::la::Vector<dim, double> &vel) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    if (_v == NULL) {
      std::cout << "ERROR setVelocity(): _v==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    for (int d = 0; d < dim; d++)
      _v[_n][d] = vel[d];
  }

  virtual tarch::la::Vector<dim, double> getPosition() const {
    tarch::la::Vector<dim, double> pos(0.0);
    for (int d = 0; d < dim; d++)
      pos[d] = _x[_n][d];
    return pos;
  }
  virtual void setPosition(const tarch::la::Vector<dim, double> &pos) {
    for (int d = 0; d < dim; d++)
      _x[_n][d] = pos[d];
  }

  virtual tarch::la::Vector<dim, double> getForce() const {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    if (_f == NULL) {
      std::cout << "ERROR getForce(): _f==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    tarch::la::Vector<dim, double> force(0.0);
    for (int d = 0; d < dim; d++)
      force[d] = _f[_n][d];
    return force;
  }
  virtual void setForce(const tarch::la::Vector<dim, double> &force) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    if (_f == NULL) {
      std::cout << "ERROR setForce(): _f==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    for (int d = 0; d < dim; d++)
      _f[_n][d] = force[d];
  }

  /** copy/paste of calculateForceAndEnergy(molecule), but ignores molecule
   * interaction, if two molecules are very close to each other. We use this
   * function
   *  from the molecule-interface to evaluate the potential energy.
   */
  virtual double getPotentialEnergy() const {
    // helper variables and molecule position
    const coupling::IndexConversion<dim> &indexConversion =
        coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell,
                                                     dim>::getInstance()
            .getMacroscopicCellService()->getIndexConversion();
    const tarch::la::Vector<dim, double> position1 = getPosition();
    const tarch::la::Vector<dim, unsigned int> linkedCellInMacroscopicCell(0);
    const tarch::la::Vector<dim, unsigned int> linkedCellsPerMacroscopicCell(1);
    coupling::interface::MDSolverInterface<LAMMPS_NS::MamicoCell, dim> *
        mdSolverInterface = coupling::interface::MamicoInterfaceProvider<
            LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface();
    const tarch::la::Vector<dim, unsigned int> cellIndex =
        mdSolverInterface->getLinkedCellIndexForMoleculePosition(position1);

    // lennard jones parameters
    const double sigma6 = pow(mdSolverInterface->getMoleculeSigma(), 6.0);
    const double epsilon = mdSolverInterface->getMoleculeEpsilon();
    const double sigma6OverCutoff6 =
        sigma6 /
        (_cutOffRadiusSquared * _cutOffRadiusSquared * _cutOffRadiusSquared);
    const double cutOffEnergy =
        4.0 * epsilon * sigma6OverCutoff6 * (sigma6OverCutoff6 - 1.0);

    // set energy to zero
    double potEnergy = 0.0;

    // init start and end coordinates for neighbour search
    tarch::la::Vector<3, unsigned int> start(0);
    tarch::la::Vector<3, unsigned int> end(1);
    tarch::la::Vector<3, unsigned int> loop(0);
    for (int d = 0; d < dim; d++) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      if (cellIndex[d] == 0) {
        std::cout << "ERROR MamicoLammpsMolecule::getPotentialEnergy(): Cannot "
                     "compute energy for molecules in ghost layer!";
        std::cout << std::endl;
        exit(EXIT_FAILURE);
      }
#endif
      start[d] = cellIndex[d] - 1;
      end[d] = start[d] + 3;
    }

    for (loop[2] = start[2]; loop[2] < end[2]; loop[2]++) {
      for (loop[1] = start[1]; loop[1] < end[1]; loop[1]++) {
        for (loop[0] = start[0]; loop[0] < end[0]; loop[0]++) {
          const tarch::la::Vector<dim, unsigned int> macroscopicCellIndex =
              coupling::initDimVector<dim>(loop);
          // we take a copy instead of a reference since we do not want to have
          // conflicts
          // in the access of the underlying molecule iterators;
          // Example VTK-plotting: the Mamico-plotter iterates over cells and
          // uses references
          //                       -> using another reference dereferences the
          // molecule that
          //                          the "outer" VTK cell iterator points to
          LAMMPS_NS::MamicoCell cell = mdSolverInterface->getLinkedCell(
              macroscopicCellIndex, linkedCellInMacroscopicCell,
              linkedCellsPerMacroscopicCell, indexConversion);
          coupling::interface::MoleculeIterator<MamicoCell, dim> *it =
              mdSolverInterface->getMoleculeIterator(cell);
          for (it->begin(); it->continueIteration(); it->next()) {
            const tarch::la::Vector<dim, double> rij =
                it->getConst().getPosition() - position1;
            const double rij2 = tarch::la::dot(rij, rij);

            // compute LJ force/energy and add the contributions, only if center
            // is within cut-off radius AND
            // if this particle is not contained in this region already. We need
            // the latter condition to re-use this function
            // for the evaluation of the potential energy in the
            // molecule-wrapper (mamico_lammps_molecule)
            if ((rij2 > _cutOffRadiusSquared) || (rij2 < _tolerance)) {
            } else {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
              std::cout << "MamicoLammpsMolecule::getPotentialEnergy(): "
                           "Compute energy with molecule at "
                        << it->getConst().getPosition() << std::endl;
#endif
              const double rij6 = rij2 * rij2 * rij2;
              potEnergy += 0.5 * (4.0 * epsilon * (sigma6 / rij6) *
                                      ((sigma6 / rij6) - 1.0) - cutOffEnergy);
            }
          }
          delete it;
        }
      }
    }

    // return potential energy
    return potEnergy;
  }

  virtual void setPotentialEnergy(const double &energy) {
    // nop, since we never explicitly store the energy
  }

private:
  double **_x;
  double **_v;
  double **_f;
  int _n;
  double _cutOffRadiusSquared;
  double _tolerance;
};

}
#endif // MAMICO_LMP_LAMMPS_MOLECULE_H
