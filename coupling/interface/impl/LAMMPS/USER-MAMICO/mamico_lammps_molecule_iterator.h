// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef MAMICO_LMP_LAMMPS_MOLECULEITERATOR_H
#define MAMICO_LMP_LAMMPS_MOLECULEITERATOR_H

#include "coupling/interface/MoleculeIterator.h"

#include "mamico_cell.h"
#include "mamico_lammps_molecule.h"

namespace LAMMPS_NS {

/** implements a molecule iterator on the mamico cells. We work with molecule
 * information x,v,f which are associated to either a ghost buffer for molecules
 * which do not belong to the current process, or to the lmp->atom->x,v,f
 * buffers provided by Lammps.
 *  @author Philipp Neumann
 */
template <unsigned int dim>
class MamicoLammpsMoleculeIterator
    : public coupling::interface::MoleculeIterator<MamicoCell, dim> {
public:
  MamicoLammpsMoleculeIterator(MoleculeInformation info, double cutoff,
                               MamicoCell &cell)
      : coupling::interface::MoleculeIterator<MamicoCell, dim>(cell),
        _info(info), _cutoff(cutoff) {
  } virtual ~MamicoLammpsMoleculeIterator() {
  }

  virtual void begin() {
    coupling::interface::MoleculeIterator<MamicoCell, dim>::_cell.begin();
  }

  virtual void next() {
    coupling::interface::MoleculeIterator<MamicoCell, dim>::_cell.next();
  }

  virtual bool continueIteration() const {
    return coupling::interface::MoleculeIterator<MamicoCell, dim>::_cell
        .continueIteration();
  }

  virtual coupling::interface::Molecule<dim> &get() {
    _molecule = MamicoLammpsMolecule<dim>(
        _info._x, _info._v, _info._f,
        coupling::interface::MoleculeIterator<MamicoCell, dim>::_cell.get(),
        _cutoff);
    return _molecule;
  }

  virtual const coupling::interface::Molecule<dim> &
  getConst() {
    _molecule = MamicoLammpsMolecule<dim>(
        _info._x, _info._v, _info._f,
        coupling::interface::MoleculeIterator<MamicoCell, dim>::_cell.get(),
        _cutoff);
    return _molecule;
  }

  private : MoleculeInformation _info;
  const double _cutoff;
  MamicoLammpsMolecule<dim> _molecule;
};

} // namespace LAMMPS_NS
#endif // MAMICO_LMP_LAMMPS_MOLECULEITERATOR_H
