// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_VELOCITYGRADIENTRELAXATIONMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_VELOCITYGRADIENTRELAXATIONMAPPING_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"
#include "tarch/la/Matrix.h"
#include <iostream>

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class VelocityGradientRelaxationMapping;
template <class LinkedCell, unsigned int dim> class VelocityGradientRelaxationTopOnlyMapping;
} // namespace cellmappings
} // namespace coupling

/** This class relaxes velocities of molecules towards an interpolated avergaged
 *velocity value. Only molecules within a three cell-wide boundary strip are
 *considered if they are inside the one cell-wide boundary strip that is spanned
 *by the outermost non-ghost coupling cell midpoints. Currently, a
 *second-order interpolation of the molecules in the outer boundary strip is
 *used (that is the first layer with three coupling cells needs to be
 *initialised with valid LB velocities).
 *	@brief This class relaxes velocities of molecules towards an
 *interpolated avergaged velocity value.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 *  @note ONLY SUPPORTS 3D AT THE MOMENT!
 *	\todo Philipp please take a look on this class
 */
template <class LinkedCell, unsigned int dim> class coupling::cellmappings::VelocityGradientRelaxationMapping {
public:
  /** obtains relaxation factor and current velocity (the velocity towards which
   *the relaxation shall be done is later obtained from the
   *microscopicmomentum-buffers).
   *	@param velocityRelaxationFactor
   *	@param currentVelocity
   *	@param cellIndex
   *	@param mdSolverInterface
   *	@param couplingCells
   */
  VelocityGradientRelaxationMapping(const double& velocityRelaxationFactor, const tarch::la::Vector<dim, double>& currentVelocity, const I02& cellIndex,
                                    coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface,
                                    const coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>* const couplingCells)
      : _mdSolverInterface(mdSolverInterface), _couplingCells(couplingCells), _velocityRelaxationFactor(velocityRelaxationFactor),
        _currentVelocity(currentVelocity), _innerLowerLeft(getInnerLowerLeftCorner()), _innerUpperRight(getInnerUpperRightCorner()),
        _outerLowerLeft(getOuterLowerLeftCorner()), _outerUpperRight(getOuterUpperRightCorner()), _cellIdx(cellIndex),
        _ignoreThisCell(ignoreThisCell(cellIndex)) {}

  /** Destructor */
  virtual ~VelocityGradientRelaxationMapping() {}

  /** empty function
   */
  void beginCellIteration() {}

  /** empty function
   */
  void endCellIteration() {}

  /** computes the current velocity directly from coupling cell and the new
   *velocity (microscopicMomentum) with second-order interpolation multiplies
   *the difference between the two with the velocity relaxation factor add it to
   *the velocity of the molecule and applies it to the molecule.
   *	@param cell
   *	@param cellIndex
   */
  void handleCell(LinkedCell& cell) {
    // if this coupling cell is not interesting, skip it immediately
    if (_ignoreThisCell) {
      return;
    }

    // std::cout << "Handle VelocityGradientRelaxation in cell " <<
    // I03{_cellIdx}.get() << "=" <<  _cellIdx <<
    // std::endl;
    coupling::interface::MoleculeIterator<LinkedCell, dim>* it = _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      coupling::interface::Molecule<dim>& wrapper(it->get());
      tarch::la::Vector<dim, double> velocity = wrapper.getVelocity();
      const tarch::la::Vector<dim, double> position(wrapper.getPosition());
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      std::cout << "VelocityGradientRelaxationMapping: Process cell " << _cellIdx << ", molecule " << position << std::endl;
#endif
      if (relaxMolecule(position)) {
        tarch::la::Vector<dim, double> normalisedPosition(0.0);
        bool secondCake = false;
        // index of coupling cell values that are involved in interpolation
        // process
        I02 couplingCellIndex[10];
        // determine cell indices and the "correct cake" for interpolation
        createCouplingCellIndex4SecondOrderInterpolation(position, normalisedPosition, secondCake, couplingCellIndex);

        tarch::la::Vector<dim, double> newVelocity(0.0);
        tarch::la::Vector<dim, double> oldVelocity(0.0);

        // get newVelocity (microscopicMomentum) with second-order interpolation
        interpolateVelocitySecondOrder(couplingCellIndex, normalisedPosition, secondCake, newVelocity);
        // get current velocity directly from coupling cell
        oldVelocity = _couplingCells[_cellIdx.get()].getCurrentVelocity();

        velocity += _velocityRelaxationFactor * (newVelocity - oldVelocity);
        wrapper.setVelocity(velocity);

#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
        for (unsigned int d = 0; d < dim; d++) {
          if (normalisedPosition[d] > 2.0) {
            std::cout << "ERROR cellmappings::VelocityGradientRelaxationMapping: "
                         "normalisedPosition>2.0! "
                      << "Local vector index=" << I03{_cellIdx} << " , global vector index=" << I01{_cellIdx} << ", currentLocalCellIndex=" << _cellIdx
                      << ", Norm. position=" << normalisedPosition << ", Actual position=" << position << std::endl;
            for (unsigned int i = 0; i < 10; i++) {
              std::cout << couplingCellIndex[i] << " ";
            }
            std::cout << std::endl;
            exit(EXIT_FAILURE);
          }
          if (normalisedPosition[d] < 0.0) {
            std::cout << "ERROR cellmappings::VelocityGradientRelaxationMapping: "
                         "normalisedPosition<0.0! "
                      << I03{_cellIdx} << " ," << normalisedPosition << std::endl;
            for (unsigned int i = 0; i < 10; i++) {
              std::cout << couplingCellIndex[i] << " ";
            }
            std::cout << std::endl;
            exit(EXIT_FAILURE);
          }
        }
        if ((tarch::la::dot(newVelocity, newVelocity) > 1000.0) || (tarch::la::dot(oldVelocity, oldVelocity) > 1000.0)) {
          for (unsigned int i = 0; i < 10; i++) {
            std::cout << couplingCellIndex[i] << " ";
          }
          std::cout << std::endl;
          std::cout << "ERROR: velocity to high! " << oldVelocity << "  " << newVelocity << "  " << I03{_cellIdx} << std::endl;
          exit(EXIT_FAILURE);
        }
#endif
      }
      it->next();
    }
    delete it;
  }

protected:
  /**
   *	@param position
   *	@returns true, if the position 'position' is within the respective
   *boundary strip that is spanned by the two outermost coupling cell
   *midpoints
   */
  virtual bool relaxMolecule(const tarch::la::Vector<dim, double>& position) const {
    // check if molecule is in the respective boundary strip
    bool isInsideOuter = true;
    bool isOutsideInner = false;
    for (unsigned int d = 0; d < dim; d++) {
      isInsideOuter = isInsideOuter && (position[d] > _outerLowerLeft[d]) && (position[d] < _outerUpperRight[d]);
      isOutsideInner = isOutsideInner || (position[d] < _innerLowerLeft[d]) || (position[d] > _innerUpperRight[d]);
    }
    return isInsideOuter && isOutsideInner;
  }

  /**
   *	@param idx
   *	@returns true if all molecules in this coupling cell do not require
   *any velocity relaxation (-> used to speed up the code)
   */
  bool ignoreThisCell(const I02& idx) const {
    const tarch::la::Vector<dim, unsigned int> globalIndex{I01{idx}.get()};
    bool innerCell = true;
    for (unsigned int d = 0; d < dim; d++) {
      innerCell = innerCell && ((globalIndex[d] > 2) && (globalIndex[d] < I01::numberCellsInDomain[d] - 3));
    }
    return innerCell;
  }

  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface;
  const coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>* const _couplingCells;
  /** relaxation factor for velocity relaxation */
  const double _velocityRelaxationFactor;
  /** current velocity in this coupling cell */
  const tarch::la::Vector<dim, double> _currentVelocity;

  /** boundary points of the boundary strip under consideration */
  const tarch::la::Vector<dim, double> _innerLowerLeft;
  const tarch::la::Vector<dim, double> _innerUpperRight;
  const tarch::la::Vector<dim, double> _outerLowerLeft;
  const tarch::la::Vector<dim, double> _outerUpperRight;

  I02 _cellIdx;
  /** true, if all molecules in this cell do not require any velocity relaxation */
  bool _ignoreThisCell;

private:
  /**
   *	@returns the inner lower left corner of the cell
   */
  tarch::la::Vector<dim, double> getInnerLowerLeftCorner() const { return IDXS.getGlobalMDDomainOffset() + 1.5 * IDXS.getCouplingCellSize(); }
  /**
   *	@returns the inner upper right corner of the cell
   */
  tarch::la::Vector<dim, double> getInnerUpperRightCorner() const {
    return IDXS.getGlobalMDDomainOffset() + IDXS.getGlobalMDDomainSize() - 1.5 * IDXS.getCouplingCellSize();
  }
  /**
   *	@returns the outer lower left corner of the cell
   */
  tarch::la::Vector<dim, double> getOuterLowerLeftCorner() const { return IDXS.getGlobalMDDomainOffset() + 0.5 * IDXS.getCouplingCellSize(); }
  /**
   *	@returns the outer upper right corner of the cell
   */
  tarch::la::Vector<dim, double> getOuterUpperRightCorner() const {
    return IDXS.getGlobalMDDomainOffset() + IDXS.getGlobalMDDomainSize() - 0.5 * IDXS.getCouplingCellSize();
  }

  /** creates a coupling cell index for the second order interpolation.
   *	@param position
   *	@param normalisedPosition
   *	@param secondCake
   *	@param res
   */
  void createCouplingCellIndex4SecondOrderInterpolation(const tarch::la::Vector<dim, double>& position, tarch::la::Vector<dim, double>& normalisedPosition,
                                                        bool& secondCake, I02* res) const {
    // compute lower left cell index and normalised position;
    // the normalised position is chosen such that the origin (0,0...,0)
    // coincides with the lower left... coupling cell's midpoint

    tarch::la::Vector<dim, double> cellMidpoint = _cellIdx.getCellMidPoint();
    auto cellSize = IDXS.getCouplingCellSize();
    tarch::la::Vector<dim, unsigned int> globalCellIndex{I01{_cellIdx}.get()};
    tarch::la::Vector<dim, unsigned int> lowerLeftCellIndex{I03{_cellIdx}.get()};

    // determine cell mid point and normalised position in this loop
    for (unsigned int d = 0; d < dim; d++) {
      // determine mid point of current coupling cell

#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      if ((position[d] < cellMidpoint[d] - 0.5 * cellSize[d]) || (position[d] > cellMidpoint[d] + 0.5 * cellSize[d])) {
        std::cout << "ERROR "
                     "VelocityGradientRelaxationMapping::"
                     "createCouplingCellIndex4SecondOrderInterpolation: "
                     "Input position "
                  << position;
        std::cout << " out of range of cell " << globalCellIndex << "!" << std::endl;
        exit(EXIT_FAILURE);
      }
#endif

      // put lowerLeftCellIndex such that it's right below the position
      if (position[d] < cellMidpoint[d]) {
        lowerLeftCellIndex[d]--;
        globalCellIndex[d]--;
      }
      // if this is on the left/ lower boundary, but still to far up, reduce
      // height
      if (lowerLeftCellIndex[d] == 1) {
        lowerLeftCellIndex[d]--;
        globalCellIndex[d]--;
      } else if (lowerLeftCellIndex[d] == I03::numberCellsInDomain[d] - 2) {
        lowerLeftCellIndex[d]--;
        globalCellIndex[d]--;
      }

      // determine normalised position
      normalisedPosition[d] = (position[d] - (IDXS.getGlobalMDDomainOffset()[d] - 0.5 * cellSize[d]) - globalCellIndex[d] * cellSize[d]) / cellSize[d];
    }
    unsigned int lowerLeftIndex = I02{I03{tarch::la::Vector<3, int>{lowerLeftCellIndex}}}.get();

    if (dim == 2) {
      std::cout << "Not implemented correctly yet!" << std::endl;
      exit(EXIT_FAILURE);
    } else if (dim == 3) {
      const unsigned int localNumberCells01 = I03::numberCellsInDomain[0] * I03::numberCellsInDomain[1];
      const unsigned int localNumber2Cells0 = 2 * I03::numberCellsInDomain[0];
      const unsigned int localNumber2Cells01Plus2Cells0 = 2 * localNumberCells01 + 2 * I03::numberCellsInDomain[0];

      tarch::la::Vector<dim, double> normal;
      normal[0] = 1.0;
      normal[1] = 1.0;
      normal[2] = 0.0;
      tarch::la::Vector<dim, double> relPos(normalisedPosition);
      relPos[0] -= 2.0;

      // THE CAKE IS A LIE

      // first cake:
      if (tarch::la::dot(relPos, normal) < 0.0) {
        secondCake = false;
        res[0] = I02{lowerLeftIndex};
        res[1] = I02{lowerLeftIndex + 1};
        res[2] = I02{lowerLeftIndex + 2};
        res[3] = I02{lowerLeftIndex + I03::numberCellsInDomain[0]};
        res[4] = I02{lowerLeftIndex + localNumber2Cells0};
        res[5] = I02{lowerLeftIndex + localNumberCells01};
        res[6] = I02{lowerLeftIndex + localNumberCells01 + I03::numberCellsInDomain[0] + 1};
        res[7] = I02{lowerLeftIndex + 2 * localNumberCells01};
        res[8] = I02{lowerLeftIndex + 2 * localNumberCells01 + 2};
        res[9] = I02{lowerLeftIndex + localNumber2Cells01Plus2Cells0};
        // second cake
      } else {
        secondCake = true;
        res[0] = I02{lowerLeftIndex + 2};
        res[1] = I02{lowerLeftIndex + I03::numberCellsInDomain[0] + 2};
        res[2] = I02{lowerLeftIndex + localNumber2Cells0};
        res[3] = I02{lowerLeftIndex + localNumber2Cells0 + 1};
        res[4] = I02{lowerLeftIndex + localNumber2Cells0 + 2};
        res[5] = I02{lowerLeftIndex + localNumberCells01 + I03::numberCellsInDomain[0] + 1};
        res[6] = I02{lowerLeftIndex + localNumberCells01 + 2 * I03::numberCellsInDomain[0] + 2};
        res[7] = I02{lowerLeftIndex + 2 * localNumberCells01 + 2};
        res[8] = I02{lowerLeftIndex + localNumber2Cells01Plus2Cells0};
        res[9] = I02{lowerLeftIndex + localNumber2Cells01Plus2Cells0 + 2};
      }
    } else {
      std::cout << "ERROR createCouplingCellIndex4Interpolation(position): "
                   "only 2D/3D supported!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  /** carries out second order interpolation of the velocity value
   *(microscopicmomentum-buffer) at position 'position'
   *	@param cellIndex
   *	@param normalisedPosition
   *	@param secondCake
   *	@param newVelocity
   */
  void interpolateVelocitySecondOrder(I02* cellIndex, const tarch::la::Vector<dim, double>& normalisedPosition, const bool& secondCake,
                                      tarch::la::Vector<dim, double>& newVelocity) const {
    if (dim == 2) {
      // TODO
    } else if (dim == 3) {
      // compute interpolation coefficients
      tarch::la::Vector<dim, double> coefficients[10];
      tarch::la::Vector<dim, double> velocity[10];

      // extract velocities from cells
      for (int i = 0; i < 10; i++) {
        velocity[i] = _couplingCells[cellIndex[i].get()].getMicroscopicMomentum();
      }

      if (secondCake) {
        const tarch::la::Vector<dim, double> v78 = velocity[7] + velocity[8];
        coefficients[0] =
            2.0 * (velocity[0] + velocity[2] + velocity[9]) + 4.0 * (velocity[5] - velocity[1] - velocity[6] - velocity[3]) + 5.0 * velocity[4] - v78;
        coefficients[1] =
            0.5 * (v78 - velocity[0]) + 2.0 * (velocity[1] - velocity[2] - velocity[5] + velocity[6]) + 4.0 * velocity[3] - 3.5 * velocity[4] - velocity[9];
        coefficients[2] =
            2.0 * (velocity[3] - velocity[0] - velocity[5] + velocity[6]) + 4.0 * velocity[1] - 3.5 * velocity[4] + 0.5 * (v78 - velocity[2]) - velocity[9];
        coefficients[3] = 0.5 * (v78 - velocity[0] - velocity[2] - velocity[4]) + 2.0 * velocity[6] - 1.5 * velocity[9];
        coefficients[4] = 0.5 * (velocity[2] + velocity[4]) - velocity[3];
        coefficients[5] = 0.5 * (velocity[0] + velocity[4]) - velocity[1];
        coefficients[6] = 0.5 * (velocity[4] + velocity[9]) - velocity[6];
        coefficients[7] =
            0.25 * (velocity[0] + velocity[2] - v78) - velocity[1] - velocity[3] + 1.5 * velocity[4] + velocity[5] - velocity[6] + 0.5 * velocity[9];
        coefficients[8] = 0.25 * (velocity[2] - velocity[4] - velocity[8] + velocity[9]);
        coefficients[9] = 0.25 * (velocity[0] - velocity[4] - velocity[7] + velocity[9]);
      } else {
        const tarch::la::Vector<dim, double> minus15V0 = -1.5 * velocity[0];
        coefficients[0] = velocity[0];
        coefficients[1] = minus15V0 + 2.0 * velocity[1] - 0.5 * velocity[2];
        coefficients[2] = minus15V0 + 2.0 * velocity[3] - 0.5 * velocity[4];
        coefficients[3] = minus15V0 + 2.0 * velocity[5] - 0.5 * velocity[7];
        coefficients[4] = 0.5 * (velocity[0] + velocity[2]) - velocity[1];
        coefficients[5] = 0.5 * (velocity[0] + velocity[4]) - velocity[3];
        coefficients[6] = 0.5 * (velocity[0] + velocity[7]) - velocity[5];
        coefficients[7] = 1.5 * velocity[0] + velocity[6] - velocity[1] - velocity[3] - velocity[5] + 0.5 * velocity[7] +
                          0.25 * (velocity[2] + velocity[4] - velocity[8] - velocity[9]);
        coefficients[8] = 0.25 * (velocity[0] - velocity[2] - velocity[7] + velocity[8]);
        coefficients[9] = 0.25 * (velocity[0] - velocity[4] - velocity[7] + velocity[9]);
      }

      // evaluate polynomial
      const tarch::la::Vector<dim, double> contribution0 =
          normalisedPosition[0] *
          (coefficients[1] + coefficients[4] * normalisedPosition[0] + coefficients[7] * normalisedPosition[1] + coefficients[8] * normalisedPosition[2]);
      const tarch::la::Vector<dim, double> contribution1 =
          normalisedPosition[1] * (coefficients[2] + coefficients[5] * normalisedPosition[1] + coefficients[9] * normalisedPosition[2]);
      const tarch::la::Vector<dim, double> contribution2 = normalisedPosition[2] * (coefficients[3] + coefficients[6] * normalisedPosition[2]);

      newVelocity = coefficients[0] + contribution0 + contribution1 + contribution2;
    }
  }
};

/** This is the same as the class
 *coupling::cellmappings::VelocityGradientRelaxationMapping, but relaxes only
 *those molecules which are located in the top boundary strip. derived from the
 *class coupling::cellmappings::VelocityGradientRelaxationMapping.
 *	@brief This is the same as the class
 *coupling::cellmappings::VelocityGradientRelaxationMapping, but relaxes only
 *those molecules which are located in the top boundary strip.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 *  @note ONLY SUPPORTS 3D AT THE MOMENT!
 *	\todo Philipp please take a look on this class
 */
template <class LinkedCell, unsigned int dim>
class coupling::cellmappings::VelocityGradientRelaxationTopOnlyMapping : public coupling::cellmappings::VelocityGradientRelaxationMapping<LinkedCell, dim> {
public:
  /** obtains relaxation factor and current velocity (the velocity towards which
   *the relaxation shall be done is later obtained from the
   *microscopicmomentum-buffers).
   *	@param velocityRelaxationFactor
   *	@param currentVelocity
   *	@param cellIndex
   *	@param mdSolverInterface
   *	@param couplingCells
   *	@note  the method ignoreThisCell() of the class
   *coupling::cellmappings::VelocityGradientRelaxationMapping is replaceed.
   *Since this function is called in the constructor of the based class, we
   *cannot overwrite it; hence, we solve it by overwriting the respective
   *variable from within this constructor (called after base object is
   *constructed)
   */
  VelocityGradientRelaxationTopOnlyMapping(const double& velocityRelaxationFactor, const tarch::la::Vector<dim, double>& currentVelocity, const I02& cellIndex,
                                           coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface,
                                           const coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>* const couplingCells)
      : coupling::cellmappings::VelocityGradientRelaxationMapping<LinkedCell, dim>(velocityRelaxationFactor, currentVelocity, cellIndex, mdSolverInterface,
                                                                                   couplingCells) {

    // the following snippet is basically a replacement of the method
    // ignoreThisCell(). Since this function is called in the constructor of the
    // based class, we cannot overwrite it; hence, we solve it by overwriting
    // the respective variable from within this constructor (called after base
    // object is constructed)
    coupling::cellmappings::VelocityGradientRelaxationMapping<LinkedCell, dim>::_ignoreThisCell =
        (I01{cellIndex}.get()[dim - 1] < (int)(I01::numberCellsInDomain[dim - 1]) - 3);
  }

protected:
  /** returns true, if the position 'position' is within the respective boundary
   *strip that is spanned by the two outermost coupling cell midpoints
   *	@param position
   *  @return true, if the position 'position' is within the respective boundary
   *strip that is spanned by the two outermost coupling cell midpoints
   */
  virtual bool relaxMolecule(const tarch::la::Vector<dim, double>& position) const {
    // check if molecule is in the respective boundary strip (only upper part is
    // considered)
    return (position[dim - 1] > coupling::cellmappings::VelocityGradientRelaxationMapping<LinkedCell, dim>::_innerUpperRight[dim - 1]) &&
           (position[dim - 1] < coupling::cellmappings::VelocityGradientRelaxationMapping<LinkedCell, dim>::_outerUpperRight[dim - 1]);
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_VELOCITYGRADIENTRELAXATIONMAPPING_H_
