#ifndef _COUPLING_NOISEREDUCTION_DATASTRUCTURES_H_
#define _COUPLING_NOISEREDUCTION_DATASTRUCTURES_H_

#include "coupling/noisereduction/NoiseReduction.h"
#include <memory>

namespace coupling {
namespace noisereduction {
/***
 * Quantities<dim>[0] := mass
 * Quantities<dim>[1..dim] := momentum
 */
template <unsigned int dim> using Quantities = tarch::la::Vector<dim + 1, double>;

/***
 *	Spacetime window of data, i.e. 4D field of T
 */
template <unsigned int dim, typename T> class Field;

/***
 *	Spacetime window of flow quantities
 */
template <unsigned int dim> using Flowfield = Field<dim, Quantities<dim>>;

/***
 *	includes local flowfield, and mean and standard deviation of quantities
 */
template <unsigned int dim> class Patch;

/***
 *	Spacetime window of patches
 */
template <unsigned int dim> using Patchfield = Field<dim, Patch<dim>>;
} // namespace noisereduction
} // namespace coupling

/***
 *	Spacetime window of data, i.e. 4D field of T
 */
template <unsigned int dim, typename T> class coupling::noisereduction::Field {
public:
  Field(const tarch::la::Vector<dim, unsigned int> &spatialSize, const unsigned int &temporalSize)
      : _spatialSize(spatialSize), _temporalSize(temporalSize), _scalarSize(computeScalarSize(spatialSize, temporalSize)),
        _data(std::allocator_traits<std::allocator<T>>::allocate(allo, _scalarSize)) {}

  template <class... Args> void construct(const tarch::la::Vector<dim, unsigned int> &pos, const unsigned int &t, Args &&...args) {
    std::allocator_traits<std::allocator<T>>::construct(allo, _data + idx(pos, t), std::forward<Args>(args)...);
  }

  void destroy(const tarch::la::Vector<dim, unsigned int> &pos, const unsigned int &t) {
    std::allocator_traits<std::allocator<T>>::destroy(allo, _data + idx(pos, t));
  }

  T &operator()(const tarch::la::Vector<dim, unsigned int> &pos, const unsigned int &t) { return _data[idx(pos, t)]; }

  const T &operator()(const tarch::la::Vector<dim, unsigned int> &pos, const unsigned int &t) const { return _data[idx(pos, t)]; }

  T &operator[](unsigned int pos) {
#ifdef NLM_DEBUG
    if (pos < 0 || pos > _scalarSize - 1) {
      std::cout << "ERROR Field T& operator[](int pos): pos out of range!" << std::endl;
      std::cout << "pos=" << pos << ", _scalarSize=" << _scalarSize << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    return _data[pos];
  }

  const T &operator[](unsigned int pos) const {
#ifdef NLM_DEBUG
    if (pos < 0 || pos > _scalarSize - 1) {
      std::cout << "ERROR Field T& operator[](int pos): pos out of range!" << std::endl;
      std::cout << "pos=" << pos << ", _scalarSize=" << _scalarSize << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    return _data[pos];
  }

  ~Field() { std::allocator_traits<std::allocator<T>>::deallocate(allo, _data, _scalarSize); }

  unsigned int getScalarSize() const { return _scalarSize; }

  const tarch::la::Vector<dim, unsigned int> &getSpatialSize() const { return _spatialSize; }

  unsigned int getTemporalSize() const { return _temporalSize; }

private:
  unsigned int computeScalarSize(const tarch::la::Vector<dim, unsigned int> &spatialSize, const unsigned int &temporalSize) const {
    unsigned int res = spatialSize[0];
    for (unsigned int d = 1; d < dim; d++) {
      res *= (spatialSize[d]);
    }
    res *= temporalSize;
    return res;
  }

  unsigned int idx(const tarch::la::Vector<dim, unsigned int> &pos, const unsigned int &t) const {
#ifdef NLM_DEBUG
    for (unsigned int d = 0; d < dim; d++) {
      if (pos[d] < 0 || pos[d] > _spatialSize[d] - 1) {
        std::cout << "ERROR Field idx(): pos out of range!" << std::endl;
        std::cout << "pos=" << pos << ", _spatialSize=" << _spatialSize << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if (t < 0 || t > _temporalSize - 1) {
      std::cout << "ERROR Field idx(): t out of range!" << std::endl;
      std::cout << "t=" << t << ", _temporalSize=" << _temporalSize << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    unsigned int idx = 0, step = 1;
    for (unsigned int d = 0; d < dim; d++) {
      idx += pos[d] * step;
      step *= _spatialSize[d];
    }
    idx += t * step;
    return idx;
  }

  static std::allocator<T> allo;
  const tarch::la::Vector<dim, unsigned int> _spatialSize;
  const unsigned int _temporalSize;
  const unsigned int _scalarSize;
  T *const _data;

  friend class coupling::noisereduction::Patch<dim>; // patches need direct access to their field's _data for memory
                                                     // efficiency
};

template <unsigned int dim, typename T> std::allocator<T> coupling::noisereduction::Field<dim, T>::allo;

/***
 *	includes local flowfield, and mean and standard deviation of quantities
 */
template <unsigned int dim> class coupling::noisereduction::Patch {
public:
  Patch(const tarch::la::Vector<dim, unsigned int> &spatialSize, const unsigned int &temporalSize, const Flowfield<dim> &basefield,
        const tarch::la::Vector<dim, unsigned int> &pos, const unsigned int &t)
      : _flowfield(spatialSize, temporalSize), _localMean(0.0), _localStandardDeviation(0.0) {
    fillFromBasefield(basefield, pos, t);
    computeLocalMean();
    computeLocalStandardDeviation();
  }

  double distance(const coupling::noisereduction::Patch<dim> &other) {
    unsigned int size = _flowfield.getScalarSize() * (dim + 1);

    double *const my_data = reinterpret_cast<double *const>(_flowfield._data);
    double *const other_data = reinterpret_cast<double *const>(other._flowfield._data);

    double res = 0;
    for (unsigned int i = 0; i < size; i++) {
      double diff(my_data[i] - other_data[i]);
      res += diff * diff;
    }
    return res;
  }

  ~Patch() {}

private:
  inline unsigned int posmod(int i, int n) { return (i % n + n) % n; }

  void fillFromBasefield(const Flowfield<dim> &basefield, const tarch::la::Vector<dim, unsigned int> &pos, const unsigned int &t) {
    static_assert(dim == 2 || dim == 3, "ERROR noisereduction::Patch::fillFromBasefield only "
                                        "implemented for 2D and 3D");

    tarch::la::Vector<dim, unsigned int> center;
    for (unsigned int d = 0; d < dim; d++) {
      center[d] = _flowfield.getSpatialSize()[d] / 2;
    }

    tarch::la::Vector<dim, unsigned int> local_pos(0);
    unsigned int local_t(0);
    tarch::la::Vector<dim, unsigned int> base_pos(0);
    unsigned int base_t(0);

    for (int offset_t = 0; offset_t > -(int)_flowfield.getTemporalSize(); offset_t--) {
      local_t = posmod(offset_t, _flowfield.getTemporalSize());
      base_t = posmod(t + offset_t, basefield.getTemporalSize());
      for (int offset_x = -(int)(_flowfield.getSpatialSize()[0]) / 2; offset_x <= (int)(_flowfield.getSpatialSize()[0]) / 2; offset_x++) {
        local_pos[0] = center[0] + offset_x;
        base_pos[0] = pos[0] + offset_x;
        for (int offset_y = -(int)(_flowfield.getSpatialSize()[1]) / 2; offset_y <= (int)(_flowfield.getSpatialSize()[1]) / 2; offset_y++) {
          local_pos[1] = center[1] + offset_y;
          base_pos[1] = pos[1] + offset_y;
          if constexpr (dim == 3) {
            for (int offset_z = -(int)(_flowfield.getSpatialSize()[2]) / 2; offset_z <= (int)(_flowfield.getSpatialSize()[2]) / 2; offset_z++) {
              local_pos[2] = center[2] + offset_z;
              base_pos[2] = pos[2] + offset_z;
              _flowfield(local_pos, local_t) = basefield(base_pos, base_t);
            }
          } else {
            _flowfield(local_pos, local_t) = basefield(base_pos, base_t);
          }
        }
      }
    }
  }

  void computeLocalMean() {
    for (unsigned int i = 0; i < _flowfield.getScalarSize(); i++)
      _localMean += _flowfield[i];
    _localMean = _localMean * (1.0 / _flowfield.getScalarSize());
  }

  void computeLocalStandardDeviation() {
    for (unsigned int i = 0; i < _flowfield.getScalarSize(); i++) {
      Quantities<dim> diff = _localMean - _flowfield[i];
      _localStandardDeviation += tarch::la::dot(diff, diff);
    }
    _localStandardDeviation = sqrt(_localStandardDeviation / _flowfield.getScalarSize());
  }

  Flowfield<dim> _flowfield;
  Quantities<dim> _localMean;
  double _localStandardDeviation;
};

#endif // _COUPLING_NOISEREDUCTION_DATASTRUCTURES_H_
