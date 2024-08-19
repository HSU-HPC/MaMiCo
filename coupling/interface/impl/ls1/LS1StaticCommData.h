#ifndef LS1_STATIC_COMM_DATA_H_
#define LS1_STATIC_COMM_DATA_H_

#include "coupling/CouplingMDDefinitions.h"
#include <string>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

namespace coupling {
namespace interface {
class LS1StaticCommData {
public:
  static LS1StaticCommData& getInstance() {
    static LS1StaticCommData singleton;
    return singleton;
  }
  LS1StaticCommData(LS1StaticCommData const&) = delete;
  void operator=(LS1StaticCommData const&) = delete;

  // data sets and gets
  void setConfigFilename(std::string name) { _ls1ConfigFilename = name; }
  const std::string getConfigFilename() { return _ls1ConfigFilename; }

  void setBoxOffsetAtDim(int dim, double offset) { _boxoffset[dim] = offset; }
  const double getBoxOffsetAtDim(int dim) { return _boxoffset[dim]; }

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  void setLocalCommunicator(MPI_Comm comm) { _localComm = comm; }
  MPI_Comm getLocalCommunicator() { return _localComm; }

  void setDomainGridDecompAtDim(int dim, int breakdown) { _domainGridDecomp[dim] = breakdown; }
  const int getDomainGridDecompAtDim(int dim) { return _domainGridDecomp[dim]; }
  std::array<int, 3> getDomainGridDecomp() { return _domainGridDecomp; }

  void setSubdomainWeights(std::array<std::vector<unsigned int>, 3> subdomainWeights) { _subdomainWeights = subdomainWeights; }
  const std::array <std::vector<unsigned int> ,3> getSubdomainWeights() { return _subdomainWeights; }
#endif

private:
  LS1StaticCommData() {}
  std::string _ls1ConfigFilename;
  std::array<double, 3> _boxoffset; // temporary till ls1 offset is natively supported
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Comm _localComm;
  std::array<int, 3> _domainGridDecomp;
  std::array<std::vector<unsigned int> ,3 >_subdomainWeights;
#endif
};
} // namespace interface
} // namespace coupling

#endif