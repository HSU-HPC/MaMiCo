#ifndef LS1_STATIC_COMM_DATA_H_
#define LS1_STATIC_COMM_DATA_H_

//possibly bad practice? will phase out
#include <string>
#include <iostream>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

namespace coupling
{
    namespace interface
    {
        class LS1StaticCommData
        {
        public:
            static LS1StaticCommData& getInstance()
            {
                static LS1StaticCommData singleton;
                return singleton;
            }
            LS1StaticCommData (LS1StaticCommData const&) = delete;
            void operator=(LS1StaticCommData const&) = delete;
            
            //data sets and gets
            void setConfigFilename(std::string name) {ls1ConfigFilename = name;}
            const std::string getConfigFilename() {return ls1ConfigFilename;}

            void setBoxOffsetAtDim(int dim, double offset) {boxoffset[dim] = offset;}
            const double getBoxOffsetAtDim(int dim) {return boxoffset[dim];}

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
            void setLocalCommunicator(MPI_Comm comm) {localComm = comm;}
            MPI_Comm getLocalCommunicator() {return localComm;}

            void setDomainGridDecompAtDim(int dim, int breakdown) {domainGridDecomp[dim] = breakdown;}
            const int getDomainGridDecompAtDim(int dim) {return domainGridDecomp[dim];}
            std::array<int,3> getDomainGridDecomp() {return domainGridDecomp;}
#endif

        private:
            LS1StaticCommData() {}
            std::string ls1ConfigFilename;
            std::array<double,3> boxoffset; //temporary till ls1 offset is natively supported
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
            MPI_Comm localComm;
            std::array<int,3> domainGridDecomp;
#endif
        };
    }
}

#endif