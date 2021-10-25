#ifndef LS1_STATIC_COMM_DATA_H
#define LS1_STATIC_COMM_DATA_H

//possibly bad practice? will phase out
#include <string>

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
            
            //data accesses
            void setConfigFilename(std::string name) {ls1ConfigFilename = name; }
            const std::string getConfigFilename() {return ls1ConfigFilename;}

        private:
            LS1StaticCommData() {}
            std::string ls1ConfigFilename;
        };
    }
}

#endif