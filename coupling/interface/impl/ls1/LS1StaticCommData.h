#ifndef LS1_STATIC_COMM_DATA_H_
#define LS1_STATIC_COMM_DATA_H_

//possibly bad practice? will phase out
#include <string>
#include <iostream>

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
            void setConfigFilename(std::string name) {ls1ConfigFilename = name; std::cout << "setConfigFilename" << name << std::endl;}
            const std::string getConfigFilename() {return ls1ConfigFilename;}

        private:
            LS1StaticCommData() {}
            std::string ls1ConfigFilename;
        };
    }
}

#endif