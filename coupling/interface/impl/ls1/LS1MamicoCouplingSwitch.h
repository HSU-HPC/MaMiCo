#ifndef LS1_MAMICO_COUPLING_SWITCH_H_
#define LS1_MAMICO_COUPLING_SWITCH_H_


namespace coupling
{
    namespace interface
    {
        class LS1MamicoCouplingSwitch
        {
        public:
            static LS1MamicoCouplingSwitch& getInstance()
            {
                static LS1MamicoCouplingSwitch singleton;
                return singleton;
            }
            LS1MamicoCouplingSwitch(LS1MamicoCouplingSwitch const&) = delete;
            void operator=(LS1MamicoCouplingSwitch const&) = delete;

            //boolmanip methods
            void setCouplingStateOn(){ _couplingState = true; }
            void setCouplingStateOff(){ _couplingState = false; }
            bool getCouplingState(){ return _couplingState; }
        private:
            LS1MamicoCouplingSwitch() {}
            bool _couplingState = false;
        };
    } //namespace interface
}//namespace coupling
#endif