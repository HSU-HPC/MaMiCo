#include"UnitTestingService.h"

//temporary entry point
int main(int argc, char *argv[]){
	//create object of (TODO: do we want OO here?) UnitTestingService	
	testing::ut::UnitTestingService utService;
	utService.runAllUnitTests();
}
