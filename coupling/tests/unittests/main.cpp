#include"UnitTestingService.h"

//temporary entry point
int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);

	//create object of UnitTestingService	
	testing::ut::UnitTestingService utService;
	utService.runAllUnitTests();

	MPI_Finalize();
}
