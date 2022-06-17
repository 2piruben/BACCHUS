#include "population.h"

int main(int argc, char* argv[]){
 	//std::cout<<"starting\n";
	population colony(10,"pars.in");
	//std::cout<<"initizalizing\n";
	colony.initialize_two();
	//std::cout<<"saving\n";
	colony.save_population();
	double totaltime = 8;
	double recordtimestep = 0.1;
	double nextrecordtime = recordtimestep;

	while(colony.currenttime()<totaltime){
		//std::cout<<"iteration "<< i << " \n";
		//std::cout<<"time: "
		//<<colony.currenttime()<<'\n';
		colony.evolve();
		if (colony.currenttime()>nextrecordtime){
			colony.save_population();
			nextrecordtime += recordtimestep;
		}

	}
	//std::cout<<"done\n";
	return 0;

}
// struct force{
// 	double parallel; // resulting force parallel to the axis of the bacterium
// 	double normal; // resulting normal force applied at extreme 1 of the bacterium (torque)
// }
