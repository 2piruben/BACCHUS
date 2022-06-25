#include "dish.h"
#include "population.h"

int main(int argc, char* argv[]){
 	//std::cout<<"starting\n";
	Population colony(10,"pars.in");
	Cytoplasm cyto;
	cyto.add_species(1,"mRNA");
	LinearReaction LR1(-1,0,0);
	cyto.add_reaction(&LR1);
	//std::cout<<"initizalizing\n";
	colony.initialize_two(cyto);
	//std::cout<<"saving\n";
	colony.save_population();
	Dish dish(colony,100);

	double totaltime = 2; //8 
	double recordtimestep = 0.1;
	double nextrecordtime = recordtimestep;

	while(colony.currenttime()<totaltime){
		dish.evolve();
		if (colony.currenttime()>nextrecordtime){
			colony.save_population();
			nextrecordtime += recordtimestep;
		}

	}
	return 0;

}




// int main(int argc, char* argv[]){
// 	HillRepReaction HR1(1,1,1,0,0);
// 	Cytoplasm cyto;
// 	cyto.add_species(1,"mRNA");
// 	cyto.add_reaction(&HR1);
// 	cyto.print();
// 	cyto.react(0.1);
// 	cyto.print();
// 	return 0;

// }
