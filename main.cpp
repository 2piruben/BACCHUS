#include "dish.h"
#include "population.h"

int main(int argc, char* argv[]){

	int nutrient_int,growth_int;

	// Defining cytoplasmic species
	Cytoplasm cyto;
	nutrient_int = cyto.add_species("nutrient_in",0.0);
	cyto.make_species_diffusible("nutrient_in",0.1); // nutrient will increase through absorption
	growth_int = cyto.add_growth_rate_modifier("growth_rate_modifier",1.0); // dummy species to control growth
	
	// Defining cytoplasmic reactions
	LinearReaction LR1(-1,nutrient_int,nutrient_int); // consumption of nutrient, it should be replace by a Hill Function at some point
	HillReaction HR1(1,1,1,nutrient_int,growth_int); // increase of growth rate with nutrient
	LinearReaction LR2(-1,growth_int,growth_int); //  natural decay of growth_rate (in the absence of nutrient) 
	
	// Adding reactions to cytoplas,
	cyto.add_reaction(&LR1);
	cyto.add_reaction(&HR1);
	cyto.add_reaction(&LR2);

	//  Initializing cells with the defined cytoplasm
	Population colony(10,"pars.in");
	colony.initialize_two(cyto);

	// Initializing dish with the created colony
	Dish dish(colony,100,10);
	dish.add_chemical("nutrient_out",10,10.0,0.001); // diffusible nutrient
	dish.set_chemical_gaussianprofile("nutrient_out",1,1,10,1); // Initial concentration of nutrient in the dish
	dish.link_chemical_bacterium("nutrient_out","nutrient_in"); // linking cytoplasmic and diffusible nutrient

	// Simulation time parameters
	double totaltime = 8; //8 
	double recordtimestep = 0.1;
	double nextrecordtime = recordtimestep;
	
	// Saving starting positions
	colony.save();
	dish.save();
	// Evolve
	while(colony.currenttime()<totaltime){
		dish.evolve();
		if (colony.currenttime()>nextrecordtime){
			colony.save();
			dish.save();
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
