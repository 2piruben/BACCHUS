#include "dish.h"
#include "population.h"

int main(int argc, char* argv[]){

	int nutrient_int,growth_int;

	// Defining cytoplasmic species cooperator and cheater
	Cytoplasm cyto_coop, cyto_cheater;
	nutrient_int = cyto_coop.add_species("shared_good",0.0);	
	nutrient_int = cyto_cheater.add_species("shared_good",0.0);	
	cyto_coop.make_species_diffusible("shared_good",10.0); // nutrient will increase through absorption
	cyto_cheater.make_species_diffusible("shared_good",10.0); // nutrient will increase through absorption
	growth_int = cyto_coop.add_growth_rate_modifier("growth_rate_modifier",1.0); // dummy species to control growth
	growth_int = cyto_cheater.add_growth_rate_modifier("growth_rate_modifier",1.0); // dummy species to control growth
	
	// Defining cytoplasmic reactions
	ConstantReaction C1(10.0,nutrient_int,nutrient_int); // production of shared good
	LinearReaction LR1(-1.0,nutrient_int,nutrient_int); // consumption of shared good
	HillReaction HR1(10,3,2,nutrient_int,growth_int); // increase of growth rate with shared_good
	LinearReaction LR2(-10.0,growth_int,growth_int); //  natural decay of growth_rate (in the absence of shared_good) 
	
	// Adding reactions to cooperator cytoplasm,
	cyto_coop.add_reaction(&C1);
	cyto_coop.add_reaction(&LR1);
	cyto_coop.add_reaction(&HR1);
	cyto_coop.add_reaction(&LR2);

	// Adding reactions to cheater cytoplasm,
	cyto_cheater.add_reaction(&LR1);
	cyto_cheater.add_reaction(&HR1);
	cyto_cheater.add_reaction(&LR2);

	//  Initializing cells with the defined cytoplasm
	Population colony(10,"pars.in");
	colony.initialize_two_coopcheat(cyto_coop, cyto_cheater, 1,2);

	// Initializing dish with the created colony
	Dish dish(colony,200,20);
	dish.add_chemical("shared_good_out",0,1000.0,0.0001); // diffusible nutrient

	// dish.set_chemical_gaussianprofile("nutrient_out",1,1,1000,0.2);  Initial concentration of nutrient in the dish
	dish.link_chemical_bacterium("shared_good_out","shared_good"); // linking cytoplasmic and diffusible nutrient

	// Simulation time parameters
	double totaltime = 11; //8 
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
