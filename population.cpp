#include "population.h"

	Population::Population(){
		id = 0 ; 
		time = 0;
		rng = gsl_rng_alloc (gsl_rng_mt19937);
		gsl_rng_set (rng,::time(NULL)*getpid());
		std::cout<<"RNG Seed used: "<<::time(NULL)*getpid()<<'\n';
	}

	Population::Population(int seed, std::string inputfile){
		id = 0 ; 
		time = 0;
		rng = gsl_rng_alloc (gsl_rng_mt19937);
		if(seed==-1){ // if seed==-1 (default) take a random number
			gsl_rng_set (rng,::time(NULL)*getpid());
			std::cout<<"RNG Seed used: "<<::time(NULL)*getpid()<<'\n';
		}
		else{
			gsl_rng_set (rng,seed);
		}

		parsfile.open("pars.in");
		std::string name, value;
		while (parsfile >> name >> value){
			if (name=="dt"){ dt = std::stod(value);};
			if (name=="r"){ r = std::stod(value);};
			if (name=="length"){ length = std::stod(value);};
			if (name=="growth_rate"){ mean_growth_rate = std::stod(value);};
			if (name=="friction"){ friction_trans = std::stod(value);};
			if (name=="springk"){ springk = std::stod(value);};
			if (name=="mem"){ mem = std::stod(value);};

		}
		std::cout<<"READING! "<<dt<<' '<<r<<' '<<length<<' '<<mean_growth_rate<<' '<<friction_trans<<' '<<'\n';

		std::cout<<"Characteristics distances:\n";
		std::cout<<"Rest Elongation: "<<mean_growth_rate*dt<<'\n';
		std::cout<<"Stress Elongation: "<<mean_growth_rate*dt*dt*springk/friction_trans<<'\n';
		std::cout<<"Membrane Repulsion allowance: "<<mem*0.01*r*r*dt/friction_trans<<'\n';
		std::cout<<"membrane Repulsion elongation: "<<mem/friction_trans*dt*pow(friction_trans*mean_growth_rate/springk+mean_growth_rate*dt,2)<<'\n';
		// std::cout<<"Automatic calculation of frictions will override input....\n";
		// friction_trans = 20*mean_growth_rate*dt*dt;
	}



	int Population::next_id(){
		id += 1;
		return id;
	}

	double Population::timestep(){
		return dt;
	}

	double Population::currenttime(){
		return time;
	}

	void Population::initialize_two(Cytoplasm cyto){ // start Population with two parallel bacteria
		std::cout<<"Initizalizing Population...\n";
		vec2d posinit1,posinit2;
		double angleinit = 3.14159/4;
		posinit1[0] = 0;
		posinit1[1] = 0;
		posinit2[0] = length*cos(angleinit);
		posinit2[1] = length*sin(angleinit);
		std::cout<<"Init bac 1 "<< posinit1[0] <<' '<< posinit1[1]<<'\n';
		cells.emplace_back(next_id(), r, posinit1, posinit2, mean_growth_rate*(1+gsl_ran_gaussian(rng,0.2)), 2.0,mem, friction_trans, springk, cyto);
		// cells.front().cyto.add_species(1.0,"GFP");
		// HillRepReaction HR1(1,1,1,0,0);
		// cells.front().cyto.add_reaction(&HR1);
		// cells.front().cyto.react(0.1);
		// cells.front().cyto.print();

		posinit1[0] = 4*0.2*cos(angleinit);
		posinit1[1] = 0;
		posinit2[0] = length*cos(angleinit)+4*0.2*cos(angleinit);
		posinit2[1] = length*sin(angleinit);

		std::cout<<"Init bac 2 "<< posinit1[0] <<' '<< posinit1[1]<<'\n';
		cells.emplace_back(next_id(), r, posinit1, posinit2, mean_growth_rate*(1+gsl_ran_gaussian(rng,0.2)), 2.0,mem, friction_trans, springk, cyto);
		for(std::list<bacterium>::iterator cell_ptr = cells.begin(); cell_ptr != cells.end(); ++cell_ptr){
		// adding pointer to the cells to the alive list
			cells_alive.push_back(&*cell_ptr);	
		}
	}


	void Population::initialize_two_coopcheat(Cytoplasm cyto_coop, Cytoplasm cyto_cheat, double growth_coop, double growth_cheat){ // start Population with two parallel bacteria
		std::cout<<"Initizalizing Population...\n";
		vec2d posinit1,posinit2;
		double angleinit = 3.14159/4;
		posinit1[0] = 0;
		posinit1[1] = 0;
		posinit2[0] = length*cos(angleinit);
		posinit2[1] = length*sin(angleinit);
		std::cout<<"Init cooperator "<< posinit1[0] <<' '<< posinit1[1]<<'\n';
		cells.emplace_back(next_id(), r, posinit1, posinit2, growth_coop*(1+gsl_ran_gaussian(rng,0.2)), 2.0,mem, friction_trans, springk, cyto_coop);
		cells.back().set_type(0); // cooperator is type 0
		// cells.front().cyto.add_species(1.0,"GFP");
		// HillRepReaction HR1(1,1,1,0,0);
		// cells.front().cyto.add_reaction(&HR1);
		// cells.front().cyto.react(0.1);
		// cells.front().cyto.print();

		posinit1[0] = 4*0.2*cos(angleinit);
		posinit1[1] = 0;
		posinit2[0] = length*cos(angleinit)+4*0.2*cos(angleinit);
		posinit2[1] = length*sin(angleinit);

		std::cout<<"Init cheater "<< posinit1[0] <<' '<< posinit1[1]<<'\n';
		cells.emplace_back(next_id(), r, posinit1, posinit2, growth_cheat*(1+gsl_ran_gaussian(rng,0.2)), 2.0,mem, friction_trans, springk, cyto_cheat);
		cells.back().set_type(1); // cheater is type 1
		for(std::list<bacterium>::iterator cell_ptr = cells.begin(); cell_ptr != cells.end(); ++cell_ptr){
		// adding pointer to the cells to the alive list
			cells_alive.push_back(&*cell_ptr);	
		}
	}




	void Population::evolve(){ // evolve the popualtion a time step dt
		/// growth and division
		// std::cout<<"################# time: "<<time<<'\n';

		for(auto cell_ptr = cells_alive.begin(); cell_ptr!=cells_alive.end();) {
			(*cell_ptr)->reset_force();
			(*cell_ptr)->cyto.react(dt);
			(*cell_ptr)->diffusecyto(dt);
			(*cell_ptr)->grow(dt);
			if ((*cell_ptr)->division_ready()) 
    		{// division
    			// adding cell to dead list
    			cells_dead.push_back((*cell_ptr));

    			// Adding daughter 1
    			cells.push_back((*cell_ptr)->get_daughter1(next_id()));
    			cells.back().set_growth_rate(mean_growth_rate*(1+gsl_ran_gaussian(rng,0.2))); // random growth_rate
    			cells_alive.push_back(&cells.back());

    			// Adding daughter 2
				cells.push_back((*cell_ptr)->get_daughter2(next_id()));   
    			cells_alive.push_back(&cells.back());
    			cells.back().set_growth_rate(mean_growth_rate*(1+gsl_ran_gaussian(rng,0.2))); // random growth_rate

    			// removing pointer to old cell
    			cell_ptr = cells_alive.erase(cell_ptr); // returned cell_ptr points to next alive cell after elimination
    		}
    		else{ // if there is no division, check next cell
    			cell_ptr++;
    		}
    	}

    	//std::cout<<"calculating forces\n";
		/// calculate forces
		for(std::list<bacterium*>::iterator cell_ptr1 = cells_alive.begin(); cell_ptr1 != cells_alive.end(); ++cell_ptr1) {
			for(std::list<bacterium*>::iterator cell_ptr2 = cells_alive.begin(); cell_ptr1 != cell_ptr2; ++cell_ptr2) {
				update_force_between(**cell_ptr1,**cell_ptr2);
    		}
    	}
		/// move cells
		//std::cout<<"applying forces\n";
		for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
			// std::cout<<"Applying force on bacterium "<<(*cell_ptr)->id_bac()<<'\n';
			if((*cell_ptr)->apply_force(dt)==1){
				// std::cout<<"Large force in Population:\n";
				// for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
    // 				std::cout<<(*cell_ptr)->id_bac()<<' '<<(*cell_ptr)->centre()[0]<<' '<<(*cell_ptr)->centre()[1]<<' ';
    // 				std::cout<<(*cell_ptr)->angle()<<' '<<(*cell_ptr)->length()<<' '<<(*cell_ptr)->length0()<<' ';
    // 				std::cout<<(*cell_ptr)->current_force_1()[0]<<' '<<(*cell_ptr)->current_force_1()[1]<<' ';	
    // 				std::cout<<(*cell_ptr)->current_force_2()[0]<<' '<<(*cell_ptr)->current_force_2()[1]<<'\n';	
    			// }
			}
		}
		time += dt;
    }

    void Population::print_population(){
    	for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
    		std::cout<<(*cell_ptr)->get_id()<<' '<<(*cell_ptr)->get_centre()[0]<<' '<<(*cell_ptr)->get_centre()[1]<<'\n';
    	}
    }

    void Population::save(){
    	std::filesystem::create_directory("output/colony"); 
    	ofilename.str("");
    	ofilename<<"output/colony/colony_"<<std::setprecision(5)<<time<<".out";
    	trajfile.open(ofilename.str()); // output trajectory
    	for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
    		trajfile<<(*cell_ptr)->get_str_physics()<<' ';
    		trajfile<<(*cell_ptr)->cyto.get_str_concentrations()<<'\n';
    	}
    	trajfile.close();
    }

	bool Population::link_diffusible_bacterium(std::string chem_in, Diffusible* diffusible_){
		bool success = false;
		for (auto& bac: cells){
			success = (bac.link_diffusible(chem_in,diffusible_)|| success);
		}
		return success;
	}