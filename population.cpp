#include "population.h"

	population::population(int seed, std::string inputfile){
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



	int population::next_id(){
		id += 1;
		return id;
	}

	double population::timestep(){
		return dt;
	}

	double population::currenttime(){
		return time;
	}

	void population::initialize_two(){ // start population with two parallel bacteria
		std::cout<<"Initizalizing population...\n";
		vec2d posinit1,posinit2;
		double angleinit = 3.14159/4;
		posinit1[0] = 0;
		posinit1[1] = 0;
		posinit2[0] = length*cos(angleinit);
		posinit2[1] = length*sin(angleinit);
		std::cout<<"Init bac 1 "<< posinit1[0] <<' '<< posinit1[1]<<'\n';
		cells.emplace_back(next_id(), r, posinit1, posinit2, mean_growth_rate*(1+gsl_ran_gaussian(rng,0.2)), 2.0,mem, friction_trans, springk);

		posinit1[0] = 4*0.2*cos(angleinit);
		posinit1[1] = 0;
		posinit2[0] = length*cos(angleinit)+4*0.2*cos(angleinit);
		posinit2[1] = length*sin(angleinit);

		std::cout<<"Init bac 2 "<< posinit1[0] <<' '<< posinit1[1]<<'\n';
		cells.emplace_back(next_id(), r, posinit1, posinit2, mean_growth_rate*(1+gsl_ran_gaussian(rng,0.2)), 2.0,mem, friction_trans, springk);
		for(std::list<bacterium>::iterator cell_ptr = cells.begin(); cell_ptr != cells.end(); ++cell_ptr){
		// adding pointer to the cells to the alive list
			cells_alive.push_back(&*cell_ptr);	
		}
	}

	void population::evolve(){ // evolve the popualtion a time step dt
		/// growth and division
		// std::cout<<"################# time: "<<time<<'\n';
		for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
			(*cell_ptr)->reset_force();
			// growth
			(*cell_ptr)->grow(dt);
			if ((*cell_ptr)->division_ready()) 
    		{// division
    			vec2d pole1,pole2;
    			// std::cout<<"Dividing cell at "<<(*cell_ptr)->centre()[0]<<' '<<(*cell_ptr)->centre()[1]<<" \n";
    			// adding cell to dead list
    			cells_dead.push_back(*cell_ptr);
    			// finding new poles of daughter 1 and adding it to cells
    			(*cell_ptr)->get_daughter1_poles(pole1,pole2);
    			cells.emplace_back(next_id(),r,pole1,pole2, mean_growth_rate*(1+gsl_ran_gaussian(rng,0.2)),2.0,mem, friction_trans, springk);
    			cells_alive.push_back(&cells.back());
    			// std::cout<<"Created cell at"<<cells.back().centre()[0]<<' '<<cells.back().centre()[1]<<" \n";

    			// finding new poles of daughter 2 and adding it to cells
    			(*cell_ptr)->get_daughter2_poles(pole1,pole2);
    			cells.emplace_back(next_id(),r,pole1,pole2, mean_growth_rate*(1+gsl_ran_gaussian(rng,0.2)),2.0,mem, friction_trans, springk);
    			// std::cout<<"Created cell at"<<cells.back().centre()[0]<<' '<<cells.back().centre()[1]<<" \n";
    			cells_alive.push_back(&cells.back());
    			// removing pointer to old cell
    			cells_alive.erase(cell_ptr);
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
				// std::cout<<"Large force in population:\n";
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

    void population::print_population(){
    	for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
    		std::cout<<(*cell_ptr)->id_bac()<<' '<<(*cell_ptr)->centre()[0]<<' '<<(*cell_ptr)->centre()[1]<<'\n';
    	}
    }

    void population::save_population(){
    	ofilename.str("");
    	ofilename<<"output/population"<<std::setprecision(5)<<time<<".out";
    	trajfile.open(ofilename.str()); // output trajectory
    	for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
    		trajfile<<(*cell_ptr)->id_bac()<<' '<<trajfile<<(*cell_ptr)->type()<<' '<<(*cell_ptr)->centre()[0]<<' '<<(*cell_ptr)->centre()[1]<<' ';
    		trajfile<<(*cell_ptr)->angle()<<' '<<(*cell_ptr)->length()<<' '<<(*cell_ptr)->length0()<<' ';
    		trajfile<<(*cell_ptr)->current_force_1()[0]<<' '<<(*cell_ptr)->current_force_1()[1]<<' ';	
    		trajfile<<(*cell_ptr)->current_force_2()[0]<<' '<<(*cell_ptr)->current_force_2()[1]<<'\n';	
    	}
    	trajfile.close();
    }