#ifndef POPULATION_H
#define POPULATION_H


#include "algebra2d.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <unistd.h>
#include "bacterium.h"
#include <list>
#include <vector>


class Population{

	protected:

		gsl_rng * rng;
		std::stringstream ofilename;
		std::ifstream parsfile;
		std::ofstream trajfile;
		double length;
		double r;
		double mean_growth_rate;
		double friction_trans;
		std::list<bacterium> cells;
		std::list<bacterium*> cells_alive;
		std::list<bacterium*> cells_dead;
		std::list<bacterium*> cells_to_die;
		double dt;
		int id;
		double time;
		double springk;
		double mem;

	public:

		Population();
		Population(int seed, std::string inputfile);
		int next_id();
		double timestep();
		double currenttime();
		void initialize_two(Cytoplasm);
		void initialize_two_coopcheat(Cytoplasm cyto_coop, Cytoplasm cyto_cheat, double growth_coop, double growth_cheat);
		void evolve();
	    void print_population();
	    void save();
		bool link_diffusible_bacterium(std::string chem_in, Diffusible* diffusible_);


};

#endif    