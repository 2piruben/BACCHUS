#ifndef POPULATION_H
#define POPULATION_H


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <list>
#include <algorithm>
#include <math.h>
#include "bacterium.h"
#include "algebra2d.h"

typedef std::unique_ptr<bacterium> p_bacterium;

class population{

	private:

	std::list<bacterium> cells; // this contains all cells alive or dead
	std::list<bacterium*> cells_alive; // points to the cells in cells_all that are alive
	std::list<bacterium*> cells_dead; // poitns to the cells in cells_all that are dead

	gsl_rng * rng; // allocator for the rng generator
	
	double time = 0;

	double mean_growth_rate = 1;
	double deviation_growth_rate = 1;
	double mean_division_length = 1;
	double timestep = 0.001;
	double length = 1;
	int id = 0;
	double r = 0.2; // bacterium radius
	double maxd = mean_growth_rate*timestep; // maximum displacement allowed 
	double friction_trans;
	double friction_rot;
	std::ofstream trajfile; // File for output traj
	std::ifstream parsfile; // File for output traj
	std::stringstream ofilename;
	std::stringstream parsline;

	public:

	population(int seed, std::string inputfile);
	int next_id();
	void initialize_two();
	void evolve();
    void print_population();
    void save_population();

};

#endif
