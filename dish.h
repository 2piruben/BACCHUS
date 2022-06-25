# ifndef DISH_H
# define DISH_H

#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>
#include "population.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


class Diffusable{

	public:

	Diffusable();
	Diffusable(std::string name_, double diff_coeff_, double decay_coeff_, double init_conc);

	double diff_coeff; // diffusion coefficient
	double decay_coeff; // decay coefficient
	double init_conc; // initial concentration
	std::string name;// name of the diffusable
};


class Dish{ // Class containing the whole dish, will contain a colony of cells plus diffusible gradients


	public:

		Dish();
		Dish(Population& colony_, int N_);
		void add_chemical(std::string name, double conc, double coeff, double decay);
		void link_chemical_bacterium(std::string chem_out, std::string chem_in, int bactype, double rate);
		void evolve();
		void set_reflective_boundary();
		void set_absorving_boundary();
		void set_constant_boundary();
		// void print();
		Population* colony;

	private:

		int N; // discretization of the dish
		int dim; // numer of diffusable species
		int L; // length of the side of the square dish
		std::vector<gsl_matrix*> diff_conc; // vector containing the different matrices of diffusible substances
		gsl_matrix* diff_conc_aux; // auxiliar vector containing the different matrices of diffusible substances
		std::vector<Diffusable> diff_pars; // vector containing the names of the different diffusers
		int boundary_condition; // 
		double dt;
		
};



#endif