# ifndef DISH_H
# define DISH_H

#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "algebra2d.h"
#include "diffusible.h"
#include "population.h"


class Dish{ // Class containing the whole dish, will contain a colony of cells plus diffusible gradients


	public:

		Dish();
		Dish(Population& colony_, int N_, double L_);
		void add_chemical(std::string name, double conc, double coeff, double decay);
		void link_chemical_bacterium(std::string chem_out, std::string chem_in);
		void set_chemical_gaussianprofile(std::string name_, double x_c, double y_c, double max, double sigma);
		void evolve();
		void set_reflective_boundary();
		void set_absorving_boundary();
		void set_constant_boundary();
		void save();
		double grid_to_dish_x(int j);// coordinate x of the centre of col j
		double grid_to_dish_y(int i);// coordinate y of the centre of row i
		int dish_to_grid_j(double x);// column of the grid corresponding to coordinate x
		int dish_to_grid_i(double y);// row  of the grid corresponding to coordinate y
		std::vector<Diffusible> diffusibles; // vector containing the different matrices of diffusible substances

		// void print();
		Population* colony;

	private:

		int N; // discretization of the dish
		int dim; // numer of diffusible species
		double L; // length of the side of the square dish
		int boundary_condition; // 
		double t;
		double dt;
		std::stringstream ofilename;
		std::ifstream parsfile;
		std::ofstream trajfile;
		
};



#endif