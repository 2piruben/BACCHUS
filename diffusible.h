# ifndef DIFFUSIBLE_H
# define DIFFUSIBLE_H

#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <gsl/gsl_matrix.h>


class Diffusible{

	private:

	int N;
	double L;

	public:

	Diffusible();
	Diffusible(std::string name_, double diff_coeff_, double decay_coeff_, double init_conc, int N, double L);

	gsl_matrix* conc; // matrix with the concentrations
	gsl_matrix* auxconc; // auxiliary matrix to perform operations
	gsl_matrix_view row0,row1,rowN,rowNm; // views to the conc matrix
	gsl_matrix_view col0,col1,colN,colNm; // views to the conc matrix
	gsl_matrix_view centre,left,right,top,bottom; // views to (N-2)x(N-2) submatrices used in the diffusion
	gsl_matrix_view auxcentre,auxleft,auxright,auxtop,auxbottom; // views to (N-2)x(N-2) submatrices used in the diffusion
	double diff_coeff; // diffusion coefficient
	double decay_coeff; // decay coefficient
	double init_conc; // initial concentration
	std::string name;// name of the diffusible
	double grid_to_dish_x(int j);
	double grid_to_dish_y(int i);
	int dish_to_grid_j(double x);
	int dish_to_grid_i(double y);
	double voxel_length();
	bool is_inside(int i, int j); // check if element i,j is inside the matrix of concentrations

};

#endif