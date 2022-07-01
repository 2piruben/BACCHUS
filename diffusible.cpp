#include "diffusible.h"
#define REFLECTIVE_BC 0
#define ABSORVING_BC 1
#define INITIAL_BC 2

Diffusible::Diffusible(){};

Diffusible::Diffusible(std::string name_, double diff_coeff_, double decay_coeff_, double init_conc, int N_, double L_){
	N = N_;
	L = L_;
	diff_coeff = diff_coeff_;
	name = name_;
	decay_coeff = decay_coeff_;
	conc = gsl_matrix_alloc(N,N);
	auxconc = gsl_matrix_alloc(N,N);
	gsl_matrix_set_all(conc,init_conc);
	std::cout<<"Initializating diffusible\n";
	// Creating views that will be used to calculate diffusion later on
	row0 = gsl_matrix_submatrix(conc,0,0,1,N);
	row1 = gsl_matrix_submatrix(conc,1,0,1,N);
	std::cout<<"Checkpoint 1\n";
	rowN = gsl_matrix_submatrix(conc,N-1,0,1,N);
	rowNm = gsl_matrix_submatrix(conc,N-2,0,1,N);
	col0 = gsl_matrix_submatrix(conc,0,0,N,1);
	col1 = gsl_matrix_submatrix(conc,0,1,N,1);
	colN = gsl_matrix_submatrix(conc,0,N-1,N,1);
	colNm = gsl_matrix_submatrix(conc,0,N-2,N,1);
	left = gsl_matrix_submatrix(conc,1,0,N-2,N-2);
	right = gsl_matrix_submatrix(conc,1,2,N-2,N-2);
	top = gsl_matrix_submatrix(conc,0,1,N-2,N-2);
	bottom = gsl_matrix_submatrix(conc,2,1,N-2,N-2);
	centre = gsl_matrix_submatrix(conc,1,1,N-2,N-2);
	auxleft = gsl_matrix_submatrix(auxconc,1,0,N-2,N-2);
	auxright = gsl_matrix_submatrix(auxconc,1,2,N-2,N-2);
	auxtop = gsl_matrix_submatrix(auxconc,0,1,N-2,N-2);
	auxbottom = gsl_matrix_submatrix(auxconc,2,1,N-2,N-2);
	auxcentre = gsl_matrix_submatrix(auxconc,1,1,N-2,N-2);
	std::cout<<"Correct Initialization of diffusible "<<name<<'\n';
}

double Diffusible::grid_to_dish_x(int j){
	return L/N*(0.5+j) - L/2;
}

double Diffusible::grid_to_dish_y(int i){
	return L/2 - L/N*(0.5+i);
}

int Diffusible::dish_to_grid_j(double x){
	return floor((x+L/2)/(L/N));//	
}
int Diffusible::dish_to_grid_i(double y){
	return floor((L/2-y)/(L/N));//	
}

double Diffusible::voxel_length(){
	return L/N;
}

bool Diffusible::is_inside(int i, int j){
	if (i>0 && i<(N-1) && j>0 && j<(N-1)){
		return true;
	}
	else{
		return false;
	}
}