#include "dish.h"
#define REFLECTIVE_BC 0
#define ABSORVING_BC 1
#define INITIAL_BC 2

Diffusable::Diffusable(){};

Diffusable::Diffusable(std::string name_, double diff_coeff_, double decay_coeff_, double init_conc){
	diff_coeff = diff_coeff_;
	name = name_;
	decay_coeff = decay_coeff_;
};


Dish::Dish(){
	N = 100;
	dim = 0;
	dt = 0.01;
}

Dish::Dish(Population& colony_, int N_){
	N = N_;
	colony = &colony_;
}

void Dish::add_chemical(std::string name, double conc, double coeff, double decay){

	diff_pars.emplace_back(name,coeff,decay); // adding diffusible constants
	diff_conc.push_back(gsl_matrix_alloc(N,N)); // creating spatial matrix 
	gsl_matrix_set_all(diff_conc.back(),conc); // initializing to specified concentration
}

void Dish::evolve(){
	colony->evolve(); // eventually some rules on timesteps should be set
	
	for(auto M: diff_conc){
		if(boundary_condition == REFLECTIVE_BC){
			for(int i;i<N;i++){
				diff_conc[i,0] = diff_conc[i,1]; 
				diff_conc[i,N-1] = diff_conc[i,N-1];
				diff_conc[0,i] = diff_conc[1,i];
				diff_conc[N-1,i] = diff_conc[N-2,i];
			}

		}

	for(int i; i<N-1; i++){
		for (int j; j<N-1; j++){
			diff_conc_aux[i,j] = gsl_matrix[i+1,j]+gsl_matrix[i-1,j]+gsl_matrix[i,j+1]+gsl_matrix[i,j-1] ; 
			diff_conc_aux[i,j] -= 4*gsl_matrix[i,j];
			diff_conc_aux[i,j] *= dt/(L/N)/(L/N);
		}
	}

	for(int i; i<N-1; i++){
		for (int j; j<N-1; j++){
			gsl_matrix[i,j] += diff_conc_aux[i,j]; 
		}
	}


	// diffusion loop

	// reaction loop 
}

void Dish::set_reflective_boundary(){
	boundary_condition = REFLECTIVE_BC;
}

void Dish::set_absorving_boundary(){
	boundary_condition = ABSORVING_BC;
}

void Dish::set_constant_boundary(){
	boundary_condition = INITIAL_BC;
}