#include "dish.h"
#define REFLECTIVE_BC 0
#define ABSORVING_BC 1
#define INITIAL_BC 2

Diffusable::Diffusable(){};

Diffusable::Diffusable(std::string name_, double diff_coeff_, double decay_coeff_, double init_conc, int N_){
	N_ = N;
	diff_coeff = diff_coeff_;
	name = name_;
	decay_coeff = decay_coeff_;
	conc = gsl_matrix_alloc(N,N);
	auxconc = gsl_matrix_alloc(N,N);

	// Creating views that will be used to calculate diffusion later on
	row0 = gsl_matrix_submatrix(conc,0,0,N,1);
	row1 = gsl_matrix_submatrix(conc,1,0,N,1);
	rowN = gsl_matrix_submatrix(conc,N-1,0,N,1);
	rowNm = gsl_matrix_submatrix(conc,N-2,0,N,1);
	col0 = gsl_matrix_submatrix(conc,0,0,N,1);
	col1 = gsl_matrix_submatrix(conc,0,1,1,N);
	colN = gsl_matrix_submatrix(conc,0,N-1,1,N);
	colNm = gsl_matrix_submatrix(conc,0,N-2,1,N);
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
	diffusables.emplace_back(name,coeff,decay,conc,N); // adding diffusible constants
}

void Dish::evolve(){
	colony->evolve(); // eventually some rules on timesteps should be set
	
	for(auto& M: diffusables){

		if(boundary_condition == REFLECTIVE_BC){ // if reflective boundaries set same vaue in two row/columns from the border
			gsl_matrix_memcpy(&M.row0.matrix,&M.row1.matrix);
			gsl_matrix_memcpy(&M.col0.matrix,&M.col1.matrix);
			gsl_matrix_memcpy(&M.rowN.matrix,&M.rowNm.matrix);
			gsl_matrix_memcpy(&M.colN.matrix,&M.colNm.matrix);
		}

	gsl_matrix_memcpy(M.auxconc,M.conc); // two matrices to avoid overwriting in incremental sum
	///// Implementation of finite difference Laplacian
	gsl_matrix_scale(M.auxconc,L*L/N/N*dt);
	gsl_matrix_add(&M.centre.matrix,&M.auxleft.matrix);
	gsl_matrix_add(&M.centre.matrix,&M.auxright.matrix);
	gsl_matrix_add(&M.centre.matrix,&M.auxtop.matrix);
	gsl_matrix_add(&M.centre.matrix,&M.auxbottom.matrix);
	gsl_matrix_scale(M.auxconc,4.0);
	gsl_matrix_sub(&M.centre.matrix,&M.auxcentre.matrix);
	// reaction loop 
	//...
	}
}

void Dish::set_reflective_boundary(){
	boundary_condition = REFLECTIVE_BC;
}

void Dish::set_absorving_boundary(){
	boundary_condition = ABSORVING_BC;
	for(auto& M: diffusables){
		gsl_matrix_set_all(&M.row0.matrix,0);
		gsl_matrix_set_all(&M.rowN.matrix,0);
		gsl_matrix_set_all(&M.col0.matrix,0);
		gsl_matrix_set_all(&M.colN.matrix,0);
	}
}

void Dish::set_constant_boundary(){
	boundary_condition = INITIAL_BC;
	for(auto& M: diffusables){
		gsl_matrix_memcpy(&M.row0.matrix,&M.row1.matrix);
		gsl_matrix_memcpy(&M.col0.matrix,&M.col1.matrix);
		gsl_matrix_memcpy(&M.rowN.matrix,&M.rowNm.matrix);
		gsl_matrix_memcpy(&M.colN.matrix,&M.colNm.matrix);
	}
}