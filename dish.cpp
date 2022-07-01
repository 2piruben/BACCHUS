#include "dish.h"
#include "diffusible.h"
#define REFLECTIVE_BC 0
#define ABSORVING_BC 1
#define INITIAL_BC 2


Dish::Dish(){
	N = 100;
	dim = 0;
	dt = 0.001;
	t = 0;
	L = 10;
}

Dish::Dish(Population& colony_, int N_, double L_){
	N = N_;
	colony = &colony_;
	t = 0;
	dt = 0.001;
	L = L_;
}

void Dish::add_chemical(std::string name, double conc, double coeff, double decay){
	diffusibles.emplace_back(name,coeff,decay,conc,N,L); // adding diffusible constants
}

void Dish::set_chemical_gaussianprofile(std::string name_, double x_c, double y_c, double max, double sigma){
	bool success = false;
	double x,y;
	for (auto& diffusible: diffusibles){
		if (diffusible.name == name_){
			success = true;
			for(int i=0; i<N; i++){
				for (int j=0; j<N; j++){
					x = grid_to_dish_x(j);
					y = grid_to_dish_y(i);
					gsl_matrix_set(diffusible.conc,i,j,max*exp(-((x-x_c)*(x-x_c)+(y-y_c)*(y-y_c))/(2*sigma*sigma)));
				}
			}
		}
	}
	if(success == false){
		std::cout<<"WARNING: Chemical species "<<name_<<" not found in Dish"<<'\n';
	}
}

void Dish::link_chemical_bacterium(std::string chem_out, std::string chem_in){

	bool success = false;

	for (auto& diffusible: diffusibles){
		if (diffusible.name == chem_out){
			success = colony->link_diffusible_bacterium(chem_in, &diffusible);
		}
	}

	if(success == false){
		std::cout<<"WARNING: Not found pair of chemicals to link "<<chem_out<<' '<<chem_in<<'\n';
	}
}



void Dish::evolve(){
	colony->evolve(); // eventually some rules on timesteps should be set
	t += dt;

	for(auto& M: diffusibles){

		if(boundary_condition == REFLECTIVE_BC){ // if reflective boundaries set same vaue in two row/columns from the border
			gsl_matrix_memcpy(&M.row0.matrix,&M.row1.matrix);
			gsl_matrix_memcpy(&M.col0.matrix,&M.col1.matrix);
			gsl_matrix_memcpy(&M.rowN.matrix,&M.rowNm.matrix);
			gsl_matrix_memcpy(&M.colN.matrix,&M.colNm.matrix);
		}

	gsl_matrix_memcpy(M.auxconc,M.conc); // two matrices to avoid overwriting in incremental sum
	///// Implementation of finite difference Laplacian
	gsl_matrix_scale(M.auxconc,M.diff_coeff*L*L/N/N*dt);
	gsl_matrix_add(&M.centre.matrix,&M.auxleft.matrix);
	gsl_matrix_add(&M.centre.matrix,&M.auxright.matrix);
	gsl_matrix_add(&M.centre.matrix,&M.auxtop.matrix);
	gsl_matrix_add(&M.centre.matrix,&M.auxbottom.matrix);
	gsl_matrix_scale(M.auxconc,4.0);
	gsl_matrix_sub(&M.centre.matrix,&M.auxcentre.matrix);
	//// Degradation of diffusible
	gsl_matrix_scale(&M.centre.matrix,(1-M.decay_coeff*dt));
	// difusible reaction loop 
	//...
	}
}

void Dish::set_reflective_boundary(){
	boundary_condition = REFLECTIVE_BC;
}

void Dish::set_absorving_boundary(){
	boundary_condition = ABSORVING_BC;
	for(auto& M: diffusibles){
		gsl_matrix_set_all(&M.row0.matrix,0);
		gsl_matrix_set_all(&M.rowN.matrix,0);
		gsl_matrix_set_all(&M.col0.matrix,0);
		gsl_matrix_set_all(&M.colN.matrix,0);
	}
}

void Dish::set_constant_boundary(){
	boundary_condition = INITIAL_BC;
	for(auto& M: diffusibles){
		gsl_matrix_memcpy(&M.row0.matrix,&M.row1.matrix);
		gsl_matrix_memcpy(&M.col0.matrix,&M.col1.matrix);
		gsl_matrix_memcpy(&M.rowN.matrix,&M.rowNm.matrix);
		gsl_matrix_memcpy(&M.colN.matrix,&M.colNm.matrix);
	}
}

double Dish::grid_to_dish_x(int j){
	return L/N*(0.5+j) - L/2;
}

double Dish::grid_to_dish_y(int i){
	return L/2 - L/N*(0.5+i);
}

int Dish::dish_to_grid_j(double x){
	return floor((x-L/2)/(L/N));//	
}
int Dish::dish_to_grid_i(double y){
	return floor((L/2-y)/(L/N));//	
}

void Dish::save(){
	for(int i = 0; i<diffusibles.size(); i++){
		std::filesystem::create_directory("output/diffusible");
		ofilename.str("");
		std::cout<<"t:"<<t<<'\n';
		ofilename<<"output/diffusible/diffusible_"<<i<<'_'<<std::setprecision(5)<<t<<".out";
    	trajfile.open(ofilename.str()); // output trajectory
		for(int j=0; j<N; j++){
			for(int k=0;k<N; k++){
				trajfile<<gsl_matrix_get(diffusibles[i].conc,j,k);
				if(k<N-1){
					trajfile<<' ';
				}
			}
			if(j<N-1){
				trajfile<<'\n';
			}
		}
		trajfile.close();
	}
}