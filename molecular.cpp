#include "molecular.h"

Cytoplasm::Cytoplasm(){
	dim_s = 0;
	dim_r = 0;
	growth_rate_idx = -1; // negative means that there is not growth rate modifier
}

// Cytoplasm::Cytoplasm(const Cytoplasm &old_cytoplasm){ // custom copy of cytoplasm
// 	// std::cout<<"Custom constuctor Cytoplasm\n";
// 	dim_s = old_cytoplasm.dim_s;
// 	dim_r = old_cytoplasm.dim_r;
// 	growth_rate_idx = old_cytoplasm.growth_rate_idx; // negative means that there is not growth rate modifier
// 	diff = old_cytoplasm.diff;

// 	for(auto & ss: old_cytoplasm.s){
// 		s.push_back(ss);
// 	}
// 	for(auto & name: old_cytoplasm.s_names){
// 		s_names.push_back(name);
// 	}
// 	for(auto reaction: old_cytoplasm.reactions){
// 		reactions.push_back(reaction);
// 	}
// 	for(auto & rname: old_cytoplasm.r_names){
// 		r_names.push_back(rname);
// 	}
// 	//std::cout<<"rsize"<<reactions.size()<<'\n';
// }

void Cytoplasm::add_reaction(Reaction* r){
	reactions.push_back(r);	
	dim_r++;
}

void Cytoplasm::set_species(int idx, double conc){
	s[idx].set_conc(conc);
}

int Cytoplasm::add_species(std::string name,double conc){
	s.emplace_back(name,conc);
	dim_s++;
	return dim_s-1;
	
}

void Cytoplasm::make_species_diffusible(std::string name,double k){
	bool success = false;
	for(auto& species: s){
		if(species.get_name() == name){
			species.make_diffusible(k);
			success = true;
		}
	}
	if(success == false){
		std::cout<<"WARNING: Species does not exist\n";
	}
}

int Cytoplasm::add_growth_rate_modifier(std::string name, double conc){
	s.emplace_back(name,conc);
	growth_rate_idx = dim_s; // this will be the index to use for the growth_rate_modifier
	dim_s++;
	return growth_rate_idx;
}

void Cytoplasm::print(){
	for(auto& species: s){
		std::cout<<species.get_name()<<' '<<species.get_conc()<<'\n';
	}
}

void Cytoplasm::print_complexity(){
	std::cout<<"Cytoplasm with "<<s.size()<<" species and "<<reactions.size()<<" reaction\n";
}

void Cytoplasm::dilute(double factor){
	for(auto & species: s){
		species.set_conc(species.get_conc()*factor); 
	}
}

double Cytoplasm::get_growth_rate_modifier(){
	if(growth_rate_idx == -1){
		return 1;
	}
	else{
		return s[growth_rate_idx].get_conc();
	}
}

std::string Cytoplasm::get_str_concentrations(){
	std::stringstream ssconc;
	for(int i=0;i<dim_s;i++){
		ssconc<<s[i].get_conc();
		if(i<(dim_s-1)){
			ssconc<<' ';
		}
	}
	return ssconc.str();
}



void Cytoplasm::react(double dt){
	// std::cout<<"Entering print\n";
	std::cout<<"Reacting in cell with content "<<get_str_concentrations()<<'\n';
	for (auto & r: reactions){
		r->react(s,dt);// for each reaction r react on the species vector s
		std::cout<<"After reacting: "<< get_str_concentrations()<<'\n';
	}
}

bool Cytoplasm::link_diffusible(std::string name_, Diffusible* diffusible_){
	bool success = false;
	for(auto& species: s){
		if(species.get_name()==name_){
			success = (species.link_diffusible(diffusible_) || success);
			diff = true;
		}
	}
	return success;
}

void Cytoplasm::diffuse(vec2d* pos, double area, double dt){
	for(auto& species: s){
		if(species.is_diffusible()){
			species.diffuse(pos,area,dt);
		}
	}
}

bool Cytoplasm::is_diffusible(){
	return diff;
}

Species::Species(){
	diff = false;
}

Species::Species(std::string name_, double conc_){
	diff = false;
	name = name_;
	conc = conc_;
}

double Species::get_conc(){
	return conc;
}

void Species::set_conc(double conc_){
	conc = conc_;
}

void Species::increase_conc(double conc_){
	conc += conc_;
}


void Species::make_diffusible(double k){
	diff = true;
	k_diff = k;
}

bool Species::is_diffusible(){
	return diff;
}


std::string Species::get_name(){
	return name;
}

bool Species::link_diffusible(Diffusible* diffusible_){
	if (diff == true){
		diffusible = diffusible_; // bacterium will interact with matrix p (be careful with ownership of p!!)
		return true;
	}
	else{
		std::cout<<"WARNING: linking to a reactant that is not set to diffuse\n";
		return false;
	}
}

void Species::diffuse(vec2d* pos, double area, double dt){
	std::vector<vec2d_int> backbone; // set of voxels to modify
	int i,j;
	int finali,finalj; // final voxel
	double tanangle;
	double currentx = pos[0][0];
	double currenty = pos[0][1]; // starting x/y
	int incrementx = 1;
	int incrementy = 1; // these will replace by -1,will depend on the orientation of the bacteria, and be constant for the full function
	double dx, dy; // proposed advanced distance 
	double centrex, centrey; // centre of the voxel studied
	double conc_out; //cout,  (c_in-c_out_segment)*length_segment 

	/// Creating backbone for the cell
	j = diffusible->dish_to_grid_j(pos[0][0]); // starting voxel
	i = diffusible->dish_to_grid_i(pos[0][1]);
	finalj = diffusible->dish_to_grid_j(pos[1][0]); // starting voxel
	finali = diffusible->dish_to_grid_i(pos[1][1]);

	// Checking the bacterium is inside the limits of the grid
	if ((diffusible->is_inside(i,j))==false or (diffusible->is_inside(finali,finalj))==false){
		std::cout<<"WARNING: Bacterium out of bounds of diffusable agent, not performing transmembrane diffusion\n";
		std::cout<<"Bacterium at position "<<pos[0]<<' '<<pos[1]<<'\n';
		std::cout<<"Corresponding to matrix ("<<i<<','<<j<<") ("<<finali<<','<<finalj<<'\n';
	}

	// setting bacterim discretized orientation
	if(pos[1][0]<pos[0][0]){incrementx = -1;};
	if(pos[1][1]<pos[0][1]){incrementy = -1;};
	centrex = diffusible->grid_to_dish_x(j); // centre of current voxel
	centrey = diffusible->grid_to_dish_y(i); 
	tanangle = tan((pos[1]-pos[0]).angle());
	backbone.emplace_back(i,j);
	while (i!=finali && j!=finalj){ // creating the backbone ...
		// std::cout<<"Diffusing on voxel "<<i<<' '<<j<<" located at a position "<<centrex<<' '<<centrey<<'\n';
		dx = (centrex - currentx) + incrementx * diffusible->voxel_length()/2;
		dy = (centrey - currenty) + incrementy * diffusible->voxel_length()/2;
		if(tanangle>dy/dx){
			i-= incrementy; // vertical advance (in matrix coordinates)
			currentx += dy/tanangle; // moving to cross point (in dish coordinates)
			currenty += dy;
			centrey += incrementy*diffusible->voxel_length(); 
		}
		else{
			j+= incrementx; // horizontal advance (in matrix coordinates)
			currentx += dx; // moving to cross point (in dish coordinates)
			currenty += dx*tanangle;
			centrex += incrementx*diffusible->voxel_length();
		}
		backbone.emplace_back(i,j);
	} // backbone finished

		// diffusing on the backbone
		for(auto& cell: backbone){
			conc_out = gsl_matrix_get(diffusible->conc,cell[0],cell[1]);
			gsl_matrix_set(diffusible->conc,cell[0],cell[1],
			    conc_out - k_diff*dt*(conc_out - conc)*area/(backbone.size()*diffusible->voxel_length()*diffusible->voxel_length())); // area correction to conserve total number of molecules
				conc += k_diff*dt*(conc_out - conc); // assumes the backbone is inscribed in the bacterium (i.e. no requirement for area correction)
			}

	conc /= backbone.size(); // average concentration change (preserving total number of molecules exchanged see notes) 
}

// Reaction functions children of the abstract class Reaction

Reaction::Reaction(){};

// Hill Repressive function
HillReaction::HillReaction(double A_, double K_, double n_, int idx_in_, int idx_out_){
	A = A_;
	K = K_;
	n = n_;
	idx_out = idx_out_;
	idx_in = idx_in_;
}
HillReaction::~HillReaction(){};
void HillReaction::react(vec_species& s, double dt){
	//std::cout<<"Calling Hill Reaction with pars"
	s[idx_out].increase_conc(dt*A/(1+pow(s[idx_in].get_conc()*K,(-1.0*n)))); 	
}

// Linear reaction
LinearReaction::LinearReaction(double k_, int idx_in_, int idx_out_){
	k = k_;
	idx_out = idx_out_;
	idx_in = idx_in_;
}
LinearReaction::~LinearReaction(){};
void LinearReaction::react(vec_species& s, double dt){
	// std::cout<<"k"<< k<<'\n';
	s[idx_out].increase_conc(dt*k*s[idx_in].get_conc()); 	
}

// Constant reaction
ConstantReaction::ConstantReaction(double k_, int idx_in_, int idx_out_){
	k = k_;
	idx_out = idx_out_;
	idx_in = idx_in_;
}
ConstantReaction::~ConstantReaction(){};
void ConstantReaction::react(vec_species& s, double dt){
	// std::cout<<"k"<< k<<'\n';
	s[idx_out].increase_conc(dt*k); 	
}



