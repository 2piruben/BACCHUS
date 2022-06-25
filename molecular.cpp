#include "molecular.h"

Cytoplasm::Cytoplasm(){
	dim_s = 0;
	dim_r = 0;
	growth_rate_idx = -1; // negative means that there is not growth rate modifier
}

void Cytoplasm::add_reaction(Reaction* r){
	reactions.push_back(r);	
	dim_r++;
}


void Cytoplasm::set_species(int idx, double conc){
	s[idx] = conc;
}

void Cytoplasm::add_species(std::string name,double conc){
	s.push_back(conc);
	s_names.push_back(name);
	dim_s++;
}

void Cytoplasm::add_growth_rate_modifier(std::string name, double conc){
	s.push_back(conc);
	s_names.push_back(name);
	growth_rate_idx = dim_s;
	dim_s++;
}

void Cytoplasm::print(){
	for(int i=0;i<dim_s;i++){
		std::cout<<s_names[i]<<' '<<s[i]<<'\n';
	}
}

void Cytoplasm::dilute(double factor){
	for(auto & ss: s){
		ss *= factor; 
	}
}

double Cytoplasm::get_growth_rate_modifier(){
	if(growth_rate_idx == -1){
		return 1;
	}
	else{
		return s[growth_rate_idx];
	}
}

std::string Cytoplasm::get_str_concentrations(){
	std::stringstream ssconc;
	for(int i=0;i<dim_s;i++){
		ssconc<<s[i];
		if(i<(dim_s-1)){
			ssconc<<' ';
		}
	}
	return ssconc.str();
}

void Cytoplasm::react(double dt){
	std::cout<<"Entering print\n";
	for (auto & r: reactions){
		r->react(s,dt);// for each reaction r react on the species vector s
	}
	std::cout<<"Exiting print\n";
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

double Species::make_diffusable(double k){
	diff = true;
	k_diff = k;
}

// Reaction functions children of the abstract class Reaction

Reaction::Reaction(){};

// Hill Repressive function
HillRepReaction::HillRepReaction(double A_, double K_, double n_, int idx_in_, int idx_out_){
	A = A_;
	K = K_;
	idx_out = idx_out_;
	idx_in = idx_in_;
}
HillRepReaction::~HillRepReaction(){};
void HillRepReaction::react(vec_species& s, double dt){
	s[idx_out] += dt*A/(1+pow(s[idx_in]/K,n)); 	
}

// Linear reaction
LinearReaction::LinearReaction(double k_, int idx_in_, int idx_out_){
	k = k_;
	idx_out = idx_out_;
	idx_in = idx_in_;
}
LinearReaction::~LinearReaction(){};
void LinearReaction::react(vec_species& s, double dt){
	s[idx_out] += dt*k*s[idx_in]; 	
}



