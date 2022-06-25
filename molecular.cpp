#include "molecular.h"

Cytoplasm::Cytoplasm(){
	dim_s = 0;
	dim_r = 0;
	growth_rate_idx = -1; // negative means that there is not growth rate modifier
}

Cytoplasm::Cytoplasm(const Cytoplasm &old_cytoplasm){ // custom copy of cytoplasm
	std::cout<<"Custom constuctor Cytoplasm\n";
	dim_s = old_cytoplasm.dim_s;
	dim_r = old_cytoplasm.dim_r;
	growth_rate_idx = old_cytoplasm.growth_rate_idx; // negative means that there is not growth rate modifier

	for(auto & ss: old_cytoplasm.s){
		s.push_back(ss);
	}
	for(auto & name: old_cytoplasm.s_names){
		s_names.push_back(name);
	}
	for(auto & reaction: old_cytoplasm.reactions){
		reactions.push_back(reaction);
	}
	for(auto & rname: old_cytoplasm.r_names){
		r_names.push_back(rname);
	}
}

void Cytoplasm::add_reaction(Reaction* r){
	reactions.push_back(r);	
	dim_r++;
}

void Cytoplasm::set_species(int idx, double conc){
	s[idx] = conc;
}

void Cytoplasm::add_species(double conc, std::string name){
	s.push_back(conc);
	s_names.push_back(name);
	dim_s++;
}

void Cytoplasm::add_growth_rate_modifier(double conc, std::string name){
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
	std::cout<<"Entering vector with "<<reactions.size()<<" elements\n";
	for (auto & r: reactions){
//		r->react(s,dt);// for each reaction r react on the species vector s
	}
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
	std::cout<<"k"<< k<<'\n';
	s[idx_out] += dt*k*s[idx_in]; 	
}



