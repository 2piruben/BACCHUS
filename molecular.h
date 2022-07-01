# ifndef MOLECULAR_H
# define MOLECULAR_H

#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include "diffusible.h"
#include "algebra2d.h"


class Species{// class to encode molecular species

	public: 

		Species();
		Species(std::string name, double conc);
		double get_conc();
		void set_conc(double conc);
		void increase_conc(double conc);
		void make_diffusible(double k);
		std::string get_name();
		bool link_diffusible(Diffusible* diffusible);
		bool is_diffusible();
		Diffusible* diffusible; // matrix of the morphogen outside the cell, 
		void diffuse(vec2d* pos, double area, double dt); // diffuse along the bacterial backbone pos for a duration dt. Perhaps this can be replace with a friend function bacterium, molecular...

	private:
		std::string name; 
		double conc; // concentration
		bool diff; // true if species is diffusable
		double k_diff; // rate of diffusion
		// pointer to diff_matrix will be created when the morphogen matrix is created
		// with the Dish class
};
typedef std::vector<Species> vec_species; // vector of molecular species

class Reaction{ // Abstract class to control different possible reactions, it will be the parent of particular reaction types with different input parameters

	public:
		Reaction();
		Reaction(int idx);
		virtual ~Reaction() = default; // abstract class destructor 
		virtual void react(vec_species& v, double dt)=0; // pure virtual, will be override by particular reactions


	private:
		int idx_out; // A reaction will operate on a vec_species, modifyinf element idx_out

};

class HillReaction : public Reaction{

	public:
		HillReaction(double A_, double K_, double n_, int idxin_, int idxout_);
		~HillReaction();
		void react(vec_species& v, double dt) override;

	private:
		double A;
		double K;
		double n;
		int idx_out;
		int idx_in;
};

class LinearReaction : public Reaction{

	public:
		LinearReaction(double k_, int idxin_, int idxout_);
		~LinearReaction();
		void react(vec_species& v, double dt) override;

	private:
		double k;
		int idx_out;
		int idx_in;
};

class ConstantReaction : public Reaction{

	public:
		ConstantReaction(double k_, int idxin_, int idxout_);
		~ConstantReaction();
		void react(vec_species& v, double dt) override;

	private:
		double k;
		int idx_out;
		int idx_in;
};

// typedef std::shared_ptr<Reaction> p_reaction; // shared pointers are required so 
typedef std::vector<Reaction*> vec_p_reaction; // vec_p_reaction contains all the reactions that occur on a cell
// they need to be pointers because Reaction is an abstract class
typedef std::vector<std::string> vec_string;
typedef std::vector<vec_string> vec_vec_string;


//////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////	
class Cytoplasm{//Class controlling the molecular content of a bacterium.
////////////
//////////

	public:

		Cytoplasm();
		// Cytoplasm(const Cytoplasm &c);
		void add_reaction(Reaction* r); // add reaction acting on species idx
		int add_growth_rate_modifier(std::string name,double conc); // add a dummy species that will track the change in growth_rate
		int add_species(std::string name,double conc); // add a chemical species to the cytoplasm
		void make_species_diffusible(std::string name, double k); // 
		void set_species(int idx, double conc); // set the concentration of a certain species
		double get_species(int idx); // return the amount of a certain species
		void react(double dt); // trigger all the reactions a window time dt
		void print(); // print the cytoplasmic information
		void print_complexity();
		void dilute(double factor); // dilute the content a certain factor, used in bacterial growth
		double get_growth_rate_modifier();
		bool link_diffusible(std::string name, Diffusible* diffusible_); // link species name to a dish diffusible matrix 
		bool is_diffusible();
		void diffuse(vec2d* pos, double area, double dt); // diffuse along the bacterium against the media the chemicals that are able to diffuse
		std::string get_str_concentrations(); //  get a string with the concentrataions of the different species

	private:

		vec_species s; // molecular species to track 
		// vec_species s_diff; // molecular species to track that are diffusable
		vec_p_reaction reactions; // reactions that modifiy species 
		int growth_rate_idx; // variable in s[] used to track the modifier on growth_rate
		int dim_s; // number of species
		int dim_r; // number of reactions
		bool diff; // Does the cytoplasm contain diffusible species?

};

#endif