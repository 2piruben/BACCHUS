# ifndef MOLECULAR_H
# define MOLECULAR_H

#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>


typedef std::vector<double> vec_species; // vector of molecular species
typedef std::vector<double> vec_species; // vector of molecular species

class Species{// class to encode molecular species

		Species();
		Species(std::string name, double conc);
		double get_conc();
		void set_conc(double conc);
		double make_diffusable(double k);

	private:
		std::string name; 
		double conc; // concentration
		bool diff; // true if species is diffusable
		double k_diff; // rate of diffusion
		gsl_matrix* diff_matrix; // matrix of the morphogen outside the cell, 
		// pointer to diff_matrix will be created when the morphogen matrix is created
		// with the Dish class
};


class Reaction{ // Abstract class to control different possible reactions, it will be the parent of particular reaction types with different input parameters

	public:
		Reaction();
		Reaction(int idx);
		virtual ~Reaction() = default; // abstract class destructor 
		virtual void react(vec_species& v, double dt)=0; // pure virtual, will be override by particular reactions


	private:
		int idx_out; // A reaction will operate on a vec_species, modifyinf element idx_out

};

class HillRepReaction : public Reaction{

	public:
		HillRepReaction(double A_, double K_, double n_, int idxin_, int idxout_);
		~HillRepReaction();
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


typedef std::vector<Reaction*> vec_p_reaction; // vec_p_reaction contains all the reactions that occur on a cell
typedef std::vector<std::string> vec_string;
typedef std::vector<vec_string> vec_vec_string;


//////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////	
class Cytoplasm{//Class controlling the molecular content of a bacterium.
////////////
//////////

	public:

		Cytoplasm();
		void add_reaction(Reaction* r); // add reaction acting on species idx
		void add_diffusible_reaction(int diff_in,int diff_out, double rate); // connect a species of the cytoplasm with a diffusible one
		void add_species(std::string name,double conc); // set the concentration of a certain species
		void add_growth_rate_modifier(std::string name,double conc); // set the concentration of a certain species
		void set_species(int idx, double conc); // set the concentration of a certain species
		double get_species(int idx); // return the amount of a certain species
		void react(double dt); // trigger all the reactions a window time dt
		void print(); // print the cytoplasmic information
		void dilute(double factor); // dilute the content a certain factor, used in bacterial growth
		double get_growth_rate_modifier(); // 
		std::string get_str_concentrations(); //  get a string with the concentrataions of the different species

	private:

		vec_species s; // molecular species to track 
		vec_species s_diff; // molecular species to track that are diffusable
		vec_string s_names; // names of molecular species to track 
		vec_p_reaction reactions; // reactions that modifiy species 
		vec_vec_string r_names; // names of the reactions
		int growth_rate_idx; // variable in s[] used to track the modifier on growth_rate
		int dim_s; // number of species
		int dim_r; // number of reactions

};

#endif