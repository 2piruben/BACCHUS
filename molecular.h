# ifndef MOLECULAR_H
# define MOLECULAR_H

#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>

typedef std::vector<double> vec_species; // vector of molecular species

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
		Cytoplasm(const Cytoplasm &c);
		void add_reaction(Reaction* r); // add reaction acting on species idx
		void add_species(double conc, std::string name); // set the concentration of a certain species
		void add_growth_rate_modifier(double conc, std::string name); // set the concentration of a certain species
		void set_species(int idx, double conc); // set the concentration of a certain species
		double get_species(int idx); // return the amount of a certain species
		void react(double dt); // trigger all the reactions a window time dt
		void print(); // print the cytoplasmic information
		void print_complexity();
		void dilute(double factor); // dilute the content a certain factor, used in bacterial growth
		double get_growth_rate_modifier(); // 
		std::string get_str_concentrations(); //  get a string with the concentrataions of the different species

	private:

		vec_species s; // molecular species to track 
		vec_string s_names; // names of molecular species to track 
		vec_p_reaction reactions; // reactions that modifiy species 
		vec_vec_string r_names; // names of the reactions
		int growth_rate_idx; // variable in s[] used to track the modifier on growth_rate
		int dim_s; // number of species
		int dim_r; // number of reactions

};

#endif