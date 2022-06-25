# ifndef BACTERIUM_H
# define BACTERIUM_H


/////////////////////////////////////////////////////////////////////////
//
// SIMULATION OF GROWING BACTERIAL COLONIES
//
//							R Perez-Carrasco
//							Created: 12 MAY 2022
//////////////////////////////////////////////////////////////////////

#include "algebra2d.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include <unistd.h>
#include <string>
#include "molecular.h"


// class chemical{
// 	std::string name;
// 	double production;
// 	double degradation;
// 	double poison;
// };



class bacterium{

protected:

	int id;
	vec2d pos[2]; // position of each pole of the bacterium
	double r; // bacterium width
	double age;
	double growth_rate; // basic growth rate for a healthy single bacterium
	double l0; // preferred length of bacterium
	double springk; // l0 is maintained with Hooke's law with stiffness springk
	double division_length;
	vec2d force[2]; // force being applied at each pole of the bacterium
	double friction_trans; // friction of bacterium
	double mem; // membrane stiffness (parameter controlling interaction forces)
	int type; // dummy label for cell type
	
//  #vector<chemical>; // vector with the chemicals that 

public:

	bacterium(int id_, double r_, vec2d pos1_, vec2d pos2_, double growth_rate_,
		double division_length_, double mem, double friction_trans_, double springk_);
	bacterium(int id_, double r_, vec2d pos1_, vec2d pos2_, double growth_rate_,
		double division_length_, double mem, double friction_trans_, double springk_, Cytoplasm cyto_);
	void get_centre(vec2d &output);
	void get_orientation(vec2d &output);
	void get_length(double &length_);
	double get_length();
	void get_angle(double &angle_);
	double get_length0();
	double get_angle();
	double get_radius();
	vec2d get_centre();
	vec2d get_current_force_1();
	vec2d get_current_force_2();
	vec2d get_pole1();
	vec2d get_pole2();
	void set_type(int type_);
	void set_growth_rate(double);
	int get_type();
	Cytoplasm cyto; // containing species and reactions through the class cytoplasm

	int get_id(){return id;};
	double radius_bac(){return r;};
	// We want functions that transform points to 3 different coordinate systems
	// glob -> with respect to the dish (centered at 0,0)
	// shift -> with respect to the center of the bacterium
	// rot  -> with respect to the center of the bacterium, rotated so extreme 1 falls in in the direction (1,0)	
	void get_glob2rot(vec2d &output);
	void get_rot2glob(vec2d &output);
	void get_shift2glob(vec2d &output);
	void get_glob2shift(vec2d &output);
	void get_rot2shift(vec2d &output);
	void get_shift2rot(vec2d &output);
	void move(vec2d Dx);
	void reset_force();
	int apply_force(double dt); // return 1 rerseved for warning flags
	void grow(double dt);
	bool division_ready();
	void get_daughter1_poles(vec2d& p1, vec2d& p2);
	void get_daughter2_poles(vec2d& p1, vec2d& p2);
	bacterium get_daughter1(int id = -1);
	bacterium get_daughter2(int id = -1);
	std::string get_str_physics();

	friend void update_force_between(bacterium &b1, bacterium &b2);

};


#endif
