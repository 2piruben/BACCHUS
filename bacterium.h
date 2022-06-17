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
#include <math.h>
#include <unistd.h>
#include <string>


class chemical{
	std::string name;
	double production;
	double degradation;
	double poison;
};



class bacterium{

protected:

	int id;
	vec2d pos[2]; // position of eadh pole of the bacterium
	double r; // bacterium width
	double age;
	double growth_rate;
	double l0; // preferred length of bacterium
	double springk; // l0 is maintained with Hooke's law with stiffness springk
	double division_length;
	vec2d force[2]; // force being applied at each pole of the bacterium
	double friction_trans; // friction of bacterium
	double mem; // membrane stiffness (parameter controlling interaction forces)
	int type; // dummy label for cell type
	vector<chemical>; // vector with the chemicals that 

public:

	bacterium(int id_, double r_, vec2d pos1_, vec2d pos2_, double growth_rate_,
		double division_length_, double mem, double friction_trans_, double springk_);

	void get_centre(vec2d &output);
	void get_orientation(vec2d &output);
	void get_length(double &length_);
	void get_angle(double &angle_);
	double length();
	double length0();
	double angle();
	double radius();
	vec2d centre();
	vec2d current_force_1();
	vec2d current_force_2();
	vec2d pole1();
	vec2d pole2();
	void set_type(int type_);
	int type_bac();


	int id_bac(){return id;};
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

	friend void update_force_between(bacterium &b1, bacterium &b2);

};

#endif
