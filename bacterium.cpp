#include "bacterium.h"
#include "algebra2d.h"

	bacterium::bacterium(int id_, double r_, vec2d pos1_, vec2d pos2_, double growth_rate_,
		double division_length_, double mem_, double friction_trans_, double springk_){
		std::cout<<"created bacteria "<<id_<<" at "<<pos1_<<' '<<pos2_<<'\n';
		id  = id_;
		pos[0][0]  = pos1_[0];
		pos[0][1]  = pos1_[1];
		pos[1][0]  = pos2_[0];
		pos[1][1]  = pos2_[1];
		r = r_;
		mem = mem_; 
		age = 0;
		l0 = dist(pos1_,pos2_); // initial target length of bacterium
		division_length = division_length_;
		growth_rate = growth_rate_;
		friction_trans = friction_trans_;
		springk = springk_; 
		reset_force();
		type = 0;
		cyto = Cytoplasm();
	}

	bacterium::bacterium(int id_, double r_, vec2d pos1_, vec2d pos2_, double growth_rate_,
		double division_length_, double mem_, double friction_trans_, double springk_, Cytoplasm cyto_){
		std::cout<<"created bacteria "<<id_<<" at "<<pos1_<<' '<<pos2_<<'\n';
		id  = id_;
		pos[0][0]  = pos1_[0];
		pos[0][1]  = pos1_[1];
		pos[1][0]  = pos2_[0];
		pos[1][1]  = pos2_[1];
		r = r_;
		mem = mem_; 
		age = 0;
		l0 = dist(pos1_,pos2_); // initial target length of bacterium
		division_length = division_length_;
		growth_rate = growth_rate_;
		friction_trans = friction_trans_;
		springk = springk_; 
		reset_force();
		type = 0;
		cyto = cyto_;
	}


/// Function returning or setting private atributes 

	void bacterium::get_centre(vec2d &output){
		output = (pos[0]+pos[1])/2;
	}


	void bacterium::get_orientation(vec2d &output){
		// vector that goes from centre of bacterium to extreme1
		output = (pos[0]-pos[1])/2;
	}

	void bacterium::get_length(double &length_){
		length_ = dist(pos[0],pos[1]);
	}

	void bacterium::get_angle(double &angle_){
		vec2d temp_vec;
		get_orientation(temp_vec);
		angle_ = temp_vec.angle();
	}

	double bacterium::get_angle(){
		vec2d temp_vec;
		get_orientation(temp_vec);
		return temp_vec.angle();
	}

	double bacterium::get_length(){
		return dist(pos[0],pos[1]);
	}

	double bacterium::get_length0(){
		return l0;
	}


	double bacterium::get_radius(){
		return r;
	}

	int bacterium::get_type(){
		return type;
	}

	vec2d bacterium::get_centre(){
		return (pos[0]+pos[1])/2.0;
	}

	vec2d bacterium::get_current_force_1(){
		return force[0];
	}

	vec2d bacterium::get_current_force_2(){
		return force[1];
	}

	vec2d bacterium::get_pole1(){
		return pos[0];
	}

	vec2d bacterium::get_pole2(){
		return pos[1];
	}

	void bacterium::set_type(int type_){
		type = type_;
	}

	void bacterium::set_growth_rate(double growth_rate_){
		growth_rate = growth_rate_;
	}

/// Functions managing reference systems of vectors with respect to bacteria

	// We want functions that transform points to 3 different coordinate systems
	// glob -> with respect to the dish (centered at 0,0)
	// shift -> with respect to the center of the bacterium
	// rot  -> with respect to the center of the bacterium, rotated so extreme 1 falls in in the direction (1,0)

	void bacterium::get_glob2rot(vec2d &output){
		// Given a vector, return the vector in the reference system of the bacterium orientated (E1 in the cartesian direction, centered at the bacterium)
		vec2d temp_vec;
		get_centre(temp_vec);
		temp_vec = output - temp_vec;
		output = rotate(temp_vec,-1.0*get_angle());
	}

	void bacterium::get_glob2shift(vec2d &output){
		// Given a global vector, return the vector in the reference system of the bacterium preserving the orientation
		output -= (pos[0]+pos[1])/2;
	}

	void bacterium::get_rot2glob(vec2d &output){
		// Opposite to vec2bac, given a vector in the reference system of the bacterium orientated, return the vector in the reference system of the dish
		output = rotate(output,get_angle());
		output += (pos[0]+pos[1])/2;
	}

	void bacterium::get_shift2glob(vec2d &output){
		// Opposite to vec2bacpos, given a vector ceneterd in the bacterium, return the vector in the global coordinate system
		output += (pos[0]+pos[1])/2;
	}

	void bacterium::get_rot2shift(vec2d &output){
		// Opposite to vec2bac, given a vector in the reference system of the bacterium orientated, return the vector in the reference system of the dish
		output = rotate(output,get_angle());
	}

	void bacterium::get_shift2rot(vec2d &output){
		output = rotate(output,-1.0*get_angle());
	}

// Mechanical functions for the bacteria

	void bacterium::move(vec2d Dx){
		pos[0] += Dx;
		pos[1] += Dx;
	}

	void bacterium::reset_force(){ // remove external forces and recalculate the pole-pole spring force
		vec2d temp_vec; 
		double l;
		get_orientation(temp_vec);
		get_length(l);
		//std::cout<<"l0 ="<<l0<<"   l1 ="<<l<<'\n';
		force[0] = -temp_vec/temp_vec.modulus()*springk*(l-l0)/2; // 2 for splitting force in two poles
		force[1] = temp_vec/temp_vec.modulus()*springk*(l-l0)/2;
		//std::cout<<"Initial force on bac: "<<id<<' '<<force[0]<<' '<<force[1]<<'\n';
	}

	int bacterium::apply_force(double dt){ // Apply the stored force and reset the force
		vec2d temp_vec;
		get_orientation(temp_vec);
		vec2d cappedforce0,cappedforce1;
		double length_before;
		length_before = get_length();
		// std::cout<<"Force on pole1 "<<force[0]<<'\n';
		// std::cout<<"Force on pole2 "<<force[1]<<'\n';
		// std::cout<<"orientation before force"<<temp_vec<<'\n';
		if (force[0].modulus()>0){
			cappedforce0 = force[0]/force[0].modulus()*fmin(force[0].modulus(),0.01);
		}
		if (force[1].modulus()>0){
			cappedforce1 = force[1]/force[1].modulus()*fmin(force[1].modulus(),0.01);
		}
		pos[0] += cappedforce0 * dt / friction_trans;
		pos[1] += cappedforce1 * dt / friction_trans;

		cyto.dilute((4.0/3.0*r+length_before)/(4.0/3.0*r+get_length()));

		get_orientation(temp_vec);
		if (((cappedforce0 * dt / friction_trans).modulus())>0.1 ||
			((cappedforce1 * dt / friction_trans).modulus())>0.1){
			vec2d temp_vec; 
			get_orientation(temp_vec);
			// std::cout<<"Large displacement detected in cell "<<id<<":\n";
			// std::cout<<"Current spring force: "<< temp_vec/temp_vec.modulus()*springk*(length()-l0)/2<<'\n';
			// std::cout<<"Total force E1:"<< cappedforce0<<'\n';
			// std::cout<<"Total force E2:"<< cappedforce1<<'\n';
			// std::cout<<"Total displacement E1:"<< (cappedforce0 * dt / friction_trans).modulus()<<'\n';
			// std::cout<<"Total displacement E2:"<< (cappedforce1 * dt / friction_trans).modulus()<<'\n';
			return 1;
		}
		// std::cout<<"orientation after force"<<temp_vec<<'\n';
		//reset_force(); // force should be reset at the start of each step, so forces
		// are accessible at after being applied
		return 0;
	}

	void bacterium::grow(double dt){ // grow natural length of the pole-pole spring
		//l0 += fmax(0,growth_rate*dt*(1-(l0-length())/l0));
		l0 += growth_rate*cyto.get_growth_rate_modifier()*dt;
	}

	bool bacterium::division_ready(){
		//std::cout<<"division condition "<<length<<' '<<division_length<<'\n';
		return (get_length()>division_length);
	}

	void bacterium::get_daughter1_poles(vec2d& p1, vec2d& p2){
		p1 = pos[0];
		get_orientation(p2);
		p2 = pos[0] - p2*(1-r/p2.modulus());
	}

	void bacterium::get_daughter2_poles(vec2d& p1, vec2d& p2){
		p1 = pos[1];
		get_orientation(p2);
		p2 = pos[1] + p2*(1-r/p2.modulus());
	}

	bacterium bacterium::get_daughter1(int id){
		vec2d pole1,pole2;
		get_daughter1_poles(pole1,pole2);
		bacterium bb(id,r,pole1,pole2, growth_rate,2.0, mem, friction_trans, springk, cyto);
		return bb;
	}

	bacterium bacterium::get_daughter2(int id){
		vec2d pole1,pole2;
		get_daughter2_poles(pole1,pole2);
		bacterium bb(id,r,pole1,pole2, growth_rate,2.0, mem, friction_trans, springk, cyto);
		return bb;
	}

	std::string bacterium::get_str_physics(){
		std::stringstream ssphys;
		ssphys<<id<<' '<<type<<' '<<get_centre()[0]<<' '<<get_centre()[1]<<' ';
    	ssphys<<get_angle()<<' '<<get_length()<<' '<<l0<<' ';
    	ssphys<<force[0][0]<<' '<<force[0][1]<<' ';	
    	ssphys<<force[1][0]<<' '<<force[1][1];
    	return ssphys.str();
	}

void update_force_between(bacterium &b1, bacterium &b2){
	vec2d v1,v2,relation;
	vec2d contact;
	vec2d forcec;
	double length,l1,l2; // length of each bacterium (from pole to pole)
	int p1, p2; // pole of bac1 and pole of bac2 that are interacting
	vec2d pos1,pos2; // position of the centre of each bacterium
	// std::cout<<"calculating forces between bac "<< b1.id << " and " << b2.id<<" \n";
	// there are several possible configurations of the interaction
	b1.get_length(l1);
	b2.get_length(l2);
	b1.get_centre(pos1);
	b2.get_centre(pos2);
	if (((l1/2.0+b1.r+l2/2.0+b2.r)*(l1/2.0+b1.r+l2/2.0+b2.r))>(dist2(pos1,pos2))) // only get into the details of the force if two cells are close enough
	{
	// std::cout<<"cells are close!\n";
	// std::cout<<"bac 1: "<<b1.pole1()<<' '<<b1.pole2()<<")\n";
	// std::cout<<"bac 2: "<<b2.pole1()<<' '<<b2.pole2()<<")\n";

	bool pole_pole = false;// extreme-extreme contact is checked until it is found or discarded 
	for(int ip1=0;ip1<2;ip1++){ // possible poles of bac1
		for(int ip2=0;ip2<2;ip2++){ // possible poles of bac2
			v1 = b1.pos[ip1];
			v2 = b2.pos[ip2];
			if(dist2(v1,v2)<(4*b1.r*b1.r)){//extremes touching
				p1 = ip1;
				p2 = ip2;
				pole_pole = true;
				goto end_pole_pole_loop;
			}
		}
	}

	end_pole_pole_loop:
	if (pole_pole){
		// std::cout<<"force between extremes \n";
		// std::cout<<"vecmean "<<v1[0]<<' '<<v1[1]<<'\n';
		// std::cout<<"vecdiff "<<v2[0]<<' '<<v2[1]<<'\n';
		//contact = (v1+v2)/2; // global coordinate system
		relation = v1-v2; // vector joining both contact E points in the direction of bac1
		// std::cout<<"vecmean "<<contact[0]<<' '<<contact[1]<<'\n';
		// std::cout<<"vecdiff "<<relation[0]<<' '<<relation[1]<<'\n';
		forcec = relation/relation.modulus()*b1.mem*pow(abs((relation.modulus()-(b1.r+b2.r)))/b1.r,1);
		// std::cout<<"force on cell 1 "<<forcec[0]<<' '<<forcec[1]<<'\n';
		// std::cout<<"at point \n"<<contact[0]<<' '<<contact[1]<<'\n';
		b1.force[p1] += forcec;
		// forcec = -1*forcec;
		// b1.get_shift2glob(contact);
		// b2.get_glob2shift(contact); // contact point with respect of bac 2 center
		// std::cout<<"force on cell 2 "<<forcec[0]<<' '<<forcec[1]<<'\n';
		// std::cout<<"at point \n"<<contact[0]<<' '<<contact[1]<<'\n';
		b2.force[p2] -= forcec;
	}

	else{ // not extreme-extreme interaction
 
		bool side_pole_interaction = false;// side_pole contact is checked until it is found or discarded 
		bacterium *b_side,*b_pole;
		for(int bac_side=0;bac_side<2;bac_side++){ // bacterium that has a side interacting
			if (bac_side == 0){
				b_side = &b1;
				b_pole = &b2;			
			}
			else{
				b_side = &b2;
				b_pole = &b1;
			}
			length = b_side->get_length();
			for(int ip=0;ip<2;ip++){ // pole of the bacterium that is interacting with its pole
				v1 = b_pole->pos[ip]; // position of interacting pole
				b_side->get_glob2rot(v1); // position of interacting pole from ref systems of b_side rotated
	 			if (abs(v1[0])<0.5*length and abs(v1[1])<(b2.r+b1.r)){ // E1 of bac2 interacting with side of bac1
		 			//std::cout<<"force between s1 and body 2 (A) \n";
		 			contact[0] = v1[0];
		 			contact[1] = v1[1]*0.5; // contact point on bac1 ref system of bac1 rotated
		 			p1 = ip;
		 			side_pole_interaction = true;
		 			// std::cout<<"applied at b1 "<< contact[0]<<' '<<contact[1]<<'\n';
		 			goto end_side_pole_loop;
		 		}		
	 		}
	 	}

	 	end_side_pole_loop:
	 	if(side_pole_interaction){
	 		// std::cout<<"force between side and body \n";
			forcec[0] = 0;
			forcec[1] = -contact[1]/(abs(contact[1]))*b1.mem*pow(abs(abs(contact[1])-b1.r)/b1.r,1); // force in ref system of bac1 rotated
			b_side->get_rot2shift(forcec);
			b_side->force[0] += forcec*(contact[0]+length/2)/(length);// force on side is split on both poles proportionally to distance to application
			b_side->force[1] += forcec*(length/2-contact[0])/(length);	
			b_pole->force[p1] -= forcec;
		}
		else{
			// std::cout<<"but not close enough....\n";
		}
	}}
}