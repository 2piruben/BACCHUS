
/////////////////////////////////////////////////////////////////////////
//
// SIMULATION OF GROWING BACTERIAL COLONIES
//
//							R Perez-Carrasco
//							Created: 12 MAY 2022
//////////////////////////////////////////////////////////////////////


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <list>
#include <algorithm>
#include <math.h>
#include <unistd.h>

//#include "bacterium.h" # eventually the class will go to a different file

typedef std::array<double,2> vec2d;

double dist2(vec2d v1, vec2d v2){
 	return ((v1[0]-v2[0])*(v1[0]-v2[0]) + (v1[1]-v2[1])*(v1[1]-v2[1])); 	
 }

double angle(vec2d v1, vec2d v2){
 	return atan((v2[1]-v1[1])/(v2[0]-v1[0]));
}

double modulus(vec2d v1){
	return sqrt(v1[0]*v1[0] + v1[1]*v1[1]);
}

double cosangle(vec2d v1, vec2d v2){
 	return (v1[0]*v2[0]+v1[1]*v2[1])/(modulus(v1)*modulus(v2));
}

double sinangle(vec2d v1, vec2d v2){
	// moving from v1 to v2 
	return (v1[0]*v2[1]-v2[0]*v1[1])/(modulus(v1)*modulus(v2));
}

vec2d dist_linepoint(vec2d line, vec2d point){ 
	vec2d distvec;
// distsance between a line (given by a vector) and a point.
// The result is two numbers: the distance perpendicular to the line, and the distance along the line. 
	distvec[0] = modulus(point)*sinangle(line,point); // positive when point is at a positive angle from line
	distvec[1] = modulus(point)*cosangle(line,point);  
	return distvec;
}

// }

vec2d vecmean(vec2d v1, vec2d v2){
	vec2d vmean;
	vmean[0] = (v2[0] + v1[0])*0.5;
	vmean[1] = (v2[1] + v1[1])*0.5;
	return vmean;
}

vec2d vecsum(vec2d v1, vec2d v2){
	vec2d vecsum;
	vecsum[0] = (v1[0] + v2[0]);
	vecsum[1] = (v1[1] + v2[1]);
	return vecsum;
}

vec2d vecdiff(vec2d v1, vec2d v2){
	vec2d vecsum;
	vecsum[0] = (v1[0] - v2[0]);
	vecsum[1] = (v1[1] - v2[1]);
	return vecsum;
}


vec2d vscaled(vec2d v1, double l){
	vec2d vscaled;
	vscaled[0] = v1[0]*l;
	vscaled[1] = v1[1]*l;
	return vscaled;
}

vec2d vnormal(vec2d v1){
	// Rotation +90C
	vec2d vn;
	vn[0] = -v1[1];
	vn[1] = v1[0];
	return vn;
}






class bacterium{

protected:

	int id;
	vec2d pos;
	double angle;
	double length; // this is the length between the extreme points of the central rectangle (the actual length of the bacterium is length + 2*radius)
	double r;
	double age;
	double growth_rate;
	double division_length;
	vec2d force; // force being applied on the centre of mass of the bacterium
	double torque; // torque applied on the bacterium
	double friction_trans;
	double friction_rot;

public:

	bacterium(int id_, double r_, vec2d pos_, double angle_, double length_, double growth_rate_,
		double division_length_, double friction_trans_, double friction_rot_){
		std::cout<<"created bacteria "<<id_<<" at "<<pos_[0]<<' '<<pos_[1]<<'\n';
		id  = id_;
		pos[0]  = pos_[0];
		pos[1]  = pos_[1];
		angle  = angle_;
		length = length_;
		r = r_;
		age = 0;
		division_length = division_length_;
		growth_rate = growth_rate_;
		friction_trans = friction_trans_;
		friction_rot = friction_rot_;
		reset_force();
	}

	double angle_bac(){
		return angle;
	}

	double length_bac(){
		return length;
	}

	vec2d pos_bac(){
		return pos;
	}

	vec2d force_bac(){
		return force;
	}

	int id_bac(){
		return id;
	}

	double radius_bac(){
		return r;
	}

	double friction_trans_bac(){
		return friction_trans;
	}

	double friction_rot_bac(){
		return friction_rot;
	}

	void get_extreme1(vec2d &output){
		output[0] = pos[0] + 0.5*length*cos(angle);
		output[1] = pos[1] + 0.5*length*sin(angle);
	}

	void get_extreme2(vec2d &output){
		output[0] = pos[0] - 0.5*length*cos(angle);
		output[1] = pos[1] - 0.5*length*sin(angle);
	}

	void get_orientation(vec2d &output){
		// vector that goes from centre of bacterium to extreme1
		output[0] = length*cos(angle)/2;
		output[1] = length*sin(angle)/2;
	}

	// We want functions that transform points to 3 different coordinate systems
	// glob -> with respect to the dish (centered at 0,0)
	// shift -> with respect to the center of the bacterium
	// rot  -> with respect to the center of the bacterium, rotated so extreme 1 falls in in the direction (1,0)

	void get_glob2rot(vec2d &output){
		// Given a vector, return the vector in the reference system of the bacterium orientated (E1 in the cartesian direction, centered at the bacterium)
		double tempx,tempy;
		tempx = output[0]-pos[0];
		tempy = output[1]-pos[1];
		output[0] = tempx*cos(angle) + tempy*sin(angle);
		output[1] = -tempx*sin(angle) + tempy*cos(angle);
	}

	void get_glob2shift(vec2d &output){
		// Given a global vector, return the vector in the reference system of the bacterium preserving the orientation
		output[0] -= pos[0];
		output[1] -= pos[1];
	}

	void get_rot2glob(vec2d &output){
		// Opposite to vec2bac, given a vector in the reference system of the bacterium orientated, return the vector in the reference system of the dish
		double tempx,tempy;
		tempx = output[0];
		tempy = output[1];
		output[0] = tempx*cos(angle) - tempy*sin(angle) + pos[0];
		output[1] = tempx*sin(angle) + tempy*cos(angle) + pos[1];
	}

	void get_shift2glob(vec2d &output){
		// Opposite to vec2bacpos, given a vector ceneterd in the bacterium, return the vector in the global coordinate system
		output[0] += pos[0];
		output[1] += pos[1];
	}

	void get_rot2shift(vec2d &output){
		double tempx,tempy;
		tempx = output[0];
		tempy = output[1];
		output[0] = tempx*cos(angle) - tempy*sin(angle);
		output[1] = tempx*sin(angle) + tempy*cos(angle)	;	
	}

	void get_shift2rot(vec2d &output){
		double tempx,tempy;
		tempx = output[0];
		tempy = output[1];
		output[0] = tempx*cos(angle) + tempy*sin(angle);
		output[1] = -tempx*sin(angle) + tempy*cos(angle);		
	}


	void move(vec2d Dx){
		pos[0] += Dx[0];
		pos[1] += Dx[1];
	}

	void update_force(vec2d pos, vec2d forc){ // update force and torque by incuding a force forc applied at position pos
		double cosforcpos = cosangle(pos,forc);
		vec2d forcetrans;
		forcetrans = vscaled(pos,modulus(forc)/modulus(pos)*cosforcpos);
		force = vecsum(force, forcetrans);// translation force applied at the center of mass
		torque -= forc[0]*pos[1]-forc[1]*pos[0];
		std::cout<<"linear force :"<<force[0]<<' '<<force[1]<<" \n";
 		std::cout<<"torque :"<<torque<<" \n";
	}

	void reset_force(){ // update force and torque by incuding a force forc applied at position pos
		force[0] = 0;
		force[1] = 0;
		torque = 0;
	}

	void apply_force(double dt, double maxd){ // Apply the stored force and move and rotate the bacterium
		// displacements bigger than maxd are not allowed
		vec2d displacement;
		std::cout<<"Applying force on cell "<<id<<" at pos:"<<pos[0]<<" "<<pos[1]<<" with force: "<<force[0]<<" "<<force[1];
		std::cout<<"and pars "<<dt<<" "<<friction_trans;
		displacement = vscaled(force,dt/friction_trans);
		displacement = vscaled(displacement, fmin(modulus(displacement),maxd));
		move(displacement);
		//std::cout<<"updated to: "<<pos[0]<<' '<<pos[1]<<'\n';
		angle += fmin(torque/friction_rot*dt,maxd/r);
		
	}

	void grow(double dt){
		length += growth_rate*dt;
	}

	bool division_ready(){
		//std::cout<<"division condition "<<length<<' '<<division_length<<'\n';
		return (length>division_length);
	}

	vec2d pos_daughter1(){
		vec2d daughter;
		daughter[0] = pos[0] + 0.25*(length+2*r)*cos(angle);
		daughter[1] = pos[1] + 0.25*(length+2*r)*sin(angle);
		return daughter;
	}

	vec2d pos_daughter2(){
		vec2d daughter;
		daughter[0] = pos[0] - 0.25*(length+2*r)*cos(angle);
		daughter[1] = pos[1] - 0.25*(length+2*r)*sin(angle);
		return daughter;
	}


	friend void update_force_between(bacterium &b1, bacterium &b2);


};



void update_force_between(bacterium &b1, bacterium &b2){
	vec2d v1,v2,relation;
	vec2d contact;
	vec2d forcec;
	std::cout<<"calculating forces between bac "<< b1.id << " and " << b2.id<<" \n";
	// there are several possible configurations of the interaction
	bool extreme_extreme = true;// extreme-extreme contact is checked until it is found or discarded


	if (((b1.length/2.0+b1.r+b2.length/2.0+b2.r)*(b1.length/2.0+b1.r+b2.length/2.0+b2.r))>(dist2(b1.pos,b2.pos))) // only get into the details of the force if two cells are close enough
	{
	std::cout<<"cells are close!\n";
	b1.get_extreme1(v1);	
	b1.get_extreme2(v2);
	std::cout<<"4 corners: "<< v1[0]<<' '<< v1[1]<<' '<< v2[0]<<' '<< v2[1]<<'\n';
	b2.get_extreme1(v1);	
	b2.get_extreme2(v2);
	std::cout<< v1[0]<<' '<< v1[1]<<' '<< v2[0]<<' '<< v2[1]<<'\n';
	b1.get_extreme2(v1);	
	b2.get_extreme1(v2);
	std::cout<< dist2(v1,v2)<<'\n';




	b1.get_extreme1(v1);	
	b2.get_extreme1(v2);
 	if (dist2(v1,v2)>(4*b1.r*b1.r)){ //extremes 1-1 not touching
		std::cout<<"checkpoint1 \n";
		b1.get_extreme2(v1);	
		b2.get_extreme2(v2);
		if(dist2(v1,v2)>(4*b1.r*b1.r)){//extremes 2-2 not touching
			b1.get_extreme1(v1);	
			b2.get_extreme2(v2);
			if(dist2(v1,v2)>(4*b1.r*b1.r)){//extremes 1-2 not touching
				b1.get_extreme2(v1);	
				b2.get_extreme1(v2);
				if(dist2(v1,v2)>(4*b1.r*b1.r)){//extremes 2-1 not touching
					extreme_extreme = false;
				}
			}
		}
	}
 	if(extreme_extreme){
 		// std::cout<<"force between extremes \n";
 		// std::cout<<"vecmean "<<v1[0]<<' '<<v1[1]<<'\n';
 		// std::cout<<"vecdiff "<<v2[0]<<' '<<v2[1]<<'\n';
 		contact = vecmean(v1,v2); // global coordinate system
 		relation = vecdiff(v1,v2); // vector joining both contact E points in the direction of bac1
 		// std::cout<<"vecmean "<<contact[0]<<' '<<contact[1]<<'\n';
 		// std::cout<<"vecdiff "<<relation[0]<<' '<<relation[1]<<'\n';
 		forcec = vscaled(relation,1.0/modulus(relation)*pow(modulus(relation)-(b1.r+b2.r),2));
 		// std::cout<<"force on cell 1 "<<forcec[0]<<' '<<forcec[1]<<'\n';
 		// std::cout<<"at point \n"<<contact[0]<<' '<<contact[1]<<'\n';
 		b1.get_glob2shift(contact); 
 		b1.update_force(contact,forcec);
 		forcec = vscaled(forcec,-1);
 		b1.get_shift2glob(contact);
 		b2.get_glob2shift(contact); // contact point with respect of bac 2 center
 		// std::cout<<"force on cell 2 "<<forcec[0]<<' '<<forcec[1]<<'\n';
 		// std::cout<<"at point \n"<<contact[0]<<' '<<contact[1]<<'\n';
 		b2.update_force(contact,forcec);
 	}

	else{ // not extreme-extreme interaction
 	
 	// Checking if the body of bac1 intereacts with an extreme of bac2	
		bool ex1_bac2 = false;
		//std::cout<<"checkpoint2 \n";
 		b2.get_extreme1(v1);
 		b2.get_extreme2(v2);
 		b1.get_glob2rot(v1); // position of E1 of bac 2 in the reference system of bac1 rotated
 		b1.get_glob2rot(v2); // position of E2 of bac 2 in the reference system of bac1 rotated
 		//std::cout<<"checkpoint3 \n";
 		if (abs(v1[0])<0.5*b1.length and abs(v1[1])<(b2.r+b1.r)){ // E1 of bac2 interacting with side of bac1
 			std::cout<<"force between side 1 and body 2 (A) \n";
 			contact[0] = v1[0];
 			contact[1] = v1[1]*0.5; // contact point on bac1 ref system of bac1 rotated
 			ex1_bac2 = true;
 		}

 		else if (abs(v2[0])<0.5*b1.length and abs(v2[1])<(b2.r+b1.r)){ // E1 of bac2 interacting with side of bac1
 			std::cout<<"force between side 1 and body 2 (B) \n";
 			contact[0] = v2[0];
 			contact[1] = v2[1]*0.5; // contact point on bac1 ref system of bac1 rotated
 			ex1_bac2 = true;
 		}

 		if(ex1_bac2){
 			std::cout<<"force between side 1 and body 2 \n";
 			forcec[0] = 0;
 			forcec[1] = -contact[1]/(abs(contact[1]))*pow(contact[1]-b1.r,2); // force in ref system of bac1 rotated

 			b1.get_rot2shift(forcec);
 			b1.get_rot2shift(contact);
 			b1.update_force(contact,forcec);

 			b1.get_shift2glob(contact);
 			b2.get_glob2shift(contact);// contact in b2 coordinate
 			forcec = vscaled(forcec,-1); // reaction force
 			b2.update_force(contact,forcec);
 		}

 		else {// not extreme1-body2 then check extreme2-body1	
			bool ex2_bac1 = false;
			//std::cout<<"checkpoint4 \n";
	 		b1.get_extreme1(v1);
	 		b1.get_extreme2(v2);
	 		b2.get_glob2rot(v1); // position of E1 of bac 1 in the reference system of bac2
	 		b2.get_glob2rot(v2); // position of E2 of bac 1 in the reference system of bac2
	 		//std::cout<<"checkpoint5 \n";
 			if (abs(v1[0])<0.5*b2.length and abs(v1[1])<(b1.r+b2.r)){ // E1 of bac1 interacting with side of bac2
 				contact[0] = v1[0];
 				contact[1] = v1[1]*0.5; // contact point on bac2 ref system of bac2
 				ex2_bac1 = true;
 				//std::cout<<"checkpoint6 \n";
 			}

 			else if (abs(v2[0])<0.5*b2.length and abs(v2[1])<(b1.r+b2.r)){ // E1 of bac1 interacting with side of bac2
 				contact[0] = v2[0];
 				contact[1] = v2[1]*0.5; // contact point on bac2 ref system of bac2
 				ex2_bac1 = true;
 			}

 			if(ex2_bac1){
 				std::cout<<"force between side 2 and body 1 \n";
 				forcec[0] = 0;
 				forcec[1] = -contact[1]/(abs(contact[1]))*pow(contact[1]-b1.r,2); // force in ref system of bac1
	 			b2.get_rot2shift(forcec);
	 			b2.get_rot2shift(contact);
	 			b2.update_force(contact,forcec);

	 			b2.get_shift2glob(contact);
	 			b1.get_glob2shift(contact);// contact in b2 coordinate
	 			forcec = vscaled(forcec,-1); // reaction force
	 			b1.update_force(contact,forcec);
		 	}

 			else{ // not extreme-side interaction. Last possibility is the rare body-body interaction with two parallel bacs}
 				if ((b1.angle - b2.angle)<0.01){ // if they are parallel 
 					//std::cout<<"parallel force \n";
 					b1.get_extreme1(v1);
 			    	b1.get_extreme2(v2);
 			    	//std::cout<< "before rot" << v1[0] <<' '<< v2[1] <<' '<< b2.pos[0]<<' '<< b2.pos[1]<< ' '<< b2.length/2 <<' '<< b2.angle << '\n';
 					b2.get_glob2rot(v1); // position of E1 of bac 1 in the reference system of bac2
 					b2.get_glob2rot(v2); // position of E2 of bac 1 in the reference system of bac2

 					// std::cout<< "after rot"<< abs(v1[0]) <<' '<< b2.length/2 <<' '<< abs(v1[1]) <<'\n';
 					if(abs(v1[1])<(b1.r+b2.r)){ //  if they are close in the normal direction
 						if ((abs(v1[0])<b2.length/2) ||  (abs(v2[0])<b2.length/2)){ // and close in the parallel one 
 							std::cout<<"parallel force\n";
							contact[0]=0;
							contact[1]=v1[1]/abs(v1[1])*b2.r; // in case of full contct we will assume a normal force in the normal axis of the bac
		 					forcec = vscaled(contact,-1.0/modulus(contact)*pow(modulus(contact)-b2.r,2));
		 					b2.get_rot2shift(contact);
		 					b2.get_rot2shift(forcec);
 							b2.update_force(contact,forcec);

		 					// change parallel force to b1 to also be applied normally?

 							b2.get_shift2glob(contact);
 							b1.get_glob2shift(contact);
 							forcec = vscaled(forcec,-1);
 							b1.update_force(contact,forcec);

 							std::cout<<"contact :"<<contact[0]<<' '<<contact[1]<<" \n";
 							std::cout<<"force :"<<forcec[0]<<' '<<forcec[1]<<" \n";
						}
 					}
 				}
 			}
 		}
 	}}
}

typedef std::unique_ptr<bacterium> p_bacterium;

class population{


	private:

	std::list<bacterium> cells; // this contains all cells alive or dead
	std::list<bacterium*> cells_alive; // points to the cells in cells_all that are alive
	std::list<bacterium*> cells_dead; // poitns to the cells in cells_all that are dead

	gsl_rng * rng; // allocator for the rng generator
	
	double time = 0;

	double mean_growth_rate = 1;
	double deviation_growth_rate = 1;
	double mean_division_length = 1;
	double timestep = 0.001;
	double length = 1;
	int id = 0;
	double r = 0.2; // bacterium radius
	double maxd = mean_growth_rate*timestep; // maximum displacement allowed 
	double friction_trans;
	double friction_rot;
	std::ofstream trajfile; // File for output traj
	std::ifstream parsfile; // File for output traj
	std::stringstream ofilename;
	std::stringstream parsline;

	public:

	population(int seed, std::string inputfile){
		rng = gsl_rng_alloc (gsl_rng_mt19937);
		if(seed==-1){ // if seed==-1 (default) take a random number
			gsl_rng_set (rng,::time(NULL)*getpid());
			std::cout<<"RNG Seed used: "<<::time(NULL)*getpid()<<'\n';
		}
		else{
			gsl_rng_set (rng,seed);
		}

		parsfile.open("pars.in");
		std::string name, value;
		while (parsfile >> name >> value){
			if (name=="timestep"){ timestep = std::stod(value);};
			if (name=="r"){ r = std::stod(value);};
			if (name=="length"){ length = std::stod(value);};
			if (name=="growth_rate"){ mean_growth_rate = std::stod(value);};
			if (name=="friction_trans"){ friction_trans = std::stod(value);};
			if (name=="friction_rot"){ friction_rot = std::stod(value);};
		}
		std::cout<<"READING! "<<timestep<<' '<<r<<' '<<length<<' '<<mean_growth_rate<<' '<<friction_trans<<' '<<friction_rot<<' '<<'\n';
	}
	int next_id(){
		id += 1;
		return id;
	}

	void initialize_two(){ // start population with two parallel bacteria
		vec2d posinit;
		double angleinit;
		posinit[0] = 0;
		posinit[1] = 0;
		angleinit = 0.78;
		std::cout<<"Init bac 1 "<< posinit[0] <<' '<< posinit[1]<<'\n';
		cells.emplace_back(next_id(), r, posinit, angleinit, length, mean_growth_rate, 2.0,friction_trans,friction_rot);
		posinit[0] = 1.1;
		posinit[1] = 1.1;
		std::cout<<"Init bac 2 "<< posinit[0] <<' '<< posinit[1]<<'\n';
		cells.emplace_back(next_id(), r, posinit, angleinit, length, mean_growth_rate, 2.0,friction_trans,friction_rot);
		for(std::list<bacterium>::iterator cell_ptr = cells.begin(); cell_ptr != cells.end(); ++cell_ptr){
		// adding pointer to the cells to the alive list
			cells_alive.push_back(&*cell_ptr);	
		}
	}

	void evolve(){ // evolve the popualtion a time step dt
		/// growth and division
		std::cout<<"################# time: "<<time<<'\n';
		for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
			(*cell_ptr)->reset_force();
			// growth
			(*cell_ptr)->grow(timestep);
			if ((*cell_ptr)->division_ready()) 
    		{// division
    			std::cout<<"Dividing cell at "<<(*cell_ptr)->pos_bac()[0]<<' '<<(*cell_ptr)->pos_bac()[1]<<" \n";
    			cells_dead.push_back(*cell_ptr);
    			cells.emplace_back(next_id(),r,(*cell_ptr)->pos_daughter1(),(*cell_ptr)->angle_bac(), ((*cell_ptr)->length_bac()-2.0*r)/2.0,1.0,2.0,
    				(*cell_ptr)->friction_trans_bac(),(*cell_ptr)->friction_rot_bac());
    			std::cout<<"Created cell at"<<cells.back().pos_bac()[0]<<' '<<cells.back().pos_bac()[1]<<" \n";
    			cells_alive.push_back(&cells.back());
    			cells.emplace_back(next_id(),r,(*cell_ptr)->pos_daughter2(),(*cell_ptr)->angle_bac(), ((*cell_ptr)->length_bac()-2.0*r)/2.0,1.0,2.0,
    				(*cell_ptr)->friction_trans_bac(),(*cell_ptr)->friction_rot_bac());
    			std::cout<<"Created cell at"<<cells.back().pos_bac()[0]<<' '<<cells.back().pos_bac()[1]<<" \n";
    			cells_alive.push_back(&cells.back());
    			cells_alive.erase(cell_ptr);
    		}
    	}
    	//std::cout<<"calculating forces\n";
		/// calculate forces
		for(std::list<bacterium*>::iterator cell_ptr1 = cells_alive.begin(); cell_ptr1 != cells_alive.end(); ++cell_ptr1) {
			for(std::list<bacterium*>::iterator cell_ptr2 = cells_alive.begin(); cell_ptr1 != cell_ptr2; ++cell_ptr2) {
				update_force_between(**cell_ptr1,**cell_ptr2);
    		}
    	}
		/// move cells
		//std::cout<<"applying forces\n";
		for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
			//std::cout<<((*cell_ptr)->force_bac())[0];
			(*cell_ptr)->apply_force(timestep,maxd);
		}
		time += timestep;
    }

    void print_population(){
    	for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
    		std::cout<<(*cell_ptr)->id_bac()<<' '<<(*cell_ptr)->pos_bac()[0]<<' '<<(*cell_ptr)->pos_bac()[1]<<'\n';
    	}
    }

    void save_population(){
    	ofilename.str("");
    	ofilename<<"population"<<std::setprecision(5)<<time<<".out";
    	trajfile.open(ofilename.str()); // output trajectory
    	for(std::list<bacterium*>::iterator cell_ptr = cells_alive.begin(); cell_ptr != cells_alive.end(); ++cell_ptr) {
    		trajfile<<(*cell_ptr)->id_bac()<<' '<<(*cell_ptr)->pos_bac()[0]<<' '<<(*cell_ptr)->pos_bac()[1]<<' ';
    		trajfile<<(*cell_ptr)->angle_bac()<<' '<<(*cell_ptr)->length_bac()<<' '<<(*cell_ptr)->radius_bac()<<' ';
    		trajfile<<(*cell_ptr)->force_bac()[0]<<' '<<(*cell_ptr)->force_bac()[1]<<'\n';	
    	}
    	trajfile.close();
    }
};


int main(int argc, char* argv[]){
 	//std::cout<<"starting\n";
	population colony(10,"pars.in");
	//std::cout<<"initizalizing\n";
	colony.initialize_two();
	//std::cout<<"saving\n";
	colony.save_population();
	//std::cout<<"integrating\n";
	for(int i=0; i<20; i++){
		//std::cout<<"iteration "<< i << " \n";
		colony.evolve();
		if (i%100==0){
			colony.save_population();
		}
	}
	//std::cout<<"done\n";
	return 0;

}
// struct force{
// 	double parallel; // resulting force parallel to the axis of the bacterium
// 	double normal; // resulting normal force applied at extreme 1 of the bacterium (torque)
// }
