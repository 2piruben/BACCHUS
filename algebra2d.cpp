#include "algebra2d.h"


std::ostream& operator<<(std::ostream& os, const vec2d& vec)
{
    os <<'('<<vec.v[0]<<','<<vec.v[1]<<')';
    return os;
}

double dist2(vec2d v1, vec2d v2){
	return ((v1-v2).pow(2)).sum();
 }

double dist(vec2d v1, vec2d v2){
	return sqrt((v1-v2).pow(2).sum());
 }


double angle(vec2d v1, vec2d v2){// angle from vector v1 to v2
 	return atan((v2[1]-v1[1])/(v2[0]-v1[0]));
}

double cosangle(vec2d v1, vec2d v2){
 	return (v1*v2)/(v1.modulus()*v2.modulus());
}

double sinangle(vec2d v1, vec2d v2){
	// moving from v1 to v2 
	return (v1[0]*v2[1]-v2[0]*v1[1])/(v1.modulus()*v2.modulus());
}

double crossprod(vec2d v1, vec2d v2){
	return (v1[0]*v2[1]-v2[0]*v1[1]);	
}

vec2d dist_linepoint(vec2d line, vec2d point){ 
	vec2d distvec;
// distsance between a line (given by a vector) and a point.
// The result is two numbers: the distance perpendicular to the line, and the distance along the line. 
	distvec[0] = point.modulus()*sinangle(line,point); // positive when point is at a positive angle from line
	distvec[1] = point.modulus()*cosangle(line,point);  
	return distvec;
}

// }

vec2d vnormal(vec2d v1){
	// Rotation +90C
	vec2d vn;
	vn[0] = -v1[1];
	vn[1] = v1[0];
	return vn;
}


vec2d rotate(vec2d v1, double angle){
	vec2d vecrot;
	vecrot[0] = v1[0]*cos(angle) - v1[1]*sin(angle);
	vecrot[1] = v1[0]*sin(angle) + v1[1]*cos(angle);
	return vecrot;
}

