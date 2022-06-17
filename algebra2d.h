# ifndef ALGEBRA2D_H
# define ALGEBRA2D_H

#include<math.h>
#include<iostream>

//typedef std::valarray<double> vec2d[2];

struct vec2d {

	double v[2];

	vec2d(double x=0, double y=0){
		v[0]=x;
		v[1]=y;
	}

	double& operator[](const int index)
	{
		return v[index];
	}

	vec2d operator+(const vec2d& vec)
	{
		return vec2d(v[0]+vec.v[0], v[1]+vec.v[1]);
	}

	vec2d operator-(const vec2d& vec)
	{
		return vec2d(v[0]-vec.v[0], v[1]-vec.v[1]);
	}

	vec2d operator-()
	{
		return vec2d(-v[0], -v[1]);
	}


	double operator*(const vec2d& vec) // inner product
	{
		return v[0]*v[0]+v[1]*v[1];
	}

	vec2d operator*(double a)
	{
		return vec2d(v[0]*a,v[1]*a);
	}

	vec2d operator/(double a)
	{
		return vec2d(v[0]/a,v[1]/a);
	}


	vec2d& operator+=(const vec2d& vec){
		v[0]+=vec.v[0];
		v[1]+=vec.v[1];
		return *this;
	}

	vec2d& operator-=(const vec2d& vec){
		v[0]-=vec.v[0];
		v[1]-=vec.v[1];
		return *this;
	}

	vec2d operator%(double mod){
		// resize a vector to have the modulus mod
		// %1 returns an unitary vector
		double modprev_ = sqrt(v[0]*v[0]+v[1]*v[1]);
		if (modprev_ <= 0){
			return vec2d(0,0);
		}
		else{
			return vec2d(v[0]*mod/modprev_,v[1]*mod/modprev_);
		}
	}

	vec2d pow(const double power){
		return vec2d(::pow(v[0],power),::pow(v[1],power));
	}

	double sum(){
		return v[0]+v[1];
	}

	double modulus(){
		return sqrt(v[0]*v[0]+v[1]*v[1]);
	}

	double angle(){
		return atan(v[1]/v[0]);
	}

};

std::ostream& operator<<(std::ostream& os, const vec2d& vec);
double dist2(vec2d v1, vec2d v2);
double dist(vec2d v1, vec2d v2);
double angle(vec2d v1, vec2d v2);
double cosangle(vec2d v1, vec2d v2);
double sinangle(vec2d v1, vec2d v2);
double crossprod(vec2d v1, vec2d v2);
vec2d dist_linepoint(vec2d line, vec2d point);
vec2d vnormal(vec2d v1);
vec2d rotate(vec2d v1, double angle);

#endif