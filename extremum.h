//****************************************************************************************
//
// SPH Code
//
// created	86/01 By Fatehi
// Last change	97/11 By Nikooei

//****************************************************************************************
#ifndef _EXTREMUM_H_
#define _EXTREMUM_H_

#include "vectors.h"

class Extremum{
private:


public:

	Vector loc;		//Location
	double vel;		//Velocity Magnitude
	double p;		//Pressure
	double h;		//Smooting length

	double acc;		//Acceleration Magnitude

	double distance;

	double vdivergence;	//Velocity divergence

};


#endif //_EXTREMUM_H_
