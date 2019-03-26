//****************************************************************************************
//
// SPH Code
//
// created	85/12 By Fatehi
// Last change	86/08 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************

#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#include "vectors.h"

enum Momentom_Type {MNONE, MINVISCID, MLAMINAR, MLAMINAR2, MTURBULENT};
enum Energy_Type {ENONE, EINVISCID, ELAMINAR, ETURBULENT};

class TOptions
{

public:

	double INITIALTIME;
	double FINALTIME;
	double MAXCFL;

	double XSPH;
	double Shifting;
	double ARTIFITIALVISCOSITY;
	double BULKVISCOSITY;

	unsigned short REPORTSTEP;
	//unsigned short SAVESTEP;
	double SAVINGTIME;

	unsigned short WRITINGMETHOD;

	unsigned short SEARCHMETHOD;

	unsigned short KERNELTYPE;


	unsigned short SOLUTIONMETHOD;


	unsigned short INTEGRATIONMETHOD;
	double FREESURFACE;


	double MAXIMUMX;
	double MAXIMUMY;

	double MAXIMUMVELOCITY;
	double MAXIMUMACCELERATION;
	double MAXIMUMPRESSURE;
	double MAXIMUMH;

	double MINIMUMX;
	double MINIMUMY;

	double MINIMUMPRESSURE;
	double MINIMUMH;
	double MINIMUMDISTANCE;



	Momentom_Type MOMENTUM;
	Energy_Type ENERGY;

	bool SUMDENSITY;
	bool NORMALDENSITY;
	bool SURFACETENSION;
	bool CORRECTGRADIENT;


	double C_0;							//(m/s)^2
	double V_0;								//m/s
	Vector BODYACC;		//Gravity Vector in m/s2

	char DATAFILE[500];
	char MATFILE[500];


	double SHIFT;


	void ReadFile();
	void Report();

};

#endif
