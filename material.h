//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change		86/09 By Fatehi
// Last change		87/11 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************


#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "vectors.h"

class TMaterial{
public:

	short type;
	
	TMaterial(): type(-1){}
	TMaterial( short i): type(i) {}
};

class TRMaterial : public TMaterial{
public:

	double rho;			//Density				(kg/m3)
	double k;			//Conductivity			(W/mK)
	double cv;			//Constant Volume Heat Capacity	(J/kgK)
	
	TRMaterial(): rho(0.0), k(0.0), cv(0.0), TMaterial() {}
	TRMaterial( short i, double r, double ka, double c): rho(r), k(ka), cv(c), TMaterial(i){}
};

class TFluid : public TRMaterial{
public:
	double mu;			//Dynamic Viscosity		(Pa.s)

	TFluid():mu(0.0), TRMaterial(){}
	TFluid( short i, double r, double ka, double c, double m /*,double t*/ ):mu(m)/*, Phi(t)*/, TRMaterial( i, r, ka, c){}
};

class TSolid : public TRMaterial{
public:
	
	double E;

	TSolid():E(0.0), TRMaterial(){}
	TSolid( short i, double r, double ka, double c, double ee):E(ee), TRMaterial(i, r, ka, c){}
};
////////////////////////////////////////////////////////////////Nikooei
class TNonNewtonian : public TRMaterial{
public:
	double mu, yieldstress, Consistency, npower, criticalstrainrate;
		TNonNewtonian() :mu(0.0), yieldstress(0.0), TRMaterial() {}
		TNonNewtonian(short i, double r, double ka, double c, double m,double y, double K, double n) :mu(m), yieldstress(y), Consistency(K), npower(n), TRMaterial(i, r, ka, c) {}
	
};
////////////////////////////////////////////////////////////////Nikooei

class TMat{
public:

	std::vector<TMaterial*> materials;
	std::vector< std::vector<double> > sigmaij;


	TMat(){}

	
	void Read( const char* MFile, bool &surfaceopt);


};

int FindType(std::vector<TMaterial*> & materials, int ptype);



#endif
