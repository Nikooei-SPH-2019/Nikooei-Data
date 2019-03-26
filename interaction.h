//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change		85/12 By Fatehi
// Last change	97/11 By Nikooei

//****************************************************************************************

#ifndef _INTERACTION_H_
#define _INTERACTION_H_

#include "vectors.h"

class Interaction{
private:
	int iindex;
	int jindex;
	double dist;
	double hmean;
	Vector vdist;
	double w;
	double dwdr;

public:
	Interaction();

	Interaction(int i, int j,double hm, double d, Vector dx, double k, double dk)
	{
		iindex=i;
		jindex=j;
		hmean=hm;
		dist=d;
		vdist=dx;
		w=k;
		dwdr=dk;
	};

	Interaction(const Interaction &inter)
	{
		iindex=inter.iindex;
		jindex=inter.jindex;
		hmean=inter.hmean;
		dist=inter.dist;
		vdist=inter.vdist;
		w=inter.w;
		dwdr=inter.dwdr;
	};

	Interaction& operator=(const Interaction &inter)
	{
		iindex=inter.iindex;
		jindex=inter.jindex;
		hmean=inter.hmean;
		dist=inter.dist;
		vdist=inter.vdist;
		w=inter.w;
		dwdr=inter.dwdr;

		return *this;
	};
	


	void SetI(int i){iindex=i;}
	void SetJ(int j){jindex=j;}
	void SetDist(double l){dist=l;}
	void SetHMean(double hm){ hmean=hm;}
	void SetVD(Vector d){vdist=d;}
	void SetW(double k){w=k;}
	void SetDW(double dk){dwdr=dk;}

	int GetI(){return iindex;}
	int GetJ(){return jindex;}
	double GetDist(){return dist;}
	double GetHMean(){return hmean;}
	Vector GetVD(){return vdist;}
	double GetW(){return w;}
	double GetDW(){return dwdr;}
};

#endif
