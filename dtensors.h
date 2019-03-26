//****************************************************************************************
//
// SPH Code
//
// created	86/01 By Fatehi
// Last change	86/05 By Fatehi
// Last change	97/11 By Nikooei

//****************************************************************************************


#ifndef _DTensor_H_
#define _DTensor_H_

#include "vectors.h"

class DTensor{
private:
	double ets[2];
public:
	DTensor(){ets[0]=ets[1]=0.0;}
	DTensor(const double nxx, const double nyy){ets[0]=nxx;ets[1]=nyy;}
	DTensor(const DTensor &p){ets[0]=p.ets[0];ets[1]=p.ets[1];}
	DTensor(const Vector &p){ets[0]=p.GetX();ets[1]=p.GetY();}
	DTensor(const double a){ets[0]=a;ets[1]=a;}

	double GetXX() const {return ets[0];}
	double GetYY() const {return ets[1];}
	double TwoNorm() const {return sqrt(ets[0]*ets[0]+ets[1]*ets[1]);}

	void SetXX(double nxx){ets[0]=nxx;}
	void SetYY(double nyy){ets[1]=nyy;}
	void SetAll(double nxx, double nyy){ets[0]=nxx;ets[1]=nyy;}

	DTensor& operator+=(const DTensor &);
	DTensor& operator-=(const DTensor &);

	friend DTensor operator+(const DTensor&);
	friend DTensor operator-(const DTensor&);
	friend DTensor operator+(const DTensor&, const DTensor&);
	friend DTensor operator-(const DTensor&, const DTensor&);
	friend DTensor operator*(const double, const DTensor&);
	friend DTensor operator*(const DTensor&, const double);
	friend Vector operator*(const DTensor&, const Vector&);
	friend Vector dot(const DTensor &, const DTensor &);
	friend double DiadProduct(const DTensor &, const DTensor &);

	double& operator[](int i) {return ets[i];}
	DTensor& operator=(const DTensor&);

	friend std::ostream& operator<<(std::ostream&, const DTensor&);
};

#endif
