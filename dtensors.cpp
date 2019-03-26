//****************************************************************************************
//
// SPH Code
//
// created	88/01 By R. Fatehi
// Last change	88/01 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************



#include "dtensors.h"

DTensor& DTensor::operator=(const DTensor& p)
{
	ets[0]=p.ets[0];
	ets[1]=p.ets[1];
	return *this;
}

DTensor & DTensor::operator+=(const DTensor& p)
{
	//for(register int i=0;i<3;i++) 
	//	ets[i] += p.ets[i];
	ets[0] += p.ets[0];
	ets[1] += p.ets[1];
	return *this;
}

DTensor & DTensor::operator-=(const DTensor& p) 
{
	//for(register int i=0;i<3;i++) 
	//	ets[i] -= p.ets[i];
	ets[0] -= p.ets[0];
	ets[1] -= p.ets[1];
	return *this;
}

DTensor operator+(const DTensor & p) 
{
	return p;
}

DTensor operator-(const DTensor& p) 
{
	return DTensor(0.0, 0.0) - p;
}

DTensor operator+(const DTensor& p1, const DTensor & p2) 
{
	DTensor sum = p1;
	sum += p2;
	return sum;
}

DTensor operator-(const DTensor& p1, const DTensor& p2) 
{
	DTensor sum = p1;
	sum -= p2;
	return sum;
}

DTensor operator*(const double scalar, const DTensor & p) 
{
	DTensor tm;
	tm[0] = scalar*p.ets[0];
	tm[1] = scalar*p.ets[1];
	return tm;
}

DTensor operator*(const DTensor & p, const double scalar) 
{
	return scalar*p; 
}

Vector operator*(const DTensor &p, const Vector &v) 
{
	Vector tm;

	tm[0]=p.ets[0]*v.GetX();
	tm[1]=p.ets[1]*v.GetY();

	return tm;
}

Vector dot(const DTensor &st1, const DTensor &st2)
{
	Vector tm;

	tm[0] = st1.ets[0]*st2.ets[0];
	tm[1] = st1.ets[1]*st2.ets[1];
	return tm;
}

std::ostream&  operator<<(std::ostream& s, const DTensor& p ) 
{
	s << p.ets[0] << '\t' << p.ets[1]  ;

	return s;
}
double DiadProduct(const DTensor &T, const DTensor &E)
{
	return T.ets[0] * E.ets[0] +  T.ets[1] * E.ets[1];

}
