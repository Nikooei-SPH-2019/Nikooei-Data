//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change		86/07 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************



#include "vectors.h"



Vector& Vector::operator=(const Vector & p)
{
	ets[0]=p.ets[0];
	ets[1]=p.ets[1];
	return *this;
}

double Vector::TwoNorm() const
{
	double norm = std::sqrt(ets[0]*ets[0]+ets[1]*ets[1]);
	return norm;
}

Vector & Vector::operator+=(const Vector & p)
{
	ets[0] += p.ets[0];
	ets[1] += p.ets[1];
	return *this;
}

Vector & Vector::operator-=(const Vector & p)
{
	ets[0] -= p.ets[0];
	ets[1] -= p.ets[1];
	return *this;
}

Vector operator*(const double scalar, const Vector & p)
{
	Vector tm;
	tm[0] = scalar*p.ets[0];
	tm[1] = scalar*p.ets[1];
	return tm;
}





Vector operator/(const Vector & p,const double scalar )
{
	Vector tm;
	tm[0] = p.ets[0]/ scalar;
	tm[1] = p.ets[1]/ scalar;
	return tm;
}


double Dot(const Vector & p1, const Vector & p2)
{
	double tm = p1.ets[0]*p2.ets[0] + p1.ets[1]*p2.ets[1];
	return tm;
}


double Vector::MaxNorm() const
{
	double norm = std::fabs(ets[0]);
	if (norm < std::fabs(ets[1]))
		norm = std::fabs(ets[1]);
	return norm;
}

std::ostream&  operator<<(std::ostream& s, const Vector& p )
{
	s << p.ets[0] << '\t' << p.ets[1];
	return s;
}

bool operator<=(const Vector& p1, const Vector& p2)
{
	return ((p1.ets[0]<=p2.ets[0]) &&
		(p1.ets[1]<=p2.ets[1])/*&&
		p1.ets[2]<=p2.ets[2]*/);
}

bool operator<(const Vector& p1, const Vector& p2)
{
	 return((p1.ets[0]<p2.ets[0]) &&
		(p1.ets[1]<p2.ets[1])/*&&
		p1.ets[2]<p2.ets[2]*/);
}

bool operator>=(const Vector& p1, const Vector& p2)
{
	 return((p1.ets[0]>=p2.ets[0]) &&
		(p1.ets[1]>=p2.ets[1])/*&&
		p1.ets[2]>=p2.ets[2]*/);
}

bool operator>(const Vector& p1, const Vector& p2)
{
	return ((p1.ets[0]>p2.ets[0]) &&
			(p1.ets[1]>p2.ets[1]) );
}

Vector Rotate90(const Vector& p)
{
	Vector tm=Vector(p.ets[1], -p.ets[0]);
	return tm;
}

Vector Rotate90unit(const Vector& p)           //unit vector of tangent           
{
	double mag = p.TwoNorm();
	Vector tm = Vector(p.ets[1]/mag, -p.ets[0]/mag);

	return tm;
}



Vector operator*(const Vector & p, const double scalar)
{
	return scalar*p;
}


Vector operator+(const Vector & p)
{
	return p;
}

Vector operator-(const Vector & p)
{
	return Vector() - p;
}

Vector operator+(const Vector & p1, const Vector & p2)
{
	Vector sum = p1;
	sum += p2;
	return sum;
}

double operator+(doublevec & p1, doublevec & p2)
{
	double sum ;
	sum = p1+p2;
	return sum;
}


Vector operator-(const Vector & p1, const Vector & p2)
{
	Vector sum = p1;
	sum -= p2;
	return sum;
}

double Cross2D(const Vector & p1, const Vector & p2)
{
	double tm = p1.ets[0]*p2.ets[1] - p1.ets[1]*p2.ets[0];
	return tm;
}


