//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change		85/12 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************

#include "stensors.h"

STensor& STensor::operator=(const STensor& p)
{
	ets[0] = p.ets[0];
	ets[1] = p.ets[1];
	ets[2] = p.ets[2];

	return *this;
}

STensor & STensor::operator+=(const STensor& p)
{
	ets[0] += p.ets[0];
	ets[1] += p.ets[1];
	ets[2] += p.ets[2];

	return *this;
}

STensor & STensor::operator-=(const STensor& p)
{
	ets[0] -= p.ets[0];
	ets[1] -= p.ets[1];
	ets[2] -= p.ets[2];

	return *this;
}

STensor operator+(const STensor & p)
{
	return p;
}

STensor operator-(const STensor& p)
{
	return STensor() - p;
}

STensor operator+(const STensor& p1, const STensor & p2)
{
	STensor sum = p1;
	sum += p2;
	return sum;
}

STensor operator-(const STensor& p1, const STensor& p2)
{
	STensor sum = p1;
	sum -= p2;
	return sum;
}

STensor operator*(const double scalar, const STensor & p)
{
	STensor tm;
	tm[0] = scalar*p.ets[0];
	tm[1] = scalar*p.ets[1];
	tm[2] = scalar*p.ets[2];

	return tm;
}

STensor operator*(const STensor & p, const double scalar)
{
	return scalar*p;
}

Vector operator*(const STensor &p, const Vector &v)
{
	Vector tm;

	tm[0] = p.ets[0] * v.GetX() + p.ets[1] * v.GetY();
	tm[1] = p.ets[1] * v.GetX() + p.ets[2] * v.GetY();

	return tm;
}

Vector operator*(const Vector &v, const STensor &p)
{
	Vector tm;

	tm[0] = p.ets[0] * v.GetX() + p.ets[1] * v.GetY();
	tm[1] = p.ets[1] * v.GetX() + p.ets[2] * v.GetY();

	return tm;
}


std::ostream&  operator<<(std::ostream& s, const STensor& p)
{
	s << p.ets[0] << '\t' << p.ets[1] << '\t' << p.ets[2]

		;

	return s;
}

double DiadProduct(const STensor &T, const STensor &E)
{

	return T.ets[0] * E.ets[0] + 2.0 * T.ets[1] * E.ets[1] + T.ets[2] * E.ets[2];

}
STensor STensor::Inverse() const
{

	const double det = ets[0] * ets[2] - ets[1] * ets[1];
	if (det)
	{
		STensor tm;
		tm.SetXX(ets[2]);
		tm.SetXY(-ets[1]);
		tm.SetYY(ets[0]);
		return 1.0 / det*tm;
	}
	else if (ets[0])
		return STensor(1.0 / ets[0], 0.0, 0.0);		//One dimensional
	else //if (ets[2])
		return STensor(0.0, 0.0, 1.0 / ets[2]);		//One dimensional

}
double STensor::Trace() const
{

	return ets[0] + ets[2];

}


double STensor::secondinvariant() const
{

	return (ets[1] * ets[1]) +0.5*(ets[0] * ets[0]) +0.5*(ets[2]*ets[2]);

}