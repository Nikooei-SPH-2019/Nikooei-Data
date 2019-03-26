//****************************************************************************************
//
// SPH Code
//
// created	89/03 By R. Fatehi
// changed		89/03 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************


#include "tensors4.h"

Tensor4 & Tensor4::operator+=(const Tensor4& p)
{
	ets[0][0][0][0]+=p.ets[0][0][0][0];
	ets[0][0][0][1]+=p.ets[0][0][0][1];
	ets[0][0][1][0]+=p.ets[0][0][1][0];
	ets[0][0][1][1]+=p.ets[0][0][1][1];
	ets[0][1][0][0]+=p.ets[0][1][0][0];
	ets[0][1][0][1]+=p.ets[0][1][0][1];
	ets[0][1][1][0]+=p.ets[0][1][1][0];
	ets[0][1][1][1]+=p.ets[0][1][1][1];
	ets[1][0][0][0]+=p.ets[1][0][0][0];
	ets[1][0][0][1]+=p.ets[1][0][0][1];
	ets[1][0][1][0]+=p.ets[1][0][1][0];
	ets[1][0][1][1]+=p.ets[1][0][1][1];
	ets[1][1][0][0]+=p.ets[1][1][0][0];
	ets[1][1][0][1]+=p.ets[1][1][0][1];
	ets[1][1][1][0]+=p.ets[1][1][1][0];
	ets[1][1][1][1]+=p.ets[1][1][1][1];

	return *this;
}

Tensor4 & Tensor4::operator-=(const Tensor4& p)
{
	ets[0][0][0][0]-=p.ets[0][0][0][0];
	ets[0][0][0][1]-=p.ets[0][0][0][1];
	ets[0][0][1][0]-=p.ets[0][0][1][0];
	ets[0][0][1][1]-=p.ets[0][0][1][1];
	ets[0][1][0][0]-=p.ets[0][1][0][0];
	ets[0][1][0][1]-=p.ets[0][1][0][1];
	ets[0][1][1][0]-=p.ets[0][1][1][0];
	ets[0][1][1][1]-=p.ets[0][1][1][1];
	ets[1][0][0][0]-=p.ets[1][0][0][0];
	ets[1][0][0][1]-=p.ets[1][0][0][1];
	ets[1][0][1][0]-=p.ets[1][0][1][0];
	ets[1][0][1][1]-=p.ets[1][0][1][1];
	ets[1][1][0][0]-=p.ets[1][1][0][0];
	ets[1][1][0][1]-=p.ets[1][1][0][1];
	ets[1][1][1][0]-=p.ets[1][1][1][0];
	ets[1][1][1][1]-=p.ets[1][1][1][1];

	return *this;
}

Tensor4 operator+(const Tensor4 & p)
{
	return p;
}

Tensor4 operator-(const Tensor4& p)
{
	return Tensor4() - p; //a.r: changed to default constructor to compatible 2D/3D versions
}

Tensor4 operator+(const Tensor4& p1, const Tensor4 & p2)
{
	Tensor4 sum = p1;
	sum += p2;
	return sum;
}

Tensor4 operator-(const Tensor4& p1, const Tensor4& p2)
{
	Tensor4 sum = p1;
	sum -= p2;
	return sum;
}

Tensor4 operator*(const double scalar, const Tensor4 & p)
{
	Tensor4 tm;
	tm.ets[0][0][0][0]=scalar*p.ets[0][0][0][0];
	tm.ets[0][0][0][1]=scalar*p.ets[0][0][0][1];
	tm.ets[0][0][1][0]=scalar*p.ets[0][0][1][0];
	tm.ets[0][0][1][1]=scalar*p.ets[0][0][1][1];
	tm.ets[0][1][0][0]=scalar*p.ets[0][1][0][0];
	tm.ets[0][1][0][1]=scalar*p.ets[0][1][0][1];
	tm.ets[0][1][1][0]=scalar*p.ets[0][1][1][0];
	tm.ets[0][1][1][1]=scalar*p.ets[0][1][1][1];
	tm.ets[1][0][0][0]=scalar*p.ets[1][0][0][0];
	tm.ets[1][0][0][1]=scalar*p.ets[1][0][0][1];
	tm.ets[1][0][1][0]=scalar*p.ets[1][0][1][0];
	tm.ets[1][0][1][1]=scalar*p.ets[1][0][1][1];
	tm.ets[1][1][0][0]=scalar*p.ets[1][1][0][0];
	tm.ets[1][1][0][1]=scalar*p.ets[1][1][0][1];
	tm.ets[1][1][1][0]=scalar*p.ets[1][1][1][0];
	tm.ets[1][1][1][1]=scalar*p.ets[1][1][1][1];

	return tm;
}

inline Tensor4 operator*(const Tensor4 & p, const double scalar)
{
	return scalar*p;
}

Tensor4 operator*(const Tensor3 & p1, const Tensor3 & p2)
{
	Tensor4 tm;
	tm.ets[0][0][0][0]= p1.ets[0][0][0]*p2.ets[0][0][0] + p1.ets[0][0][1]*p2.ets[1][0][0];
	tm.ets[0][0][0][1]= p1.ets[0][0][0]*p2.ets[0][0][1] + p1.ets[0][0][1]*p2.ets[1][0][1];
	tm.ets[0][0][1][0]= p1.ets[0][0][0]*p2.ets[0][1][0] + p1.ets[0][0][1]*p2.ets[1][1][0];
	tm.ets[0][0][1][1]= p1.ets[0][0][0]*p2.ets[0][1][1] + p1.ets[0][0][1]*p2.ets[1][1][1];
	tm.ets[0][1][0][0]= p1.ets[0][1][0]*p2.ets[0][0][0] + p1.ets[0][1][1]*p2.ets[1][0][0];
	tm.ets[0][1][0][1]= p1.ets[0][1][0]*p2.ets[0][0][1] + p1.ets[0][1][1]*p2.ets[1][0][1];
	tm.ets[0][1][1][0]= p1.ets[0][1][0]*p2.ets[0][1][0] + p1.ets[0][1][1]*p2.ets[1][1][0];
	tm.ets[0][1][1][1]= p1.ets[0][1][0]*p2.ets[0][1][1] + p1.ets[0][1][1]*p2.ets[1][1][1];
	tm.ets[1][0][0][0]= p1.ets[1][0][0]*p2.ets[0][0][0] + p1.ets[1][0][1]*p2.ets[1][0][0];
	tm.ets[1][0][0][1]= p1.ets[1][0][0]*p2.ets[0][0][1] + p1.ets[1][0][1]*p2.ets[1][0][1];
	tm.ets[1][0][1][0]= p1.ets[1][0][0]*p2.ets[0][1][0] + p1.ets[1][0][1]*p2.ets[1][1][0];
	tm.ets[1][0][1][1]= p1.ets[1][0][0]*p2.ets[0][1][1] + p1.ets[1][0][1]*p2.ets[1][1][1];
	tm.ets[1][1][0][0]= p1.ets[1][1][0]*p2.ets[0][0][0] + p1.ets[1][1][1]*p2.ets[1][0][0];
	tm.ets[1][1][0][1]= p1.ets[1][1][0]*p2.ets[0][0][1] + p1.ets[1][1][1]*p2.ets[1][0][1];
	tm.ets[1][1][1][0]= p1.ets[1][1][0]*p2.ets[0][1][0] + p1.ets[1][1][1]*p2.ets[1][1][0];
	tm.ets[1][1][1][1]= p1.ets[1][1][0]*p2.ets[0][1][1] + p1.ets[1][1][1]*p2.ets[1][1][1];

	return tm;
}


