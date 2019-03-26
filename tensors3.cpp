//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// changed		85/12 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************

#include "tensors3.h"

Tensor3 & Tensor3::operator+=(const Tensor3& p)
{
	ets[0][0][0]+=p.ets[0][0][0];
	ets[0][0][1]+=p.ets[0][0][1];
	ets[0][1][0]+=p.ets[0][1][0];
	ets[0][1][1]+=p.ets[0][1][1];
	ets[1][0][0]+=p.ets[1][0][0];
	ets[1][0][1]+=p.ets[1][0][1];
	ets[1][1][0]+=p.ets[1][1][0];
	ets[1][1][1]+=p.ets[1][1][1];

	return *this;
}

Tensor3 & Tensor3::operator-=(const Tensor3& p)
{
	ets[0][0][0]-=p.ets[0][0][0];
	ets[0][0][1]-=p.ets[0][0][1];
	ets[0][1][0]-=p.ets[0][1][0];
	ets[0][1][1]-=p.ets[0][1][1];
	ets[1][0][0]-=p.ets[1][0][0];
	ets[1][0][1]-=p.ets[1][0][1];
	ets[1][1][0]-=p.ets[1][1][0];
	ets[1][1][1]-=p.ets[1][1][1];

	return *this;
}

Tensor3 operator+(const Tensor3 & p)
{
	return p;
}

Tensor3 operator-(const Tensor3& p)
{
	return Tensor3() - p; //a.r: changed to default constructor to compatible 2D/3D versions
}

Tensor3 operator+(const Tensor3& p1, const Tensor3 & p2)
{
	Tensor3 sum = p1;
	sum += p2;
	return sum;
}

Tensor3 operator-(const Tensor3& p1, const Tensor3& p2)
{
	Tensor3 sum = p1;
	sum -= p2;
	return sum;
}

Tensor3 operator*(const double scalar, const Tensor3 & p)
{
	Tensor3 tm;
	tm.ets[0][0][0]= scalar*p.ets[0][0][0];
	tm.ets[0][0][1]= scalar*p.ets[0][0][1];
	tm.ets[0][1][0]= scalar*p.ets[0][1][0];
	tm.ets[0][1][1]= scalar*p.ets[0][1][1];
	tm.ets[1][0][0]= scalar*p.ets[1][0][0];
	tm.ets[1][0][1]= scalar*p.ets[1][0][1];
	tm.ets[1][1][0]= scalar*p.ets[1][1][0];
	tm.ets[1][1][1]= scalar*p.ets[1][1][1];

	return tm;
}

inline Tensor3 operator*(const Tensor3 & p, const double scalar)
{
	return scalar*p;
}

Tensor3 operator*(const Tensor3 & p1, const Tensor & p2)
{
	Tensor3 tm;

	tm.ets[0][0][0]= p1.ets[0][0][0]*p2.GetXX()+p1.ets[0][0][1]*p2.GetYX();
	tm.ets[0][0][1]= p1.ets[0][0][0]*p2.GetXY()+p1.ets[0][0][1]*p2.GetYY();
	tm.ets[0][1][0]= p1.ets[0][1][0]*p2.GetXX()+p1.ets[0][1][1]*p2.GetYX();
	tm.ets[0][1][1]= p1.ets[0][1][0]*p2.GetXY()+p1.ets[0][1][1]*p2.GetYY();
	tm.ets[1][0][0]= p1.ets[1][0][0]*p2.GetXX()+p1.ets[1][0][1]*p2.GetYX();
	tm.ets[1][0][1]= p1.ets[1][0][0]*p2.GetXY()+p1.ets[1][0][1]*p2.GetYY();
	tm.ets[1][1][0]= p1.ets[1][1][0]*p2.GetXX()+p1.ets[1][1][1]*p2.GetYX();
	tm.ets[1][1][1]= p1.ets[1][1][0]*p2.GetXY()+p1.ets[1][1][1]*p2.GetYY();

	return tm;
}

Tensor3 operator*(const Tensor3 & p1, const STensor & p2)
{
	Tensor3 tm;

	tm.ets[0][0][0]= p1.ets[0][0][0]*p2.GetXX()+p1.ets[0][0][1]*p2.GetXY();
	tm.ets[0][0][1]= p1.ets[0][0][0]*p2.GetXY()+p1.ets[0][0][1]*p2.GetYY();
	tm.ets[0][1][0]= p1.ets[0][1][0]*p2.GetXX()+p1.ets[0][1][1]*p2.GetXY();
	tm.ets[0][1][1]= p1.ets[0][1][0]*p2.GetXY()+p1.ets[0][1][1]*p2.GetYY();
	tm.ets[1][0][0]= p1.ets[1][0][0]*p2.GetXX()+p1.ets[1][0][1]*p2.GetXY();
	tm.ets[1][0][1]= p1.ets[1][0][0]*p2.GetXY()+p1.ets[1][0][1]*p2.GetYY();
	tm.ets[1][1][0]= p1.ets[1][1][0]*p2.GetXX()+p1.ets[1][1][1]*p2.GetXY();
	tm.ets[1][1][1]= p1.ets[1][1][0]*p2.GetXY()+p1.ets[1][1][1]*p2.GetYY();


	return tm;
}
Tensor3 operator*(const STensor & p2, const Tensor3 & p1)
{
	Tensor3 tm;
	tm.ets[0][0][0]= p1.ets[0][0][0]*p2.GetXX()+p1.ets[1][0][0]*p2.GetXY();
	tm.ets[0][0][1]= p1.ets[0][0][1]*p2.GetXX()+p1.ets[1][0][1]*p2.GetXY();
	tm.ets[0][1][0]= p1.ets[0][1][0]*p2.GetXX()+p1.ets[1][1][0]*p2.GetXY();
	tm.ets[0][1][1]= p1.ets[0][1][1]*p2.GetXX()+p1.ets[1][1][1]*p2.GetXY();
	tm.ets[1][0][0]= p1.ets[0][0][0]*p2.GetXY()+p1.ets[1][0][0]*p2.GetYY();
	tm.ets[1][0][1]= p1.ets[0][0][1]*p2.GetXY()+p1.ets[1][0][1]*p2.GetYY();
	tm.ets[1][1][0]= p1.ets[0][1][0]*p2.GetXY()+p1.ets[1][1][0]*p2.GetYY();
	tm.ets[1][1][1]= p1.ets[0][1][1]*p2.GetXY()+p1.ets[1][1][1]*p2.GetYY();


	return tm;
}
Tensor3 operator*(const Tensor & p2, const Tensor3 & p1)
{
	Tensor3 tm;
	tm.ets[0][0][0]= p1.ets[0][0][0]*p2.GetXX()+p1.ets[1][0][0]*p2.GetXY();
	tm.ets[0][0][1]= p1.ets[0][0][1]*p2.GetXX()+p1.ets[1][0][1]*p2.GetXY();
	tm.ets[0][1][0]= p1.ets[0][1][0]*p2.GetXX()+p1.ets[1][1][0]*p2.GetXY();
	tm.ets[0][1][1]= p1.ets[0][1][1]*p2.GetXX()+p1.ets[1][1][1]*p2.GetXY();
	tm.ets[1][0][0]= p1.ets[0][0][0]*p2.GetYX()+p1.ets[1][0][0]*p2.GetYY();
	tm.ets[1][0][1]= p1.ets[0][0][1]*p2.GetYX()+p1.ets[1][0][1]*p2.GetYY();
	tm.ets[1][1][0]= p1.ets[0][1][0]*p2.GetYX()+p1.ets[1][1][0]*p2.GetYY();
	tm.ets[1][1][1]= p1.ets[0][1][1]*p2.GetYX()+p1.ets[1][1][1]*p2.GetYY();


	return tm;
}

Tensor operator*(const Tensor3 &p, const Vector &v)
{
	Tensor tm;

	tm.SetXX( p.ets[0][0][0]*v.GetX()+p.ets[0][0][1]*v.GetY() );
	tm.SetXY( p.ets[0][1][0]*v.GetX()+p.ets[0][1][1]*v.GetY() );
	tm.SetYX( p.ets[1][0][0]*v.GetX()+p.ets[1][0][1]*v.GetY() );
	tm.SetYY( p.ets[1][1][0]*v.GetX()+p.ets[1][1][1]*v.GetY() );

	return tm;
}

Tensor operator*(const Vector &v, const Tensor3 &p)
{
	Tensor tm;

	tm.SetXX( p.ets[0][0][0]*v.GetX()+p.ets[1][0][0]*v.GetY() );
	tm.SetXY( p.ets[0][1][0]*v.GetX()+p.ets[1][1][0]*v.GetY() );
	tm.SetYX( p.ets[0][0][1]*v.GetX()+p.ets[1][0][1]*v.GetY() );
	tm.SetYY( p.ets[0][1][1]*v.GetX()+p.ets[1][1][1]*v.GetY() );

	return tm;
}

Tensor3 OProduct(const Tensor &p1, const Vector &v)
{
	Tensor3 tm;

	tm.ets[0][0][0]= p1.GetXX()*v.GetX();
	tm.ets[0][0][1]= p1.GetXX()*v.GetY();
	tm.ets[0][1][0]= p1.GetXY()*v.GetX();
	tm.ets[0][1][1]= p1.GetXY()*v.GetY();
	tm.ets[1][0][0]= p1.GetYX()*v.GetX();
	tm.ets[1][0][1]= p1.GetYX()*v.GetY();
	tm.ets[1][1][0]= p1.GetYY()*v.GetX();
	tm.ets[1][1][1]= p1.GetYY()*v.GetY();


	return tm;
}
Tensor3 OProduct(const Vector &v, const Tensor &p1)
{
	Tensor3 tm;

	tm.ets[0][0][0]= p1.GetXX()*v.GetX();
	tm.ets[0][0][1]= p1.GetXY()*v.GetX();
	tm.ets[0][1][0]= p1.GetYX()*v.GetX();
	tm.ets[0][1][1]= p1.GetYY()*v.GetX();
	tm.ets[1][0][0]= p1.GetXX()*v.GetY();
	tm.ets[1][0][1]= p1.GetXY()*v.GetY();
	tm.ets[1][1][0]= p1.GetYX()*v.GetY();
	tm.ets[1][1][1]= p1.GetYY()*v.GetY();

	return tm;
}
	
Vector DiadProduct(const Tensor &p1, const Tensor3 &p2)
{
	Vector tm;


	tm.SetX( p1.GetXX()*p2.ets[0][0][0]+p1.GetXY()*p2.ets[0][1][0]+p1.GetYX()*p2.ets[1][0][0]+p1.GetYY()*p2.ets[1][1][0] );
	tm.SetY( p1.GetXX()*p2.ets[0][0][1]+p1.GetXY()*p2.ets[0][1][1]+p1.GetYX()*p2.ets[1][0][1]+p1.GetYY()*p2.ets[1][1][1] );

	return tm;

}
Tensor DiadProduct(const Tensor3 &p1, const Tensor3 &p2)
{
	Tensor tm;


	tm.SetXX( p1.ets[0][0][0]*p2.ets[0][0][0]+p1.ets[0][1][0]*p2.ets[0][1][0]+p1.ets[0][0][1]*p2.ets[1][0][0]+p1.ets[0][1][1]*p2.ets[1][1][0] );
	tm.SetXY( p1.ets[0][0][0]*p2.ets[0][0][1]+p1.ets[0][1][0]*p2.ets[0][1][1]+p1.ets[0][0][1]*p2.ets[1][0][1]+p1.ets[0][1][1]*p2.ets[1][1][1] );
	tm.SetYX( p1.ets[1][0][0]*p2.ets[0][0][0]+p1.ets[1][1][0]*p2.ets[0][1][0]+p1.ets[1][0][1]*p2.ets[1][0][0]+p1.ets[1][1][1]*p2.ets[1][1][0] );
	tm.SetYY( p1.ets[1][0][0]*p2.ets[0][0][1]+p1.ets[1][1][0]*p2.ets[0][1][1]+p1.ets[1][0][1]*p2.ets[1][0][1]+p1.ets[1][1][1]*p2.ets[1][1][1] );


	return tm;

}


