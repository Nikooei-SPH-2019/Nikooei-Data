//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change		85/12 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************

#include "tensors.h"

STensor Tensor::GetSymPart() const
{
	STensor tm;

	tm.SetXX(ets[0]);
	tm.SetXY(0.5*(ets[1]+ets[2]));
	tm.SetYY(ets[3]);

	return tm;
}

Tensor Tensor::Inverse() const
{

	const double det = ets[0]*ets[3]-ets[1]*ets[2];
	if (det)
	{
		Tensor tm;
		tm.SetXX(ets[3]);
		tm.SetXY(-ets[1]);
		tm.SetYX(-ets[2]);
		tm.SetYY(ets[0]);
		return 1.0/det*tm;
	}
	else
		return Tensor(1/ets[0], 0.0, 0.0, 0.0);		//One dimensional

}


Tensor& Tensor::operator=(const Tensor& p)
{

	ets[0]=p.ets[0];
	ets[1]=p.ets[1];
	ets[2]=p.ets[2];
	ets[3]=p.ets[3];

	return *this;
}

Tensor & Tensor::operator+=(const Tensor& p)
{
	ets[0] +=p.ets[0];
	ets[1] +=p.ets[1];
	ets[2] +=p.ets[2];
	ets[3] +=p.ets[3];

	return *this;
}

Tensor & Tensor::operator-=(const Tensor& p)
{
	ets[0] -=p.ets[0];
	ets[1] -=p.ets[1];
	ets[2] -=p.ets[2];
	ets[3] -=p.ets[3];

	return *this;
}

Tensor operator+(const Tensor & p)
{
	return p;
}

Tensor operator-(const Tensor& p)
{
	return Tensor() - p; //a.r: changed to default constructor to compatible 2D/3D versions
}

Tensor operator+(const Tensor& p1, const Tensor & p2)
{
	Tensor sum = p1;
	sum += p2;
	return sum;
}

Tensor operator-(const Tensor& p1, const Tensor& p2)
{
	Tensor sum = p1;
	sum -= p2;
	return sum;
}

Tensor operator*(const double scalar, const Tensor & p)
{
	Tensor tm;
	tm[0] = scalar*p.ets[0];
	tm[1] = scalar*p.ets[1];
	tm[2] = scalar*p.ets[2];

	return tm;
}

inline Tensor operator*(const Tensor & p, const double scalar)
{
	return scalar*p;
}

Tensor operator*(const Tensor & p1, const Tensor & p2)
{
	Tensor tm;
	tm[0]=p1.ets[0]*p2.ets[0]+p1.ets[1]*p2.ets[2];
	tm[1]=p1.ets[0]*p2.ets[1]+p1.ets[1]*p2.ets[3];
	tm[2]=p1.ets[2]*p2.ets[0]+p1.ets[3]*p2.ets[2];
	tm[3]=p1.ets[2]*p2.ets[1]+p1.ets[3]*p2.ets[3];

	return tm;
}

Tensor operator*(const STensor & p1, const STensor & p2)
{
	Tensor tm;

	tm[0]=p1.GetXX()*p2.GetXX()+p1.GetXY()*p2.GetXY();
	tm[1]=p1.GetXX()*p2.GetXY()+p1.GetXY()*p2.GetYY();
	tm[2]=p1.GetXY()*p2.GetXX()+p1.GetYY()*p2.GetXY();
	tm[3]=p1.GetXY()*p2.GetXY()+p1.GetYY()*p2.GetYY();


	return tm;
}
Tensor operator*(const STensor & p1, const Tensor & p2)
{
	Tensor tm;

	tm[0]=p1.GetXX()*p2.GetXX()+p1.GetXY()*p2.GetYX();
	tm[1]=p1.GetXX()*p2.GetXY()+p1.GetXY()*p2.GetYY();
	tm[2]=p1.GetXY()*p2.GetXX()+p1.GetYY()*p2.GetYX();
	tm[3]=p1.GetXY()*p2.GetXY()+p1.GetYY()*p2.GetYY();

	return tm;
}

Vector operator*(const Tensor &p, const Vector &v)
{
	Vector tm;

	tm[0]=p.ets[0]*v.GetX()+p.ets[1]*v.GetY();
	tm[1]=p.ets[2]*v.GetX()+p.ets[3]*v.GetY();

	return tm;
}

Vector operator*(const Vector &v, const Tensor &p )
{
	Vector tm;
	tm[0]=p.ets[0]*v.GetX()+p.ets[2]*v.GetY();
	tm[1]=p.ets[1]*v.GetX()+p.ets[3]*v.GetY();

	return tm;
}

Tensor OProduct(const Vector &v1, const Vector &v2)
{
	Tensor tm;

	tm[0]=v1.GetX()*v2.GetX();
	tm[1]=v1.GetX()*v2.GetY();
	tm[2]=v1.GetY()*v2.GetX();
	tm[3]=v1.GetY()*v2.GetY();

	return tm;
}



STensor OProduct_same_vec(const Vector &v1, const Vector &v2)
{
	STensor tm;

	tm.SetXX(v1.GetX()*v2.GetX());
	tm.SetXY(v1.GetX()*v2.GetY());
	tm.SetYY(v1.GetY()*v2.GetY());

	return tm;
}


Tensor ToTensor(STensor st)
{
	Tensor t;

	t[0]=st.GetXX();
	t[1]=st.GetXY();
	t[2]=st.GetXY();
	t[3]=st.GetYY();

	return t;
}

std::ostream&  operator<<(std::ostream& s, const Tensor& p )
{
	s << p.ets[0] << '\t' << p.ets[1] << '\t' << p.ets[2] << '\t' << p.ets[3]

	  ;
	return s;
}

void Tensor::Transpose()
{

	const double xy=ets[1];
	ets[1]=ets[2];
	ets[2]=xy;

}
double DiadProduct(const Tensor &T, const Tensor &E)
{

	return T.ets[0] * E.ets[0] + T.ets[1] * E.ets[1] + T.ets[2] * E.ets[2] + T.ets[3] * E.ets[3];

}

