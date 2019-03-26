//****************************************************************************************
//
// SPH Code
//
// created	88/05 By  Fatehi
// changed		89/03 By Fatehi
////// Last change	97/11 By Nikooei
//****************************************************************************************

#ifndef _TENSOR3_H_
#define _TENSOR3_H_

#include "vectors.h"
#include "stensors.h"
#include "tensors.h"

class Tensor3{
	public:
		double ets[DIM][DIM][DIM];

		Tensor3()
		{
			ets[0][0][0]=ets[0][0][1]=ets[0][1][0]=ets[0][1][1]=ets[1][0][0]=ets[1][0][1]=ets[1][1][0]=ets[1][1][1]=0.0;
		}
		Tensor3(Vector v1, Vector v2, Vector v3)
		{
			ets[0][0][0]=v1[0]*v2[0]*v3[0];
			ets[0][0][1]=v1[0]*v2[0]*v3[1];
			ets[0][1][0]=v1[0]*v2[1]*v3[0];
			ets[0][1][1]=v1[0]*v2[1]*v3[1];
			ets[1][0][0]=v1[1]*v2[0]*v3[0];
			ets[1][0][1]=v1[1]*v2[0]*v3[1];
			ets[1][1][0]=v1[1]*v2[1]*v3[0];
			ets[1][1][1]=v1[1]*v2[1]*v3[1];
		}
		Tensor3(Vector v1)
		{
			ets[0][0][0]=v1[0]*v1[0]*v1[0];
			ets[0][0][1]=v1[0]*v1[0]*v1[1];
			ets[0][1][0]=v1[0]*v1[1]*v1[0];
			ets[0][1][1]=v1[0]*v1[1]*v1[1];
			ets[1][0][0]=v1[1]*v1[0]*v1[0];
			ets[1][0][1]=v1[1]*v1[0]*v1[1];
			ets[1][1][0]=v1[1]*v1[1]*v1[0];
			ets[1][1][1]=v1[1]*v1[1]*v1[1];
		}
		Tensor3(const Tensor3 &p)
		{
			ets[0][0][0]=p.ets[0][0][0];
			ets[0][0][1]=p.ets[0][0][1];
			ets[0][1][0]=p.ets[0][1][0];
			ets[0][1][1]=p.ets[0][1][1];
			ets[1][0][0]=p.ets[1][0][0];
			ets[1][0][1]=p.ets[1][0][1];
			ets[1][1][0]=p.ets[1][1][0];
			ets[1][1][1]=p.ets[1][1][1];
		}

	double GetXXX() const {return ets[0][0][0];}
	double GetXYX() const {return ets[0][1][0];}
	double GetYXX() const {return ets[1][0][0];}
	double GetYYX() const {return ets[1][1][0];}
	double GetXXY() const {return ets[0][0][1];}
	double GetXYY() const {return ets[0][1][1];}
	double GetYXY() const {return ets[1][0][1];}
	double GetYYY() const {return ets[1][1][1];}

	Tensor3& operator+=(const Tensor3 &);
	Tensor3& operator-=(const Tensor3 &);

	friend Tensor3 operator+(const Tensor3&);
	friend Tensor3 operator-(const Tensor3&);
	friend Tensor3 operator+(const Tensor3&, const Tensor3&);
	friend Tensor3 operator-(const Tensor3&, const Tensor3&);
	friend Tensor3 operator*(const double, const Tensor3&);
	friend Tensor3 operator*(const Tensor3&, const double);
	friend Tensor3 operator*(const Tensor3&, const Tensor&);
	friend Tensor3 operator*(const Tensor3&, const STensor&);
	friend Tensor3 operator*(const STensor&, const Tensor3&);
	friend Tensor3 operator*(const Tensor&, const Tensor3&);
	friend Tensor operator*(const Tensor3 &, const Vector &);
	friend Tensor operator*(const Vector &, const Tensor3 &);
	friend Tensor3 OProduct(const Tensor &, const Vector &);
	friend Tensor3 OProduct(const Vector &, const Tensor &);
	friend Vector DiadProduct(const Tensor &, const Tensor3 &);
	friend Tensor DiadProduct(const Tensor3 &, const Tensor3 &);
	//friend Tensor DiadProduct(const Tensor3 &, const Tensor3 &);

	double& operator()(const int i,const int j,const int k) {return ets[i][j][k];}

	//friend std::ostream& operator<<(std::ostream&, const Tensor3&);
};

typedef std::vector<Tensor3> Tensor3vec;


#endif

