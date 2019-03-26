//****************************************************************************************
//
// SPH Code
//
// created	89/03 By R. Fatehi
// changed		89/03 By Fatehi
////// Last change	97/11 By Nikooei
//****************************************************************************************

#ifndef _TENSOR4_H_
#define _TENSOR4_H_

#include "vectors.h"
#include "stensors.h"
#include "tensors.h"
#include "tensors3.h"

class Tensor4{
	public:
		double ets[DIM][DIM][DIM][DIM];

		Tensor4()
		{

			ets[0][0][0][0]=ets[0][0][0][1]=ets[0][0][1][0]=ets[0][0][1][1]=ets[0][1][0][0]=ets[0][1][0][1]=ets[0][1][1][0]=ets[0][1][1][1]=
			ets[1][0][0][0]=ets[1][0][0][1]=ets[1][0][1][0]=ets[1][0][1][1]=ets[1][1][0][0]=ets[1][1][0][1]=ets[1][1][1][0]=ets[1][1][1][1]=0.0;

		}

		Tensor4(Vector v1, Vector v2, Vector v3, Vector v4)
		{
			ets[0][0][0][0]=v1[0]*v2[0]*v3[0]*v4[0];
			ets[0][0][1][0]=v1[0]*v2[0]*v3[1]*v4[0];
			ets[0][1][0][0]=v1[0]*v2[1]*v3[0]*v4[0];
			ets[0][1][1][0]=v1[0]*v2[1]*v3[1]*v4[0];
			ets[1][0][0][0]=v1[1]*v2[0]*v3[0]*v4[0];
			ets[1][0][1][0]=v1[1]*v2[0]*v3[1]*v4[0];
			ets[1][1][0][0]=v1[1]*v2[1]*v3[0]*v4[0];
			ets[1][1][1][0]=v1[1]*v2[1]*v3[1]*v4[0];
			ets[0][0][0][1]=v1[0]*v2[0]*v3[0]*v4[1];
			ets[0][0][1][1]=v1[0]*v2[0]*v3[1]*v4[1];
			ets[0][1][0][1]=v1[0]*v2[1]*v3[0]*v4[1];
			ets[0][1][1][1]=v1[0]*v2[1]*v3[1]*v4[1];
			ets[1][0][0][1]=v1[1]*v2[0]*v3[0]*v4[1];
			ets[1][0][1][1]=v1[1]*v2[0]*v3[1]*v4[1];
			ets[1][1][0][1]=v1[1]*v2[1]*v3[0]*v4[1];
			ets[1][1][1][1]=v1[1]*v2[1]*v3[1]*v4[1];
		}
		Tensor4(Vector v1)
		{
			ets[0][0][0][0]=v1[0]*v1[0]*v1[0]*v1[0];
			ets[0][0][1][0]=v1[0]*v1[0]*v1[1]*v1[0];
			ets[0][1][0][0]=v1[0]*v1[1]*v1[0]*v1[0];
			ets[0][1][1][0]=v1[0]*v1[1]*v1[1]*v1[0];
			ets[1][0][0][0]=v1[1]*v1[0]*v1[0]*v1[0];
			ets[1][0][1][0]=v1[1]*v1[0]*v1[1]*v1[0];
			ets[1][1][0][0]=v1[1]*v1[1]*v1[0]*v1[0];
			ets[1][1][1][0]=v1[1]*v1[1]*v1[1]*v1[0];
			ets[0][0][0][1]=v1[0]*v1[0]*v1[0]*v1[1];
			ets[0][0][1][1]=v1[0]*v1[0]*v1[1]*v1[1];
			ets[0][1][0][1]=v1[0]*v1[1]*v1[0]*v1[1];
			ets[0][1][1][1]=v1[0]*v1[1]*v1[1]*v1[1];
			ets[1][0][0][1]=v1[1]*v1[0]*v1[0]*v1[1];
			ets[1][0][1][1]=v1[1]*v1[0]*v1[1]*v1[1];
			ets[1][1][0][1]=v1[1]*v1[1]*v1[0]*v1[1];
			ets[1][1][1][1]=v1[1]*v1[1]*v1[1]*v1[1];
		}
		
		Tensor4(const Tensor4 &p)
		{
		ets[0][0][0][0]=p.ets[0][0][0][0];
		ets[0][0][0][1]=p.ets[0][0][0][1];
		ets[0][0][1][0]=p.ets[0][0][1][0];
		ets[0][0][1][1]=p.ets[0][0][1][1];
		ets[0][1][0][0]=p.ets[0][1][0][0];
		ets[0][1][0][1]=p.ets[0][1][0][1];
		ets[0][1][1][0]=p.ets[0][1][1][0];
		ets[0][1][1][1]=p.ets[0][1][1][1];
		ets[1][0][0][0]=p.ets[1][0][0][0];
		ets[1][0][0][1]=p.ets[1][0][0][1];
		ets[1][0][1][0]=p.ets[1][0][1][0];
		ets[1][0][1][1]=p.ets[1][0][1][1];
		ets[1][1][0][0]=p.ets[1][1][0][0];
		ets[1][1][0][1]=p.ets[1][1][0][1];
		ets[1][1][1][0]=p.ets[1][1][1][0];
		ets[1][1][1][1]=p.ets[1][1][1][1];
		}


	Tensor4& operator+=(const Tensor4 &);
	Tensor4& operator-=(const Tensor4 &);

	friend Tensor4 operator+(const Tensor4&);
	friend Tensor4 operator-(const Tensor4&);
	friend Tensor4 operator+(const Tensor4&, const Tensor4&);
	friend Tensor4 operator-(const Tensor4&, const Tensor4&);
	friend Tensor4 operator*(const double, const Tensor4&);
	friend Tensor4 operator*(const Tensor4&, const double);
	friend Tensor4 operator*(const Tensor3&, const Tensor3&);

};

typedef std::vector<Tensor4> Tensor4vec;


#endif

