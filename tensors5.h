//****************************************************************************************
//
// SPH Code
//
// created	89/03 By R. Fatehi
// changed		89/03 By Fatehi
////// Last change	97/11 By Nikooei
//****************************************************************************************

#ifndef _TENSOR5_H_
#define _TENSOR5_H_

#include "vectors.h"
#include "stensors.h"
#include "tensors.h"
#include "tensors3.h"
#include "tensors4.h"

class Tensor5{
	public:
		double ets[DIM][DIM][DIM][DIM][DIM];
		Tensor5()
		{

			ets[0][0][0][0][0]=ets[0][0][0][0][1]=ets[0][0][0][1][0]=ets[0][0][0][1][1]=ets[0][0][1][0][0]=ets[0][0][1][0][1]=ets[0][0][1][1][0]=ets[0][0][1][1][1]=
			ets[0][1][0][0][0]=ets[0][1][0][0][1]=ets[0][1][0][1][0]=ets[0][1][0][1][1]=ets[0][1][1][0][0]=ets[0][1][1][0][1]=ets[0][1][1][1][0]=ets[0][1][1][1][1]=			
			ets[1][0][0][0][0]=ets[1][0][0][0][1]=ets[1][0][0][1][0]=ets[1][0][0][1][1]=ets[1][0][1][0][0]=ets[1][0][1][0][1]=ets[1][0][1][1][0]=ets[1][0][1][1][1]=
			ets[1][1][0][0][0]=ets[1][1][0][0][1]=ets[1][1][0][1][0]=ets[1][1][0][1][1]=ets[1][1][1][0][0]=ets[1][1][1][0][1]=ets[1][1][1][1][0]=ets[1][1][1][1][1]=0.0;

		}

		Tensor5(Vector v1, Vector v2, Vector v3, Vector v4, Vector v5)
		{
			for(short i=0; i<DIM; i++)
			for(short j=0; j<DIM; j++)
			for(short k=0; k<DIM; k++)
			for(short l=0; l<DIM; l++)
			for(short m=0; m<DIM; m++)
			ets[i][j][k][l][m]=v1[i]*v2[j]*v3[k]*v4[l]*v5[m];
		}
		Tensor5(Vector v1)
		{
			for(short i=0; i<DIM; i++)
			for(short j=0; j<DIM; j++)
			for(short k=0; k<DIM; k++)
			for(short l=0; l<DIM; l++)
			for(short m=0; m<DIM; m++)
			ets[i][j][k][l][m]=v1[i]*v1[j]*v1[k]*v1[l]*v1[m];
		}
		
		Tensor5(const Tensor5 &p)
		{
		ets[0][0][0][0][0]=p.ets[0][0][0][0][0];
		ets[0][0][0][0][1]=p.ets[0][0][0][0][1];
		ets[0][0][0][1][0]=p.ets[0][0][0][1][0];
		ets[0][0][0][1][1]=p.ets[0][0][0][1][1];
		ets[0][0][1][0][0]=p.ets[0][0][1][0][0];
		ets[0][0][1][0][1]=p.ets[0][0][1][0][1];
		ets[0][0][1][1][0]=p.ets[0][0][1][1][0];
		ets[0][0][1][1][1]=p.ets[0][0][1][1][1];
		ets[0][1][0][0][0]=p.ets[0][1][0][0][0];
		ets[0][1][0][0][1]=p.ets[0][1][0][0][1];
		ets[0][1][0][1][0]=p.ets[0][1][0][1][0];
		ets[0][1][0][1][1]=p.ets[0][1][0][1][1];
		ets[0][1][1][0][0]=p.ets[0][1][1][0][0];
		ets[0][1][1][0][1]=p.ets[0][1][1][0][1];
		ets[0][1][1][1][0]=p.ets[0][1][1][1][0];
		ets[0][1][1][1][1]=p.ets[0][1][1][1][1];
		ets[1][0][0][0][0]=p.ets[1][0][0][0][0];
		ets[1][0][0][0][1]=p.ets[1][0][0][0][1];
		ets[1][0][0][1][0]=p.ets[1][0][0][1][0];
		ets[1][0][0][1][1]=p.ets[1][0][0][1][1];
		ets[1][0][1][0][0]=p.ets[1][0][1][0][0];
		ets[1][0][1][0][1]=p.ets[1][0][1][0][1];
		ets[1][0][1][1][0]=p.ets[1][0][1][1][0];
		ets[1][0][1][1][1]=p.ets[1][0][1][1][1];
		ets[1][1][0][0][0]=p.ets[1][1][0][0][0];
		ets[1][1][0][0][1]=p.ets[1][1][0][0][1];
		ets[1][1][0][1][0]=p.ets[1][1][0][1][0];
		ets[1][1][0][1][1]=p.ets[1][1][0][1][1];
		ets[1][1][1][0][0]=p.ets[1][1][1][0][0];
		ets[1][1][1][0][1]=p.ets[1][1][1][0][1];
		ets[1][1][1][1][0]=p.ets[1][1][1][1][0];
		ets[1][1][1][1][1]=p.ets[1][1][1][1][1];
		}


	Tensor5& operator+=(const Tensor5 &);
	Tensor5& operator-=(const Tensor5 &);

	friend Tensor5 operator+(const Tensor5&);
	friend Tensor5 operator-(const Tensor5&);
	friend Tensor5 operator+(const Tensor5&, const Tensor5&);
	friend Tensor5 operator-(const Tensor5&, const Tensor5&);
	friend Tensor5 operator*(const double, const Tensor5&);
	friend Tensor5 operator*(const Tensor5&, const double);
	friend Tensor5 operator*(const Tensor5&, const Tensor&);
	friend Tensor5 operator*(const Tensor4&, const Tensor3&);

	//friend std::ostream& operator<<(std::ostream&, const Tensor5&);
};

typedef std::vector<Tensor5> Tensor5vec;


#endif

