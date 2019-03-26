//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change	85/12 By Fatehi
// Last change	87/07 By Fatehi
// changed		89/03 By Fatehi
//Last change	97/11 By Nikooei
//****************************************************************************************


#ifndef _VECTORS_H_
#define _VECTORS_H_

#include "general.h"

class Vector{
	private:
		double ets[DIM];	
	public:

		Vector(){ets[0]=ets[1]=0.0;}
		Vector(double n){ets[0]=n;ets[1]=n;}
		Vector(double nx, double ny){ets[0]=nx;ets[1]=ny;}
		Vector(const Vector &p){ets[0]=p.ets[0];ets[1]=p.ets[1];}


	double GetX() const {return ets[0];}
	double GetY() const {return ets[1];}

	void SetX(double nx){ets[0]=nx;}
	void SetY(double ny){ets[1]=ny;}
	void SetXY(double nx, double ny){ets[0]=nx;ets[1]=ny;}

	Vector& operator+=(const Vector &);
	Vector& operator-=(const Vector &);
	double MaxNorm() const;
	double TwoNorm() const;

	double& operator[](int i) {return ets[i];}
	Vector& operator=(const Vector&);

	friend bool operator<=(const Vector&, const Vector&);
	friend bool operator<(const Vector&, const Vector&);
	friend bool operator>=(const Vector&, const Vector&);
	friend bool operator>(const Vector&, const Vector&);

	friend Vector operator+(const Vector&);
	friend Vector operator-(const Vector&);
	friend Vector operator+(const Vector&, const Vector&);
	friend Vector operator-(const Vector&, const Vector&);
	friend double Dot(const Vector&, const Vector&);
	friend Vector Cross(const Vector&, const Vector&);
	friend Vector operator*(const double, const Vector&);
	friend Vector operator*(const Vector&, const double);
	friend Vector operator/(const Vector&, const double);                   
	friend double Cross2D(const Vector &, const Vector &);
	friend Vector Rotate90(const Vector&);

	friend Vector Rotate90unit(const Vector&);                      


	friend std::ostream& operator<<(std::ostream&, const Vector&);
};


typedef std::vector<Vector> Vectorvec;

#endif
