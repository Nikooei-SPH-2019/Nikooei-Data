//****************************************************************************************
//
// SPH Code
//
// created	86/01 By R. Ramazani
// Last change	86/01 By Fatehi
////// Last change	97/11 By Nikooei
//****************************************************************************************

#ifndef _TENSOR_H_
#define _TENSOR_H_

#include "vectors.h"
#include "stensors.h"

class Tensor{

	private:
		double ets[4];
	public:
		Tensor(){ets[0]=ets[1]=ets[2]=ets[3]=0.0;}
		Tensor(double nxx, double nxy, double nyx, double nyy){ets[0]=nxx;ets[1]=nxy;ets[2]=nyx;ets[3]=nyy;}
		Tensor(const Tensor &p){ets[0]=p.ets[0];ets[1]=p.ets[1];ets[2]=p.ets[2];ets[3]=p.ets[3];}
		Tensor(double a){ets[0]=a;ets[1]=0.0;ets[2]=0.0;ets[3]=a;}
		double TwoNorm() const {return sqrt(ets[0]*ets[0]+ets[1]*ets[1]+ets[2]*ets[2]+ets[3]*ets[3]);}
		void SetAll(double nxx, double nxy, double nyx, double nyy){ets[0]=nxx;ets[1]=nxy;ets[2]=nyx;ets[3]=nyy;}


	double GetXX() const {return ets[0];}
	double GetXY() const {return ets[1];}
	double GetYX() const {return ets[2];}
	double GetYY() const {return ets[3];}

	STensor GetSymPart() const;
	Tensor Inverse() const;
	void Transpose();
	double Trace() const{

	 return ets[0]+ets[3];

	}
		

	void SetXX(double nxx){ets[0]=nxx;}
	void SetXY(double nxy){ets[1]=nxy;}
	void SetYX(double nyx){ets[2]=nyx;}
	void SetYY(double nyy){ets[3]=nyy;}

	Tensor& operator+=(const Tensor &);
	Tensor& operator-=(const Tensor &);

	friend Tensor operator+(const Tensor&);
	friend Tensor operator-(const Tensor&);
	friend Tensor operator+(const Tensor&, const Tensor&);
	friend Tensor operator-(const Tensor&, const Tensor&);
	friend Tensor operator*(const double, const Tensor&);
	friend Tensor operator*(const Tensor&, const double);
	friend Tensor operator*(const Tensor&, const Tensor&);
	friend Tensor operator*(const STensor&, const Tensor&);
	friend Tensor operator*(const STensor&, const STensor&);
	friend Vector operator*(const Tensor&, const Vector&);
	friend Vector operator*(const Vector &, const Tensor &);
	friend Tensor OProduct(const Vector&, const Vector&);
	friend Tensor ToTensor(STensor);
	friend double DiadProduct(const Tensor &, const Tensor &);
	friend STensor OProduct_same_vec(const Vector&, const Vector&);


	double& operator[](int i) {return ets[i];}
	Tensor& operator=(const Tensor&);

	friend std::ostream& operator<<(std::ostream&, const Tensor&);
};

typedef std::vector<Tensor> Tensorvec;

Tensor operator+(const Tensor&);
Tensor operator-(const Tensor&);
Tensor operator+(const Tensor&, const Tensor&);
Tensor operator-(const Tensor&, const Tensor&);
Tensor operator*(const double, const Tensor&);
Tensor operator*(const Tensor&, const double);
Tensor operator*(const Tensor&, const Tensor&);
Tensor operator*(const STensor&, const STensor&);
Vector operator*(const Tensor&, const Vector&);
Vector operator*(const Vector &, const Tensor &);
Tensor OProduct(const Vector&, const Vector&);

Tensor ToTensor(STensor);
double DiadProduct(const Tensor &, const Tensor &);


 inline Tensor Rotate(double t)		//counter-clockwise rotation by angle t
 {
	 return Tensor(cos(t), -sin(t), sin(t), cos(t));
 }




#endif

