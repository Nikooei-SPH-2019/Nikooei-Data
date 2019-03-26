//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change		85/12 By Fatehi
////// Last change	97/11 By Nikooei
//****************************************************************************************


#ifndef _STENSOR_H_
#define _STENSOR_H_

#include "vectors.h"

class STensor {

private:
	double ets[3];
public:
	STensor() { ets[0] = ets[1] = ets[2] = 0.0; }
	STensor(double d) { ets[0] = d; ets[1] = 0.0; ets[2] = d; }
	STensor(double nxx, double nxy, double nyy) { ets[0] = nxx; ets[1] = nxy; ets[2] = nyy; }
	STensor(const STensor &p) { ets[0] = p.ets[0]; ets[1] = p.ets[1]; ets[2] = p.ets[2]; }
	void SetAll(double nxx, double nxy, double nyy) { ets[0] = nxx; ets[1] = nxy; ets[2] = nyy; }


	double GetXX() const { return ets[0]; }
	double GetXY() const { return ets[1]; }
	double GetYY() const { return ets[2]; }

	void SetXX(double nxx) { ets[0] = nxx; }
	void SetXY(double nxy) { ets[1] = nxy; }
	void SetYY(double nyy) { ets[2] = nyy; }

	STensor& operator+=(const STensor &);
	STensor& operator-=(const STensor &);

	STensor Inverse() const;
	double Trace() const;
	double secondinvariant() const;

	friend STensor operator+(const STensor&);
	friend STensor operator-(const STensor&);
	friend STensor operator+(const STensor&, const STensor&);
	friend STensor operator-(const STensor&, const STensor&);
	friend STensor operator*(const double, const STensor&);
	friend STensor operator*(const STensor&, const double);
	friend Vector operator*(const STensor&, const Vector&);
	friend Vector operator*(const Vector&, const STensor&);
	friend double DiadProduct(const STensor &, const STensor &);

	double& operator[](int i) { return ets[i]; }
	STensor& operator=(const STensor&);

	friend std::ostream& operator<<(std::ostream&, const STensor&);
};
typedef std::vector<STensor> STensorvec;

#endif
