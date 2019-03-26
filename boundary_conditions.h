//****************************************************************************************
//
// SPH Code
//
// created	86/09 By Fatehi
// Last change	86/09 By Fatehi
// Last change	87/08 By Fatehi
// Last change	87/11 By Fatehi
// changed	88/10 By Fatehi
// changed	89/02 By Fatehi
// changed	89/03 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
// boundary conditions file
//	Decleration of all boundary conditions.
//****************************************************************************************

#ifndef _BOUNDARY_CONDITIONS_H_
#define _BOUNDARY_CONDITIONS_H_

#include "particle.h"
#include "interaction.h"
#include "options.h"
class World;		


class TPointBC
{
protected:
	int index;
	std::vector<TIParticle> points;

public:
	static World* w;

	TPointBC(void) {}
	TPointBC(int i) : index(i) {}
	~TPointBC(void) {}

	int GetSize() { return points.size(); }
	int GetIndex() { return index; }
	void AddPoint(const TIParticle &point);
	void Write(std::ofstream &dfile, int pNo);
	virtual void Update(std::vector<TFParticle> &, double) {}

};

//////////////////////////////////////////////////////////////////////////////
// C  L  A  S  S     D  E  C  L  E  R  A  T  I  O  N
// Class:       TFlowBC
// Purpose:     Base for in- and outflow boundary conditions.
// Author:      R. Fatehi
// Revision:    1.0
// Orig. Date:  Dec. 2007
// Rev. Date:   March. 2019
//////////////////////////////////////////////////////////////////////////////
class TFlowBC : public TPointBC
{
protected:
	std::vector<std::vector<TFParticle>::iterator> particles_iter;

public:

	TFlowBC(void) {}
	TFlowBC(int i) :TPointBC(i) {}
	~TFlowBC(void) {}
		virtual void InitBCValue() {}
	virtual void InitBCValue(double) {}
	virtual void SetV() {}
	virtual void SetA() {}
	virtual void SetP() {}

};

//***********************************************************BY US
//////////////////////////////////////////////////////////////////////////////
// C  L  A  S  S     D  E  C  L  E  R  A  T  I  O  N
// Class:       TWallBC
// Purpose:     Wall boundary condition. 
// Author:      US
// Revision:    1.0
// Orig. Date:  Oct. 2008
// Rev. Date:   March. 2019
//////////////////////////////////////////////////////////////////////////////
class TWallBC : public TPointBC
{
public:
	TWallBC(void) {}
	~TWallBC(void) {}

	void EvalNormals();

};

//////////////////////////////////////////////////////////////////////////////
// C  L  A  S  S     D  E  C  L  E  R  A  T  I  O  N
// Class:       TwallBCVector
// Purpose:     Vector of wall boundary conditions. 
// Author:      US
// Revision:    1.0
// Orig. Date:  Oct. 2008
// Rev. Date:   March. 2019
//////////////////////////////////////////////////////////////////////////////
class TwallBCVector : public std::vector <TWallBC>
{
private:

public:
	std::vector <std::vector <TFParticle>::iterator> imagesIter;
	//std::vector <TIParticle> imageparticles;

	TwallBCVector(void) {
		assign(20, TWallBC());
	}
	~TwallBCVector(void) {}

	void Update(std::vector<TFParticle> &particles, bool SLIP);
	void UpdateValues(std::vector<TFParticle> &particles, bool SLIP);
	void EvalNormals();
	void clearImages(std::vector<TFParticle> &particles);
};
//***********************************************************END BY US

#endif //_BOUNDARY_CONDITIONS_H_
