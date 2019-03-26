//****************************************************************************************
//
// SPH Code
//
// created	86/05 By Fatehi
// Last change	86/08 By Fatehi
// Last change	87/11 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
//  linkedlist file
//	Finds interaction through a grid-based search.
//****************************************************************************************

#ifndef _LINKEDLIST_H_
#define _LINKEDLIST_H_

#include "vectors.h"
#include "particle.h"

class TCell{
public:
	Vector lbloc;		//Location of left-bottom point of the cell
	std::vector<TCell*> neighbours;
	std::vector<int> iparticles;

	TCell()  {}
	~TCell() {}
};

class TGrid{
private:
	int xNo, yNo;			//Number of cells in X and Y direction

	Vector lbloc,rtloc;		//Location of left-bottom and right-top of the cell
	double delta;			//increment in X and Y direction

public:

	std::vector< std::vector<TCell> > cells;

	TGrid(Vector p1, Vector p2, double hmax);
	~TGrid(){}
	void Initialize();
	void AssignParticles(World &mw);
	void FindNeighbourCells();
};



#endif //_LINKEDLIST_H_

