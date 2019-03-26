//****************************************************************************************
//
// SPH Code
//
// created	86/07 By Fatehi
// Last change	86/09 By Fatehi
// Last change	87/11 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
// Periodic boundary file
//	Handles periodic boundary conditions.
//****************************************************************************************

#include "boundary_conditions.h"

void TPointBC::AddPoint(const TIParticle &point)
{
	points.push_back(point);
}
void TPointBC::Write(std::ofstream &dfile, int pNo)
{
	int ppNo = points.size();
	for (register int i = 0; i<ppNo; i++)
	{
		dfile << std::endl;
		dfile << i + pNo << "  ";
		dfile << points[i];
	}
}



