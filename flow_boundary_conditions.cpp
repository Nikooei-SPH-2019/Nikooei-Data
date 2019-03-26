//****************************************************************************************
//
// SPH Code
//
// created	86/09 By Fatehi
// Last change	86/09 By Fatehi
// Last change	87/11 By Fatehi
// changed	88/02 By Fatehi
// changed	88/09 By Fatehi
// changed	88/10 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
// Flow boundary file
//	Handles flow boundary conditions.
//****************************************************************************************

#include "boundary_conditions.h"
#include "world.h"

static const double trshout = 1.91;
static const double trshin = 1.91;
World* (TPointBC::w) = 0;
