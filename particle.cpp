//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change	86/07 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
// Particle code file
//	Defines Class Particle & enumerator Particle_Type
//****************************************************************************************

#include "particle.h"


std::ostream &operator<<(std::ostream &stream, TFParticle &par)
{
//	if (par.GetBCindex()==wallimage)  
	//	stream << par.GetBCindex() << '\t';
	//else 
		stream << par.GetMaterial()->type << '\t';

	stream << par.GetLoc() << '\t';
	stream << par.GetVel() << '\t';
	stream << par.GetP() << '\t';
	stream << par.GetTau() << '\t';
	stream << par.GetM() << '\t';
	stream << par.GetH() << '\t';
	stream << par.GetMaterial()->rho<< '\t';		//Estimated Density
	stream << par.GetGammaDot() << '\t';
	stream << par.GetVorticity() << '\t';
	stream << par.Get_Okubo_Weiss();
	return stream;  // return the stream
}

std::ostream &operator<<(std::ostream &stream, TIParticle &par)
{
	stream << par.index << '\t';
	stream << par.GetLoc() << '\t';
	stream << par.normal << '\t';
	stream << 0.0 << '\t';
	stream << STensor() << '\t';
	stream << 0.0 << '\t';
	stream << 0.0 << '\t';
	stream << 0.0;	

	return stream;  // return the stream
}

