//****************************************************************************************
//
// SPH Code
//
// created	85/06 By Safdari
// Last change	86/07 By Fatehi
// Last change	87/08 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
// Read ASCII file
//
//****************************************************************************************
#include "world.h"
#include <cstring>


bool ReadASCIIFile(World & mw, const char *DFile,unsigned long int loadbyte){

	bool newFile=true;
	std::ifstream dfile(DFile);
	if(!dfile)
	{
		FatalError("Cannot open data file.");
	}

	//Go to Read bytes Position

		int pNo=0;
	short ptype;
	double ptype_r, x, y, vx, vy, p, sxx, sxy, syy, m, h, rho, GammaDot,Vorticity,Okubo_Weiss;

	std::string temp;

	std::cout << "Reading ASCII data file...";

	//Seeking to the beginning of the last zone
	if (loadbyte>0)
		dfile.seekg( loadbyte );
	else
	{
		dfile.ignore(200,'\n');
		dfile >> temp;
	}

	dfile >> temp;
	int IB_Number = -1; /// By S.J
	bool yes_no = false;
	while(!dfile.eof())
	{

		dfile >> ptype_r >> x >> y >>
						 vx >> vy >>
						 p >>
						 sxx >> sxy >> syy >>
						 m >> h >> rho>>GammaDot>> Vorticity>> Okubo_Weiss;


		ptype=(short)ptype_r;

			mw.particles.push_back(TFParticle(
					(Particle_Type)ptype,
					Vector( x, y ),
					m, h,
					Vector(vx, vy ),
					p,
					STensor(sxx, sxy, syy ) ) );

			const int id=FindType(mw.mat.materials, ptype);
			mw.particles[pNo].SetMaterial((TFluid*)mw.mat.materials[id]);

			mw.particles[pNo].SetSi( rho/mw.particles[pNo].GetM() );
		

		pNo++;
		dfile >> temp;
		if (temp=="Zone")
		{
			newFile = false;
			break;
		}
	}

dfile.clear();
	dfile.close();


	for (int i = 19; i >= 0; i--)
		if (!mw.wallBCVector[i].GetSize())
			mw.wallBCVector.erase(mw.wallBCVector.begin() + i);
		else break;

	if (!newFile)
	{
		std::ofstream dfile(DFile);
		dfile << "Variables= \"Num\" \"PType\" \"X\" \"Y\" \"VX\" \"VY\" \"P\" \"sxx\" \"sxy\" \"syy\" \"m\" \"h\" \"rho\" \"GammaDot\"\"Vorticity\"\"Okubo_Weiss\" ";
		return false;
	}


	return true;


}

