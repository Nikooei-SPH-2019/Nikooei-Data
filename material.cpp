//****************************************************************************************
//
// SPH Code
//
// created	86/09 By Fatehi
// Last change	86/09 By Fatehi
// Last change	87/11 By Fatehi
// Last change	88/02 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************

#include "material.h"
#include "particle.h"

void TMat::Read( const char* MFile, bool &surfaceopt)
{
	short phases=0;

	std::ifstream mfile(MFile);
	if(!mfile) 
	{
		FatalError("Cannot open material file.");
	}

	int matNo;
	short type;
	//double dtemp;
	std::string temp;
	
	std::cout << "Reading material file..." << std::endl;
	

	mfile >> temp >> matNo;
	
	for(register short i=0; i<matNo; i++)
	{

		mfile >> temp >> type;
		if (type==-1)
		{
			materials.push_back(new TMaterial(type));
			std::cout << " A dummy material is loaded." << std::endl;
		}
		else if (type>=wall1 && type<solid1)
		{
			double rho, k, cv;
			mfile >> temp >> rho;
			mfile >> temp >> k;
			mfile >> temp >> cv;

			materials.push_back(new TRMaterial(type, rho, k, cv));
			std::cout << " Wall"<<type<<" is loaded." << std::endl;
		}
		else if (type>=solid1 && type<fluid1)
		{
			double rho, k, cv, E;
			mfile >> temp >> rho;
			mfile >> temp >> k;
			mfile >> temp >> cv;
			mfile >> temp >> E;

			materials.push_back(new TSolid(type, rho, k, cv, E));
			std::cout << " Solid"<<type<<" is loaded." << std::endl;
		}
		else if (type>=fluid1 && type<20)
		{
			double rho, k, cv, mu;
			mfile >> temp >> rho;
			mfile >> temp >> k;
			mfile >> temp >> cv;
			mfile >> temp >> mu;
		
			materials.push_back(new TFluid(type, rho, k, cv, mu ));
			std::cout << " Fluid"<<type-9<<" is loaded." << std::endl;

			phases= max( type, phases);
		}


	}

	doublevec st;
	st.assign(phases+1, 0.0);
	sigmaij.assign(phases+1, st);

	short t1, t2;
	double stij;
	mfile >> temp >> t1 >> t2 >> stij;
	while (!mfile.eof())
	{
		if (stij)
			surfaceopt=true;
		sigmaij[min(t1, t2)][max(t1, t2)]= stij;
		sigmaij[max(t1, t2)][min(t1, t2)]=-stij;
		mfile >> temp >> t1 >> t2 >> stij;
	}


	mfile.clear();
	mfile.close();

	std::cout << "Number of materials= " << materials.size() << std::endl;
	std::cout << "###############################################################################"<< std::endl;

}


int FindType(std::vector<TMaterial*> & materials, int ptype)
{
	register unsigned int i;
	for (i=0; i<materials.size(); i++)
		if (materials[i]->type==ptype)
			break;
	
	return i;
}
