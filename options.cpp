//****************************************************************************************
//
// SPH Code
//
// created	85/12 By Fatehi
// Last change	86/08 By Fatehi
// Last change	87/07 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************

#include "options.h"
#include <cstring>
//#include "parallelDebug.h"

const char* const OPTIONFILE="options.opt";  //Add const a.r

void TOptions::ReadFile()
{

	std::string temp;
	int tempint;
	double tempdouble1,tempdouble2;

	//Check option file.

	std::ifstream ofile(OPTIONFILE);
	if(!ofile)
	{
		FatalError("Cannot open option file.");
	}

	//Read Options from option file.

	std::cout << "Reading option file...";

	ofile >> temp >> INITIALTIME;
	ofile >> temp >> FINALTIME;
	ofile >> temp >> MAXCFL;

	ofile >> temp >> XSPH;
	ofile >> temp >> Shifting;
	ofile >> temp >> ARTIFITIALVISCOSITY;
	ofile >> temp >> BULKVISCOSITY;

	ofile >> temp >> REPORTSTEP;
	ofile >> temp >> SAVINGTIME;

	ofile >> temp >> WRITINGMETHOD;

	ofile >> temp >> SEARCHMETHOD;
	ofile >> temp >> KERNELTYPE;
	ofile >> temp >> SOLUTIONMETHOD;
	ofile >> temp >> FREESURFACE;

	ofile >> temp >> MAXIMUMX;
	ofile >> temp >> MAXIMUMY;

	ofile >> temp >> MAXIMUMVELOCITY;
	ofile >> temp >> MAXIMUMACCELERATION;
	ofile >> temp >> MAXIMUMPRESSURE;
	ofile >> temp >> MAXIMUMH;

	ofile >> temp >> MINIMUMX;
	ofile >> temp >> MINIMUMY;

	ofile >> temp >> MINIMUMPRESSURE;
	ofile >> temp >> MINIMUMH;
	ofile >> temp >> MINIMUMDISTANCE;


	ofile >> temp >> tempint;
	MOMENTUM=(Momentom_Type)tempint;
	ofile >> temp >> tempint;
	ENERGY=(Energy_Type)tempint;

	ofile >> temp >> NORMALDENSITY;


	ofile >> temp >> C_0;
	ofile >> temp >> V_0;
	ofile >> temp >> tempdouble1 >> tempdouble2;

	BODYACC=Vector(tempdouble1,tempdouble2);

	ofile >> temp >> temp;
	strcpy(DATAFILE, temp.c_str());

	ofile >> temp >> temp;
	strcpy(MATFILE, temp.c_str());


	ofile.close();
	std::cout << " OK." << std::endl;


	SUMDENSITY=(SOLUTIONMETHOD==1);
	SURFACETENSION=false;

}



void TOptions::Report()
{
	// Report input data

	std::cout << "Data file name= " << DATAFILE << std::endl;
	std::cout << "Material file name= " << MATFILE << std::endl;
	std::cout << "###############################################################################"<< std::endl;
	std::cout << "Initial time (s)= " << INITIALTIME << std::endl;
	std::cout << "Final time (s)= " << FINALTIME << std::endl;
	std::cout << "Maximum CFL Number= " << MAXCFL << std::endl;
	std::cout << "Saving time (s)= " << SAVINGTIME << std::endl;
	std::cout << "###############################################################################"<< std::endl;
	std::cout << "Search method =	";
	switch(SEARCHMETHOD)
	{

	case 1:
		std::cout << "Linked List";
		break;
	
	}
	std::cout << std::endl;


	std::cout << "Kernel type =	";
	switch(KERNELTYPE)
	{

	case 1:
		std::cout << "Quintic Function (Wendland)";
		break;
		}
	std::cout << std::endl;

	std::cout << "###############################################################################"<< std::endl;
	if (MOMENTUM)
	{
		if (XSPH)
			std::cout << "XSPH= " << XSPH << std::endl;
		if (ARTIFITIALVISCOSITY)
			std::cout << "ARTIFITIAL VISCOSITY= " << ARTIFITIALVISCOSITY << std::endl;
		if (BULKVISCOSITY)
			std::cout << "BULK VISCOSITY= " << BULKVISCOSITY << std::endl;
		if (NORMALDENSITY)
			std::cout << "NORMAL DENSITY= " << NORMALDENSITY << std::endl;
		std::cout << "###############################################################################"<< std::endl;
	}

	switch(MOMENTUM)
	{
	
	case MLAMINAR2:
		std::cout << "Momentum is solved as LAMINAR2";
		break;
	default:
		FatalError("Wrong Momentum Solver!");
	}
	std::cout << std::endl;


	if (MOMENTUM)
	{
		std::cout << "Solution method = ";
		switch(SOLUTIONMETHOD)
		{
		
		case 1:
			std::cout << "Predictor-Corrector (Summation)";
			break;
	
		default:
			FatalError("Wrong solution method!");
		}
		std::cout << std::endl;
	}

	if (FREESURFACE && FREESURFACE<1.0)
	{
		std::cout << "Free surface option is activated";
		std::cout << std::endl;

	}

	switch(ENERGY)
	{
	case ENONE:
		break;
		default:
		FatalError("Wrong Energy Solver!");
	}
	std::cout << std::endl;
	std::cout << "###############################################################################"<< std::endl;

}



