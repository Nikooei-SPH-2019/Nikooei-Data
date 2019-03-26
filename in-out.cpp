//****************************************************************************************
//
// SPH Code
//
// created	85/12 By Fatehi
// Last change 86/08 By Fatehi
// Last change 87/08 By Fatehi
// Last change 88/04 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************

#include"general.h"

void Wellcome()
{
	std::cout << "###############################################################################"<< std::endl;
	std::cout << "#	                   This is SePeHr v 1.8 an SPH code                   #"<< std::endl;
	std::cout << "###############################################################################"<< std::endl;
}


void WriteHeader()
{
	std::cout << "###############################################################################"<< std::endl;
	std::cout << "  Step       Time      Max.Acc     Max.V       C     dx/h     Err.  CPU.T/step" << std::endl;
	std::cout << "-------------------------------------------------------------------------------" << std::endl;
}

