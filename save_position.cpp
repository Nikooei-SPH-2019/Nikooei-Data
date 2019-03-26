//****************************************************************************************
//
// SPH Code
//
// created	85/06 By Safdari
// Last change	86/07 By R. Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
// Writes last position of saving point in "lastsave" file
//
//****************************************************************************************
#include "general.h"

void SavePosition(std::ofstream &ofs)
{
	std::streampos y = ofs.tellp( );
    unsigned long int savebyte;
    savebyte=y;

	std::ofstream save("lastsave", std::ios::out | std::ios::binary );
    save.write((char *) &savebyte, sizeof savebyte);

	std::cout << std::endl <<"Last saving byte= "<< savebyte << "...";
    save.close();
}