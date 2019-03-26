//****************************************************************************************
//
// SPH Code
//
// created	86/05 By Fatehi
// Last change	86/07 By Fatehi
// Last change	87/08 By Fatehi
// Last change	87/11 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
//  writeASCII file
//	writes ASCII.
//****************************************************************************************
#include "world.h"
void ExactSolution(World & mw)
{
	int pNo=mw.particles.size();		///Conduction problem
	int minf=50000;
	if (mw.GetTime())
		minf=minf/10;
	const double sinhpi=sinh(pi);
	for(register int i=0; i<pNo; i++)
	{
		const double x=mw.particles[i].GetLoc().GetX();
		const double y=mw.particles[i].GetLoc().GetY();
		double sum=0.0;
		for(int m=minf; m>0; m--)
		{
			const double beta=pi*pi*(1.0+(double)m*(double)m);
			const int minustom= 1-2*(m % 2);
			double ebeta=1.0;
			if (mw.GetTime())
				ebeta=exp(-beta*mw.GetTime());
			sum += (2.0*m)/pi/(1.0+(double)m*(double)m)*minustom*ebeta*sin(m*pi*y);
		}

		const double Tex=(sum+sinh(pi*y)/sinhpi)*sin(pi*x);
		Vector u;
		u.SetX(Tex);
		mw.particles[i].SetVel(u);
	}
}

void WriteASCIIFile(World & mw, const char* DFile,const char* szTitle)
{

	std::cout << "Saving ASCII data file...";

	std::ofstream dfile;
	dfile.open( DFile, std::ios::app );

	int pNo=mw.particles.size();
	
	dfile << std::endl << szTitle;

	dfile << std::setprecision(16);
	dfile << std::scientific;

	SavePosition(dfile);

	for(register int i=0; i<pNo; i++)
	{
		dfile << std::endl;
		dfile << i << '\t';
		dfile << mw.particles[i];
	}


	
	
	//***********************************************************BY MAFB & MSS
	for (register unsigned short i=0; i<mw.wallBCVector.size(); i++)			// Wall Boundary Conditions
	{
		mw.wallBCVector[i].Write(dfile, pNo);
		pNo+=mw.wallBCVector[i].GetSize();
	}
	//***********************************************************End BY MAFB & MSS

		
	dfile.close();


}
void ClearASCIIFile( const char* DFile)
{

	std::cout << "Clearing ASCII data file...";

	std::ofstream dfile;
	dfile.open( DFile, std::ios::out | std::ios::trunc);

	dfile.close();

}



