//****************************************************************************************
//                                 In the name of GOD
//****************************************************************************************
// main file
//	Main body of the program
//****************************************************************************************




#include "timer.h"
#include "world.h"
#include "inout.h"
#include <conio.h>
#include <cmath>



int main()
{
//----------------------------------------------------------------------------------------
// Decleration Part
	
	World world;
	Timer timer;
	std::time_t tstart, tfinal;
	int	step=0;
	TPointBC::w=&world;

//----------------------------------------------------------------------------------------
//Body of the code
	
	

 	world.options.ReadFile();

	world.mat.Read(world.options.MATFILE, world.options.SURFACETENSION);

	world.ReadDataFile();
	
	std::time(&tstart);
	world.Initialize();
	
	//Time-stepping loop
	//----------------------------------------------------------------------------------------
	timer.Start();
	timer.Pause();

	world.SetOMEGA(0.0);
	

	if (world.GetDT() < 1e-10 || world.GetDT() > 1e1)   /// It's default value is -9.25e+61  ... this is useful when we resume.
		world.EvalTimeStep();

	while (world.GetTime()< world.options.FINALTIME)
	{
		timer.Resume();
		
		world.SetSTEP(step+1);
		world.DoTimeIntegration();			

		world.UpdateBoundary();

		world.March();

		timer.Pause();
		step++;

		
 		if (!(step%world.options.REPORTSTEP)) 
		{
		world.EvalMaximum();

			world.EvalMinimum();
			world.CheckLimits();
			world.EvalTimeStep();
			timer.Resume();

			timer.Stop();
			world.Report(step,timer.CPUTime());
			timer.Start();
			timer.Pause();
		}

		double ts=world.GetTime()/world.options.SAVINGTIME;
		if ((ts - (int)ts) > (1.0-(world.GetDT()/world.options.SAVINGTIME)))
			world.SaveDataFile();

} 
	//----------------------------------------------------------------------------------------


	std::time(&tfinal);
	PrintTotalTime(tstart, tfinal);
	std::cout.width(8);
	std::cout.precision(5);
	std::cout << "Pmax= "<<world.maximum.p<< '\t'+ "Pmin= "<< world.minimum.p <<'\t'+"DP= "<< world.maximum.p-world.minimum.p<<'\n';


//----------------------------------------------------------------------------------------

	//_getch();

	return 0;

}