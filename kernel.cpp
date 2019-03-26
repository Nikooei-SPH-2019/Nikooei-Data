//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change		86/07 By Fatehi
// Last change		87/07 By Fatehi
// Last change		87/12 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************


#include "kernel.h"

double Kernel(double d, double h, short KERNELTYPE)
{
	const double  r = d / h;
	double  w;

	switch (KERNELTYPE)
	{

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 1:                         //Lucy Function [Lucy 1977]

		w = 5.0 / (pi*h*h);


		if (r<1.0)
			return w*((1.0 + 3.0*r)*(1.0 - r)*(1.0 - r)*(1.0 - r));

		return 0.0;

		
	}


	FatalError("Wrong Kernel type");
}



double DKernel(double d, double h, short KERNELTYPE)
{
	const double  r = d / h;
	double  w;

	switch (KERNELTYPE)
	{

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 1:                         //Lucy Function [Lucy 1977]

		w = -60.0 / (pi*h*h*h);

		if (r<1.0)
			return w*(r*(1.0 - r)*(1.0 - r));

		return 0.0;
		
	}
	

	FatalError("Wrong Kernel type");
}

double D2Kernel(double d, double h, short KERNELTYPE)
{
	const double  r = d / h;
	double  w;

	switch (KERNELTYPE)
	{

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	case 1:                         //Lucy Function [Lucy 1977]

		w = -60.0 / (pi*h*h*h*h);


		if (r<1.0)
			return w*((1.0 - 3.0*r)*(1.0 - r));

		return 0.0;
		
	}

	FatalError("Wrong Kernel type");
}



