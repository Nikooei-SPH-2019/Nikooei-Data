
#include "kernel.h"
#include "inout.h"
#include <stack>

#include "world.h"



static double rhomax = 0.0;						//(kg/m3)Density
static double rhomin = 1e10;						//(kg/m3)Density
static double numax = 0.0;						//(m2/s)Momentum diffusivity
static double alphamax = 0.0;					//(m2/s)Energy diffusivity
static double delta = 0.005;				//Maximum allowable relative density variation



void World::DoTimeIntegration()
{
	if (options.MOMENTUM)
	{
		switch (options.SOLUTIONMETHOD)
		{
		
		case 1:   //// By Mohammd Nikooei

			Projection3(*this);

			break;

	
		}
	}
}



void World::EvalMomentum()
{

	BodyForce(*this);


	switch (options.MOMENTUM)
	{

	case MLAMINAR2:



		if (options.SOLUTIONMETHOD == 1)
		{
			Granularmaterial1(*this);

		}

				
		break;


	default:
		FatalError("Wrong momentum type!");

	}




}



void World::EvalVStar() //////////////////////Nikooei
{


	switch (options.MOMENTUM)
	{

	case MLAMINAR2:


		break;


	default:
		FatalError("Wrong momentum type!");

	}





}







void World::UpdateInteractions()
{
	int iNo = interactions.size();

	for (register int i = 0; i<iNo; i++)
	{
		const int ii = interactions[i].GetI();
		const int jj = interactions[i].GetJ();

		interactions[i].SetVD(particles[ii].GetLoc() - particles[jj].GetLoc());
		interactions[i].SetDist(interactions[i].GetVD().TwoNorm());
		interactions[i].SetHMean(0.5*(particles[ii].GetH() + particles[jj].GetH()));
		interactions[i].SetW(Kernel(interactions[i].GetDist(), interactions[i].GetHMean(), options.KERNELTYPE));
		interactions[i].SetDW(DKernel(interactions[i].GetDist(), interactions[i].GetHMean(), options.KERNELTYPE));
	}


}



void World::EvalDivergence3()
//[Bonet 2004] [Randles, Libersky 1996]
{
	const int pNo = particles.size();
	std::vector<double> DelDotr(pNo);

	for (register int i = 0; i<pNo; i++)
		particles[i].SetVDivergence(0.0);


	const int iNo = interactions.size();

	for (register int i = 0; i<iNo; i++)				//Compute Divergence of Velocity
	{
		const int ii = interactions[i].GetI();
		const int jj = interactions[i].GetJ();
		if ((particles[ii].GetMaterial()->type<fluid1) && (particles[jj].GetMaterial()->type<fluid1))
			continue;
		const double d = interactions[i].GetDist();
		Vector dx = interactions[i].GetVD();
		const Vector dw = interactions[i].GetDW() / d*dx;

		const Vector dv = (particles[jj].GetVel() - particles[ii].GetVel());
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction

		const double rhoij = 0.5*(particles[ii].GetRho() + particles[jj].GetRho());

		particles[ii].SetVDivergence(particles[ii].GetVDivergence() + particles[jj].GetM() / rhoij*Dot(dv, particles[ii].B*dw));
		particles[jj].SetVDivergence(particles[jj].GetVDivergence() + particles[ii].GetM() / rhoij*Dot(dv, particles[jj].B*dw));


		DelDotr[ii] = DelDotr[ii] + particles[jj].GetM() / rhoij*Dot(dx, dw);
		DelDotr[jj] = DelDotr[jj] + particles[ii].GetM() / rhoij*Dot(dx, dw);

		int fluidneighbors = 0;
		double hm;
		double dr;

	}
	for (register int i = 0; i < pNo; i++)
		particles[i].SetDelDotr(DelDotr[i]);
}




void World::UpdateBoundary()
{
	//***********************************************************BY MAFB & MSS & SQ
	wallBCVector.clearImages(particles);   // Clear Wall Boundary Conditions

	wallBCVector.Update(particles, (options.MOMENTUM == 1));   // Wall Boundary Conditions

															   //***********************************************************END BY MAFB & MSS & SQ



}



void World::Report(int stepc, double t)   //OK
{
	std::cout.width(7);
	std::cout << stepc;
	std::cout.precision(3);
	std::cout.width(10);
	std::cout << time;
	std::cout.precision(2);
	std::cout.width(11);
	std::cout << maximum.acc;
	std::cout.precision(4);
	std::cout.width(11);
	std::cout << maximum.vel;
	std::cout.width(11);
	if ((options.MOMENTUM) && !(options.SOLUTIONMETHOD == 30) && !(options.SOLUTIONMETHOD == 40))
	{
		std::cout.precision(3);
		std::cout << sqrt(c2);
	}
	else

		std::cout << options.BODYACC.GetX();// by S.J

	std::cout.precision(2);
	std::cout.width(9);
	std::cout << minimum.distance / minimum.h;
	std::cout.precision(2);
	std::cout.width(9);
	std::cout << Error();
	std::cout.precision(2);
	std::cout.width(8);
	std::cout << t / double(options.REPORTSTEP);
	std::cout << std::endl;

}





void World::EvalMinimum()  //OK
{

	minimum.loc = Vector(options.MAXIMUMX, options.MAXIMUMY);

	minimum.p = options.MAXIMUMPRESSURE;
	minimum.h = options.MAXIMUMH;
	minimum.distance = options.MAXIMUMH;


	int pNo = particles.size();

	for (register int i = 0; i<pNo; i++)
	{
		if (particles[i].GetLoc().GetX()<minimum.loc.GetX())
			minimum.loc.SetX(particles[i].GetLoc().GetX());

		if (particles[i].GetLoc().GetY()<minimum.loc.GetY())
			minimum.loc.SetY(particles[i].GetLoc().GetY());

		if (particles[i].GetMaterial()->type < fluid1)
			continue;

		if (particles[i].GetP()<minimum.p)
			minimum.p = particles[i].GetP();

		if (particles[i].GetH()<minimum.h)
			minimum.h = particles[i].GetH();
	}

	int iNo = interactions.size();

	for (register int i = 0; i<iNo; i++)
	{
		if (interactions[i].GetDist()<minimum.distance)
			minimum.distance = interactions[i].GetDist();
	}

}

void World::EvalMaximum()  //OK
{

	maximum.loc = Vector(options.MINIMUMX, options.MINIMUMY);

	maximum.vel = 0.0;

	maximum.acc = 0.0;
	maximum.p = options.MINIMUMPRESSURE;
	maximum.h = 0.0;
	maximum.vdivergence = 0.0;


	int pNo = particles.size();

	for (register int i = 0; i<pNo; i++)
	{
		if (particles[i].GetLoc().GetX() > maximum.loc.GetX())
			maximum.loc.SetX(particles[i].GetLoc().GetX());

		if (particles[i].GetLoc().GetY() > maximum.loc.GetY())
			maximum.loc.SetY(particles[i].GetLoc().GetY());


		if (particles[i].GetVel().TwoNorm() > maximum.vel && particles[i].GetLoc().GetX()>0.02)   ////////For ignoring the velocity of particles near the vertical wall
			maximum.vel = particles[i].GetVel().TwoNorm();

		if (particles[i].GetMaterial()->type < fluid1)
			continue;

		if (particles[i].GetAcc().TwoNorm() > maximum.acc)
			maximum.acc = particles[i].GetAcc().TwoNorm();

		if (particles[i].GetP() > maximum.p)
			maximum.p = particles[i].GetP();

		if (particles[i].GetH() > maximum.h)
			maximum.h = particles[i].GetH();

		if (std::abs(particles[i].GetVDivergence()) > maximum.vdivergence)
			maximum.vdivergence = std::abs(particles[i].GetVDivergence());
	}

}


void World::CheckLimits()  //OK
{

	if (minimum.loc.GetX()< options.MINIMUMX
		|| minimum.loc.GetY()< options.MINIMUMY)
		FatalError("Minimum limit for (x or y) is exceeded!");
	if (minimum.p < options.MINIMUMPRESSURE)
		FatalError("Minimum limit for prssure is exceeded!");

	if (maximum.loc.GetX()> options.MAXIMUMX
		|| maximum.loc.GetY()> options.MAXIMUMY)
		FatalError("Maximum limit for (x or y) is exceeded!");
	if (maximum.vel > options.MAXIMUMVELOCITY)
		FatalError("Maximum limit for velocity magnitude is exceeded!");
	if (maximum.acc > options.MAXIMUMACCELERATION)
		FatalError("Maximum limit for acceleration magnitude is exceeded!");
	if (maximum.p > options.MAXIMUMPRESSURE)
		FatalError("Maximum limit for prssure is exceeded!");
	if (maximum.h > options.MAXIMUMH)
		FatalError("Maximum limit for h is exceeded!");

}



void World::EvalTimeStep()   //OK
{

	deltat = options.FINALTIME - options.INITIALTIME;

	double dx = max(minimum.distance, 0.1*minimum.h);
	const double amax = sqrt(rhomax / rhomin)*options.BODYACC.TwoNorm();

	if (options.MOMENTUM)
	{
		double vmax;

		const double L = (maximum.loc - minimum.loc).TwoNorm();
		const double aL = Dot(options.BODYACC, (maximum.loc - minimum.loc));

		c2 = max(max(max(options.C_0*options.C_0					//computing speed of sound
			, numax*maximum.vel / L / delta)
			, sqrt(rhomax / rhomin)*aL / delta)
			, maximum.vel*maximum.vel / delta);

		vmax = sqrt(c2) + maximum.vel;



		if (options.MOMENTUM > MINVISCID)
			deltat = min(deltat
				, dx*dx / numax);

		int pNo = particles.size();

		double Max_Mu = 0;
		for (int i = 0; i < pNo; i++)

		{
			double MULTIPLICATION = 10000;
			double Stress = particles[i].GetTau().secondinvariant();
			double GammaDot = particles[i].GetGammaDot();
			double Mu;
			if (Stress == 0)
				Mu = MULTIPLICATION*(particles[i].GetMaterial()->mu);
			else
				Mu = Stress / GammaDot;

			if (Mu > 1000)
				Mu = 1000;

			if (Mu > Max_Mu)
				Max_Mu = Mu;

		}


		if (options.FREESURFACE && options.FREESURFACE < 1.0)

		{
			double vmax = maximum.vel;
			double L = 1;
			double h = 2.6*dx;

			double Density = particles[0].GetRho();
			double nu = Max_Mu / Density;
			deltat = min(min(dx / vmax, dx / sqrt(2 * 9.81*L)), dx*dx / nu);

		}

	}

	if (options.ENERGY)
		deltat = min(deltat
			, dx*dx / alphamax);


	deltat = options.MAXCFL*deltat;

}




void World::Initialize()  //OK
{

	time = options.INITIALTIME;

	//***********************************************************BY  MAFB & MSS & SQ
	wallBCVector.EvalNormals();

	wallBCVector.Update(particles, (options.MOMENTUM == 1));
	//***********************************************************END BY  MAFB & MSS & SQ

	EvalMinimum();  //OK
	EvalMaximum();  //OK
	CheckLimits();  //OK



	std::cout << "Searching for interactions...";
	FindInteractions();
	std::cout << " OK" << std::endl;
	std::cout << "Initial number of interactions= " << interactions.size() << std::endl;
	WriteHeader();

	if (options.SUMDENSITY)
	{
		EvalParticleNumDensity();
	}
	int pNo = particles.size();
	for (int i = 0; i<pNo; i++)
	{
		particles[i].SetRho(particles[i].GetM()*particles[i].GetSi());
		particles[i].p1 = 0.0;
	}

	EvalMinimum();
	dx0 = (particles[0].GetLoc() - particles[1].GetLoc()).TwoNorm(); //minimum.distance;


	TPointBC::w = this;

	int matNo = mat.materials.size();

	for (register int i = 0; i<matNo; i++)
	{
		//if (mat.materials[i]->type == dummy)
			//continue;

		TRMaterial* m = (TRMaterial*)mat.materials[i];

		if (options.MOMENTUM)
		{
			double rho = m->rho;
			if (rho > rhomax)
				rhomax = rho;
			if (rho < rhomin)
				rhomin = rho;
		}

		if (mat.materials[i]->type<10)
			continue;

		TFluid* fluid = (TFluid*)mat.materials[i];
		TNonNewtonian* Nonnewtonfluid = (TNonNewtonian*)mat.materials[i];          

		if (options.MOMENTUM > MINVISCID)
		{
			double nu = fluid->mu / fluid->rho;
			if (nu > numax)
				numax = nu;
		}

		}


}



void World::EvalInitialParticleNumDensity()

{
	int pNo = particles.size();

	int iNo = interactions.size();
	for (register int i = 0; i<pNo; i++)
	{
		particles[i].SetSi_0(0);
	}
	for (register int i = 0; i < iNo; i++)
	{
		const int ii = interactions[i].GetI();
		const int jj = interactions[i].GetJ();
		if (STEP == 1)
		{
			particles[ii].SetSi_0(particles[ii].GetSi_0() + interactions[i].GetW());
			particles[jj].SetSi_0(particles[jj].GetSi_0() + interactions[i].GetW());
		}
	}

}

void World::EvalParticleNumDensity()
{

	int pNo = particles.size();

	int iNo = interactions.size();

	for (register int i = 0; i<pNo; i++)
	{
		particles[i].SetSi(Kernel(0.0, particles[i].GetH(), options.KERNELTYPE));	//Self action
		particles[i].SetSi_0(0);
	}

	for (register int i = 0; i<iNo; i++)
	{
		const int ii = interactions[i].GetI();
		const int jj = interactions[i].GetJ();

		particles[ii].SetSi(particles[ii].GetSi() + interactions[i].GetW());
		particles[jj].SetSi(particles[jj].GetSi() + interactions[i].GetW());


	}


	//***********************************************************BY MAFB & MSS & SQ
	if (wallBCVector.size())
		for (register int i = 0; i<pNo; i++)
		{
			if (particles[i].mainPIndex>-1)
				particles[i].SetSi(particles[particles[i].mainPIndex].GetSi());
		}
	//***********************************************************BY MAFB & MSS & SQ
}


void World::FindInteractions()
{
	switch (options.SEARCHMETHOD)
	{

	case 1:
		LinkedList(*this);
		break;


	}


}




void World::ReadDataFile()
{
	unsigned long int loadbyte;
	std::ifstream load("lastsave", std::ios::in | std::ios::binary);
	load.read((char *)&loadbyte, sizeof loadbyte);
	load.close();

	bool readfromfirst = (loadbyte == 0 || !load);
	char ask = ' ';
	if (!readfromfirst)
	{
		std::cout << "Data file is loading from byte " << loadbyte << std::endl;
		std::cout << "Do you want to contineu (y/n)? ";
		std::cin >> ask;
	}


	if (readfromfirst || ask == 'n')
	{
		std::cout << "Data file is loading from the beginning. " << std::endl;
		bool newFile = ReadASCIIFile(*this, options.DATAFILE, 0);

		if (!newFile)
			WriteASCIIFile(*this, options.DATAFILE, "Zone");

	}
	else if (ask == 'y')
		switch (options.WRITINGMETHOD)
		{
		case 0:
			ReadASCIIFile(*this, options.DATAFILE, loadbyte);

			break;


		}
	else
		FatalError("Improper input!");


	std::cout << " OK." << std::endl;
	std::cout << "Number of particles= " << particles.size() << std::endl;
	std::cout << "###############################################################################" << std::endl;


}




void World::SaveDataFile()
{
	switch (options.WRITINGMETHOD)
	{
	case 0:
		WriteASCIIFile(*this, options.DATAFILE, "Zone");

		break;

	
	}


	std::cout << "OK." << std::endl;
	std::cout << "Number of particles= " << particles.size() << std::endl;
	std::cout << "Number of interactions= " << interactions.size() << std::endl;
	WriteHeader();

}

double World::Error()
{

	return deltat*maximum.vdivergence;
}
