//****************************************************************************************
//
// SPH Code
//
// created	88/02 By Fatehi
// changed		89/03 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
// Correction file
//	Difinition of functions used for renormalization of kernel gradient.
//****************************************************************************************


#include "world.h"
#include "tensors3.h"
#include "tensors4.h"


void World::UpdateCorrectors()		////[JP Vila 1999, 2005] [Randles, Libersky 1996]
{
	const int iNo = interactions.size();
	const int pNo = particles.size();
	intvec Bed_neighbor(pNo, 0), Flow_neighbor(pNo, 0);
	intvec Interfaceneighbors(pNo, 0);


	gradSum_r.assign(pNo, STensor());
	gradSum.assign(pNo, Vector());

	for (int i = 0; i<iNo; i++)
	{
		const int ii = interactions[i].GetI();
		const int jj = interactions[i].GetJ();
		const Vector dx = interactions[i].GetVD();
		const double d = interactions[i].GetDist();
		const Vector dw = interactions[i].GetDW() / d*dx;
		const double rhoij = 0.5*(particles[ii].GetRho() + particles[jj].GetRho());
		const STensor dwRbar = interactions[i].GetDW() / d*OProduct(dx, -dx).GetSymPart();
		/////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction
		int typei = particles[ii].GetMaterial()->type;
		int typej = particles[jj].GetMaterial()->type;


		gradSum_r[ii] += particles[jj].GetM() / rhoij*dwRbar;
		gradSum_r[jj] += particles[ii].GetM() / rhoij*dwRbar;

		gradSum[ii] += particles[jj].GetM() / rhoij *dw;
		gradSum[jj] -= particles[ii].GetM() / rhoij *dw;
	

		if (typei < fluid1)
		{
			if (typej >= 15)
				Bed_neighbor[ii] = Bed_neighbor[ii] + 1;

			if (typej >= fluid1 && typej < 15)
				Flow_neighbor[ii] = Flow_neighbor[ii] + 1;

		}

		if (typej < fluid1)
		{
			if (typei >= 15)
				Bed_neighbor[jj] = Bed_neighbor[jj] + 1;

			if (typei >= fluid1 && typei < 15)
				Flow_neighbor[jj] = Flow_neighbor[jj] + 1;

		}


		if ((typei == 11 && typej == 15) || (typej == 11 && typei == 15))
		{
			Interfaceneighbors[ii] = Interfaceneighbors[ii] + 1;
			Interfaceneighbors[jj] = Interfaceneighbors[jj] + 1;
		}



	}
	for (int i = 0; i < pNo; i++)
	{
		particles[i].B = gradSum_r[i].Inverse();
		particles[i].Set_bed_neigh(Bed_neighbor[i]);
		particles[i].Set_overlaying_flow_neigh(Flow_neighbor[i]);
		particles[i].SetInterfaceneighbor(Interfaceneighbors[i]);

	}


	}

void World::UpdateCorrectors(const short m, const short n)		////[JP Vila 1999, 2005] [Randles, Libersky 1996]
{
	const int iNo = interactions.size();
	const int pNo = particles.size();

	STensorvec gradSum_r(pNo, STensor());

	for (int i = 0; i<iNo; i++)
	{
		const int ii = interactions[i].GetI();
		const int jj = interactions[i].GetJ();
		const short typei = max(particles[ii].GetMaterial()->type, (const short)0);
		const short typej = max(particles[jj].GetMaterial()->type, (const short)0);
		if ((min(typei, typej) == m && max(typei, typej) == n)		// m < n
			|| (typei == m && typej == m)
			|| (typei == n && typej == n))
		{
			const Vector dx = particles[ii].GetLoc() - particles[jj].GetLoc();
			const double d = interactions[i].GetDist();

			const STensor dwRbar = interactions[i].GetDW() / d*OProduct(dx, -dx).GetSymPart();
			gradSum_r[ii] += 1.0 / (particles[jj].GetSi())*dwRbar;
			gradSum_r[jj] += 1.0 / (particles[ii].GetSi())*dwRbar;
		}
	}

	for (int i = 0; i<pNo; i++)
	{
		const short typei = max(particles[i].GetMaterial()->type, (const short)0);
		if (typei == m || typei == n)		// m < n
			particles[i].B = gradSum_r[i].Inverse();
	}
}
