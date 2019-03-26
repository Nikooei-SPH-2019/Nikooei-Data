#include "world.h"



void AddXSPH(World &mw)            /////Monaghan 1989
{
	const int pNo = mw.particles.size();
	Vectorvec dvel(pNo, Vector());


	const int iNo = mw.interactions.size();

	for (register int i = 0; i<iNo; i++)				//Compute Average deviation for velocity
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();
		int typei = mw.particles[ii].GetMaterial()->type;
		int typej = mw.particles[jj].GetMaterial()->type;
			const Vector vv = (mw.particles[jj].GetVel() - mw.particles[ii].GetVel())*mw.interactions[i].GetW();

		dvel[ii] += mw.particles[ii].GetM() / (mw.particles[ii].GetRho() + mw.particles[jj].GetRho())*vv;
		dvel[jj] -= mw.particles[jj].GetM() / (mw.particles[ii].GetRho() + mw.particles[jj].GetRho())*vv;
	}

	for (register int i = 0; i<pNo; i++)			//Add Average deviation to velocity
	{
		const short type = mw.particles[i].GetMaterial()->type;

		if ((mw.particles[i].GetMaterial()->type<fluid1)|| (type >= 15 && mw.particles[i].GetMu_effective() == 1000))
			continue;
	
		else if (mw.particles[i].GetMaterial()->type >= 15 && mw.Walldistance[i]<0.01)
				mw.particles[i].SetVel(mw.particles[i].GetVel() + 40 * mw.options.XSPH*dvel[i]);

		else

			mw.particles[i].SetVel(mw.particles[i].GetVel() + mw.options.XSPH*dvel[i]);

	}


}


void ShiftParticles3(World &mw)		           ///////Nikooei                 
{
	
	const int pNo = mw.particles.size();
	Vectorvec DP(pNo, Vector());                  ///////Nikooei
	doublevec neighboringwallparticle(pNo);       /////Nikooei-For Shifting
	doublevec rbar(pNo, 0.0);			//mean distance to neighbors
	intvec neigh(pNo, 0);				//Number of neighbors
	Vectorvec R(pNo, Vector());		//Shifting vector
	Tensorvec DV(pNo, Tensor());
	
	const int iNo = mw.interactions.size();
	for (register int i = 0; i < iNo; i++)				//Compute rbar
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		
		const double d = mw.interactions[i].GetDist();

		rbar[ii] += d;
		rbar[jj] += d;

	
			neigh[ii]++;
			neigh[jj]++;

		const Vector dx = mw.interactions[i].GetVD();
		const Vector dw = mw.interactions[i].GetDW() / d*dx;
		const Vector dv = (mw.particles[jj].GetVel() - mw.particles[ii].GetVel());

		const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());
		const double dp = (mw.particles[jj].GetP() - mw.particles[ii].GetP());

		const Vector dwdp = dw*dp;
		const Tensor dwdv = OProduct(dw, dv);


		//////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction

		DP[ii] = DP[ii] + mw.particles[jj].GetM() / rhoij*dwdp;
		DP[jj] = DP[jj] + mw.particles[ii].GetM() / rhoij*dwdp;

		DV[ii] = DV[ii] + mw.particles[jj].GetM() / rhoij*dwdv;
		DV[jj] = DV[jj] + mw.particles[ii].GetM() / rhoij*dwdv;

		}


	for (register int i = 0; i<pNo; i++)
	{
		rbar[i] = 1.0 / neigh[i] * rbar[i];
	}

	for (register int i = 0; i < iNo; i++)				//Compute rbar
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const double d = mw.interactions[i].GetDist();
		const Vector eij = 1.0 / d*mw.interactions[i].GetVD();

		//if (d == 0)
		//	continue;

		const Vector Rij = 1.0 / d / d / d*eij;

		R[ii] += rbar[ii] * rbar[ii] * rbar[ii] * Rij;
		R[jj] -= rbar[jj] * rbar[jj] * rbar[jj] * Rij;
	
	}
	

	 double alpha = mw.options.SHIFT*mw.GetDT();



	for (register int i = 0; i < pNo; i++)
	{
		int type = mw.particles[i].GetMaterial()->type;
		if (mw.particles[i].GetMaterial()->type < fluid1 || (type >= 15 && mw.particles[i].GetMu_effective() == 1000))
			continue;


			Vector dr;
		

			dr = alpha*mw.particles[i].GetVel().TwoNorm()*R[i];



			if (mw.particles[i].GetP() == 0 ) 
			{
				Vector normaltofreesurface = mw.gradSum[i] / mw.gradSum[i].TwoNorm();

				Vector Tangent;

				Tangent.SetX(normaltofreesurface.GetY());
				Tangent.SetY(-normaltofreesurface.GetX());

				double Costangentdr = Dot(dr, Tangent) / (Tangent.TwoNorm()*dr.TwoNorm());

				dr = (dr.TwoNorm())*Costangentdr*Tangent;
			
			}
			
			const Vector drnormal = dr / (dr.TwoNorm());


			if (dr.TwoNorm() > 0.5*mw.particles[i].GetH())                   /////For violant shifting distances (Stansby 2012)
			{
				dr.SetX(drnormal.GetX()*0.2*mw.particles[i].GetH());
				dr.SetY(drnormal.GetY()*0.2*mw.particles[i].GetH());
			}

			mw.particles[i].SetLoc(mw.particles[i].GetLoc() + dr);
			
			mw.particles[i].SetVel(mw.particles[i].GetVel() + dr*DV[i]);

			const double dp = Dot(dr, DP[i]);

			if(mw.particles[i].GetP()!=0)
			mw.particles[i].SetP(mw.particles[i].GetP() + dp);

		}

	}





