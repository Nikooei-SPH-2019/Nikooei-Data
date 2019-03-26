#include "world.h"


void BodyForce(World &mw)
{
	
	int pNo = mw.particles.size();

	for (register int i = 0; i<pNo; i++)
	{


			mw.particles[i].SetAcc(mw.options.BODYACC);

	}

}



void PressureForceSymmetric(World &mw)                             //Nikooei
{

	int pNo = mw.particles.size();
	Vectorvec DelPOverRho(pNo, Vector());

	int iNo = mw.interactions.size();

	for (register int i = 0; i<iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();
		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d*dx;

		const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

		double RhoAsquare = 1 / (mw.particles[ii].GetRho()*mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho()*mw.particles[jj].GetRho());
						
		double h = mw.particles[ii].GetH();

		
		DelPOverRho[ii] = DelPOverRho[ii] + mw.particles[jj].GetM() *(mw.particles[ii].GetP()* RhoAsquare +
			mw.particles[jj].GetP()*RhoBsquare)*dw;
		DelPOverRho[jj] = DelPOverRho[jj] - mw.particles[ii].GetM() *(mw.particles[ii].GetP()* RhoAsquare +
			mw.particles[jj].GetP()*RhoBsquare)*dw;





	}
	for (int i = 0; i < pNo; i++)
		mw.particles[i].SetAcc(mw.particles[i].GetAcc() - DelPOverRho[i]);


}


