//****************************************************************************************
//
// SPH Code
//
// created	85/12 By Fatehi
// Last change	86/02 By Fatehi
// changed		89/03 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************

#include "world.h"
#include "tensors3.h"


Tensorvec EvalVelocityGradient(World& mw)			//Evaluates < grad V>
{
	const int iNo = mw.interactions.size();
	const int pNo = mw.particles.size();
	mw.Wallneighbor.assign(pNo, double());                      ////////////////Important

	Tensorvec DV(pNo, Tensor());

	for (register int i = 0; i<iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		if (typei<fluid1 && typej<fluid1)
			continue;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d*dx;


		const Vector dv = (mw.particles[jj].GetVel() - mw.particles[ii].GetVel());
		const Tensor dwdv = OProduct(dw, dv);

		const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction

		DV[ii] = DV[ii] + mw.particles[jj].GetM() / rhoij*dwdv;
		DV[jj] = DV[jj] + mw.particles[ii].GetM() / rhoij*dwdv;
		

		double hm = (mw.particles[jj].GetH() + mw.particles[jj].GetH()) / 2;

		if (mw.particles[ii].GetMaterial()->type >= fluid1 && mw.particles[jj].GetMaterial()->type < fluid1 && d < hm)
		{
			mw.Wallneighbor[ii] = mw.Wallneighbor[ii] + 1;
					}
		if (mw.particles[jj].GetMaterial()->type >= fluid1 && mw.particles[ii].GetMaterial()->type < fluid1 && d < hm)
		{
			mw.Wallneighbor[jj] = mw.Wallneighbor[jj] + 1;
		}

			
	}
		
	return DV;
}



void Granularmaterial1(World &mw)			//Pressure-dependent Rheology           
{
	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	double GammaDotMagi = 0, GammaDotMagj = 0, muAppij, muBingij = 0, TauMagi, TauMagj, YieldStress = 0, mueffj, TurbulentViscosity;
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor());;
	mw.DelV = EvalVelocityGradient(mw);
	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());
	mw.Walldistance.assign(pNo, double());
	double Phi_f = 30;
	double Phi_Bed = 10;

	double moving_mass_internalangle = Phi_f*0.01745329251994329576923690768489;
	double Substrate_internalangle = Phi_Bed*0.01745329251994329576923690768489;
	
	double Substrate_Mus = tan(Substrate_internalangle);
	double moving_mass_Mus = tan(moving_mass_internalangle);
	double  Substrate_YieldStress, moving_mass_YieldStress;
	

	double Substrate_BingViscosity = 0.1;   //Pa.s
	double moving_mass_BingViscosity = 0.1; //Pa.s
	double effectivP;
	double High_Viscosity_Factor = 10000;



	for (register int i = 0; i < pNo; i++)

	{
		mw.Walldistance[i] = 1;        //Initializing

	}

	for (int i = 0; i < pNo; i++)

	{
		effectivP = mw.particles[i].GetP();

		const short type = mw.particles[i].GetMaterial()->type;

		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    
		
		double mueffi = 0;
		double Betha = 0.000001;
		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();
		double Mus;		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double Tausecondinv = sqrt(mw.particles[i].GetTau().secondinvariant());
		
	if (type >= 15 || (type == 0 && mw.particles[i].Get_bed_neigh()>mw.particles[i].Get_overlaying_flow_neigh()))
		{
		
			Substrate_YieldStress = effectivP*Substrate_Mus;
			
			double Subst__Crit_GammaDot = Substrate_YieldStress / (High_Viscosity_Factor*Substrate_BingViscosity) + 1;


			if (mw.particles[i].GetGammaDot() > Subst__Crit_GammaDot)
			{
				mueffi = Substrate_BingViscosity + (Substrate_YieldStress - Substrate_BingViscosity*Subst__Crit_GammaDot) / (secondinvariant);

			}
			else
			{
				mueffi = High_Viscosity_Factor * Substrate_BingViscosity;
			}
		}
			

		if ((type >= fluid1 && type < 15) || (type == 0 && mw.particles[i].Get_bed_neigh()< mw.particles[i].Get_overlaying_flow_neigh()))
		{

			moving_mass_YieldStress = effectivP*moving_mass_Mus;
		
			double Mov_Mass_Crit_GammaDot = moving_mass_YieldStress / (High_Viscosity_Factor*moving_mass_BingViscosity);

			if (mw.particles[i].GetGammaDot() > Mov_Mass_Crit_GammaDot )
			{
				mueffi = moving_mass_BingViscosity + (moving_mass_YieldStress - moving_mass_BingViscosity*Mov_Mass_Crit_GammaDot) / (secondinvariant);

			}
			else
			{
				mueffi = (High_Viscosity_Factor*moving_mass_BingViscosity);
			}
		}
		
		double currStep = mw.GetSTEP();

		mw.particles[i].SetTau(mueffi*GammaDot[i]);
		mw.particles[i].SetMu_effective(mueffi);

	}
	std::vector<double> DelDotr_Erosion(pNo);


	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d*dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d*dx;
		double hm = (mw.particles[ii].GetH() + mw.particles[jj].GetH()) / 2;


		if (typei < fluid1 && typej < fluid1)
			continue;

		
		double RhoAsquare = 1 / (mw.particles[ii].GetRho()*mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho()*mw.particles[jj].GetRho());
	
		const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

			
	
		double DeltaX = mw.particles[ii].GetH() / 2.6;

			TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() *((dw)*(mw.particles[ii].GetTau()* RhoAsquare +
				mw.particles[jj].GetTau()*RhoBsquare));
			TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() *((dw)*(mw.particles[ii].GetTau()* RhoAsquare +
				mw.particles[jj].GetTau()*RhoBsquare));
		
			
			///////////////////////////////////////////


			if (mw.particles[ii].GetMaterial()->type >= fluid1 && mw.particles[jj].GetMaterial()->type < fluid1 && d < hm)
			{
				mw.Wallneighbor[ii] = mw.Wallneighbor[ii] + 1;
				mw.Walldistance[ii] = min(d, mw.Walldistance[ii]);

			}
			if (mw.particles[jj].GetMaterial()->type >= fluid1 && mw.particles[ii].GetMaterial()->type < fluid1 && d < hm)
			{
				mw.Wallneighbor[jj] = mw.Wallneighbor[jj] + 1;
				mw.Walldistance[jj] = min(d, mw.Walldistance[jj]);

			}


	}
	
	for (int i = 0; i < pNo; i++)
	{
		
		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]); 
			

	}


}









