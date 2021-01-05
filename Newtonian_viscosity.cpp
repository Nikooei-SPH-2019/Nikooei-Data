//****************************************************************************************
//
// SPH Code
//
// created	85/12 By Fatehi
// Last change	86/02 By Fatehi
// changed		89/03 By Fatehi
//
//****************************************************************************************


#include "world.h"
#include "tensors3.h"


Tensorvec EvalMuVelocityGradient(World& mw)			//Evaluates <mu grad V>i2
{
	const int iNo = mw.interactions.size();
	const int pNo = mw.particles.size();

	Tensorvec DV(pNo, Tensor());

	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		if (typei < fluid1 && typej < fluid1)
			continue;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx;

		double muij;
		if (typei < fluid1)
			muij = mw.particles[jj].GetMaterial()->mu;
		else if (typej < fluid1)
			muij = mw.particles[ii].GetMaterial()->mu;
		else
		{
			const double mui = mw.particles[ii].GetMaterial()->mu;
			const double muj = mw.particles[jj].GetMaterial()->mu;
			muij = 2.0 / (1.0 / mui + 1.0 / muj);
		}

		const Vector dv = muij * (mw.particles[jj].GetVel() - mw.particles[ii].GetVel());
		const Tensor dwdv = OProduct(dw, dv);

		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction


		//DV[ii] = DV[ii] + 1.0 / mw.particles[jj].GetSi()*dwdv;
		//DV[jj] = DV[jj] + 1.0 / mw.particles[ii].GetSi()*dwdv;



		DV[ii] = DV[ii] + mw.particles[jj].GetM() / rhoij * dwdv;
		DV[jj] = DV[jj] + mw.particles[ii].GetM() / rhoij * dwdv;


	}

	for (register int i = 0; i < pNo; i++)
	{
		DV[i] = mw.particles[i].B * DV[i];
		//	mw.particles[i].SetTau(2.0*DV[i].GetSymPart());
	}

	return DV;
}

Tensorvec EvalVelocityGradient(World& mw)			//Evaluates < grad V>
{
	const int iNo = mw.interactions.size();
	const int pNo = mw.particles.size();
	doublevec Flow_Neighbor_volume(pNo, 0), Bed_Neighbor_volume(pNo, 0), Total_Neighbor_volume(pNo, 0),
		Flow_Neighbor_count(pNo, 0), Bed_Neighbor_count(pNo, 0), Total_Neighbor_count(pNo, 0);
	mw.gradSum.assign(pNo, Vector());
	mw.gradSum_Erosion.assign(pNo, Vector());
	mw.Walldistance.assign(pNo, double());
	mw.Wallneighbor.assign(pNo, double());                      ////////////////Important-Nikooei
//	mw.Moving_neighbor.assign(pNo, double());                      ////////////////Important-Nikooei
	Flow_Neighbor_volume.assign(pNo, double());                      ////////////////Important-Nikooei
	Bed_Neighbor_volume.assign(pNo, double());                      ////////////////Important-Nikooei
	Total_Neighbor_volume.assign(pNo, double());                      ////////////////Important-Nikooei
	intvec Bed_neighbor_of_wall(pNo, 0), Flow_neighbor_of_wall(pNo, 0) /*Interfaceneighbors(pNo, 0)*/;
	std::vector<double> DelDotr(pNo);
	std::vector<double> DelDotr_Erosion(pNo);


	Vector Normal_erosion;  //Normal Vector on Entrainment interface



	for (int i = 0; i < pNo; i++)

	{
		mw.Walldistance[i] = 1000;        //Initializing


		//Flow_neighbor_of_wall[i] = 0;
		//Flow_neighbor_of_wall[i] = 0;
	}



	Tensorvec DV(pNo, Tensor());

	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		if (typei < fluid1 && typej < fluid1)
			continue;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx;

		const Vector dv = (mw.particles[jj].GetVel() - mw.particles[ii].GetVel());
		const Tensor dwdv = OProduct(dw, dv);

		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());
		const STensor dwRbar = mw.interactions[i].GetDW() / d * OProduct(dx, -dx).GetSymPart();

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction


		//DV[ii] = DV[ii] + 1.0 / mw.particles[jj].GetSi()*dwdv;
		//DV[jj] = DV[jj] + 1.0 / mw.particles[ii].GetSi()*dwdv;



		DV[ii] = DV[ii] + mw.particles[jj].GetM() / rhoij * dwdv;
		DV[jj] = DV[jj] + mw.particles[ii].GetM() / rhoij * dwdv;

		
		DelDotr[ii] = DelDotr[ii] + mw.particles[jj].GetM() / rhoij * Dot(dx, dw);
	
		DelDotr[jj] = DelDotr[jj] + mw.particles[ii].GetM() / rhoij * Dot(dx, dw);


		if (typei == 15 && typej == 15 && mw.particles[ii].GetVel().TwoNorm() == 0 &&
			mw.particles[jj].GetVel().TwoNorm() == 0 && ii != jj && mw.particles[ii].GetP() != 0 && mw.particles[jj].GetP() != 0)
		{
			DelDotr_Erosion[ii] = DelDotr_Erosion[ii] + mw.particles[jj].GetM() / rhoij * Dot(dx, dw);
			DelDotr_Erosion[jj] = DelDotr_Erosion[jj] + mw.particles[ii].GetM() / rhoij * Dot(dx, dw);


			mw.gradSum_Erosion[ii] += mw.particles[jj].GetM() / rhoij * dw;
			mw.gradSum_Erosion[jj] -= mw.particles[ii].GetM() / rhoij * dw;




		}


		double hm = (mw.particles[jj].GetH() + mw.particles[jj].GetH()) / 2;

		if (typei >= fluid1 && typej < fluid1 && d < hm)
		{
			//	mw.Wallneighbor[ii] = mw.Wallneighbor[ii] + 1;
			mw.Walldistance[ii] = min(d, mw.Walldistance[ii]);

			//	mw.wallNormal[ii] = mw.wallNormal[ii] + mw.gradSum[jj] / mw.Wallneighbor[ii];

		}
		if (typej >= fluid1 && typei < fluid1 && d < hm)
		{
			//	mw.Wallneighbor[jj] = mw.Wallneighbor[jj] + 1;
			mw.Walldistance[jj] = min(d, mw.Walldistance[jj]);

			//	mw.wallNormal[jj] = mw.wallNormal[jj] + mw.gradSum[ii] / mw.Wallneighbor[jj];
		}







		///////////////////////////////////////////////////////////////////
#if 0
		if (typei < fluid1)
		{
			if (typej >= 15)
				Bed_neighbor_of_wall[ii] = Bed_neighbor_of_wall[ii] + 1;

			if (typej >= fluid1 && typej < 15)
				Flow_neighbor_of_wall[ii] = Flow_neighbor_of_wall[ii] + 1;

		}

		if (typej < fluid1)
		{
			if (typei >= 15)
				Bed_neighbor_of_wall[jj] = Bed_neighbor_of_wall[jj] + 1;

			if (typei >= fluid1 && typei < 15)
				Flow_neighbor_of_wall[jj] = Flow_neighbor_of_wall[jj] + 1;

		}



		if ((typei == 11 && typej == 15) || (typej == 11 && typei == 15))
		{
			Interfaceneighbors[ii] = Interfaceneighbors[ii] + 1;
			Interfaceneighbors[jj] = Interfaceneighbors[jj] + 1;

		}

#endif // 0

		////////////////////////////////////////////////////////////Volume_fraction


		if (typej == 15)
		{
			Bed_Neighbor_volume[ii] = Bed_Neighbor_volume[ii] + mw.particles[jj].GetM() / mw.particles[jj].GetRho();
			Bed_Neighbor_count[ii] = Bed_Neighbor_count[ii] + 1;

		}
		if (typej == 11)
		{
			Flow_Neighbor_volume[ii] = Flow_Neighbor_volume[ii] + mw.particles[jj].GetM() / mw.particles[jj].GetRho();
			Flow_Neighbor_count[ii] = Flow_Neighbor_count[ii] + 1;

		}
		if (typej >= fluid1)
		{
			Total_Neighbor_volume[ii] = Total_Neighbor_volume[ii] + mw.particles[jj].GetM() / mw.particles[jj].GetRho();
			Total_Neighbor_count[ii] = Total_Neighbor_count[ii] + 1;
		}
		if (typei == 15)
		{
			Bed_Neighbor_volume[jj] = Bed_Neighbor_volume[jj] + mw.particles[ii].GetM() / mw.particles[ii].GetRho();
			Bed_Neighbor_count[jj] = Bed_Neighbor_count[jj] + 1;

		}
		if (typei == 11)
		{
			Flow_Neighbor_volume[jj] = Flow_Neighbor_volume[jj] + mw.particles[ii].GetM() / mw.particles[ii].GetRho();
			Flow_Neighbor_count[jj] = Flow_Neighbor_count[jj] + 1;

		}
		if (typei >= fluid1)
		{
			Total_Neighbor_volume[jj] = Total_Neighbor_volume[jj] + mw.particles[ii].GetM() / mw.particles[ii].GetRho();
			Total_Neighbor_count[jj] = Total_Neighbor_count[jj] + 1;

		}
		////////////////////////////////////////////////////////////Volume_fraction

		if (typei < fluid1)
		{
			if (typej >= 15)
				Bed_neighbor_of_wall[ii] = Bed_neighbor_of_wall[ii] + 1;

			if (typej >= fluid1 && typej < 15)
				Flow_neighbor_of_wall[ii] = Flow_neighbor_of_wall[ii] + 1;

		}

		if (typej < fluid1)
		{
			if (typei >= 15)
				Bed_neighbor_of_wall[jj] = Bed_neighbor_of_wall[jj] + 1;

			if (typei >= fluid1 && typei < 15)
				Flow_neighbor_of_wall[jj] = Flow_neighbor_of_wall[jj] + 1;

		}


#if 0
		if (/*particles[ii].GetMu_effective()==1000*/mw.particles[ii].GetVel().TwoNorm() < 0.07 * sqrt(9.81 * 1.32) && mw.particles[jj].GetVel().TwoNorm() > 0.07 * sqrt(9.81 * 0.14)/*particles[jj].GetMu_effective() != 1000*/)
			mw.Moving_neighbor[ii] = mw.Moving_neighbor[ii] + 1;

		if (/*particles[jj].GetMu_effective() == 1000 */  mw.particles[jj].GetVel().TwoNorm() < 0.07 * sqrt(9.81 * 1.32) && mw.particles[ii].GetVel().TwoNorm() > 0.07 * sqrt(9.81 * 0.14) /*particles[ii].GetMu_effective() != 1000*/)
			mw.Moving_neighbor[jj] = mw.Moving_neighbor[jj] + 1;

#endif // 0

		mw.gradSum[ii] += mw.particles[jj].GetM() / rhoij * dw;
		mw.gradSum[jj] -= mw.particles[ii].GetM() / rhoij * dw;

	}
	//////////////////////////////////////////////////
	for (int i = 0; i < pNo; i++)
	{
		mw.particles[i].Set_bed_neigh(Bed_neighbor_of_wall[i]);
		mw.particles[i].Set_overlaying_flow_neigh(Flow_neighbor_of_wall[i]);
		//mw.particles[i].SetInterfaceneighbor(Interfaceneighbors[i]);
	//	mw.particles[i].SetMoving_neighbor(mw.Moving_neighbor[i]);

		double denumnormal = sqrt(mw.gradSum_Erosion[i].GetX() * mw.gradSum_Erosion[i].GetX() + mw.gradSum_Erosion[i].GetY() * mw.gradSum_Erosion[i].GetY());
		Normal_erosion.SetX(mw.gradSum_Erosion[i].GetX() / denumnormal);
		Normal_erosion.SetY(mw.gradSum_Erosion[i].GetY() / denumnormal);

		if (abs(mw.gradSum_Erosion[i].TwoNorm()) > 100000 || denumnormal == 0)
		{
			Normal_erosion.SetX(0);
			Normal_erosion.SetY(0);

		}
		mw.particles[i].SetNormal_erosion(Normal_erosion);

		mw.particles[i].Set_Volume_fraction_flow(Flow_Neighbor_volume[i] / Total_Neighbor_volume[i]);
		mw.particles[i].Set_Volume_fraction_bed(Bed_Neighbor_volume[i] / Total_Neighbor_volume[i]);
		mw.particles[i].SetDelDotr_Erosion(DelDotr_Erosion[i]);
		if (mw.GetSTEP() == 1)
			mw.particles[i].SetDelDotr_Erosion(2);

		mw.particles[i].SetDelDotr(DelDotr[i]);
		mw.particles[i].SetFluidneighbor(Flow_Neighbor_count[i]);
		mw.particles[i].Set_bed_neigh(Bed_Neighbor_count[i]);
		mw.particles[i].Set_Total_neigh(Total_Neighbor_count[i]);


	}

	   
	/*
	for (register int i = 0; i<pNo; i++)
	{
	DV[i] = mw.particles[i].B*DV[i];
	//	mw.particles[i].SetTau(2.0*DV[i].GetSymPart());
	}
	*/
	return DV;
}


void NewtonianViscosity(World & mw)			//Nikooei
{

	int pNo = mw.particles.size();

	int iNo = mw.interactions.size();
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor());;
	mw.muDV = EvalMuVelocityGradient(mw);                                        //Not suitable for Non-Newtonian
	Vectorvec TaudivergenceoverRho(pNo, Vector());
	Tensorvec muDelV(pNo, Tensor());
	Vectorvec nuDel2U(pNo, Vector());



	for (register int i = 0; i < iNo; i++)
	{
		/*
		const int ii=mw.interactions[i].GetI();
		const int jj=mw.interactions[i].GetJ();

		const short typei=mw.particles[ii].GetMaterial()->type;
		const short typej=mw.particles[jj].GetMaterial()->type;

		if (typei<fluid1 && typej<fluid1)
		continue;

		const Vector dx=mw.interactions[i].GetVD();
		const double d=mw.interactions[i].GetDist();
		const double dw=mw.interactions[i].GetDW();
		const double hm=mw.interactions[i].GetHMean();

		const Vector dv=(mw.particles[ii].GetVel() - mw.particles[jj].GetVel());

		double muij;
		if (typei<fluid1)
		muij=mw.particles[jj].GetMaterial()->mu;
		else if (typej<fluid1)
		muij=mw.particles[ii].GetMaterial()->mu;
		else
		{
		const double mui=mw.particles[ii].GetMaterial()->mu;
		const double muj=mw.particles[jj].GetMaterial()->mu;
		muij=2.0/(1.0/mui+1.0/muj);
		}

		const Vector force=2.0*muij*dw/(d+1e-4*hm) /(mw.particles[ii].GetSi()*mw.particles[jj].GetSi())* dv ;		//A force in the direction of velocity
		//difference and proportinal to inverse of distance

		mw.particles[ii].SetAcc(mw.particles[ii].GetAcc()+ 1.0/mw.particles[ii].GetM()*force);		//Action
		mw.particles[jj].SetAcc(mw.particles[jj].GetAcc()- 1.0/mw.particles[jj].GetM()*force);		//Reaction
		*/

		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		if (typei < fluid1 && typej < fluid1)
			continue;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx;

		/*


		double muij;
		if (typei < fluid1)
		muij = mw.particles[jj].GetMaterial()->mu;
		else if (typej < fluid1)
		muij = mw.particles[ii].GetMaterial()->mu;
		else
		{
		const double mui = mw.particles[ii].GetMaterial()->mu;
		const double muj = mw.particles[jj].GetMaterial()->mu;
		muij = 2.0 / (1.0 / mui + 1.0 / muj);
		}
		const Vector mudv = muij*(mw.particles[jj].GetVel() - mw.particles[ii].GetVel());
		const Tensor muDelwdv = OProduct(dw, mudv);
		const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

		//////////////////////////////////Nikooei-omega correction

		muDelV[ii] = muDelV[ii] + mw.particles[jj].GetM() / rhoij*muDelwdv;          // m/rho=omega
		muDelV[jj] = muDelV[jj] + mw.particles[ii].GetM() / rhoij*muDelwdv;






		mw.particles[ii].SetTau(2 * muDelV[ii].GetSymPart());

		mw.particles[jj].SetTau(2 * muDelV[jj].GetSymPart());

		const STensor tauij = 1.0 / rhoij* (mw.particles[jj].GetTau() - mw.particles[ii].GetTau());
		//	TaudivergenceoverRho[ii] += 1.0 / mw.particles[jj].GetSi()*((mw.particles[ii].B*dw)*tauij);
		//TaudivergenceoverRho[jj] += 1.0 / mw.particles[ii].GetSi()*((mw.particles[jj].B*dw)*tauij);           //Doubt    Negative or positive

		TaudivergenceoverRho[ii] += mw.particles[jj].GetM() / rhoij*((mw.particles[ii].B*dw)*tauij);
		TaudivergenceoverRho[jj] += mw.particles[ii].GetM() / rhoij*((mw.particles[jj].B*dw)*tauij);           //Doubt    Negative or positive
		*/

		////////////////////////////////////////////////////////////Bui 2016 ,....
		const Vector dv = (mw.particles[ii].GetVel() - mw.particles[jj].GetVel());
		double h = mw.particles[ii].GetH();
		double nu = 7.9554494828957836117740652346858e-4;
		//double nu = 1e-6;   //mu=0.01
		double Etha = 0.001 * h;   //////////////Bui 2016


								 //double Etha = 0.1*h;
								 //Viscous flows
								 /*   /////////////////////////////////////////////////////         Morris 1997
								 nuDel2U[ii] = nuDel2U[ii] + 2 * mw.particles[jj].GetM() / mw.particles[ii].GetRho()*nu*(Dot(dx, dw)*(dv) / (d*d + Etha*Etha));
								 nuDel2U[jj] = nuDel2U[jj] - 2 * mw.particles[ii].GetM() / mw.particles[ii].GetRho()*nu*(Dot(dx, dw)*(dv) / (d*d + Etha*Etha));
								 */

								 ////////////////////////////////////////////////////////////Monaghan & Gingold 1983
		nuDel2U[ii] = nuDel2U[ii] + 8 * mw.particles[jj].GetM() / mw.particles[ii].GetRho() * nu * (Dot(dx, dv) * (dw) / (d * d + Etha * Etha));
		nuDel2U[jj] = nuDel2U[jj] - 8 * mw.particles[ii].GetM() / mw.particles[ii].GetRho() * nu * (Dot(dx, dv) * (dw) / (d * d + Etha * Etha));

		//////////////////////Barcarolo Thesis mentioned two above procedures  The second operator  does not suffer from the issue mentioned by Colagrossi, hence, whenever
		//	dealing with viscous free - surface flow such an operator is preferable.

	}


	for (register int i = 0; i < pNo; i++)
	{

		/*
		double DeltaX = mw.particles[i].GetH() / 2.6;
		GammaDot[i] = mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)
		double mueffi;
		double secondinvariant = GammaDot[i].secondinvariant();   //GammaDot[i].secondinvariant();

		double GammaDotMagi = 2 * sqrt(secondinvariant);


		STensor Turbulency;
		double Cs = 0.1;
		double TurbulentViscosity = (Cs*DeltaX)*(Cs*DeltaX)*GammaDotMagi;

		double Kinetic = pow(TurbulentViscosity / (0.08*DeltaX), 2);
		Turbulency = 2 * mw.particles[i].GetRho()*TurbulentViscosity*GammaDot[i];

		Turbulency.SetXX(Turbulency.GetXX() - (2 / 3)* mw.particles[i].GetRho()* Kinetic);
		Turbulency.SetYY(Turbulency.GetYY() - (2 / 3) * mw.particles[i].GetRho()* Kinetic);


		mw.particles[i].SetTurbulTau(Turbulency);
		*/
		//mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TaudivergenceoverRho[i]);

		//if (mw.particles[i].GetP() == 0 && mw.particles[i].GetFluidneighbor() < 5)   ///////////////////For free particles which are in leading viscous force not computed
		//	continue;

		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + nuDel2U[i] /*+ TaudivergenceoverRho[i]*/);


		//DV[i] = mw.particles[i].B*DV[i];
	}



}





/////////////////////////////////////////////////////////////////Nikooe-Nonnewtonian

void Bingham(World & mw)

{

	double Cs = 0.1;
	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	double GammaDotMagi = 0, GammaDotMagj, muAppij, muBingij = 0, TauMagi, TauMagj, YieldStress = 25, mueffi, mueffj, TurbulentViscosity;
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor());;
	mw.DelV = EvalVelocityGradient(mw);

	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());


	/*

	for (register int i = 0; i < iNo; i++)
	{
	const int ii = mw.interactions[i].GetI(); //creates index
	//interactions is a matrix
	const int jj = mw.interactions[i].GetJ();
	const short typei = mw.particles[ii].GetMaterial()->type; // az andise bala estefade mikonad
	const short typej = mw.particles[jj].GetMaterial()->type;



	if (typei < fluid1 && typej < fluid1)
	continue;
	const Vector dx = mw.interactions[i].GetVD();
	const double d = mw.interactions[i].GetDist();
	const Vector dw = mw.interactions[i].GetDW() / d*dx; // Del W= W_ij/|r_ij|*r_ij
	const Vector eij = 1.0 / d*dx;

	///Assume one type fluid ( bingham) is adjacent to wall

	double mubingij;
	if (typei < fluid1) // i particles are boundary
	mubingij = mw.particles[jj].GetMaterial()->mu;
	else if (typej < fluid1) // j particles are boundary
	mubingij = mw.particles[ii].GetMaterial()->mu;
	else // both of i and j particles are Nonnewtonfluid
	{
	const double mubingi = mw.particles[ii].GetMaterial()->mu;
	const double mubingj = mw.particles[jj].GetMaterial()->mu;
	mubingij = 2.0 / (1.0 / mubingi + 1.0 / mubingj);
	}
	/////////////////////////////////////////////////////////////////Modified Bingham Model(Bi-Viscosity)-Tanner and Milthorpe 1983  - E.Mitsoulis 2007 -  Dent & Lang 1983

	const Vector dv = (mw.particles[jj].GetVel() - mw.particles[ii].GetVel());             //for first iteration??????????????????????
	const Tensor dwdv = OProduct(dw, dv); //dwdv=Del W *dV= W_ij/|r_ij|*r_ij *dV
	//const STensor dwe = OProduct(dw, eij).GetSymPart();

	///////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction

	const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

	//	DV[ii] = DV[ii] + 1.0 / mw.particles[jj].GetSi()*dwdv;

	//	DV[jj] = DV[jj] + 1.0 / mw.particles[ii].GetSi()*dwdv;

	DV[ii] = DV[ii] + mw.particles[jj].GetM() / rhoij*dwdv;
	DV[jj] = DV[jj] + mw.particles[ii].GetM() / rhoij*dwdv;



	STensor a = DV[ii].GetSymPart();
	double b = DV[ii].GetSymPart().secondinvariant();

	Gammadot[ii] = Gammadot[ii] + DV[ii].GetSymPart();
	Gammadot[jj] = Gammadot[jj] + DV[jj].GetSymPart();


	//Gammadotmagj = sqrt((0.5)*(DV[jj].GetSymPart().secondinvariant()));

	//mw.particles[ii].SetTau(mw.particles[ii].GetTau()+mubingij*Gammadot[ii]);
	//mw.particles[jj].SetTau(mw.particles[jj].GetTau() + mubingij*Gammadot[jj]);


	}

	*/


	for (int i = 0; i < pNo; i++)
	{
		double DeltaX = 0.01;
		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();
		double muBingij = 0.07;
		double a = muBingij / YieldStress;
		GammaDotMagi = sqrt((0.5) * (GammaDot[i].secondinvariant()));
		double b = a * GammaDotMagi;
		mueffi = (1000 * muBingij) * ((1 + b) / (1 + 1000 * b));
		mw.particles[i].SetTau(mueffi * GammaDot[i]);
		TauMagi = sqrt((0.5) * (mw.particles[i].GetTau().secondinvariant()));       //     |τ|         second invarient of stress tensor

		TurbulentViscosity = (Cs * DeltaX) * (Cs * DeltaX) * GammaDotMagi;
		mw.particles[i].SetTurbulTau(2 * mw.particles[i].GetRho() * TurbulentViscosity * GammaDot[i]);




	}

	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;


		if (typei < fluid1 && typej < fluid1)
			continue;


		/*
		const Vector dv = (mw.particles[jj].GetVel() - mw.particles[ii].GetVel());
		const Tensor dwdv = OProduct(dw, dv);

		const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());


		DV[ii] = DV[ii] + mw.particles[jj].GetM() / rhoij*dwdv;
		DV[jj] = DV[jj] + mw.particles[ii].GetM() / rhoij*dwdv;

		*/
		/*
		GammaDot[ii] = 2 * mw.DelV[ii].GetSymPart();          //  γdot= ∇u+transpose∇u=2D=2 (∇u+transpose∇u)/2               DV[ii].GetSymPart()=D
		GammaDot[jj] = 2 * DV[jj].GetSymPart();

		GammaDotMagi = sqrt((0.5)*(GammaDot[ii].secondinvariant()));
		GammaDotMagj = sqrt((0.5)*(GammaDot[jj].secondinvariant()));


		double muBingij = 0.07;
		double a = muBingij / YieldStress;
		double b = a*GammaDotMagi;
		double c = a*GammaDotMagi;

		mueffi = (1000 * muBingij)*((1 + b) / (1 + 1000 * b));
		mueffj = (1000 * muBingij)*((1 + c) / (1 + 1000 * c));


		mw.particles[ii].SetTau(mueffi*GammaDot[ii]);
		mw.particles[jj].SetTau(mueffj*GammaDot[jj]);


		TauMagi = sqrt((0.5)*(mw.particles[ii].GetTau().secondinvariant()));       //     |τ|         second invarient of stress tensor
		TauMagj = sqrt((0.5)*(mw.particles[jj].GetTau().secondinvariant()));       //     |τ|         second invarient of stress tensor

		*/

		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());


		//	const STensor tauij = 1.0 / rhoij* (mw.particles[jj].GetTau() - mw.particles[ii].GetTau());
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction
		//TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() / rhoij*((mw.particles[ii].B*dw)*tauij);
		//	TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] + mw.particles[ii].GetM() / rhoij*((mw.particles[jj].B*dw)*tauij);
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((mw.particles[ii].B * dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((mw.particles[jj].B * dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));


		TurbulTauDivergenceOverRho[ii] = TurbulTauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((mw.particles[ii].B * dw) * (mw.particles[ii].GetTurbulTau() * RhoAsquare +
			mw.particles[jj].GetTurbulTau() * RhoBsquare));

		TurbulTauDivergenceOverRho[jj] = TurbulTauDivergenceOverRho[jj] + mw.particles[ii].GetM() * ((mw.particles[jj].B * dw) * (mw.particles[ii].GetTurbulTau() * RhoAsquare +
			mw.particles[jj].GetTurbulTau() * RhoBsquare));





	}



	for (int i = 0; i < pNo; i++)
		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i] + TurbulTauDivergenceOverRho[i]);
}
/*
void ModifiedBingham1phase(World &mw)			//coulomb yield stress
{
int pNo = mw.particles.size();
int iNo = mw.interactions.size();
double gammadotmagi, gammadotmagj, muappij, mubingij = 0, taumag, yieldstress=0;
STensor Tau;
Tensorvec muDelV(pNo, Tensor());

//mw.DelV = EvalVelocityGradient(mw);
Tensorvec DV(pNo, Tensor());


Vectorvec TaudivergenceoverRho(pNo, Vector());

for (register int i = 0; i < iNo; i++)
{
const int ii = mw.interactions[i].GetI(); //creates index
//interactions is a matrix
const int jj = mw.interactions[i].GetJ();
const short typei = mw.particles[ii].GetMaterial()->type; // az andise bala estefade mikonad
const short typej = mw.particles[jj].GetMaterial()->type;

if (typei < fluid1 && typej < fluid1)
continue;
const Vector dx = mw.interactions[i].GetVD();
const double d = mw.interactions[i].GetDist();
const Vector dw = mw.interactions[i].GetDW() / d*dx; // Del W= W_ij/|r_ij|*r_ij
const Vector eij = 1.0 / d*dx;

///Assume one type fluid ( bingham) is adjacent to wall

double mubingij;
if (typei < fluid1) // i particles are boundary
mubingij = mw.particles[jj].GetMaterial()->mu;
else if (typej < fluid1) // j particles are boundary
mubingij = mw.particles[ii].GetMaterial()->mu;
else // both of i and j particles are Nonnewtonfluid
{
const double mubingi = mw.particles[ii].GetMaterial()->mu;
const double mubingj = mw.particles[jj].GetMaterial()->mu;
mubingij = 2.0 / (1.0 / mubingi + 1.0 / mubingj);
}
/////////////////////////////////////////////////////////////////Modified Bingham Model(Bi-Viscosity)-Tanner and Milthorpe 1983  - E.Mitsoulis 2007 -  Dent & Lang 1983

const Vector dv = (mw.particles[jj].GetVel() - mw.particles[ii].GetVel());
const Tensor dwdv = OProduct(dw, dv); //dwdv=Del W *dV= W_ij/|r_ij|*r_ij *dV
//const STensor dwe = OProduct(dw, eij).GetSymPart();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction

const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

//	DV[ii] = DV[ii] + 1.0 / mw.particles[jj].GetSi()*dwdv;

//	DV[jj] = DV[jj] + 1.0 / mw.particles[ii].GetSi()*dwdv;

DV[ii] = DV[ii] + mw.particles[jj].GetM() / rhoij*dwdv;
DV[jj] = DV[jj] + mw.particles[ii].GetM() / rhoij*dwdv;



STensor a = DV[ii].GetSymPart();
double b = DV[ii].GetSymPart().secondinvariant();

gammadotmagi = sqrt((0.5)*(DV[ii].GetSymPart().secondinvariant()));
gammadotmagj = sqrt((0.5)*(DV[jj].GetSymPart().secondinvariant()));



double criticalstrainrate = 0.3571;

double mueffective = 0.0;

double mu0 = 1000 * mubingij;
//double yieldstress = mu0*criticalstrainrate;


if (mw.particles[ii].GetMaterial()->type == fluid1)
double yieldstress = 0.477*mw.particles[ii].GetP();


if (mw.particles[ii].GetMaterial()->type == fluid1)
double yieldstress = 0.477*mw.particles[ii].GetP();


///Stress of ii
if (gammadotmagi <= criticalstrainrate)
{
mueffective = mu0;
mw.particles[ii].SetTau(mueffective*DV[ii].GetSymPart());

}
else if (gammadotmagi > criticalstrainrate)
{
mw.particles[ii].SetTau((DV[ii].GetSymPart() - criticalstrainrate)*mubingij + yieldstress);
}



/////////stress of jj
if (gammadotmagj <= criticalstrainrate)
{
mueffective = mu0;
mw.particles[jj].SetTau(mueffective*DV[jj].GetSymPart());

}
else if (gammadotmagj > criticalstrainrate)
{
mw.particles[jj].SetTau((DV[jj].GetSymPart() - criticalstrainrate)*mubingij + yieldstress);
}




const STensor tauij = 1.0 / rhoij* (mw.particles[jj].GetTau() - mw.particles[ii].GetTau());
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction
TaudivergenceoverRho[ii] += mw.particles[jj].GetM() / rhoij*((mw.particles[ii].B*dw)*tauij);
TaudivergenceoverRho[jj] -= mw.particles[ii].GetM() / rhoij*((mw.particles[jj].B*dw)*tauij);



}

for (int i = 0; i<pNo; i++)
mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TaudivergenceoverRho[i]);
}
*/

void ModifiedBingham2phase(World & mw)			//Nonnewtonian Rheology
{

	double Cs = 0.1;
	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	double GammaDotMagi = 0, GammaDotMagj = 0, muAppij, muBingij = 0, TauMagi, TauMagj, YieldStress = 0, mueffj, TurbulentViscosity;
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor());;
	mw.DelV = EvalVelocityGradient(mw);

	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());


	for (int i = 0; i < pNo; i++)

	{
		double DeltaX = mw.particles[i].GetH() / 3;
		GammaDot[i] = mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    
		double mueffi;
		double secondinvariant = GammaDot[i].secondinvariant();   //GammaDot[i].secondinvariant();

		GammaDotMagi = 2 * sqrt(secondinvariant);


		if (mw.particles[i].GetMaterial()->type == fluid1)
		{

			mueffi = mw.particles[i].GetMaterial()->mu;

		}


		if ((mw.particles[i].GetMaterial()->type == BinghamCoulomb) || (mw.particles[i].GetMaterial()->type < fluid1))      //Scouring dam break
		{
			double Cohesion = 0.0;                   // Fraccarollo -PVC
			double internalangle = 0.471; //27 Degree        // Fraccarollo -PVC
										  //	double muBingij = mw.particles[i].GetMaterial()->mu;

			double muBingij = 0.071;       //Ionescu (2015)



										   //YieldStress = 23.5;
			YieldStress = 60;
			/////////////////////////////////////////////////////////Cross
			/*
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////artificialviscosity
			double artificialviscosity = 0.2*(mw.particles[i].GetVel().TwoNorm())*(mw.particles[i].GetH() / 3);

			double modifiedmubingij = muBingij + artificialviscosity;

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			*/

			double modifiedmubingij = muBingij;


			double b = modifiedmubingij * GammaDotMagi;

			if (mw.GetSTEP() != 1)
				mueffi = (1000 * modifiedmubingij) * ((YieldStress + b) / (YieldStress + 1000 * b));
			else
				mueffi = 0;




			/*

			double muBingij = 0.071;       //Ionescu (2015)
			double mu0 = 10;
			double Theta = 0.92;  //////Xie 2014  Cv=15%
			double K = 8.4;     //////Xie 2014  Cv=15%



			double KGammadot = K*GammaDotMagi;

			double KGammadotPowTheta = pow(KGammadot, Theta);
			mueffi = (mu0 + 0.001*mu0*KGammadotPowTheta) / (1 + KGammadotPowTheta);

			*/



















			//////////////////////////////////////////////Fatehi

			/*

			if (GammaDotMagi >= YieldStress / (1000 * modifiedmubingij))
			mueffi = modifiedmubingij + YieldStress / (GammaDotMagi);

			else if (GammaDotMagi < YieldStress / (1000 * modifiedmubingij))
			mueffi = 1001*modifiedmubingij;

			*/



			//mueffi = ((YieldStress + b) / (1 + GammaDotMagi));
			//mueffi = 1+(YieldStress / (0.001+GammaDotMagi));

			//	if (mueffi > 80)
			//	mueffi = 120;


			mw.particles[i].SetTau(2 * mueffi * GammaDot[i]);
			/*
			//////////////////////////////////////////////////////////////////////////Turbulency

			STensor Turbulency;
			Cs = 0.1;
			TurbulentViscosity = (Cs*DeltaX)*(Cs*DeltaX)*GammaDotMagi;

			double Kinetic = pow(TurbulentViscosity / (0.08*DeltaX), 2);
			Turbulency = 2 * mw.particles[i].GetRho()*TurbulentViscosity*GammaDot[i];

			Turbulency.SetXX(Turbulency.GetXX() - (2 / 3)* mw.particles[i].GetRho()* Kinetic);
			Turbulency.SetYY(Turbulency.GetYY() - (2 / 3) * mw.particles[i].GetRho()* Kinetic);


			mw.particles[i].SetTurbulTau(Turbulency);
			///////////////////////////////////////////////////////////////////////////////////////
			*/
		}

		//////////////////////////////////Bingham-Bi viscous
		/*

		double GammaCrit = YieldStress / (1000 * muBingij);
		if (GammaDotMagi < GammaCrit)
		mw.particles[i].SetTau((1000 * muBingij)*GammaDot[i]);
		else
		mw.particles[i].SetTau((muBingij+(YieldStress- muBingij*GammaCrit)/ GammaDotMagi)*GammaDot[i]);
		*/

		/////////////


		////////////////////////////////////////////////////////////////////////////////// //Papanastasiou (1987), Zho 2010
		/*
		double m=0.1;
		if (mw.GetSTEP() != 1)
		mueffi = (muBingij + YieldStress*(1 - exp(-m*GammaDotMagi)) / GammaDotMagi);
		else
		mueffi = 0;
		*/
		///////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////// //Herschel-Bulkley model
		/*
		double n = 0;
		double K =0 ;
		if (mw.GetSTEP() != 1)            ///If Tau=psin(Phi) it should be applied
		mueffi = K* pow(GammaDotMagi,(n - 1)) + (YieldStress / GammaDotMagi);
		else
		mw.particles[i].SetTau(0);
		*/
		////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////Casson model
		/*
		if (mw.GetSTEP() != 1)
		mueffi = (sqrt(muBingij) + sqrt(YieldStress / GammaDotMagi))* (sqrt(muBingij) + sqrt(YieldStress / GammaDotMagi));
		else
		mw.particles[i].SetTau(0);

		///////////////////////////////////
		double currStep = mw.GetSTEP();
		//if (mw.GetSTEP() != 1)            ///If Tau=psin(Phi) it should be applied
		mw.particles[i].SetTau(2 * mueffi*GammaDot[i]);
		//else
		//	mw.particles[i].SetTau(0);

		}


		/*
		////////////////////////////////////////////Bi-viscous
		double r = sin(internalangle);
		YieldStress = Cohesion*cos(internalangle) + mw.particles[i].GetP()*sin(internalangle);
		double criticalstrainrate = YieldStress / muBingij;
		if (GammaDotMagi >= criticalstrainrate)
		mueffi = muBingij + ((YieldStress - muBingij*criticalstrainrate) / GammaDotMagi);
		if (GammaDotMagi < criticalstrainrate)
		mueffi =1000* muBingij ;

		/////////////////////////////////////////
		}

		*/




		///////////////////////////////////////////////////////////////////////////////



	}

	doublevec neighwallcounter;
	mw.Wallneighbor.assign(pNo, double());                      ////////////////Important
	mw.wallNormal.assign(pNo, Vector());
	mw.Walldistance.assign(pNo, double());



	for (register int i = 0; i < pNo; i++)
		mw.Walldistance[i] = 1;        //Initializing



	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;


		if (typei < fluid1 && typej < fluid1)
			continue;
		/*
		if (typei < fluid1)
		{

		TauMagj = sqrt((0.5)*(mw.particles[jj].GetTau().secondinvariant()));       //
		GammaDotMagj = sqrt((0.5)*(GammaDot[jj].secondinvariant()));



		mw.particles[ii].SetTau((TauMagj/ GammaDotMagj)*GammaDot[ii]);
		}



		if (typej < fluid1)
		{
		TauMagi = sqrt((0.5)*(mw.particles[ii].GetTau().secondinvariant()));       //
		GammaDotMagi = sqrt((0.5)*(GammaDot[ii].secondinvariant()));

		mw.particles[jj].SetTau((TauMagi / GammaDotMagi)*GammaDot[jj]);
		}

		*/




		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());

		/*
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() *((mw.particles[ii].B*dw)*(mw.particles[ii].GetTau()* RhoAsquare +
		mw.particles[jj].GetTau()*RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() *((mw.particles[jj].B*dw)*(mw.particles[ii].GetTau()* RhoAsquare +
		mw.particles[jj].GetTau()*RhoBsquare));

		*/


		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * (mw.particles[ii].GetTau() * RhoAsquare + mw.particles[jj].GetTau() * RhoBsquare) * dw;
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * (mw.particles[ii].GetTau() * RhoAsquare + mw.particles[jj].GetTau() * RhoBsquare) * dw;







		/*
		TurbulTauDivergenceOverRho[ii] = TurbulTauDivergenceOverRho[ii] + mw.particles[jj].GetM() *((mw.particles[ii].B*dw)*(mw.particles[ii].GetTurbulTau()* RhoAsquare +
		mw.particles[jj].GetTurbulTau()*RhoBsquare));

		TurbulTauDivergenceOverRho[jj] = TurbulTauDivergenceOverRho[jj] + mw.particles[ii].GetM() *((mw.particles[jj].B*dw)*(mw.particles[ii].GetTurbulTau()* RhoAsquare +
		mw.particles[jj].GetTurbulTau()*RhoBsquare));

		*/



		/*
		const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());
		const STensor tauij = 1.0 / rhoij* (mw.particles[jj].GetTau() - mw.particles[ii].GetTau());
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() / rhoij*((mw.particles[ii].B*dw)*tauij);
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] + mw.particles[ii].GetM() / rhoij*((mw.particles[jj].B*dw)*tauij);
		*/
		/*
		TurbulTauDivergenceOverRho[ii] = TurbulTauDivergenceOverRho[ii] + mw.particles[jj].GetM() *((mw.particles[ii].B*dw)*(mw.particles[ii].GetTurbulTau()* RhoAsquare +
		mw.particles[jj].GetTurbulTau()*RhoBsquare));

		TurbulTauDivergenceOverRho[jj] = TurbulTauDivergenceOverRho[jj] + mw.particles[ii].GetM() *((mw.particles[jj].B*dw)*(mw.particles[ii].GetTurbulTau()* RhoAsquare +
		mw.particles[jj].GetTurbulTau()*RhoBsquare));
		*/

		double hm = (mw.particles[ii].GetH() + mw.particles[jj].GetH()) / 2;

		/////////////

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
		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]);  //+ TurbulTauDivergenceOverRho[i]);

																					 //mw.Walldistance[i] = 1;                            //For finding Normal distance to wall
	}

}


//////////////////////////////////////////////pressure dependent Rheologies




void Granularmaterial_Interface_continuity(World & mw)			//Pressure-dependent Rheology           constant viscosity model (ETHA=cst)            Ionescu 2015   
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
	double Phi_f = 25.5;
	double Phi_Bed = 10;

	double moving_mass_internalangle = Phi_f * 0.01745329251994329576923690768489;
	double Substrate_internalangle = Phi_Bed * 0.01745329251994329576923690768489;


	double Wall_Friction_angle = 0.34906585039886591538473815369772;    //20   Crosta

	double Substrate_Mus = tan(Substrate_internalangle);
	double moving_mass_Mus = tan(moving_mass_internalangle);
	double Rhos, d, Substrate_YieldStress, moving_mass_YieldStress;
	double InertialNum;
	double Substrate_Critical_inertialnumber;
	double moving_mass_Critical_inertialnumber;

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


		double mueffi, MuI;

		double Betha = 0.000001;


		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();


																				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Mu(I) Rheology

																				//	double InertialNum =0.001*(2*sqrt(secondinvariant))/(sqrt(abs(mw.particles[i].GetP())/mw.particles[i].GetRho()));
																				//	MuI = 0.38 + 0.26/(0.279 / InertialNum + 1);
																				//mueffi =MuI/(2 * sqrt(secondinvariant))* mw.particles[i].GetP(); //max((MuI / (2 * sqrt(secondinvariant)))* mw.particles[i].GetP(), 0.0);     //
		double Mus;		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		double Tausecondinv = sqrt(mw.particles[i].GetTau().secondinvariant());


		//	if (mw.GetTime() > 0.1)
		//{
		if (type == 15 || (type == 0 && mw.particles[i].GetLoc().GetY() < 0))
		{
			//Rhos = 1200;      //////////////////////////Grain diameter    Crosta -Rice
			//d = 0.002;  /////////////////particle diameter      Crosta
			Substrate_YieldStress = effectivP * Substrate_Mus;
			//InertialNum = secondinvariant*d / sqrt(effectivP / Rhos);
			//Substrate_Critical_inertialnumber = sqrt(abs(mw.particles[i].GetP())*Rhos)*tan(Substrate_internalangle)*d / (High_Viscosity_Factor * Substrate_BingViscosity);
			double Subst__Crit_GammaDot = Substrate_YieldStress / (High_Viscosity_Factor * Substrate_BingViscosity);


			if (/*InertialNum>Substrate_Critical_inertialnumber*/mw.particles[i].GetGammaDot() > Subst__Crit_GammaDot
				/*Tausecondinv > Substrate_YieldStress*/)
			{
				mueffi = Substrate_BingViscosity + (Substrate_YieldStress - Substrate_BingViscosity * Subst__Crit_GammaDot) / (secondinvariant);

			}
			else
			{
				mueffi = High_Viscosity_Factor * Substrate_BingViscosity;
			}
		}
		if (type == 11 || (type == 0 && mw.particles[i].GetLoc().GetY() > -0.005))
		{

			//Rhos = 1460;      //////////////////////////Grain diameter    Crosta -Rice
			//d = 0.003;  /////////////////particle diameter      Crosta
			moving_mass_YieldStress = effectivP * moving_mass_Mus;
			//InertialNum = secondinvariant*d / sqrt(effectivP / Rhos);
			//moving_mass_Critical_inertialnumber = sqrt(abs(mw.particles[i].GetP())*Rhos)*tan(moving_mass_internalangle)*d / (3000 * 0.1);
			double Mov_Mass_Crit_GammaDot = moving_mass_YieldStress / (High_Viscosity_Factor * moving_mass_BingViscosity);

			if (mw.particles[i].GetGammaDot() > Mov_Mass_Crit_GammaDot /*Tausecondinv > Substrate_YieldStress*/)
			{
				mueffi = moving_mass_BingViscosity + (moving_mass_YieldStress - moving_mass_BingViscosity * Mov_Mass_Crit_GammaDot) / (secondinvariant);

			}
			else
			{
				mueffi = (High_Viscosity_Factor * moving_mass_BingViscosity);
			}
		}
		//	}
		//else
		//	mueffi = 1;


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////



		double currStep = mw.GetSTEP();
		//	if (mw.GetSTEP() != 1)
		mw.particles[i].SetTau(mueffi * GammaDot[i]);

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
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;
		double hm = (mw.particles[ii].GetH() + mw.particles[jj].GetH()) / 2;


		if (typei < fluid1 && typej < fluid1)
			continue;



		/*if ((mw.particles[ii].GetMaterial()->type == 0 && mw.particles[ii].GetLoc().GetY() > 0)
		|| (mw.particles[jj].GetMaterial()->type == 0 && mw.particles[jj].GetLoc().GetY() > 0))
		continue;*/


		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());
		/*
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() *((mw.particles[ii].B*dw)*(mw.particles[ii].GetTau()* RhoAsquare +
		mw.particles[jj].GetTau()*RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() *((mw.particles[jj].B*dw)*(mw.particles[ii].GetTau()* RhoAsquare +
		mw.particles[jj].GetTau()*RhoBsquare));
		*/
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));

		/*	if(typei==15 && typej==11)
		YieldStress_average[ii]= YieldStress_average[ii]+ abs(mw.particles[jj].GetP())*tan(moving_mass_internalangle);


		if (typei == 11 && typej == 15)
		YieldStress_average[jj] = YieldStress_average[jj] + abs(mw.particles[ii].GetP())*tan(moving_mass_internalangle);

		*/
		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());
		//	const STensor tauij = 1.0 / rhoij* (mw.particles[jj].GetTau() - mw.particles[ii].GetTau());
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction
		//TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() / rhoij*((mw.particles[ii].B*dw)*tauij);
		//TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] + mw.particles[ii].GetM() / rhoij*((mw.particles[jj].B*dw)*tauij);

		/*
		TurbulTauDivergenceOverRho[ii] = TurbulTauDivergenceOverRho[ii] + mw.particles[jj].GetM() *((mw.particles[ii].B*dw)*(mw.particles[ii].GetTurbulTau()* RhoAsquare +
		mw.particles[jj].GetTurbulTau()*RhoBsquare));

		TurbulTauDivergenceOverRho[jj] = TurbulTauDivergenceOverRho[jj] + mw.particles[ii].GetM() *((mw.particles[jj].B*dw)*(mw.particles[ii].GetTurbulTau()* RhoAsquare +
		mw.particles[jj].GetTurbulTau()*RhoBsquare));
		*/
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

		///////////////






	}



	for (int i = 0; i < pNo; i++)
	{



		if (/*mw.GetTime()<0.05&&*/ mw.particles[i].GetMaterial()->type == 11 && mw.Walldistance[i]<0.01 && mw.particles[i].GetLoc().GetY()>0.02 /*&& mw.particles[i].GetLoc().GetY()<0.1*/)

			mw.particles[i].SetAcc(mw.particles[i].GetAcc() + 0.2 * TauDivergenceOverRho[i]);



		else

			mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);




	}


}

/*
Tensorvec EvalDisplacementVector(World& mw)			//Evaluates Displacement Vector for elastic constitutive relation
{

	const int iNo = mw.interactions.size();
	const int pNo = mw.particles.size();
	Vectorvec Loc_0;
	std::ifstream Initial_Loc;
	double L_I, H_Max;
	int num, ptype;
	double X0, Y0;
	Tensorvec DX(pNo, Tensor());

	Initial_Loc.open("Initial_Loc.txt");

	if (!Initial_Loc)
	{
		std::cout << "Unable to open file Initial_Loc.txt";

	}



	while (!Initial_Loc.eof())
	{
		Initial_Loc >> num >> ptype >> X0 >> Y0;        //Initializing


		Loc_0.push_back(Vector(X0, Y0));
		//Loc_0 = (X0, Y0);


	}

	Initial_Loc.clear();
	Initial_Loc.close();


	for (register int i = 0; i < pNo; i++)
	{

		Vector Delta_X_Elastic = mw.particles[i].GetLoc() - Loc_0[i];
		DX[i] = oproduct();

	}
	return DX;





}

*/


void Granularmaterial1(World& mw)			//Pressure-dependent Rheology           constant viscosity model (ETHA=cst)            Ionescu 2015   
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
	//mw.Walldistance.assign(pNo, double());
	double Phi_f = 34;   ///Manzella
	double Phi_Bed = 28;  ///Manzella
	//	double Phi_Basal = 10;

	double moving_mass_internalangle = Phi_f * 0.017453;
	double Substrate_internalangle = Phi_Bed * 0.0174532;
	//double Wall_Friction_angle = Phi_Basal*0.01745329251994329576923690768489;    //10

	double Substrate_Mus = tan(Substrate_internalangle);
	double moving_mass_Mus = tan(moving_mass_internalangle);
	//	double Wall_Mus = tan(Wall_Friction_angle);

	double Rhos, d, Substrate_YieldStress, moving_mass_YieldStress, Nearwall_YieldStress;
	double InertialNum;
	double Substrate_Critical_inertialnumber;
	double moving_mass_Critical_inertialnumber;

	double Substrate_BingViscosity = 0.1;   //Pa.s
	double moving_mass_BingViscosity = 0.1; //Pa.s
	double effectivP;
	double High_Viscosity_Factor = 10000;



#if 0
	for (register int i = 0; i < pNo; i++)

	{
		mw.Walldistance[i] = 1;        //Initializing

	}
#endif // 0


	for (int i = 0; i < pNo; i++)

	{
		effectivP = abs(mw.particles[i].GetP());
		const short type = mw.particles[i].GetMaterial()->type;
		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    
		double mueffi = 0, MuI;
		double Betha = 0.000001;
		double R = GammaDot[i].secondinvariant();
		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();
		double Mus;		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double Tausecondinv = sqrt(mw.particles[i].GetTau().secondinvariant());


		//	if (mw.GetTime() > 0.1)
		//{
#ifdef D2
		if (type == 15 || (type == 0 && mw.particles[i].GetLoc().GetY() < 0)) ///Horizontal
#else
		if (type >= 15 || (type == 0 && mw.particles[i].Get_bed_neigh() > mw.particles[i].Get_overlaying_flow_neigh()))

#endif
		{

			Substrate_YieldStress = effectivP * Substrate_Mus;


			//if (mw.particles[i].GetMaterial()->type == 15 && (effectivP == 0 || abs(mw.particles[i].GetDelDotr()) < 1.9))
				//Substrate_YieldStress = effectivP * Substrate_Mus + 10;

			double Subst__Crit_GammaDot = Substrate_YieldStress / (High_Viscosity_Factor * Substrate_BingViscosity);
			if (secondinvariant > Subst__Crit_GammaDot)
			{
				mueffi = Substrate_BingViscosity + (Substrate_YieldStress - Substrate_BingViscosity * Subst__Crit_GammaDot) / (secondinvariant);
			}
			else
			{
				mueffi = High_Viscosity_Factor * Substrate_BingViscosity;
			}
		}
#ifdef D2
		if (type == 11 || (type == 0 && mw.particles[i].GetLoc().GetY() > -0.005))
		{
#else
		if (type >= fluid1 && type < 15 || (type == 0 && mw.particles[i].Get_bed_neigh() < mw.particles[i].Get_overlaying_flow_neigh())) /*&& mw.Walldistance[i] >= mw.particles[i].GetH() / 2)*/
		{
#endif

			moving_mass_YieldStress = effectivP * moving_mass_Mus;
			double Mov_Mass_Crit_GammaDot = moving_mass_YieldStress / (High_Viscosity_Factor * moving_mass_BingViscosity);

			if (secondinvariant > Mov_Mass_Crit_GammaDot)
			{
				mueffi = moving_mass_BingViscosity + (moving_mass_YieldStress - moving_mass_BingViscosity * Mov_Mass_Crit_GammaDot) / (secondinvariant);

			}
			else
			{
				mueffi = (High_Viscosity_Factor * moving_mass_BingViscosity);
			}


		}
		/*
			if ((type >= fluid1 && type < 15 && mw.Walldistance[i] < mw.particles[i].GetH() / 2) || (type == 0 && mw.particles[i].Get_bed_neigh()< mw.particles[i].Get_overlaying_flow_neigh()))
			{

				Nearwall_YieldStress = effectivP*Wall_Mus;

				double Mov_Mass_Crit_GammaDot = Nearwall_YieldStress / (High_Viscosity_Factor*moving_mass_BingViscosity);

				if (mw.particles[i].GetGammaDot() > Mov_Mass_Crit_GammaDot )
				{
					mueffi = moving_mass_BingViscosity + (Nearwall_YieldStress - moving_mass_BingViscosity*Mov_Mass_Crit_GammaDot) / (secondinvariant);

				}
				else
				{
					mueffi = (High_Viscosity_Factor*moving_mass_BingViscosity);
				}
				}
			*/
#if 0
		if ((mw.Walldistance[i] <= 0.01 && type == 11))      /////////////////////Slip BC
		{
			double Boundary_Mus = tan(Wall_Friction_angle);
			moving_mass_YieldStress = effectivP * Boundary_Mus;
			double Mov_Mass_Crit_GammaDot = moving_mass_YieldStress / (High_Viscosity_Factor * moving_mass_BingViscosity);

			if (mw.particles[i].GetGammaDot() > Mov_Mass_Crit_GammaDot)
			{
				mueffi = 1 * (moving_mass_BingViscosity + (moving_mass_YieldStress -
					moving_mass_BingViscosity * Mov_Mass_Crit_GammaDot) / (secondinvariant));
				mw.particles[i].Set_Static_state(0);

			}
			else
			{
				mueffi = 1 * (High_Viscosity_Factor * moving_mass_BingViscosity);
				mw.particles[i].Set_Static_state(1);
			}

		}
#endif // 0

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		mw.particles[i].SetTau(mueffi * GammaDot[i]);
		mw.particles[i].SetMu_effective(mueffi);


		//	else
		//	mw.particles[i].SetTau(0);



		/*
		double Cohesion = 0.0;                   // Fraccarollo -PVC
		double internalangle = 0.59341194567807205615405486128613; //34 Degree         // Babaie -Comsol
		double mueffi;						  //	double muBingij = mw.particles[i].GetMaterial()->mu;
		//	double muBingij = 3000; //0.04;
		double muBingij = 1;
		////////////////////////////Cross
		double I = GammaDot[i].secondinvariant();   //GammaDot[i].secondinvariant();

		double GammaDotMagi = 2 * sqrt(0.075+I);
		double r = sin(internalangle);

		YieldStress = 0.001; //mw.particles[i].GetP()*sin(internalangle);
		/////////////////////////////////////////////////////////Cross
		double b = muBingij*GammaDotMagi;
		//	if (GammaDotMagi != 0)
		mueffi =(1000 * muBingij)*((YieldStress + b) / (YieldStress + 1000 * b));
		double presenttimestep = mw.GetSTEP();
		//if(mw.GetSTEP()!=1)
		mw.particles[i].SetTau(2 * mueffi*GammaDot[i]);
		//else
		//	mw.particles[i].SetTau(0);
		//	else
		//	mueffi = 0.0;

		*/
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
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;
		double hm = (mw.particles[ii].GetH() + mw.particles[jj].GetH()) / 2;


		if (typei < fluid1 && typej < fluid1)
			continue;



		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());

		/*
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
				*/


				///////////////
			/*	double Mu_i = mw.particles[ii].GetMu_effective();
				double Mu_j = mw.particles[jj].GetMu_effective();
				double Mu_interparticle_harmonic = 2 * Mu_i*Mu_j / (Mu_i + Mu_j);   //Hu & Adams (2006) , NodoushanShan & Shakibaeinia (2018)
				double Mu_interparticle_arithmetic = (Mu_i + Mu_j) / 2;   //Hu & Adams (2006) , NodoushanShan & Shakibaeinia (2018)

			double DeltaX = mw.particles[ii].GetH() / 2.6;

				// NodoushanShan & Shakibaeinia(2018)
			#ifdef D2
				double D = 2; //
				double Landa = 1;  //////////???????????????????????
				double n_0 = 1/(DeltaX*DeltaX);
			#else
				double D = 3; //
				double Landa =1 ;  //////////???????????????????????
				double n_0 =1/(DeltaX*DeltaX*DeltaX) ;  // initial particle number density
			#endif

		#if 0
				if (mw.particles[ii].GetInterfaceneighbor() != 0) ////// close to Flow-Bed interface
					TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + 2 * D / (Landa*n_0* mw.particles[ii].GetRho())*Mu_interparticle_harmonic*(mw.particles[jj].GetVel() - mw.particles[ii].GetVel())*(mw.interactions[i].GetW());
				else if (mw.particles[jj].GetInterfaceneighbor() != 0) ////// close to Flow-Bed interface
					TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] + 2 * D / (Landa*n_0* mw.particles[jj].GetRho())*Mu_interparticle_harmonic*(mw.particles[ii].GetVel() - mw.particles[jj].GetVel())*(mw.interactions[i].GetW());


		#endif // 0

				//}
				////////////////////////
				*/
				//else
		{
			TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
				mw.particles[jj].GetTau() * RhoBsquare));
			TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
				mw.particles[jj].GetTau() * RhoBsquare));
		}







	}



	for (int i = 0; i < pNo; i++)
	{






		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);




	}


	}




void Granularmaterial3(World & mw)			//Pressure-dependent Rheology           Variable viscosity model (Mu(I) Rheology)            Ionescu 2015   
{
	double Cs = 0.1;
	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	double GammaDotMagi = 0, GammaDotMagj = 0, YieldStress = 0, YieldStress_Bed = 0, YieldStress_close_to_wall=0, Mu2, Mus, Mus_Bed, I0, k, InertialNum = 0;
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor());;
	mw.DelV = EvalVelocityGradient(mw);
	mw.Wallneighbor.assign(pNo, double());                      ////////////////Important
	// EvalDisplacementVector(mw);

	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());

#if 0

	std::ofstream Inertialnum;
	Inertialnum.open("Inertialnum.dat", std::ofstream::out | std::ofstream::app);
	Inertialnum << "Variables = Num  X Y I" << '\n';
	Inertialnum << "Zone" << '\n';

#endif // 0

	for (int i = 0; i < pNo; i++)

	{
		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    
		double  MuI = 1000;

		double Betha = 0.075;

		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();


		double Phi_flow = 20;
		double Phi_Bed = 30;
		double Phi_Basal = 28;


		double internalangle = Phi_flow * 0.017453;          //25.5 Degree        // Ionescu 2015
		double internalangle_Bed = Phi_Bed * 0.017453;          //25.5 Degree        // Ionescu 2015
		double Wall_Friction_angle = Phi_Basal * 0.017453;   

		int type = mw.particles[i].GetMaterial()->type;



		double effectivP = abs(mw.particles[i].GetP());

		Mus = sin(internalangle);
		Mus_Bed = sin(internalangle_Bed);



		YieldStress = effectivP * Mus;
		YieldStress_Bed = effectivP * Mus_Bed;
		if (effectivP == 0 && type == 15)
			YieldStress_Bed = effectivP * Mus_Bed + 10000;

		YieldStress_close_to_wall = effectivP * sin(Wall_Friction_angle);


		I0 = 0.279;
		k = 0.035;
		Mu2 = 0.73;
		double d = 0.0007;  /////////////////particle diameter
		double Rhos = 2500;      //////////////////////////Grain density

		//double Critical_I_Flow = sqrt(effectivP * Rhos) * tan(internalangle) * d / 1000;
		//double Critical_I_Bed = sqrt(effectivP * Rhos) * tan(internalangle_Bed) * d / 1;

		double Mov_Mass_Crit_GammaDot = YieldStress / (3000 * 0.1);
		//	double Bed_Crit_GammaDot = YieldStress_Bed / (10000 * 0.1)/* + 1*/;
		
		double Mov_Mass_Crit_GammaDot_close_to_wall = YieldStress_close_to_wall / (3000 * 0.1);

		double Bed_Crit_GammaDot = YieldStress_Bed / (3000 * 0.1);
		//double Bed_Crit_GammaDot = 0.01*sqrt(effectivP/ Rhos) / 0.007;

		if (type == 11 || (type == 0 && mw.particles[i].Get_bed_neigh() < mw.particles[i].Get_overlaying_flow_neigh()))
		{
			InertialNum = 2 * secondinvariant * d / sqrt(effectivP / Rhos);     ////////////I

			if (/*InertialNum > 0.022 /*Critical_I_Flow /*&& mw.particles[i].GetP()!=0*/secondinvariant > Mov_Mass_Crit_GammaDot)

				MuI = 2 * (Mu2 - Mus) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt(effectivP + 0.000001) / k) + (YieldStress) / (secondinvariant);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}

		/*if ((type == 11 && mw.Walldistance[i] < mw.particles[i].GetH() / 2.2) || (type == 0 && mw.particles[i].Get_bed_neigh() < mw.particles[i].Get_overlaying_flow_neigh()))
		{
			InertialNum = 2 * secondinvariant * d / sqrt(effectivP / Rhos);     ////////////I

			if (secondinvariant > Mov_Mass_Crit_GammaDot_close_to_wall)

				MuI = 2 * (Mu2 - Mus) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt(effectivP + 0.000001) / k) + (YieldStress_close_to_wall) / (secondinvariant);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}
		*/
			   		 	  	  	   

		if (type == 15 || (type == 0 && mw.particles[i].Get_bed_neigh() > mw.particles[i].Get_overlaying_flow_neigh()))
		{
			InertialNum = 2 * secondinvariant * d / sqrt(effectivP / Rhos);     ////////////I

			if (/*InertialNum > 0.003*/ /*Critical_I_Bed/* && mw.particles[i].GetP() != 0*/secondinvariant > Bed_Crit_GammaDot)
				//	if (mw.particles[i].GetGammaDot() > Bed_Crit_GammaDot /*Tausecondinv > Substrate_YieldStress*/)

				MuI = 2 * (Mu2 - Mus) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt((effectivP + 0.000001)) / k) + (YieldStress_Bed) / (secondinvariant);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}

		mw.particles[i].Set_I(InertialNum);


		if (mw.GetSTEP() > 1)
		{
			mw.particles[i].Set_I(InertialNum);
			//Inertialnum << i << '\t';
			//Inertialnum << mw.particles[i].GetLoc() << '\t';
			//Inertialnum << InertialNum << '\n';
		}


		mw.particles[i].SetMu_effective(MuI);

		mw.particles[i].SetTau(MuI * GammaDot[i]);


		STensor SigmaTangent;


	}
	//Inertialnum.close();
	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;


		if (typei < fluid1 && typej < fluid1)
			continue;


		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());
		/*
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() *((mw.particles[ii].B*dw)*(mw.particles[ii].GetTau()* RhoAsquare +
		mw.particles[jj].GetTau()*RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() *((mw.particles[jj].B*dw)*(mw.particles[ii].GetTau()* RhoAsquare +
		mw.particles[jj].GetTau()*RhoBsquare));
		*/
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));


		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

		double hm = (mw.particles[jj].GetH() + mw.particles[jj].GetH()) / 2;

		if (mw.particles[ii].GetMaterial()->type >= fluid1 && mw.particles[jj].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[ii] = mw.Wallneighbor[ii] + 1;



			//	mw.wallNormal[ii] = mw.wallNormal[ii] + mw.gradSum[jj] / mw.Wallneighbor[ii];

		}
		if (mw.particles[jj].GetMaterial()->type >= fluid1 && mw.particles[ii].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[jj] = mw.Wallneighbor[jj] + 1;



			//	mw.wallNormal[jj] = mw.wallNormal[jj] + mw.gradSum[ii] / mw.Wallneighbor[jj];
		}




	}



	for (int i = 0; i < pNo; i++)
	{


		/*
		double Betha = 0.075;
		double effectivP = abs(mw.particles[i].GetP());

		double secondinvariant = 2 * sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();

		double internalangle = 0.4450589592585540421155411459646;          //25.5 Degree        // Ionescu 2015
		double Mus = sin(internalangle);
		double d = 0.0007;  /////////////////particle diameter
		double Rhos = 2500;      //////////////////////////Grain Density


		double InertialNum = 2 * secondinvariant*d / sqrt(effectivP / Rhos);     ////////////I



		double	YieldStress = effectivP*Mus;

		double Tausecondinv = 2 * sqrt(mw.particles[i].GetTau().secondinvariant());

		if (InertialNum > 0.01 && mw.GetSTEP() != 1 && mw.particles[i].GetMaterial()->type != 0)

		*/
		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);

		

	}



}   //Orig



/////////////////Backup


void Granularmaterial4(World & mw)			//Pressure-dependent Rheology           Variable viscosity model (Mu(I) Rheology)            Ionescu 2015   
{
	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	double GammaDotMagi = 0, GammaDotMagj = 0, YieldStress = 0, YieldStress_Bed = 0, Mu2, Mus, Mus_Bed, I0, k, InertialNum = 0;
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor());;
	mw.DelV = EvalVelocityGradient(mw);
	mw.Wallneighbor.assign(pNo, double());                      ////////////////Important
																// EvalDisplacementVector(mw);
	doublevec Vol_Frac_Flow(pNo, double()), Vol_Frac_Bed(pNo, double());

	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());


	double thresh_volume_conc = 0.5;


	std::ofstream Inertialnum;
	Inertialnum.open("Inertialnum.dat", std::ofstream::out | std::ofstream::app);
	Inertialnum << "Variables = Num  X Y I" << '\n';
	Inertialnum << "Zone" << '\n';

	for (int i = 0; i < pNo; i++)

	{
		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    
		double  MuI = 1000;

		double Betha = 0.075;

		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();


		double Phi_flow = 30;
		double Phi_Bed = 30;

		double internalangle = Phi_flow * 0.01745329251994329576923690768489;          //25.5 Degree        // Ionescu 2015
		double internalangle_Bed = Phi_Bed * 0.01745329251994329576923690768489;          //25.5 Degree        // Ionescu 2015

		double effectivP = abs(mw.particles[i].GetP());





		Mus = sin(internalangle);
		Mus_Bed = sin(internalangle_Bed);



		YieldStress = effectivP * Mus;
		YieldStress_Bed = effectivP * Mus_Bed;
		if (mw.particles[i].GetMaterial()->type == 15 && (effectivP == 0 || abs(mw.particles[i].GetDelDotr()) < 1.9))
			YieldStress_Bed = effectivP * Mus_Bed + 100;


		I0 = 0.279;
		Mu2 = 0.73;
		//double d = 0.0008;  /////////////////particle diameter


		//double Rho_Flow = 4100; //////////////////////////Grain density-Moberely
		//double Rho_Bed = 4100; //////////////////////////Grain density-Moberely

		double Rho_Flow = 2500; //////////////////////////Grain density-Mangeney
		double Rho_Bed = 2500; //////////////////////////Grain density-Mangeney

							   //	double Rho_Flow = 0.98; //////////////////////////Grain density-Larcher 2018
							   //double Rho_Bed = 0.98; //////////////////////////Grain density-Larcher 2018

		double D_Flow = 0.0007;  /////////////////particle diameter
		double D_Bed = 0.0007;   /////////////////particle diameter



		Vol_Frac_Flow[i] = mw.particles[i].Get_Volume_fraction_flow();
		Vol_Frac_Bed[i] = 1 - Vol_Frac_Flow[i];




		int type = mw.particles[i].GetMaterial()->type;
		//double Critical_I_Flow = sqrt(effectivP*Rhos)*tan(internalangle)*d / 1000;
		//double Critical_I_Bed = sqrt(effectivP*Rhos)*tan(internalangle_Bed)*d / 1;

		double Mov_Mass_Crit_GammaDot = YieldStress / (3000 * 0.1);

		double Bed_Crit_GammaDot = YieldStress_Bed / (3000 * 0.1);



#if 0

		if ((type == 11 && Vol_Frac_Flow[i] < 0.6))
		{
			double d_mixture = Vol_Frac_Flow[i] * D_Flow + Vol_Frac_Bed[i] * D_Bed;
			double Rho_mixture = Vol_Frac_Flow[i] * Rho_Flow + Vol_Frac_Bed[i] * Rho_Bed;
			k = d_mixture * sqrt(Rho_mixture);
			//InertialNum = 2 * secondinvariant*d_mixture / sqrt(effectivP / Rho_mixture);     ////////////I

			MuI = 2 * (Mu2 - Mus) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt(effectivP + 0.000001) / k) + (YieldStress) / (secondinvariant);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}

#endif // 0

		if ((type == 11 /*&& Vol_Frac_Flow[i] > 0.5*/) || (type == 0 && mw.particles[i].Get_Volume_fraction_flow() > mw.particles[i].Get_Volume_fraction_bed()))
		{
			InertialNum = 2 * secondinvariant * D_Flow / sqrt(effectivP / Rho_Flow);     ////////////I
			k = D_Flow * sqrt(Rho_Flow);
			if (/*InertialNum > 0.022 /*Critical_I_Flow /*&& mw.particles[i].GetP()!=0*/secondinvariant > Mov_Mass_Crit_GammaDot)

				MuI = 2 * (Mu2 - Mus) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt(effectivP + 0.000001) / k) + (YieldStress) / (secondinvariant);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}

		if (type == 15 || (type == 0 && mw.particles[i].Get_Volume_fraction_flow() < mw.particles[i].Get_Volume_fraction_bed()))
		{
			InertialNum = 2 * secondinvariant * D_Bed / sqrt(effectivP / Rho_Bed);     ////////////I
			k = D_Bed * sqrt(Rho_Bed);

			if (/*InertialNum > 0.003*/ /*Critical_I_Bed/* && mw.particles[i].GetP() != 0*/secondinvariant > Bed_Crit_GammaDot)
				//	if (mw.particles[i].GetGammaDot() > Bed_Crit_GammaDot /*Tausecondinv > Substrate_YieldStress*/)

				MuI = 2 * (Mu2 - Mus) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt((effectivP + 0.000001)) / k) + (YieldStress_Bed) / (secondinvariant);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}

		mw.particles[i].Set_I(InertialNum);


		if (mw.GetSTEP() > 1)
		{
			mw.particles[i].Set_I(InertialNum);
			Inertialnum << i << '\t';
			Inertialnum << mw.particles[i].GetLoc() << '\t';
			Inertialnum << InertialNum << '\n';
		}


		mw.particles[i].SetMu_effective(MuI);

		mw.particles[i].SetTau(MuI * GammaDot[i]);


		STensor SigmaTangent;


	}
	Inertialnum.close();
	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;


		if (typei < fluid1 && typej < fluid1)
			continue;


		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());
		/*
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() *((mw.particles[ii].B*dw)*(mw.particles[ii].GetTau()* RhoAsquare +
		mw.particles[jj].GetTau()*RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() *((mw.particles[jj].B*dw)*(mw.particles[ii].GetTau()* RhoAsquare +
		mw.particles[jj].GetTau()*RhoBsquare));
		*/
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));


		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

		double hm = (mw.particles[jj].GetH() + mw.particles[jj].GetH()) / 2;

		if (mw.particles[ii].GetMaterial()->type >= fluid1 && mw.particles[jj].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[ii] = mw.Wallneighbor[ii] + 1;



			//	mw.wallNormal[ii] = mw.wallNormal[ii] + mw.gradSum[jj] / mw.Wallneighbor[ii];

		}
		if (mw.particles[jj].GetMaterial()->type >= fluid1 && mw.particles[ii].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[jj] = mw.Wallneighbor[jj] + 1;



			//	mw.wallNormal[jj] = mw.wallNormal[jj] + mw.gradSum[ii] / mw.Wallneighbor[jj];
		}




	}



	for (int i = 0; i < pNo; i++)
	{


		/*
		double Betha = 0.075;
		double effectivP = abs(mw.particles[i].GetP());

		double secondinvariant = 2 * sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();

		double internalangle = 0.4450589592585540421155411459646;          //25.5 Degree        // Ionescu 2015
		double Mus = sin(internalangle);
		double d = 0.0007;  /////////////////particle diameter
		double Rhos = 2500;      //////////////////////////Grain Density


		double InertialNum = 2 * secondinvariant*d / sqrt(effectivP / Rhos);     ////////////I



		double	YieldStress = effectivP*Mus;

		double Tausecondinv = 2 * sqrt(mw.particles[i].GetTau().secondinvariant());

		if (InertialNum > 0.01 && mw.GetSTEP() != 1 && mw.particles[i].GetMaterial()->type != 0)

		*/
		//if (mw.particles[i].GetMaterial()->type==11 && mw.Walldistance[i] < 0.02)
		//	mw.particles[i].SetAcc(mw.particles[i].GetAcc() + 0.1*TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);


		//else
		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);




	}



}   //Orig



#if 0

	//Mixing definition - modification of mixong layer
void Granularmaterial4(World& mw)			//Mixing definition - modification of mixing layer
{
	double GammaDotMagi = 0, GammaDotMagj = 0, YieldStress_Flow = 0, YieldStress_Bed = 0, Mu2_Flow, Mu2_Bed, Mus_Flow, Mus_Bed, I0_Flow, I0_Bed, k_Flow, k_Bed, InertialNum = 0;

	double thresh_volume_conc = 0.9;

	//Same density for Flow and Bed
/*	double Rho_Flow = 2500; //////////////////////////Grain density-Mangeney
	double Rho_Bed = 2500; //////////////////////////Grain density-Mangeney
	double D_Flow = 0.0007;  /////////////////particle diameter
	double D_Bed = 0.0007;   /////////////////particle diameter
	*/



	//Higher density for Flow
	/*double Rho_Flow = 4100; //////////////////////////Grain density-Mangeney
	double Rho_Bed = 2500; //////////////////////////Grain density-Mangeney
	double D_Flow = 0.002;  /////////////////particle diameter
	double D_Bed = 0.0007;   /////////////////particle diameter
	*/


	//Higher density for Bed
	double Rho_Flow = 2500; //////////////////////////Grain density-Mangeney
	double Rho_Bed = 4100; //////////////////////////Grain density-Mangeney
	double D_Flow = 0.0007;  /////////////////particle diameter
	double D_Bed = 0.002;  /////////////////particle diameter


	double Phi_flow = 30;
	double Phi_Bed = 30;

	I0_Flow = 0.279;
	I0_Bed = 0.279;
	Mu2_Flow = 0.73;
	Mu2_Bed = 0.73;
	//double d = 0.0008;  /////////////////particle diameter
	double Theta = 0.5235987756;    //Angle of bed =30 degree


	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor());
	mw.DelV = EvalVelocityGradient(mw);
	mw.Wallneighbor.assign(pNo, double());
	doublevec Vol_Frac_Flow(pNo, double()), Vol_Frac_Bed(pNo, double());

	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());
	double ts = mw.GetTime() / mw.options.SAVINGTIME;
	std::ofstream Inertialnum;

	if ((ts - (int)ts) > (1.0 - (mw.GetDT() / mw.options.SAVINGTIME)))
	{
		Inertialnum.open("Inertialnum.dat", std::ofstream::out | std::ofstream::app);
		Inertialnum << "Variables = Num  X Y GammaDot Mu_eff I" << '\n';
		Inertialnum << "Zone" << '\n';
	}


	for (int i = 0; i < pNo; i++)

	{
		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    
		double  MuI = 1000;

		double Betha = 0.075;
		int type = mw.particles[i].GetMaterial()->type;

		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();
		double internalangle_Flow = Phi_flow * 0.01745329;          //25.5 Degree        // Ionescu 2015
		double internalangle_Bed = Phi_Bed * 0.01745329;          //25.5 Degree        // Ionescu 2015
		double effectivP = abs(mw.particles[i].GetP());

		Mus_Flow = sin(internalangle_Flow);
		Mus_Bed = sin(internalangle_Bed);

		YieldStress_Flow = effectivP * Mus_Flow;
		YieldStress_Bed = effectivP * Mus_Bed;
		if (mw.particles[i].GetMaterial()->type == 15 && (effectivP == 0 || abs(mw.particles[i].GetDelDotr()) < 1.9))
			YieldStress_Bed = effectivP * Mus_Bed + 100;

		k_Flow = D_Flow * sqrt(Rho_Flow);
		k_Bed = D_Bed * sqrt(Rho_Bed);

		Vol_Frac_Flow[i] = mw.particles[i].Get_Volume_fraction_flow();
		Vol_Frac_Bed[i] = 1 - Vol_Frac_Flow[i];

		double Mov_Mass_Crit_GammaDot = YieldStress_Flow / (3000 * 0.1);
		double Bed_Crit_GammaDot = YieldStress_Bed / (3000 * 0.1);

		double secondinvariant_New = sqrt(Rotation_to_tetha(GammaDot[i], -Theta).secondinvariant());  //Total Stress

		///////////////////////Entrainment criteria-Start
		bool Entrained = false;
		if (type == 15 /*&& mw.particles[i].GetVel().TwoNorm()==0*/)
		{
			/*
						double W_Bed = 0.5236 * Rho_Bed * 9.81 * (D_Bed * D_Bed * D_Bed);   /////////////////particle weight
						double A_Bed = 0.7584 * (D_Bed * D_Bed);  /////////////////particle area
						STensor Identity_tensor;
						Identity_tensor.SetAll(1.0, 0.0, 1.0);

						Vector Normal_to_interaface, Tang_to_interaface;

						Normal_to_interaface.SetX(sin(Theta));
						Normal_to_interaface.SetY(cos(Theta));
						Tang_to_interaface.SetX(cos(Theta));
						Tang_to_interaface.SetY(-sin(Theta));

						STensor Tang_proj_op = Identity_tensor - OProduct_same_vec(Normal_to_interaface, Normal_to_interaface);

						//	Normal_interaface << i << '\t';
						//	Normal_interaface << type << '\t';
						//	Normal_interaface << mw.particles[i].GetLoc() << '\t';
						//	Normal_interaface << Normal_to_interaface << '\n';

							//double Theta = atan(mw.particles[i].GetNormal_erosion().GetX() / mw.particles[i].GetNormal_erosion().GetY());

						//STensor Stress_on_interface = MuI * GammaDot[i]
						==+ mw.particles[i].GetP()*Identity_tensor;   //فشار رو وارد کن

					//	STensor Active_Total_Stress = Rotation_to_tetha(Stress_on_interface, -Theta);  //Total Stress

					//	Vector Active_Normal_Stress = Active_Total_Stress * Normal_to_interaface; // Stress vector Normal to interface

						//Vector Active_Shear_Stress = Tang_proj_op*Active_Normal_Stress;  // Stress vector tangent to interface   درست نیست
					//	Vector Active_Shear_Stress = Active_Total_Stress * Tang_to_interaface;  // Stress vector tangent to interface

						Vector Active_Shear_Stress = Active_Total_Stress - Active_Normal_Stress;
						Vector Passive_Shear_Stress = -Active_Normal_Stress.TwoNorm() * tan(internalangle_Bed) * Tang_to_interaface;
						double s = Active_Shear_Stress.TwoNorm();
						double n = Passive_Shear_Stress.TwoNorm();
						Entrained =(((Active_Shear_Stress.TwoNorm()*A_Bed + W_Bed*sin(Theta)) > Passive_Shear_Stress.TwoNorm()*A_Bed)
						 || (Active_Normal_Stress.TwoNorm()*A_Bed > W_Bed*cos(Theta))) ;

						*/
			Entrained = secondinvariant > Bed_Crit_GammaDot;
			//Entrained = Active_Shear_Stress.TwoNorm() > Passive_Shear_Stress.TwoNorm();
		}
		//////////////////////////Entrainment criteria-End

		if ((type == 11 && Vol_Frac_Flow[i] > thresh_volume_conc)
			||
			(type == 0 && mw.particles[i].Get_Volume_fraction_flow() > mw.particles[i].Get_Volume_fraction_bed()))///////////Wall particles neighbor to Flow
		{
			InertialNum = 2 * secondinvariant * D_Flow / sqrt(effectivP / Rho_Flow);     ////////////I

			if (secondinvariant > Mov_Mass_Crit_GammaDot)

				MuI = 2 * (Mu2_Flow - Mus_Flow) * (effectivP + 0.000001) / (2 * secondinvariant + I0_Flow * sqrt(effectivP + 0.000001) / k_Flow) + (YieldStress_Flow) / (secondinvariant);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}


		/////////////////////////////////Mixing layer particle

		if ((type == 11 && Vol_Frac_Flow[i] < thresh_volume_conc) || (type == 15 && Entrained))
		{
			if (type == 15)
			{
				InertialNum = 2 * secondinvariant * D_Bed / sqrt(effectivP / Rho_Bed);     ////////////I

				double Phi_Max_Bed = 0.5;   ////Maximum volume fraction of the bed
				double Phi_Min_Bed = 0.4;   ///Minimum volume fraction of the bed
				if (InertialNum < 1)
					Vol_Frac_Bed[i] = Phi_Max_Bed + (Phi_Min_Bed - Phi_Max_Bed) * InertialNum;
				Vol_Frac_Flow[i] = 1 - Vol_Frac_Bed[i];
			}

			double d_mixture = Vol_Frac_Flow[i] * D_Flow + Vol_Frac_Bed[i] * D_Bed;
			double Rho_mixture = Vol_Frac_Flow[i] * Rho_Flow + Vol_Frac_Bed[i] * Rho_Bed;
			double k_mixture = d_mixture * sqrt(Rho_mixture);
			double Phi_Mixture = Vol_Frac_Flow[i] * internalangle_Flow + Vol_Frac_Bed[i] * internalangle_Bed;
			double Mus_mixture = sin(Phi_Mixture);
			double I0_mixture = Vol_Frac_Flow[i] * I0_Flow + Vol_Frac_Bed[i] * I0_Bed;
			double Mu2_mixture = Vol_Frac_Flow[i] * Mu2_Flow + Vol_Frac_Bed[i] * Mu2_Bed;

			MuI = 2 * (Mu2_mixture - Mus_mixture) * (effectivP + 0.000001) / (2 * secondinvariant + I0_mixture * sqrt(effectivP + 0.000001) / k_mixture) + (effectivP * Mus_mixture) / (secondinvariant);

		}

		if (type == 0 && mw.particles[i].Get_Volume_fraction_flow() < mw.particles[i].Get_Volume_fraction_bed())  ///////////Wall particles neighbor to bed
		{
			if (secondinvariant > Bed_Crit_GammaDot)
				MuI = 2 * (Mu2_Bed - Mus_Bed) * (effectivP + 0.000001) / (2 * secondinvariant + I0_Bed * sqrt((effectivP + 0.000001)) / k_Bed) + (YieldStress_Bed) / (secondinvariant);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}

		if (mw.GetSTEP() > 1)
		{
			mw.particles[i].Set_I(InertialNum);
			Inertialnum << i << '\t';
			Inertialnum << mw.particles[i].GetLoc() << '\t';
			Inertialnum << secondinvariant << '\t';
			Inertialnum << MuI << '\t';
			Inertialnum << InertialNum << '\n';
		}

		mw.particles[i].Set_I(InertialNum);
		mw.particles[i].SetMu_effective(MuI);
		mw.particles[i].SetTau(MuI * GammaDot[i]);

		STensor SigmaTangent;

	}
	Inertialnum.close();

	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;


		if (typei < fluid1 && typej < fluid1)
			continue;


		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());

		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));


		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

		double hm = (mw.particles[jj].GetH() + mw.particles[jj].GetH()) / 2;

		if (mw.particles[ii].GetMaterial()->type >= fluid1 && mw.particles[jj].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[ii] = mw.Wallneighbor[ii] + 1;



			//	mw.wallNormal[ii] = mw.wallNormal[ii] + mw.gradSum[jj] / mw.Wallneighbor[ii];

		}
		if (mw.particles[jj].GetMaterial()->type >= fluid1 && mw.particles[ii].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[jj] = mw.Wallneighbor[jj] + 1;



			//	mw.wallNormal[jj] = mw.wallNormal[jj] + mw.gradSum[ii] / mw.Wallneighbor[jj];
		}




	}



	for (int i = 0; i < pNo; i++)
	{


		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);


	}


}   //Orig  


#endif // 0

	//Without  Mixing definition 
void Granularmaterial5(World& mw)			//Without  Mixing definition 
{
	double GammaDotMagi = 0, GammaDotMagj = 0, YieldStress = 0, YieldStress_Bed = 0, Mu2_Flow, Mu2_Bed, Mus_Flow, Mus_Bed, I0, k, InertialNum = 0;

	double Phi_flow = 30;
	double Phi_Bed = 30;

	I0 = 0.279;
	Mu2_Flow = 0.73;
	Mu2_Bed = 0.73;

	//Higher density for Flow
	double Rho_Flow = 4100; //////////////////////////Grain density-Mangeney
	double Rho_Bed = 2500; //////////////////////////Grain density-Mangeney
	double D_Flow = 0.002;  /////////////////particle diameter
	double D_Bed = 0.0007;   /////////////////particle diameter
	
	
	//Higher density for Bed
	/*double Rho_Flow = 2500; //////////////////////////Grain density-Mangeney
	double Rho_Bed = 4100; //////////////////////////Grain density-Mangeney
	double D_Flow = 0.0007;  /////////////////particle diameter
	double D_Bed = 0.002;  /////////////////particle diameter
	*/
	
	
	double Theta = 0.52359877559829887307710723054658;    //Angle of bed

	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor());;
	mw.DelV = EvalVelocityGradient(mw);
	mw.Wallneighbor.assign(pNo, double());                      ////////////////Important
																// EvalDisplacementVector(mw);
	doublevec Vol_Frac_Flow(pNo, double()), Vol_Frac_Bed(pNo, double());

	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());

	//double thresh_volume_conc = 0.4;

	for (int i = 0; i < pNo; i++)

	{
		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    
		double  MuI = 1000;

		double Betha = 0.075;


		double internalangle = Phi_flow * 0.0174533;          //25.5 Degree        // Ionescu 2015
		double internalangle_Bed = Phi_Bed * 0.0174533;          //25.5 Degree        // Ionescu 2015

		double effectivP = abs(mw.particles[i].GetP())+0.00001;

		Mus_Flow = sin(internalangle);
		Mus_Bed = sin(internalangle_Bed);
		int type = mw.particles[i].GetMaterial()->type;

		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();
		//double secondinvariant_New = sqrt(Rotation_to_tetha(GammaDot[i], -Theta).secondinvariant());  //Total Stress

		YieldStress = effectivP * Mus_Flow;
		YieldStress_Bed = effectivP * Mus_Bed;
		if (mw.particles[i].GetMaterial()->type == 15 && (effectivP == 0 || abs(mw.particles[i].GetDelDotr()) < 1.9))
			YieldStress_Bed = effectivP * Mus_Bed + 100;


		double Mov_Mass_Crit_GammaDot = YieldStress / (3000 * 0.1);

		double Bed_Crit_GammaDot = YieldStress_Bed / (3000 * 0.1);

		///////////////////////Entrainment criteria-Start


		if (type == 11 || (type == 0 && mw.particles[i].Get_bed_neigh() < mw.particles[i].Get_overlaying_flow_neigh()))
		{
			InertialNum = 2 * secondinvariant * D_Flow / sqrt(effectivP / Rho_Flow);     ////////////I
			k = D_Flow * sqrt(Rho_Flow);
			if (secondinvariant > Mov_Mass_Crit_GammaDot)

				MuI = 2 * (Mu2_Flow - Mus_Flow) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt(effectivP + 0.000001) / k) + (YieldStress) / (secondinvariant);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}


		if (type == 15 || (type == 0 && mw.particles[i].Get_bed_neigh() > mw.particles[i].Get_overlaying_flow_neigh()))
		{
			InertialNum = 2 * secondinvariant * D_Bed / sqrt(effectivP / Rho_Bed);     ////////////I
			k = D_Bed * sqrt(Rho_Bed);
			if (secondinvariant > Bed_Crit_GammaDot)
				MuI = 2 * (Mu2_Bed - Mus_Bed) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt((effectivP + 0.000001)) / k) + (YieldStress_Bed) / (secondinvariant);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}

		mw.particles[i].Set_I(InertialNum);
		mw.particles[i].SetMu_effective(MuI);
		mw.particles[i].SetTau(MuI * GammaDot[i]);


		STensor SigmaTangent;

	}

	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;


		if (typei < fluid1 && typej < fluid1)
			continue;


		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());

		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));


		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

		double hm = (mw.particles[jj].GetH() + mw.particles[jj].GetH()) / 2;

		if (mw.particles[ii].GetMaterial()->type >= fluid1 && mw.particles[jj].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[ii] = mw.Wallneighbor[ii] + 1;


		}
		if (mw.particles[jj].GetMaterial()->type >= fluid1 && mw.particles[ii].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[jj] = mw.Wallneighbor[jj] + 1;

		}


	}



	for (int i = 0; i < pNo; i++)
	{


		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);


	}


}   //Orig


#if 0

void Granularmaterial4(World & mw)			//Pressure-dependent Rheology           Variable viscosity model (Mu(I) Rheology)            Ionescu 2015   
{
	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	double GammaDotMagi = 0, GammaDotMagj = 0, YieldStress = 0, YieldStress_Bed = 0, Mu2, Mus, Mus_Bed, I0, k, InertialNum = 0;
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor());;
	mw.DelV = EvalVelocityGradient(mw);
	mw.Wallneighbor.assign(pNo, double());                      ////////////////Important
	// EvalDisplacementVector(mw);
	doublevec Vol_Frac_Flow(pNo, double()), Vol_Frac_Bed(pNo, double());

	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());


	double thresh_volume_conc = 0.5;


#if 0
	std::ofstream Inertialnum;
	Inertialnum.open("Inertialnum.dat", std::ofstream::out | std::ofstream::app);
	Inertialnum << "Variables = Num  X Y I" << '\n';
	Inertialnum << "Zone" << '\n';
#endif // 0

	///////////////////////
	std::ofstream Normal_interaface;
	Normal_interaface.open(("Normal_interaface.plt"), std::ios::app);
	Normal_interaface << "Variables = Num PType X Y Normal_to_interaface_x  Normal_to_interaface_y " << '\n';
	Normal_interaface << "Zone" << '\n';


	std::ifstream Particle_tag;
	int num, tag;
	intvec TAG(pNo, 0);
	Tensorvec DX(pNo, Tensor());

	Particle_tag.open("Particle_tag.plt");

	if (mw.GetSTEP() != 1)
	{
		if (!Particle_tag)
		{
			std::cout << "Unable to open file Initial_Loc.txt";

		}

		//while (!Particle_tag.eof())
		//{
			//Particle_tag >> num >> tag;        //Initializing


			//TAG[i]=tag;

		//}






	}

	/////////////////////////////////////////
	for (int i = 0; i < pNo; i++)

	{

		Particle_tag >> num >> tag;        //Initializing
		TAG[i] = tag;

		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    
		double  MuI = 1000;

		double Betha = 0.075;

		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();
		int type = mw.particles[i].GetMaterial()->type;


		double Phi_flow = 30;
		double Phi_Bed = 10;

		double internalangle = Phi_flow * 0.01745329251994329576923690768489;          //25.5 Degree        // Ionescu 2015
		double internalangle_Bed = Phi_Bed * 0.01745329251994329576923690768489;          //25.5 Degree        // Ionescu 2015

  ///////////////////////////////////////////////////////		Yield Stress

		double effectivP = abs(mw.particles[i].GetP());
		Mus = sin(internalangle);
		Mus_Bed = sin(internalangle_Bed);

		YieldStress = effectivP * Mus;
		YieldStress_Bed = effectivP * Mus_Bed;
		if (type == 15 && (effectivP == 0 || abs(mw.particles[i].GetDelDotr()) < 1.9))
			YieldStress_Bed = effectivP * Mus_Bed * +100;
		///////////////////////////////////////////////////////		Yield Stress-End

		I0 = 0.279;
		Mu2 = 0.73;
		//double d = 0.0008;  /////////////////particle diameter


		//double Rho_Flow = 4100; //////////////////////////Grain density-Moberely
		//double Rho_Bed = 4100; //////////////////////////Grain density-Moberely

		double Rho_Flow = 2500; //////////////////////////Grain density-Mangeney
		double Rho_Bed = 2500; //////////////////////////Grain density-Mangeney

	//	double Rho_Flow = 0.98; //////////////////////////Grain density-Larcher 2018
		//double Rho_Bed = 0.98; //////////////////////////Grain density-Larcher 2018

		double D_Flow = 0.0007;  /////////////////particle diameter
		double D_Bed = 0.0007;   /////////////////particle diameter

		double W_Bed = 0.5236 * Rho_Bed * 9.81 * (D_Bed * D_Bed * D_Bed);   /////////////////particle weight
		double A_Bed = 0.7584 * (D_Bed * D_Bed);  /////////////////particle area

		Vol_Frac_Flow[i] = mw.particles[i].Get_Volume_fraction_flow();
		Vol_Frac_Bed[i] = 1 - Vol_Frac_Flow[i];

		TAG[i] = mw.particles[i].Get_tag();

		double Mov_Mass_Crit_GammaDot = YieldStress / (3000 * 0.1);

		double Bed_Crit_GammaDot = YieldStress_Bed / (3000 * 0.1);


		if (mw.GetSTEP() == 1)
		{
			if (type == 11)
				TAG[i] = 3;

		}



		if (type == 15 && mw.particles[i].GetVel().TwoNorm() == 0)    ///Interface particle
		{

			////////////////////
			STensor Identity_tensor;
			Identity_tensor.SetAll(1.0, 0.0, 1.0);

			Vector Normal_to_interaface = -mw.particles[i].GetNormal_erosion();
			STensor Tang_proj_op = Identity_tensor - OProduct_same_vec(Normal_to_interaface, Normal_to_interaface);

			Normal_interaface << i << '\t';
			Normal_interaface << type << '\t';
			Normal_interaface << mw.particles[i].GetLoc() << '\t';
			Normal_interaface << Normal_to_interaface << '\n';

			//double Theta = atan(mw.particles[i].GetNormal_erosion().GetX() / mw.particles[i].GetNormal_erosion().GetY());
			double Theta = 0.52359877559829887307710723054658;

			STensor Stress_on_interface = MuI * GammaDot[i]/*+ mw.particles[i].GetP()*Identity_tensor*/;   //فشار رو وارد کن


			STensor Active_Total_Stress = Rotation_to_tetha(Stress_on_interface, Theta);  //Total Stress


			Vector Active_Normal_Stress = Normal_to_interaface * Active_Total_Stress; // Stress vector Normal to interface

			Vector Active_Shear_Stress = Tang_proj_op * Active_Total_Stress * Normal_to_interaface;  // Stress vector tangent to interface


			Vector Passive_Shear_Stress = Active_Normal_Stress * tan(internalangle_Bed);

			////////////////////

						/////////////////////////////////Entrainment criteria
			if ((Active_Shear_Stress.TwoNorm() != 0) && (((Active_Shear_Stress.TwoNorm() * A_Bed + W_Bed * sin(Theta)) > Passive_Shear_Stress.TwoNorm() * A_Bed)
				|| (Active_Normal_Stress.TwoNorm() > 0 && (Active_Normal_Stress.TwoNorm() * A_Bed > W_Bed * cos(Theta)))))
			{

				TAG[i] = 2;

			}

		}



		if (TAG[i] == 2)    //Mixing layer particle
		{
			double d_mixture = Vol_Frac_Flow[i] * D_Flow + Vol_Frac_Bed[i] * D_Bed;
			double Rho_mixture = Vol_Frac_Flow[i] * Rho_Flow + Vol_Frac_Bed[i] * Rho_Bed;
			k = d_mixture * sqrt(Rho_mixture);
			if (type == 11)
				MuI = 2 * (Mu2 - Mus) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt(effectivP + 0.000001) / k) + (YieldStress) / (secondinvariant);
			if (type == 15)
				MuI = 2 * (Mu2 - Mus_Bed) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt(effectivP + 0.000001) / k) + (YieldStress_Bed) / (secondinvariant);
			//////////////////////یک تنش بحرانی معادل هم  میشود گذاشت
		}

		if (TAG[i] == 3)   //Overlaying flow particle
		{
			InertialNum = 2 * secondinvariant * D_Flow / sqrt(effectivP / Rho_Flow);     ////////////I
			k = D_Flow * sqrt(Rho_Flow);
			if (secondinvariant > Mov_Mass_Crit_GammaDot)
				MuI = 2 * (Mu2 - Mus) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt(effectivP + 0.000001) / k) + (YieldStress) / (secondinvariant);

		}


#if 0

		if (mw.GetSTEP() > 1)
		{
			mw.particles[i].Set_I(InertialNum);
			Inertialnum << i << '\t';
			Inertialnum << mw.particles[i].GetLoc() << '\t';
			Inertialnum << InertialNum << '\n';
		}

#endif // 0


		mw.particles[i].SetMu_effective(MuI);

		mw.particles[i].SetTau(MuI * GammaDot[i]);


		STensor SigmaTangent;


	}

	Normal_interaface.close();
	Particle_tag.close();



	//Inertialnum.close();
	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;


		if (typei < fluid1 && typej < fluid1)
			continue;


		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());

		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));


		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

		double hm = (mw.particles[jj].GetH() + mw.particles[jj].GetH()) / 2;

		if (mw.particles[ii].GetMaterial()->type >= fluid1 && mw.particles[jj].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[ii] = mw.Wallneighbor[ii] + 1;

		}
		if (mw.particles[jj].GetMaterial()->type >= fluid1 && mw.particles[ii].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[jj] = mw.Wallneighbor[jj] + 1;

		}


	}



	for (int i = 0; i < pNo; i++)
	{


		/*
		double Betha = 0.075;
		double effectivP = abs(mw.particles[i].GetP());

		double secondinvariant = 2 * sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();

		double internalangle = 0.4450589592585540421155411459646;          //25.5 Degree        // Ionescu 2015
		double Mus = sin(internalangle);
		double d = 0.0007;  /////////////////particle diameter
		double Rhos = 2500;      //////////////////////////Grain Density


		double InertialNum = 2 * secondinvariant*d / sqrt(effectivP / Rhos);     ////////////I



		double	YieldStress = effectivP*Mus;

		double Tausecondinv = 2 * sqrt(mw.particles[i].GetTau().secondinvariant());

		if (InertialNum > 0.01 && mw.GetSTEP() != 1 && mw.particles[i].GetMaterial()->type != 0)

		*/
		//if (mw.particles[i].GetMaterial()->type==11 && mw.Walldistance[i] < 0.02)
		//	mw.particles[i].SetAcc(mw.particles[i].GetAcc() + 0.1*TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);


		//else
		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);

	}


}


/*
void Granularmaterial2(World &mw)			//pressure dependent Rheology ////////Implicit Viscosity
{
double Cs = 0.1;
int pNo = mw.particles.size();
int iNo = mw.interactions.size();
double GammaDotMagi = 0, GammaDotMagj = 0, muAppij, muBingij = 0, TauMagi, TauMagj, YieldStress = 0, mueffj, TurbulentViscosity;
STensor Tau;
Tensorvec DV(pNo, Tensor());
STensorvec GammaDot(pNo, STensor());;
mw.DelV = EvalVelocityGradient(mw);

Vectorvec TauDivergenceOverRho(pNo, Vector());
Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());


for (int i = 0; i < pNo; i++)

{


GammaDot[i] = mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)


double mueffi, MuI;

double secondinvariant = GammaDot[i].secondinvariant();   //GammaDot[i].secondinvariant();



double internalangle = 0.44680428851054837169246483673308; //25.6 Degree        // Ionescu 2015
double effectivP = mw.particles[i].GetP();
double alpha = 1;  //Pa.s
//	double betha = 0.075;   //1/(s^-2)
YieldStress = abs(effectivP) *tan(internalangle);

mw.particles[i].SetEffectiveviscosity(alpha + (YieldStress) / (2 * sqrt(betha + secondinvariant)));


}

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


if (typei < fluid1 && typej < fluid1)
continue;

const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());
Vector Vstar[ii] = mw.particles[ii].GetVel();
Vector Vstar[jj] = mw.particles[jj].GetVel();


DelVstar[ii] = DelVstar[ii] + (mw.particles[jj].GetM() / rhoij)*(mw.particles[ii].B*dw)*(Vstar[jj] - Vstar[ii]);
DelVstar[jj] = DelVstar[jj] +or - (mw.particles[ii].GetM() / rhoij)*(mw.particles[jj].B*dw)*(Vstar[jj] - Vstar[ii]);



COEF[2 * i] = (mw.particles[jj].GetM() / (rhoij*rhoij)) * ((2.0*(1.0 / d)* DiadProduct(mw.particles[ii].Bhat, dwe)) + (Dot(T[ii], dw)));
denum[ii] = denum[ii] + COEF[2 * i];


COEF[2 * i + 1] = (mw.particles[ii].GetM() / (rhoij*rhoij)) * ((2.0*(1.0 / d)* DiadProduct(mw.particles[jj].Bhat, dwe)) + (Dot(T[jj], dw)));
denum[jj] = denum[jj] + COEF[2 * i + 1];



num[ii]=mw.particles[ii].GetVel()+





}



for (int Iter = 0; Iter < 100; Iter++)
{




OldVelo[i] = mw.particles[i].GetVel();


Vector Vstarnew=(OldVelo[i]+num[i]) / denum[i];

mw.particles[i].SetVel(Vstarnew);


}




















}
*/


#endif // 0




//Orig

void Elastic_Viscoplastic(World & mw)			//Ghaitanellis (2018) Modelling bed-load sediment transport through a granular approach in SPH   
{
	double Cs = 0.1;
	double E = 5000;  //Young’s modulus
	double ν = 0.3;   //Poisson’s ratio
	double G = E / (2 * (1 + ν));       //Lame’s constant=1923
	double K;
	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	double  YieldStress = 0, YieldStress_Bed = 0, Mu2, Mus, Mus_Bed, I0, k, InertialNum = 0;
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor()), Gamma(pNo, STensor()), Old_Elastic_Tau(pNo, STensor());
	mw.DelV = EvalVelocityGradient(mw);
	//mw.DelX = EvalDisplacementVector(mw);

	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());

	double E_S_XX = 0, E_S_XY = 0, E_S_YY = 0;

#if 0
	std::ofstream Inertialnum;
	Inertialnum.open("Inertialnum.dat", std::ofstream::out);
	Inertialnum << "Variables = Num  X Y I" << '\n';
	Inertialnum << "Zone" << '\n';

#endif // 0
	double d = 0.007;   //Grain diameter
	double Rho_g = 2600;   //Grain density
	double yield_viscosity = 3.14 * d * sqrt(Rho_g * G) / (0.163 * ν + 0.877);  //=53



	std::ifstream Elastic_stress_in;
	Elastic_stress_in.open("Elastic_stress.dat", std::ios::in);



	if (!Elastic_stress_in)	// if in1 not attached to a valid input source, abort
	{
		std::cout << "Sorry, bad file.";
		exit(0);	// special system call to abort program
					//  may require <cstdlib> to be included
	}

	//while (!Elastic_stress_in.eof())
	for (int i = 0; i < pNo; i++)
	{

		Elastic_stress_in >> E_S_XX >> E_S_XY >> E_S_YY;

		Old_Elastic_Tau[i].SetXX(E_S_XX);
		Old_Elastic_Tau[i].SetXY(E_S_XY);
		Old_Elastic_Tau[i].SetYY(E_S_YY);
	}


	Elastic_stress_in.close();


	std::ofstream Elastic_stress_out;
	Elastic_stress_out.open("Elastic_stress.dat", std::ofstream::out);
	Elastic_stress_out << "Variables = EL_st_XX EL_st_XY EL_st_YY" << '\n';
	Elastic_stress_out << "Zone" << '\n';



	for (int i = 0; i < pNo; i++)

	{


		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    
		//Gamma[i] =  mw.DelX[i].GetSymPart();

		double  MuI = 1000;

		double Betha = 0.075;

		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();


		double Phi_flow = 25.5;
		double Phi_Bed = 10;

		double internalangle = Phi_flow * 0.01745329251994329576923690768489;          //25.5 Degree        // Ionescu 2015
		double internalangle_Bed = Phi_Bed * 0.01745329251994329576923690768489;          //25.5 Degree        // Ionescu 2015

		double effectivP = abs(mw.particles[i].GetP());

		Mus = sin(internalangle);
		Mus_Bed = sin(internalangle_Bed);

		YieldStress = effectivP * Mus;
		YieldStress_Bed = effectivP * Mus_Bed + 0.0001;



		I0 = 0.279;
		k = 0.035;
		Mu2 = 0.73;
		double d = 0.0007;  /////////////////particle diameter
		double Rhos = 2500;      //////////////////////////Grain density
		int q = 2;
		int type = mw.particles[i].GetMaterial()->type;
		double Critical_I_Flow = sqrt(effectivP * Rhos) * tan(internalangle) * d / 1000;
		double Critical_I_Bed = sqrt(effectivP * Rhos) * tan(internalangle_Bed) * d / 1;

		//double Mov_Mass_Crit_GammaDot = YieldStress / (10000 * 0.1);


		double Mov_Mass_Crit_GammaDot = YieldStress / yield_viscosity;   //Ghaitanellis (2018)



		double Bed_Crit_GammaDot = YieldStress_Bed / (10000 * 0.1)/* + 1*/;


		if (type == 11 || (type == 0 && mw.particles[i].Get_bed_neigh() < mw.particles[i].Get_overlaying_flow_neigh()))
		{
			InertialNum = 2 * secondinvariant * d / sqrt(effectivP / Rhos);     ////////////I

			if (secondinvariant < Mov_Mass_Crit_GammaDot)

				K = (q + 1) * pow((secondinvariant / Mov_Mass_Crit_GammaDot), q) - q * pow((secondinvariant / Mov_Mass_Crit_GammaDot), (q + 1));
			else
			{
				K = 1;
				MuI = 2 * (Mu2 - Mus) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt(effectivP + 0.000001) / k) + (YieldStress) / (secondinvariant);

			}
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}

		else if (type == 15 || (type == 0 && mw.particles[i].Get_bed_neigh() > mw.particles[i].Get_overlaying_flow_neigh()))
		{
			InertialNum = 2 * secondinvariant * d / sqrt(effectivP / Rhos);     ////////////I

			if (secondinvariant < Bed_Crit_GammaDot)

				K = (q + 1) * pow((secondinvariant / Bed_Crit_GammaDot), q) - q * pow((secondinvariant / Bed_Crit_GammaDot), (q + 1));
			else
			{
				K = 1;
				MuI = 2 * (Mu2 - Mus) * (effectivP + 0.000001) / (2 * secondinvariant + I0 * sqrt((effectivP + 0.000001)) / k) + (YieldStress_Bed) / (secondinvariant);

			}
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		}
		else

			K = 0;




		mw.particles[i].Set_I(InertialNum);





		mw.particles[i].SetMu_effective(MuI);
		STensor Viscoplastic_Stress = MuI * GammaDot[i];

		Tensor Rot_Rate = mw.DelV[i].Rotation_Rate_Tensor();
		STensor Elastic_Stress = Old_Elastic_Tau[i] + (2 * G * GammaDot[i] + Jaumann_rate_tensor(Old_Elastic_Tau[i], Rot_Rate)) * mw.GetDT();   //Ghaitanellis (2018) Modelling bed-load sediment transport through a granular approach in SPH
		mw.particles[i].Set_Elastic_Tau(Elastic_Stress);
		mw.particles[i].SetTau(K * Viscoplastic_Stress + (1 - K) * Elastic_Stress);

		Elastic_stress_out << Elastic_Stress << '\n';


		STensor SigmaTangent;


	}


	Elastic_stress_out.close();

	//Inertialnum.close();

	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;


		if (typei < fluid1 && typej < fluid1)
			continue;


		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());

		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));


		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());

		double hm = (mw.particles[jj].GetH() + mw.particles[jj].GetH()) / 2;

		if (mw.particles[ii].GetMaterial()->type >= fluid1 && mw.particles[jj].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[ii] = mw.Wallneighbor[ii] + 1;


		}
		if (mw.particles[jj].GetMaterial()->type >= fluid1 && mw.particles[ii].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[jj] = mw.Wallneighbor[jj] + 1;

		}




	}



	for (int i = 0; i < pNo; i++)
	{


		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);




	}



}   //Orig







#if 0

////////////////////فرعی
void Granularmaterial3(World & mw)			//Pressure-dependent Rheology           Variable viscosity model (Mu(I) Rheology)            Ionescu 2015   
{
	double Cs = 0.1;
	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	double GammaDotMagi = 0, GammaDotMagj = 0, YieldStress = 0, Mu2, Mus, I0, k;
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	STensorvec GammaDot(pNo, STensor());;
	mw.DelV = EvalVelocityGradient(mw);
	//mw.Wallneighbor.assign(pNo, double());                      ////////////////Important

	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());




	for (int i = 0; i < pNo; i++)

	{


		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    


		double  MuI;

		double Betha = 0.075;

		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();




		double internalangle = 0.4450589592585540421155411459646;          //25.5 Degree        // Ionescu 2015
		double effectivP = abs(mw.particles[i].GetP());

		Mus = sin(internalangle);

		YieldStress = effectivP * Mus;

		I0 = 0.279;
		k = 0.035;
		Mu2 = 0.73;
		double d = 0.0007;  /////////////////particle diameter
		double Rhos = 2500;      //////////////////////////Grain diameter


		double InertialNum = 2 * secondinvariant * d / sqrt(effectivP / Rhos);     ////////////I

		if (InertialNum > 0.001)
			MuI = 2 * (Mu2 - Mus) * effectivP / (2 * secondinvariant + I0 * sqrt(effectivP) / k) + (YieldStress) / (secondinvariant);
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////



		//	double currStep = mw.GetSTEP();
		//if (mw.GetSTEP() != 1)
		mw.particles[i].SetTau(MuI * GammaDot[i]);
		//else
		//	mw.particles[i].SetTau(0);

		STensor SigmaTangent;
		/*
		if (mw.Wallneighbor[i] != 0 && mw.particles[i].GetLoc().GetX() > 0.005 && mw.particles[i].GetLoc().GetY() <0.008) (mw.particles[i].GetH())) /////Liu 1989-2006
		{
		SigmaTangent.SetXX(mw.particles[i].GetTau().GetXX());
		SigmaTangent.SetXY(-Mus*effectivP);
		SigmaTangent.SetYY(mw.particles[i].GetTau().GetYY());

		mw.particles[i].SetTau(SigmaTangent);


		}


		*/





	}

	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;


		if (typei < fluid1 && typej < fluid1)
			continue;


		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());
		/*
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() *((mw.particles[ii].B*dw)*(mw.particles[ii].GetTau()* RhoAsquare +
		mw.particles[jj].GetTau()*RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() *((mw.particles[jj].B*dw)*(mw.particles[ii].GetTau()* RhoAsquare +
		mw.particles[jj].GetTau()*RhoBsquare));
		*/
		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));


		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());
		//	const STensor tauij = 1.0 / rhoij* (mw.particles[jj].GetTau() - mw.particles[ii].GetTau());
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Nikooei-omega correction
		//TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() / rhoij*((mw.particles[ii].B*dw)*tauij);
		//TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] + mw.particles[ii].GetM() / rhoij*((mw.particles[jj].B*dw)*tauij);

		/*
		TurbulTauDivergenceOverRho[ii] = TurbulTauDivergenceOverRho[ii] + mw.particles[jj].GetM() *((mw.particles[ii].B*dw)*(mw.particles[ii].GetTurbulTau()* RhoAsquare +
		mw.particles[jj].GetTurbulTau()*RhoBsquare));

		TurbulTauDivergenceOverRho[jj] = TurbulTauDivergenceOverRho[jj] + mw.particles[ii].GetM() *((mw.particles[jj].B*dw)*(mw.particles[ii].GetTurbulTau()* RhoAsquare +
		mw.particles[jj].GetTurbulTau()*RhoBsquare));
		*/





		double hm = (mw.particles[jj].GetH() + mw.particles[jj].GetH()) / 2;

		if (mw.particles[ii].GetMaterial()->type >= fluid1 && mw.particles[jj].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[ii] = mw.Wallneighbor[ii] + 1;



			//	mw.wallNormal[ii] = mw.wallNormal[ii] + mw.gradSum[jj] / mw.Wallneighbor[ii];

		}
		if (mw.particles[jj].GetMaterial()->type >= fluid1 && mw.particles[ii].GetMaterial()->type < fluid1 && d < 0.8 * hm)
		{
			mw.Wallneighbor[jj] = mw.Wallneighbor[jj] + 1;



			//	mw.wallNormal[jj] = mw.wallNormal[jj] + mw.gradSum[ii] / mw.Wallneighbor[jj];
		}

















	}



	for (int i = 0; i < pNo; i++)
	{


		/*
		double Betha = 0.075;
		double effectivP = abs(mw.particles[i].GetP());

		double secondinvariant = 2 * sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();

		double internalangle = 0.4450589592585540421155411459646;          //25.5 Degree        // Ionescu 2015
		double Mus = sin(internalangle);
		double d = 0.0007;  /////////////////particle diameter
		double Rhos = 2500;      //////////////////////////Grain Density


		double InertialNum = 2 * secondinvariant*d / sqrt(effectivP / Rhos);     ////////////I



		double	YieldStress = effectivP*Mus;

		double Tausecondinv = 2 * sqrt(mw.particles[i].GetTau().secondinvariant());

		if (InertialNum > 0.01 && mw.GetSTEP() != 1 && mw.particles[i].GetMaterial()->type != 0)

		*/
		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]); //+ TurbulTauDivergenceOverRho[i]);




	}



}


#endif // 0




void NonNewtonianViscosityHerschelBulkley(World & mw)			//Symmetric stress -simple difference  (Papanastasiu  of Bi-viscous Modification)
{
	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	double GammaDotMagi = 0, GammaDotMagj = 0, muAppij, muBingij = 0, TauMagi, TauMagj, YieldStress = 0, mueffj, consistencyindex, Powerlawindex;
	STensor Tau;
	Tensorvec DV(pNo, Tensor());
	Vectorvec TaudivergenceoverRho(pNo, Vector());
	Vectorvec nuDel2U(pNo, Vector());
	STensorvec GammaDot(pNo, STensor());;
	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());
	double mueffi;
	double Betha = 0.000001;
	double secondinvarianti, MuEffective;
	doublevec neighwallcounter;
	mw.Wallneighbor.assign(pNo, double());                      ////////////////Important
	mw.wallNormal.assign(pNo, Vector());
	mw.Walldistance.assign(pNo, double());

	mw.muDV = EvalMuVelocityGradient(mw);
	mw.DelV = EvalVelocityGradient(mw);

	double Max_Viscosity = 0;

	for (int i = 0; i < pNo; i++)

	{

		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)    


		double secondinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();

																				///////////////////////////////Rhological properties
																				/*	YieldStress = 89;
																				double n = 0.415;
																				double Powerlawindex = -n + 1;
																				double K = 47.68;*/   ////////////////Bates, Ancey 2016
																				///////////////////////////////
		YieldStress = 58;
		double n = 0.33;
		double Powerlawindex = -n + 1;
		double K = 35;

		//	double Mu = 1;





		///////////////////////////////////////////////////////////Modified
		//	if (secondinvarianti < YieldStress / (100 * K))
		//	MuEffective = 100 * K / secondinvarianti;
		//else
		//	MuEffective = (YieldStress / (secondinvarianti)+K / pow((secondinvarianti), Powerlawindex));
		///////////////////////////////////////////////////////////


		///////////////////////////////////////////Papanastasiu
		double m = 100;

		double Critical_Mu = m * YieldStress;


		MuEffective = K / pow(secondinvariant, Powerlawindex) + (YieldStress / secondinvariant) * (1 - exp(-m * secondinvariant));

		if (MuEffective > Critical_Mu)
			MuEffective = Critical_Mu;


		//	if(MuEffective>Max_Viscosity)
		//	mw.

		mw.particles[i].SetTau(MuEffective * GammaDot[i]);                /////Herschel-Bulkley-Papanastasiou model


																		//////////////////////////////////////////////////////////////////////////////////////////Ezafi
																		//double Taumag = 2 * sqrt((MuEffective*GammaDot[i]).secondinvariant());
																		//	TurbulentViscosity = (Cs*DeltaX)*(Cs*DeltaX)*GammaDotMagi;
																		//mw.particles[i].SetTurbulTau(2 * mw.particles[i].GetRho()*TurbulentViscosity*GammaDot[i]);
																		//////////////////////////////////////////////////////////////////////////////////////////
	}


	for (register int i = 0; i < pNo; i++)
		mw.Walldistance[i] = 1;        //Initializing



	for (register int i = 0; i < iNo; i++)
	{

		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();

		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;

		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;
		const Vector dv = (mw.particles[ii].GetVel() - mw.particles[jj].GetVel());
		double h = mw.particles[ii].GetH();

		double Etha = 0.001 * h;   //////////////Bui 2016

		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());


		if (typei < fluid1 && typej < fluid1)
			continue;




		if (typei == wall2 || typej == wall2)
			continue;


		double RhoAsquare = 1 / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho());
		double RhoBsquare = 1 / (mw.particles[jj].GetRho() * mw.particles[jj].GetRho());


		TauDivergenceOverRho[ii] = TauDivergenceOverRho[ii] + mw.particles[jj].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));
		TauDivergenceOverRho[jj] = TauDivergenceOverRho[jj] - mw.particles[ii].GetM() * ((dw) * (mw.particles[ii].GetTau() * RhoAsquare +
			mw.particles[jj].GetTau() * RhoBsquare));

		///////////////////////////////////Simple stress
		/*
		const STensor tauij = 1.0 / rhoij* (mw.particles[jj].GetTau() - mw.particles[ii].GetTau());
		TaudivergenceoverRho[ii] += mw.particles[jj].GetM() / mw.particles[ii].GetRho()*((dw)*tauij);
		TaudivergenceoverRho[jj] += mw.particles[ii].GetM() / mw.particles[jj].GetRho()*((dw)*tauij);
		*/

		///////////////////////////////////Simple stress



	}

	for (int i = 0; i < pNo; i++)
	{
		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + TauDivergenceOverRho[i]);
	}




}








#if 0

void NonNewtonianViscosityHerschelBulkley(World& mw)			//Muijstress
{
	int pNo = mw.particles.size();
	int iNo = mw.interactions.size();
	double GammaDotMagi = 0, GammaDotMagj = 0, muAppij, muBingij = 0, TauMagi, TauMagj, YieldStress = 0, mueffj, consistencyindex, Powerlawindex;
	STensor Tau;
	mw.muDV = EvalMuVelocityGradient(mw);
	Tensorvec DV(pNo, Tensor());
	Vectorvec TaudivergenceoverRho(pNo, Vector());
	Vectorvec Muijstress(pNo, Vector());
	STensorvec GammaDot(pNo, STensor());
	Vectorvec TauDivergenceOverRho(pNo, Vector());
	Vectorvec TurbulTauDivergenceOverRho(pNo, Vector());
	doublevec neighwallcounter;


	mw.DelV = EvalVelocityGradient(mw);
	mw.Wallneighbor.assign(pNo, double());                      ////////////////Important
	mw.wallNormal.assign(pNo, Vector());
	mw.Walldistance.assign(pNo, double());



	for (register int i = 0; i < pNo; i++)

		mw.Walldistance[i] = 1;        //Initializing




	for (register int i = 0; i < iNo; i++)
	{
		const int ii = mw.interactions[i].GetI();
		const int jj = mw.interactions[i].GetJ();
		const short typei = mw.particles[ii].GetMaterial()->type;
		const short typej = mw.particles[jj].GetMaterial()->type;
		const Vector dx = mw.interactions[i].GetVD();
		const double d = mw.interactions[i].GetDist();
		const Vector dw = mw.interactions[i].GetDW() / d * dx; // Del W= W_ij/|r_ij|*r_ij	
		const Vector eij = 1.0 / d * dx;
		const double rhoij = 0.5 * (mw.particles[ii].GetRho() + mw.particles[jj].GetRho());
		const Vector dv = (mw.particles[ii].GetVel() - mw.particles[jj].GetVel());
		double mueffi;
		double Betha = 0.0000001;
		double MuEffectivei, MuEffectivej;
		double Etha = 0.1 * mw.particles[1].GetH();
		double hm = (mw.particles[ii].GetH() + mw.particles[jj].GetH()) / 2;

		if (typei < fluid1 && typej < fluid1)
			continue;

		/////////////////////////////////////////////////////Rheological properties		
		/*	YieldStress = 58;
		double n = 0.33;
		double Powerlawindex = -n + 1;
		double K = 35;   ////////////////Bates, Ancey 2016*/
		YieldStress = 58;
		double n = 0.33;
		double Powerlawindex = -n + 1;
		double K = 35;   ////////////////Bates, Ancey 2016
						 /////////////////////////////////////////////////////	
		double m = 10;


		GammaDot[ii] = mw.DelV[ii].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)
		GammaDot[jj] = mw.DelV[jj].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)


		double secondinvarianti = 2 * sqrt(Betha + GammaDot[ii].secondinvariant());   //GammaDot[i].secondinvariant();
		double secondinvariantj = 2 * sqrt(Betha + GammaDot[jj].secondinvariant());   //GammaDot[i].secondinvariant();
		double Critical_Mu = m * YieldStress;

		//	if (secondinvarianti < YieldStress / (m * K))
		//	MuEffectivei = m * YieldStress;
		//else
		//	MuEffectivei = YieldStress / (secondinvarianti)+K / pow(secondinvarianti, Powerlawindex);
		MuEffectivei = K / pow(secondinvarianti, Powerlawindex) + (YieldStress / secondinvarianti) * (1 - exp(-m * secondinvarianti));

		if (MuEffectivei > Critical_Mu)
			MuEffectivei = Critical_Mu;

		//	if (secondinvariantj < YieldStress / (m * K))
		//	MuEffectivej = m *YieldStress;
		//	else
		//	MuEffectivej = YieldStress / (secondinvariantj)+K / pow(secondinvariantj, Powerlawindex);
		MuEffectivej = K / pow(secondinvariantj, Powerlawindex) + (YieldStress / secondinvariantj) * (1 - exp(-m * secondinvariantj));
		if (MuEffectivej > Critical_Mu)
			MuEffectivej = Critical_Mu;


		if ((mw.particles[ii].GetMaterial()->type == 1 /*&& mw.particles[ii].GetLoc().GetY() >= 0.55*/) || (mw.particles[jj].GetMaterial()->type == 1 /*&& mw.particles[jj].GetLoc().GetY() >= 0.55*/))
			continue;



		///////////////////////////////////First form 
		//	Muijstress[ii] = Muijstress[ii] +2* mw.particles[jj].GetM() / (mw.particles[ii].GetRho()*mw.particles[ii].GetRho())*(MuEffectivei + MuEffectivej)*(Dot(dx, dw)*(dv) / (d*d /*+ Etha*Etha*/));
		//Muijstress[jj] = Muijstress[jj] -2* mw.particles[ii].GetM() / (mw.particles[ii].GetRho()*mw.particles[ii].GetRho())*(MuEffectivei + MuEffectivej)*(Dot(dx, dw)*(dv) / (d*d/* + Etha*Etha*/));
		///////////////////////////////////


		///////////////////////////////////Second form    Pasculii 2013
		double Muij = MuEffectivei * MuEffectivej / (MuEffectivei + MuEffectivej);
		Muijstress[ii] = Muijstress[ii] + 2 * Muij * mw.particles[jj].GetM() / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho()) * (Dot(dx, dv) * (dw) / (d * d + Etha * Etha));
		Muijstress[jj] = Muijstress[jj] - 2 * Muij * mw.particles[ii].GetM() / (mw.particles[ii].GetRho() * mw.particles[ii].GetRho()) * (Dot(dx, dv) * (dw) / (d * d + Etha * Etha));
		////////////////////////////////////////////////////////////





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
		///////////////////////////////////////////
	}



	for (int i = 0; i < pNo; i++)
	{
		mw.particles[i].SetAcc(mw.particles[i].GetAcc() + Muijstress[i]);



		/////////////
		YieldStress = 58;
		double n = 0.33;
		double Powerlawindex = -n + 1;
		double K = 35;   ////////////////Bates, Ancey 2016
						 /////////////////////////////////////////////////////	
		double m = 10000;
		double Betha = 0.0000001;



		GammaDot[i] = mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)
		double secondinvariant = 2 * sqrt(Betha + GammaDot[i].secondinvariant());   //GammaDot[i].secondinvariant();

		double MuEffective = K / pow(secondinvariant, Powerlawindex) + (YieldStress / secondinvariant) * (1 - exp(-m * secondinvariant));

		mw.particles[i].SetTau(MuEffective * GammaDot[i]);


	}

}



#endif // 0




