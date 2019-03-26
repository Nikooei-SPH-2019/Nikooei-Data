
#include "world.h"
#include "kernel.h"
#include "tensors3.h"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/itl/smoother/jacobi.hpp>"
using namespace mtl;
typedef matrix::compressed2D<double> SparseMatrix;  //Nikooei
typedef vector::dense_vector<double> DenseVector;   //Nikooei
 


void PressureSolver(World &mw)
{
	//Definitions/////////////////////////////////////////////
	int pNo = mw.particles.size();
	std::vector<double> sumpressure(pNo), sumnumber(pNo);   //Nikooei-Free surface
															//double si0 = 40000;            //Shao-non-newtonian dam
	bool Freesurfacecondition, Freesurfaceconditioni, Freesurfaceconditionj;
	double FlowDelDotr = 1.6, BedDelDotr = 1.6;


	std::vector<double> RHS2(pNo);
	std::vector<Vector> omegadelW(pNo);
	int ii, jj;
	double d;
	double siij, coef, sij, sii;
	int iNo = mw.interactions.size();
	double CriticalFluidneigh = 3;
	int K = 0;

	///////////////////////////////////////////////Definitions//
	
	//Finding Dry wall particles//////////////////////////////////

	for (register int i = 0; i < pNo; i++)
	{
		double hm;
		Vector dx;
		double dr;
		Vector normal;
		Vector Tangent;


		if (mw.particles[i].GetMaterial()->type < fluid1)
		{

			double denumnormal = sqrt(mw.gradSum[i].GetX()*mw.gradSum[i].GetX() + mw.gradSum[i].GetY()*mw.gradSum[i].GetY());
		
			normal.SetX(mw.gradSum[i].GetX() / denumnormal);
			normal.SetY(mw.gradSum[i].GetY() / denumnormal);
		
			
			int fluidneighbors = 0;


			if (mw.particles[i].Get_overlaying_flow_neigh() != 0)
				fluidneighbors = mw.particles[i].Get_overlaying_flow_neigh();
			else if (mw.particles[i].Get_bed_neigh() != 0)
				fluidneighbors = mw.particles[i].Get_bed_neigh();


			mw.particles[i].SetFluidneighbor(fluidneighbors);

		}

		//////////////////////////////////////////////////////////   Finding normal vector towards the fluid particle close to the wall
		double min_distance = 1;
		int index = -1;
		Vector NormalVec_fluid;
		if (mw.particles[i].GetMaterial()->type >= fluid1)
		{
			for (register int j = 0; j < pNo; j++)
			{
				if (mw.particles[j].GetMaterial()->type >= fluid1)
					continue;


				dx = mw.particles[i].GetLoc() - mw.particles[j].GetLoc();		//Distance Vector from j to i: ri-rj
				dr = dx.TwoNorm();

				if (dr < min_distance)
				{

					min_distance = dr;
					index = j;
				}
			}

			for (register int j = 0; j < pNo; j++)
			{
				if (j == index)
				{
					double denumnormal = sqrt(mw.gradSum[j].GetX()*mw.gradSum[j].GetX() + mw.gradSum[j].GetY()*mw.gradSum[j].GetY());
					NormalVec_fluid.SetX(mw.gradSum[j].GetX() / denumnormal);
					NormalVec_fluid.SetY(mw.gradSum[j].GetY() / denumnormal);
					mw.particles[i].SetNormalVec_fluid(NormalVec_fluid);

				}
			}
			//////////////////////////////////////////////////////////
		}
		
	}
	
	int T = 0;
	for (register int i = 0; i < pNo; i++)
	{
		Freesurfacecondition = (mw.particles[i].GetMaterial()->type != 15
			&& mw.particles[i].GetMaterial()->type >fluid1 && abs(mw.particles[i].GetDelDotr()) <= FlowDelDotr)
			|| (mw.particles[i].GetMaterial()->type == 15 && abs(mw.particles[i].GetDelDotr()) <= BedDelDotr);

	if ((mw.particles[i].GetMaterial()->type<fluid1  && mw.particles[i].GetFluidneighbor() < CriticalFluidneigh) || Freesurfacecondition)
			K++;
		else
		{
			mw.particles[i].SetNewIndex(T);
			mw.particles[i].SetOldIndex(i);
			T++;
		}
	}


	int NewpNo = pNo - K;

	SparseMatrix Smatrix(NewpNo, NewpNo);// = SparseMatrixGen(mw);
	DenseVector RHS(NewpNo), R0(NewpNo);// = DenseVectorGen(mw);
	DenseVector X(NewpNo, 0.0);
	
	//////////////////////////////////Finding Dry wall particles//


	//Construction of MATRIX////////////////////////////////////////////////
	{
		mtl::matrix::inserter<SparseMatrix> ins(Smatrix);
		std::vector<double> Ki(pNo);
		std::vector<double> DelDotr_Erosion(pNo);
		STensorvec GammaDot(pNo, STensor());;
		mw.gradSum_Erosion.assign(pNo, Vector());
		
		for (register int i = 0; i < iNo; i++)
		{

			ii = mw.interactions[i].GetI();
			jj = mw.interactions[i].GetJ();
			d = mw.interactions[i].GetDist();
			const Vector dx = mw.interactions[i].GetVD();
			const Vector eij = 1.0 / d*dx;
			const Vector dw = mw.interactions[i].GetDW()*eij;
			const double rhoij = 0.5*(mw.particles[ii].GetRho() + mw.particles[jj].GetRho());
			const STensor dwe = OProduct(dw, eij).GetSymPart();
			Vector normali; //= mw.gradSum[ii];
			Vector normalj; //= mw.gradSum[jj];
			Vectorvec DeldotTau(pNo, 0);


			double Fij, Fji, Sij, Sji, Deraieij, Deraieji;
			Vector Dij, Dji;
			double hm = (mw.particles[ii].GetH() + mw.particles[jj].GetH()) / 2;
			int typei = mw.particles[ii].GetMaterial()->type;
			int typej = mw.particles[jj].GetMaterial()->type;
			//////////Normals

			if (typei < fluid1)
			{
				double denumnormali = sqrt(mw.gradSum[ii].GetX()*mw.gradSum[ii].GetX() + mw.gradSum[ii].GetY()*mw.gradSum[ii].GetY());
			
				
					normali.SetX(mw.gradSum[ii].GetX() / denumnormali);
					normali.SetY(mw.gradSum[ii].GetY() / denumnormali);
			


			}



			////FSI    Fluid-Solid interaction   DelDotTau   Fatehi-2012
			double RhoAsquare = 1 / (mw.particles[ii].GetRho()*mw.particles[ii].GetRho());
			double RhoBsquare = 1 / (mw.particles[jj].GetRho()*mw.particles[jj].GetRho());
			if (typei < fluid1 )
			{
				DeldotTau[ii] = DeldotTau[ii] + mw.particles[jj].GetM() *((dw)*(mw.particles[ii].GetTau()* RhoAsquare +
					mw.particles[jj].GetTau()*RhoBsquare));


			}
			if (typej < fluid1 )
			{
				DeldotTau[jj] = DeldotTau[jj] - mw.particles[ii].GetM() *((dw)*(mw.particles[ii].GetTau()* RhoAsquare +
					mw.particles[jj].GetTau()*RhoBsquare));

			}
			

			if (typej < fluid1)
			{
				double denumnormalj = sqrt(mw.gradSum[jj].GetX()*mw.gradSum[jj].GetX() + mw.gradSum[jj].GetY()*mw.gradSum[jj].GetY());
			
			
					normalj.SetX(mw.gradSum[jj].GetX() / denumnormalj);
					normalj.SetY(mw.gradSum[jj].GetY() / denumnormalj);
			
				}
			//////////Normals
			
			Fij = (mw.particles[jj].GetM() / (rhoij*rhoij)) * (2.0*(1.0 / (d*d + 0.001*hm*hm))*Dot(dx, dw));
			Fji = (mw.particles[ii].GetM() / (rhoij*rhoij)) * (2.0*(1.0 / (d*d + 0.001*hm*hm))*Dot(dx, dw));

			Sij = Dot((mw.particles[jj].GetM() / (rhoij*rhoij)*dw), normali);
			Sji = Dot((mw.particles[ii].GetM() / (rhoij*rhoij)*dw), normalj);
			

			if (typei < fluid1 && mw.particles[ii].GetFluidneighbor() >= CriticalFluidneigh)
			{

				Deraieij = -Sij;
				Ki[ii] = Ki[ii] + Sij;
			}

			if (typej < fluid1 && mw.particles[jj].GetFluidneighbor() >= CriticalFluidneigh)
			{

				Deraieji = Sji;
				Ki[jj] = Ki[jj] - Sji;

			}

			Freesurfaceconditioni = (mw.particles[ii].GetMaterial()->type != 15 && mw.particles[ii].GetMaterial()->type >fluid1
				&& abs(mw.particles[ii].GetDelDotr()) <= FlowDelDotr)
				|| (mw.particles[ii].GetMaterial()->type == 15 && abs(mw.particles[ii].GetDelDotr()) <= BedDelDotr);

			Freesurfaceconditionj = (mw.particles[jj].GetMaterial()->type != 15 && mw.particles[jj].GetMaterial()->type >fluid1
				&& abs(mw.particles[jj].GetDelDotr()) <= FlowDelDotr)
				|| (mw.particles[jj].GetMaterial()->type == 15 && abs(mw.particles[jj].GetDelDotr()) <= BedDelDotr);

			if (typei >= fluid1 && !Freesurfaceconditioni)
			{
				Deraieij = -Fij;

				Ki[ii] = Ki[ii] + Fij;


			}
			if (typej >= fluid1 && !Freesurfaceconditionj)
			{
				Deraieji = -Fji;

				Ki[jj] = Ki[jj] + Fji;

			}

			//Insertion/////////////////////////////////////
		

			if (Freesurfaceconditioni || (mw.particles[ii].GetMaterial()->type < fluid1 &&  mw.particles[ii].GetFluidneighbor() < CriticalFluidneigh)
				|| Freesurfaceconditionj || (mw.particles[jj].GetMaterial()->type < fluid1 &&  mw.particles[jj].GetFluidneighbor() < CriticalFluidneigh))
				continue;

						
			ins(mw.particles[ii].GetNewIndex(), mw.particles[jj].GetNewIndex()) << Deraieij;
			ins(mw.particles[jj].GetNewIndex(), mw.particles[ii].GetNewIndex()) << Deraieji;
	
			omegadelW[ii] += mw.particles[jj].GetM() / (rhoij)*dw;
			omegadelW[jj] -= mw.particles[ii].GetM() / (rhoij)*dw;

			sumpressure[ii] = sumpressure[ii] + mw.particles[jj].GetP();                                                     //Nikooei-Free surface
			sumpressure[jj] = sumpressure[jj] + mw.particles[ii].GetP();

			sumnumber[ii] = sumnumber[ii] + 1;
			sumnumber[jj] = sumnumber[jj] + 1;


			////////////////////////////////////////////////////
			
			double Betha = 0.0000001;
			double Rhos = 1300;
			double effectivPi = abs(mw.particles[ii].GetP());
			double effectivPj = abs(mw.particles[jj].GetP());

			GammaDot[ii] = 2 * mw.DelV[ii].GetSymPart();   ///Strain rate tensor (Viscous fluid white 2006)    
			GammaDot[jj] = 2 * mw.DelV[jj].GetSymPart();   ///Strain rate tensor (Viscous fluid white 2006)    

			double Substrate_internalangle = 0.17453292519943295769236907684886;    //10   Crosta


			double secondinvarianti = sqrt(Betha + GammaDot[ii].secondinvariant());
			double secondinvariantj = sqrt(Betha + GammaDot[jj].secondinvariant());

			double InertialNumi = secondinvarianti*d / sqrt(effectivPi / Rhos);
			double InertialNumj = secondinvariantj*d / sqrt(effectivPj / Rhos);
			double Substrate_Critical_inertialnumberi = sqrt(abs(mw.particles[ii].GetP())*Rhos)*tan(Substrate_internalangle)*d / (1000 * 0.1);
			double Substrate_Critical_inertialnumberj = sqrt(abs(mw.particles[jj].GetP())*Rhos)*tan(Substrate_internalangle)*d / (1000 * 0.1);

			if (/*typei == 15 && typej == 15 && */mw.particles[ii].GetVel().TwoNorm()<0.1 &&
				mw.particles[jj].GetVel().TwoNorm()<0.1&& ii != jj && mw.particles[ii].GetP() != 0 && mw.particles[jj].GetP() != 0)
			{
				DelDotr_Erosion[ii] = DelDotr_Erosion[ii] + mw.particles[jj].GetM() / rhoij*Dot(dx, dw);
				DelDotr_Erosion[jj] = DelDotr_Erosion[jj] + mw.particles[ii].GetM() / rhoij*Dot(dx, dw);


				if (mw.particles[ii].GetP() != 0 && mw.particles[jj].GetP() != 0)
				{
					mw.particles[ii].Set_Static_state(1);
					mw.particles[jj].Set_Static_state(1);

				}

				
				mw.gradSum_Erosion[ii] += mw.particles[jj].GetM() / rhoij *dw;
				mw.gradSum_Erosion[jj] -= mw.particles[ii].GetM() / rhoij *dw;

			}





			}

		////////////////////////////////////////////////////////////RHS and Diagonal entries//////////////////////////////////////////////////////////
		int counter = 0;
		for (int i = 0; i < pNo; i++)
		{

			mw.particles[i].SetDelDotr_Erosion(DelDotr_Erosion[i]);


			Freesurfacecondition = (mw.particles[i].GetMaterial()->type != 15
				&& mw.particles[i].GetMaterial()->type>fluid1 && abs(mw.particles[i].GetDelDotr()) <= FlowDelDotr)
				|| (mw.particles[i].GetMaterial()->type == 15 && abs(mw.particles[i].GetDelDotr()) <= BedDelDotr);

			if ((mw.particles[i].GetMaterial()->type < fluid1  && mw.particles[i].GetFluidneighbor() < CriticalFluidneigh) || Freesurfacecondition)
				continue;

			
			int typei = mw.particles[i].GetMaterial()->type;

			ins(counter, counter) << Ki[i];

			double dt = mw.GetDT();

			/////////////////////////////////////////////////////////////////////////////////////RHS
			if (typei >= fluid1)
				RHS[counter] = mw.particles[i].GetVDivergence() / dt;

			else if (typei < fluid1)
			{
				Vector normali;
				double denumnormali = sqrt(mw.gradSum[i].GetX()*mw.gradSum[i].GetX() + mw.gradSum[i].GetY()*mw.gradSum[i].GetY());
			
				if (typei != wall2)
				{
					normali.SetX(mw.gradSum[ii].GetX() / denumnormali);
					normali.SetY(mw.gradSum[ii].GetY() / denumnormali);
				}
				else
				{
					double Tetha = 0 * 3.1415926535897932384626433832795 / 180;
					normali.SetX(-cos(Tetha));
					normali.SetY(sin(Tetha));
				}
				
				Vector RHS_Wall = mw.options.BODYACC ;
				RHS[counter] = Dot(-RHS_Wall, normali);

			}

			counter++;
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////
		}
	////////////////////////////////////////////////Construction of MATRIX//


	itl::basic_iteration<double> iter(RHS, 100, 1.e-10, 0);
	
	/////////////////////////////////////////////////////////////////////////////////////////Solver
	itl::pc::ilu_0<SparseMatrix> PCILU0(Smatrix);
	itl::bicgstab(Smatrix, X, RHS, PCILU0, iter);


	
	for (register int i = 0; i < pNo; i++)
	{


		Freesurfaceconditioni = (mw.particles[i].GetMaterial()->type != 15 && mw.particles[i].GetMaterial()->type>fluid1 && abs(mw.particles[i].GetDelDotr()) <= FlowDelDotr)
			|| (mw.particles[i].GetMaterial()->type == 15 && abs(mw.particles[i].GetDelDotr()) <= BedDelDotr);

		if ((mw.particles[i].GetMaterial()->type < fluid1 && mw.particles[i].GetFluidneighbor() < CriticalFluidneigh) || Freesurfaceconditioni)
			mw.particles[i].SetP(0);
		else
			mw.particles[i].SetP(X[mw.particles[i].GetNewIndex()]);


	}
	


}


void Projection3(World &mw)     /// Projection- New Wall boundary condition - moving boundary
{

	int pNo = mw.particles.size();

	Vector Loc_all;

	int count_par = 0;


	int m = 0;


	std::vector<Vector> vold(pNo), locold(pNo), v_bar(pNo), vnew(pNo);     //Nikooei vnew(pNo)

	double dt = mw.GetDT();

	mw.FindInteractions();   /// FINDS INTERACTIONS
	mw.UpdateCorrectors();   /// CALCULATE MATRIX B

	int iNo = mw.interactions.size();
	double d;
	Vector dx, dw;
	
	
	
	mw.EvalMomentum();   /// CALCULATE THE NEWTONIAN VISCOSITY 
	for (register int i = 0; i < pNo; i++)
	{
		vold[i] = mw.particles[i].GetVel();				//Save old time level velocity and positions
		locold[i] = mw.particles[i].GetLoc();
		if (mw.particles[i].GetMaterial()->type < fluid1 || (mw.particles[i].GetMaterial()->type >= 15 && mw.particles[i].GetMu_effective() == 1000))
			continue;
		mw.particles[i].SetVel(mw.particles[i].GetVel() + mw.particles[i].GetAcc()*dt);   /// V(star)
		
	}
	
	mw.UpdateInteractions();					// Update Values for intermediate positions
	mw.EvalDivergence3();						//Evaluate velocity divergence of each particle


	PressureSolver(mw);

	PressureForceSymmetric(mw);

	
	STensorvec GammaDot(pNo, STensor());
	
	for (register int i = 0; i < pNo; i++)								//New Velocity n+1
	{

		GammaDot[i] = 2 * mw.DelV[i].GetSymPart();                              ///Strain rate tensor (Viscous fluid white 2006)
		double Betha = 0.0000001;

		double secinvariant = sqrt(Betha + GammaDot[i].secondinvariant());   //Ga

		mw.particles[i].SetGammaDot(secinvariant);
		
			
		if (mw.particles[i].GetMaterial()->type < fluid1||(mw.particles[i].GetMaterial()->type >= 15 && mw.particles[i].GetMu_effective() == 1000))
			continue;

		////////////////////////////////////////////////Partial-Slip BC   -Start


		if (mw.particles[i].GetMaterial()->type == 11 && mw.Walldistance[i] < mw.particles[i].GetH() / 2.2 && mw.GetSTEP()>1 && mw.particles[i].GetLoc().GetX()<0.07)
		{
			double Ls = 0.018;                ///Slip length(m)=2*particle diameter
			

			Vector Tang_Traction, Slip_Vel;
			STensor Identity_tensor;

			Identity_tensor.SetAll(1.0, 0.0, 1.0);


			STensor Stress = mw.particles[i].GetTau();
			Vector Normal_i = mw.particles[i].GetNormalVec_fluid();

		STensor Tang_proj_op = Identity_tensor - OProduct_same_vec(Normal_i, Normal_i);
			Tang_Traction = Tang_proj_op*(Stress)*Normal_i;
			Slip_Vel = (Ls / mw.particles[i].GetMu_effective()) *Tang_Traction;

			mw.particles[i].SetVel(Slip_Vel);
		}

		else
			
		mw.particles[i].SetVel(vold[i] + mw.particles[i].GetAcc()*dt);

	}
	////////////////////////////////////////////////Partial-Slip BC -End  

	AddXSPH(mw);

	Vector vel_old = mw.GetVELOCITY();
	

	double L_0 = 0, H_0 = 0;/////////////Initial position of Gate in collapse
	std::ofstream Data_Ini;

	Vector Radi, Norm_here;
	for (int i = 0; i<pNo; i++)
	{
	
		if ((mw.particles[i].GetMaterial()->type == 11) && (mw.particles[i].GetLoc().GetX() > L_0) && mw.GetSTEP() == 1)
		{
			L_0 = mw.particles[i].GetLoc().GetX();
		}
		if ((mw.particles[i].GetMaterial()->type == 11) && (mw.particles[i].GetLoc().GetY() > H_0) && mw.GetSTEP() == 1)
		{
			H_0 = mw.particles[i].GetLoc().GetY();
		}
	}
	if (mw.GetSTEP() == 1)
	{
		Data_Ini.open("Initial_Data.txt", std::ofstream::out);

		Data_Ini << L_0 << '\n';
		Data_Ini << H_0;
		Data_Ini.close();

	}

	std::ifstream A;
	double initial_Height, initial_length;
		
		A.open("Initial_Data.txt");
	while (!A.eof())
	{
		A >> initial_length>> initial_Height ;
	}
	

	double Threshold_V = 0.07*sqrt(9.81*initial_Height);

	////// CALCULATE FORCES AND TORSION on the SOLID BODY ... and NEW LOCATION OF SOLID BODY - END ///

	double r = mw.minimum.loc.GetY();
	

	std::ifstream Infile2;
	double L_I, H_Max;
	Infile2.open("Initial_Data.txt");
	if (!Infile2)
	{
		std::cout << "Unable to open file datafile.txt";
	}


	while (!Infile2.eof())
	{
		Infile2 >> L_I >> H_Max;
	}
	
	Infile2.close();
	
	double ts = mw.GetTime() / mw.options.SAVINGTIME;

		///////////////////////////////////////////////////////////////////////////////////////Extraction of Local Data from flow within probing sections -Start
	if ((ts - (int)ts) > (1.0 - (mw.GetDT() / mw.options.SAVINGTIME)))
	{

		const int M = 25;           //////////////////Number of Vertical probes

		const int N = 20;           //////////////////Number of Horizontal probes
		double DeltaX = L_I;              ///////////probe distances
		double DeltaY;
		

		if (mw.minimum.loc.GetY() == 0)
			DeltaY = H_Max / M;              ///////////probe distances within the initial mass (non-erodible case)
		else
			DeltaY = abs(mw.minimum.loc.GetY()) / 3;              ///////////probe distances within the erodible mass (erodible case)




		double H;                     //////Maximum H in a position (Free surface profile)
		doublevec KE(N + 1, 0);
		doublevec XMom(N + 1, 0);
		doublevec Mass(N + 1, 0);
		doublevec MaxY(N + 1, 0);   ///////////////////////Maximum height in a fixed location
		doublevec VX_total(N + 1, 0);//////////////////Total & Averaged X_Velocity in a location
		doublevec VY_total(N + 1, 0);//////////////////Total & Averaged  Y_Velocity in a location
		doublevec V_total(N + 1, 0);//////////////////Total & Averaged Velocity in a location
		doublevec P_total(N + 1, 0);//////////////////Total & Averaged Pressure in a location
		doublevec V_max(N + 1, 0);
		double E_X = mw.minimum.h / 2.6;
		double E_Y = mw.minimum.h / 2.6;
		intvec O(N + 1, 0);
		double KE2[N + 1][M], XMom2[N + 1][M], Mass2[N + 1][M];
		int X_sec_index;
		int Y_sec_index;
		double  X_sec, Y_sec;
		double X_f_flow = 0;
		double X_f_bed = 0;


		/////////////////initialize
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				KE2[i][j] = 0;
				XMom2[i][j] = 0;
				Mass2[i][j] = 0;

			}

		for (int j = 0; j < pNo; j++)
		{


			if (mw.particles[j].GetMaterial()->type != 0)
			{
				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				Vector NewLoc = mw.particles[j].GetLoc() + (0.5*vold[j] + 0.5*mw.particles[j].GetVel())*dt;

				double NewX = NewLoc.GetX();
				if (NewX >= L_I)
				{
					int e = int((NewX - L_I) / DeltaX);
					X_sec_index = min(e, N);
				}
				else
					X_sec_index = 0;

				X_sec = X_sec_index*DeltaX + L_I;  //i=0 ----> X_sec=L_I , i=1 ----> X_sec=L_I+DeltaX 

				Y_sec_index = min(int((mw.particles[j].GetLoc().GetY() - mw.minimum.loc.GetY()) / DeltaY), M);

				Y_sec = Y_sec_index*DeltaY + mw.minimum.loc.GetY();  //i=0 ----> X_sec=L_I , i=1 ----> X_sec=L_I+DeltaX 
																	 ////////////////////////////////////////////////////////////////////////////////////////////////////////

				if (NewX <= (X_sec + E_X) && NewX >= (X_sec - E_X))
				{
					KE[X_sec_index] = KE[X_sec_index] + 0.5*mw.particles[j].GetM()*(mw.particles[j].GetVel().TwoNorm())*(mw.particles[j].GetVel().TwoNorm());
					XMom[X_sec_index] = XMom[X_sec_index] + mw.particles[j].GetM()*(mw.particles[j].GetVel().GetX());
					Mass[X_sec_index] = Mass[X_sec_index] + mw.particles[j].GetM();

					if (mw.particles[j].GetLoc().GetY() < (Y_sec + E_Y) && mw.particles[j].GetLoc().GetY() > (Y_sec - E_Y))
					{
						KE2[X_sec_index][Y_sec_index] = KE2[X_sec_index][Y_sec_index] + 0.5*mw.particles[j].GetM()*(mw.particles[j].GetVel().TwoNorm())*(mw.particles[j].GetVel().TwoNorm());
						XMom2[X_sec_index][Y_sec_index] = XMom2[X_sec_index][Y_sec_index] + mw.particles[j].GetM()*(mw.particles[j].GetVel().GetX());
						Mass2[X_sec_index][Y_sec_index] = Mass2[X_sec_index][Y_sec_index] + mw.particles[j].GetM();
					}
					/////////////////////////////////////////////////////////
					if (mw.particles[j].GetLoc().GetY() > MaxY[X_sec_index])
					{
						MaxY[X_sec_index] = mw.particles[j].GetLoc().GetY();
					}

					VX_total[X_sec_index] = VX_total[X_sec_index] + mw.particles[j].GetVel().GetX();
					VY_total[X_sec_index] = VY_total[X_sec_index] + mw.particles[j].GetVel().GetY();
					V_total[X_sec_index] = VX_total[X_sec_index] + mw.particles[j].GetVel().TwoNorm();
					P_total[X_sec_index] = P_total[X_sec_index] + mw.particles[j].GetP();
					O[X_sec_index] = O[X_sec_index] + 1;

					if (mw.particles[j].GetVel().TwoNorm() > V_max[X_sec_index])
						V_max[X_sec_index] = mw.particles[j].GetVel().TwoNorm();
					/////////////////////////////////////////////////////////////////

				}

				std::ifstream X;
				double initial_Height, initial_length;

				X.open("Initial_Data.txt");
				while (!X.eof())
				{
					X >> initial_length >> initial_Height;
				}

				double Threshold_V = 0.07*sqrt(9.81*initial_Height);



				////////////////////////////////////////////////////////////////Runout
				double type = mw.particles[j].GetMaterial()->type;
				double X_loc = mw.particles[j].GetLoc().GetX();
				if (type == 11 && X_loc >= X_f_flow)
				{
					X_f_flow = X_loc;

				}
				double V_m = mw.particles[j].GetVel().TwoNorm();
				if (type == 15 && X_loc >= X_f_bed && mw.particles[j].GetVel().TwoNorm() >= Threshold_V && X_loc<(0.95*mw.maximum.loc.GetX()))
				{
					X_f_bed = X_loc;
				}
			////////////////////////////////////////////////////////////////Runout


			}
		}

		double KE_Flux2[N][M];
		double Mom_Flux2[N][M];
		double Mass_Flux2[N][M];
		double VX_Ave;
		double VY_Ave;
		double V_Ave;
		double P_Ave;
		double Max_H;
		double Froude;


		std::string name = "type1_data_X_";     //base pattern of file name
		std::ofstream outfstr[N];      //creating array of N output file streams
	//	std::string name2 = "type2_data_X_";     //base pattern of file name
		std::ofstream outfstr2[N];      //creating array of N output file streams*/
		for (int i = 0; i < (N); i++)
		{   //open all file streams 

			outfstr[i].open(name + char('0' + i) + ".dat", std::ofstream::out | std::ofstream::app);
			//outfstr2[i].open(name2 + char('0' + i) + ".dat", std::ofstream::out | std::ofstream::app);

			if (mw.GetTime() < (1.1* mw.options.SAVINGTIME))			 ///////////////File Header
			{
				outfstr[i] << "Variables= t   Probe_X   KE  XMom Mass Maximum_H  VX_Ave  VY_Ave  V_Ave Froude  P_Ave  V_max" << '\n';
				outfstr[i] << "Zone " << '\n';

			}														///////////////////////////
			outfstr2[i] << "Variables= t   Probe_X  Probe_Y  KE  XMom  Mass" << '\n';
			outfstr2[i] << "Zone " << '\n';


			VX_Ave = VX_total[i] / O[i];
			VY_Ave = VY_total[i] / O[i];
			V_Ave = V_total[i] / O[i];
			P_Ave = P_total[i] / O[i];
			Max_H = MaxY[i];
			if (Max_H != 0)
				Froude = V_Ave / sqrt(9.81*Max_H);
			else
				Froude = 0;

			if (O[i] == 0)           ////////////////These parameters are nil until the flow front reaches the position
			{
				VX_Ave = 0;
				VY_Ave = 0;
				V_Ave = 0;
				P_Ave = 0;

			}

			
			outfstr[i] << mw.GetTime() << '\t';
			outfstr[i] << (i*DeltaX + L_I) << '\t';
			outfstr[i] << KE[i] << '\t';
			outfstr[i] << XMom[i] << '\t';
			outfstr[i] << Mass[i] << '\t';
			outfstr[i] << Max_H << '\t';
			outfstr[i] << VX_Ave << '\t';
			outfstr[i] << VY_Ave << '\t';
			outfstr[i] << V_Ave << '\t';
			outfstr[i] << Froude << '\t';
			outfstr[i] << P_Ave << '\t';
			outfstr[i] << V_max[i] << '\n';


		
	}
		
		std::ofstream X_f;     //  
		X_f.open("Runout.dat", std::ofstream::out | std::ofstream::app);

		if (mw.GetTime() < (1.1* mw.options.SAVINGTIME))			 ///////////////File Header
		{
			X_f << "Variables = t    X_initial_mass  X_bed   t_Norm   X_initial_mass_Norm  X_bed_Norm   " << '\n';
			X_f << "Zone " << '\n';
		}

		X_f << mw.GetTime() << '\t';
		X_f << X_f_flow << '\t';
		X_f << X_f_bed << '\t';
		X_f << mw.GetTime() / sqrt(initial_Height / 9.81) << '\t';
		X_f << (X_f_flow - initial_length) / initial_length << '\t';
		X_f << (X_f_bed - initial_length) / initial_length << '\n';

		X_f.close();


	}

	////////////////////////////////////////////////////////////////////////////Extraction of Local Data from flow -End
	




	//// CALCULATE THE NEW LOCATION OF THE PARTICLES
	for (register int i = 0; i<pNo; i++)							//New Location n+1
	{

		mw.particles[i].SetLoc(mw.particles[i].GetLoc() + (0.5*vold[i] + 0.5*mw.particles[i].GetVel())*dt);

	}

	//// SHIFTING THE PARTICLES

	mw.options.SHIFT = mw.options.Shifting;
	if (mw.options.SHIFT)				//Shift particles slightly to avoid interference
		ShiftParticles3(mw);
	
}
