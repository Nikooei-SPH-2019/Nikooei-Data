//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change	86/05 By Fatehi
// Last change	86/09 By Fatehi
// Last change	87/11 By Fatehi
// changed		89/03 By Fatehi
////// Last change	97/11 By Nikooei
//****************************************************************************************
// Particle header file
//	Defines Class Particle & enumerator Particle_Type
//****************************************************************************************


#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "vectors.h"
#include "stensors.h"
#include "tensors.h"
#include "tensors4.h"
#include "dtensors.h"
#include "material.h"

enum Particle_Type
{	
	//dummy=-1,
	wall1=				0,	wall2,
	solid1=				5,	solid2,
	fluid1=				10, BinghamCoulomb,Binghamcross, HerschelBulkley, Thixotropy ,      //Nikooei-Non-newtonian

};


//////////////////////////////////////////////////////////////////////////////
// C  L  A  S  S     D  E  C  L  E  R  A  T  I  O  N
// Class:       TGParticle
// Purpose:     Ultimate general particle including position only
// Author:      R. Fatehi
// Revision:    1.0
// Orig. Date:  Sept. 2007
// Rev. Date:   March. 2019
//////////////////////////////////////////////////////////////////////////////
class TGParticle
{
protected:

	Vector loc;		//Location

public:
	//***********************************************************BY MAFB & MSS
	int mainPIndex;   //A holder of the main particle number
					  //Could be added to TwallBCVector 
	TGParticle(): mainPIndex(-1) {}

	TGParticle(Vector r ):loc(r), mainPIndex(-1) {}
	//***********************************************************End BY MAFB & MSS

	void SetLoc(Vector l){loc=l;}

	double Distance(TGParticle &pa){Vector d=loc-pa.loc;return d.TwoNorm();}

	Vector GetLoc(){return loc;}

};


//////////////////////////////////
/////////// BY SINA JAHANGIRI
///////////////////////////////////


class IbInteraction: public TGParticle
{
	
public:
	
	int index;
	Vector DX;
	double DR;
		
	IbInteraction() {}
		
	//int GetSize(){return IbInteraction.size();}
	int GetIndex(){return index;}
	
	void SetIndex (int PN){index=PN;}
	
	void SetDR(double DRR){DR=DRR;}
	double GetDR() {return DR;}

	void SetDX(Vector DXX){DX=DXX;}
	Vector GetDX() {return DX;}
	
	
};

class ParInteraction: public TGParticle
{


public:
	
	int index;
	double Aij;
	
	ParInteraction() {}
	
	//int GetSize(){return IbInteraction.size();}
	int GetIndex(){return index;}
	void SetIndex (int PN){index=PN;}

	void SetAij(double DRR){Aij=DRR;}
	double GetAij() {return Aij;}
	
};
/////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// C  L  A  S  S     D  E  C  L  E  R  A  T  I  O  N
// Class:       TIParticle
// Purpose:     General image particle including position and an index refers
//				to the real particle
// Author:      R. Fatehi
// Revision:    1.0
// Orig. Date:  Sept. 2007
// Rev. Date:   March. 2019
//////////////////////////////////////////////////////////////////////////////
class TIParticle: public TGParticle
{
public:
	int index ;
	
	Vector Force; // By S.J
	Vector Velocity; // Velocity --- for IB : By S.Jahangiri
	double Omega; // By S.J
	Vector normal;		//normal
	double Mass; // By S.J
	Vector Arm; // By S.J : Arm to Center particle!
	Vector TmpVel; // By S.J : Arm to Center particle!
	Vector TmpLoc; // By S.J : Temp Location of particles when cross periodic boundaries
	double P; /// By S.J
	STensor Tau; // By S.J
	
	std::vector<IbInteraction> Interactions; // by Sina Jahangiri

	TIParticle() {}
	TIParticle( int i, Vector r , Vector V, double omega): index(i), Velocity(V),Omega(omega), TGParticle(r){}   // By S.Jahangiri
	TIParticle( int i, Vector r ): index(i), normal(), TGParticle(r){}
	TIParticle( int i, Vector r , Vector n ): index(i), normal(n), TGParticle(r){}
	TIParticle( int i, Vector r , Vector n , Vector V): index(i), normal(n), Velocity(V), TGParticle(r){}   // By S.Jahangiri
	
	void SetForce(Vector f){Force=f;}
	Vector GetForce(){return Force;}
	void SetVel(Vector v){Velocity=v;}
	Vector GetVel(){return Velocity;}
	Vector GetNorm(){return normal;}
	void SetNorm(Vector n){normal=n;}
	double GetNormalSize(){return normal.TwoNorm() ;}
	void SetOmega(double om){Omega=om;}
	double GetOmega(){return Omega;}
	void SetMass(double m){Mass=m;}
	double GetMass(){return Mass;}
	Vector GetArm() {return Arm;}
	void SetArm(Vector ARM){Arm=ARM;}
	Vector GetTmpVel() {return TmpVel;}
	void SetTmpVel(Vector TMPARM){TmpVel=TMPARM;}
	Vector GetTmpLoc() {return TmpLoc;}
	void SetTmpLoc(Vector TMPLOC){TmpLoc=TMPLOC;}
	double GetP(){return P;}
	void SetP (double Pres){P=Pres;}
	STensor GetTau() {return Tau;}
	void SetTau(STensor TAU){Tau=TAU;}


	double Distance(TIParticle &pa){Vector d=loc-pa.loc;return d.TwoNorm();}

};
//////////////////////////////////////////////////////////////////////////////
// C  L  A  S  S     D  E  C  L  E  R  A  T  I  O  N
// Class:       TRParticle
// Purpose:     General real particle including position, velocity, mass, h, 
//				number density, stress, acc and material.
// Author:      R. Fatehi
// Revision:    1.0
// Orig. Date:  Sept. 2007
// Rev. Date:   March. 2019
//////////////////////////////////////////////////////////////////////////////
class TRParticle: public TGParticle
{
private:
	TRMaterial* mat;	//Material
protected:
	int BCindex;

	Vector vel;		//Velocity
	double m;		//Mass
	double h;		//Smooting length

	double si;		//Particle number density
	double si_0;    //Initial Particle number density -Nikooei
	Vector acc;		//Acceleration

	int Interface_num;    //Interface number for Flow-bed interface


	STensor tau;	//Stress
	STensor Turbultau;      //Turbulent viscosity                            // Nikooei
	double rho;		//(variable) density 
	double Head;
	int fluidneighbors;             ///////////////////////////////////////////Nikooei
	int Index;
	int OldIndex;
	double DelDotr;
	Vector tmp_vel; // Temp Velocity  : S.Jahangiri
	Vector tmp_acc; // Temp Acceleration : S.Jahangiri
	double Mu_effective;
	Vector NormalVec_fluid;
	
public:

	TRParticle(){}

	TRParticle(Particle_Type , Vector r, double mass, double hw, Vector v, STensor s)
				:BCindex(0), vel(v), m(mass), h(hw), si(0.0), acc(), tau(s), TGParticle(r){}

	void SetBCindex(int i){BCindex=i;}

	void SetVel(Vector v){vel=v;}
	void SetM(double mass){m=mass;}
	void SetH(double sl){h=sl;}

	void SetSi(double s){si=s;}
	void SetSi_0(double si0) { si_0 = si0; }  //Si in time=0
	void SetAcc(Vector ss){acc=ss;}

	void SetMaterial(TRMaterial* m){mat=m;}
	void SetTau(STensor t){tau=t;}
	void SetTurbulTau(STensor t) { Turbultau = t; }

	void SetNormalVec_fluid(Vector  NVF) { NormalVec_fluid = NVF; }
	Vector GetNormalVec_fluid() { return NormalVec_fluid; }
	
	
	void SetRho(double r){rho=r;}

	double Distance(TRParticle &pa){Vector d=loc-pa.loc;return d.TwoNorm();}

	int GetBCindex(){return BCindex;}

	Vector GetVel(){return vel;}
	double GetM(){return m;}
	double GetH(){return h;}

	void SetMu_effective(double mueff) { Mu_effective = mueff; }
	double GetMu_effective() { return Mu_effective; }

	void SetInterfaceneighbor(int a) { Interface_num = a; }
	int GetInterfaceneighbor() { return Interface_num; }

	double GetSi(){return si;}
	double GetSi_0() { return si_0; }
	

	double GetRho(){return rho;}

	double GetHead(){ return Head; }
	void SetHead(double H) { Head = H; }
	
	
	void SetNewIndex(int ind) { Index = ind; }   /////////////////////////////Nikooei
	int GetNewIndex() { return Index; }   /////////////////////////////Nikooei

	void SetOldIndex(int a) { OldIndex = a; }   /////////////////////////////Nikooei
	int GetOldIndex() { return OldIndex; }   /////////////////////////////Nikooei
	
	void SetDelDotr(double y) { DelDotr = y; }   /////////////////////////////Nikooei
	float GetDelDotr() { return DelDotr; }   /////////////////////////////Nikooei
	
	Vector GetAcc(){return acc;}

	TRMaterial* GetMaterial(){return mat;}
	STensor GetTau(){return tau;}
	STensor GetTurbulTau(){ return Turbultau; }                            //By Nikooei

	friend std::ostream &operator<<(std::ostream &stream, TIParticle &par);


};
	std::ostream &operator<<(std::ostream &stream, TIParticle &par);

//////////////////////////////////////////////////////////////////////////////
// C  L  A  S  S     D  E  C  L  E  R  A  T  I  O  N
// Class:       TFParticle
// Purpose:     General real particle including position, velocity, mass, h, 
//				number density, stress, acc and material.
// Author:      R. Fatehi
// Revision:    1.0
// Orig. Date:  Sept. 2007
// Rev. Date:   March. 2019
//////////////////////////////////////////////////////////////////////////////
class TFParticle: public TRParticle
{
protected:
	TFluid* mat;	//Material
	TNonNewtonian* Nonnewtonmat;            //////////////////////Nikooei-Non-Newtonian

	double p;		//Pressure
	double TotalP;   //Total P

	double RHS2;

	double u;		//Internal energy
	double du;		//Time differential for Internal energy

	double vdivergence;	//Velocity divergence
	int fluidneighbors;             ///////////////////////////////////////////Nikooei
	int Index;
	int OldIndex;
public:
	double p1;		//additive pressure for bulk viscosity effect.

	STensor B;		//Corrector1 for kernel gradient.
	STensor Bhat;		//Corrector2 for kernel gradient in second derivatives.
	//Tensor4 Bhat4;		//Corrector2 for kernel gradient in second derivatives.
	double KE;
	double GammaDot;
	double Vorticity;
	double Okubo_Weiss;
	int Static_state=0;   ///For erosible beds

	double freesurfacetracking;
	double DelDotr;                          /////////////////nikooei-Stansby 2012
	double DelDotr_Erosion;
	std::vector<int> AllInter; // By S.J
	std::vector<ParInteraction> Interactions;   /// By Sina Jahangiri
	int overlaying_flow_neigh, bed_neigh;
	TFParticle() {}

	TFParticle(Particle_Type p_t, Vector r, double mass, double hw, Vector v, double pres=0.0, STensor s=STensor(), double ui=0 )
				: p(pres), u(ui), du(0.0), vdivergence(0.0), TRParticle( p_t, r, mass, hw, v, s ){}
	
	void SetMaterial(TFluid* m){mat=m;}
	void SetMaterial(TNonNewtonian* m) { Nonnewtonmat = m; }               //////////////////////Nikooei-Non-Newtonian

	void SetP(double press){p=press;}
	void SetU(double inte){u=inte;}
	void SetDU(double dinte){du=dinte;}

	void SetVDivergence(double vdiv){vdivergence=vdiv;}

	void SetFluidneighbor(int FN) { fluidneighbors = FN; }   /////////////////////////////Nikooei
	int GetFluidneighbor() { return fluidneighbors; }   /////////////////////////////Nikooei

	void SetGammaDot(double Ga) { GammaDot = Ga; }
	double GetGammaDot() { return GammaDot; }


	void Set_overlaying_flow_neigh(int fn) { overlaying_flow_neigh = fn; }/////////////////////////////Nikooei
	int Get_overlaying_flow_neigh() { return overlaying_flow_neigh; }


	void Set_bed_neigh(int bn) { bed_neigh = bn; }
	int Get_bed_neigh() { return bed_neigh; }/////////////////////////////Nikooei


	void SetNewIndex(int ind) { Index = ind; }   /////////////////////////////Nikooei
	int GetNewIndex() { return Index; }   /////////////////////////////Nikooei

	void SetOldIndex(int a) { OldIndex = a; }   /////////////////////////////Nikooei
	int GetOldIndex() { return OldIndex; }   /////////////////////////////Nikooei

	void SetDelDotr(double y) { DelDotr = y; }   /////////////////////////////Nikooei
	float GetDelDotr() { return DelDotr; }   /////////////////////////////Nikooei

	void SetDelDotr_Erosion(double n) { DelDotr_Erosion = n; }   /////////////////////////////Nikooei
	float GetDelDotr_Erosion() { return DelDotr_Erosion; } 

	void SetVorticity(double w)
	{
		Vorticity = w;
	};
	double GetVorticity() { return Vorticity; };
	
	void Set_Okubo_Weiss(double r)
	{
		Okubo_Weiss = r;
	};
	double Get_Okubo_Weiss() { return Okubo_Weiss; };

	void Set_Static_state(int l)
	{
		Static_state = l;
	};
	double Get_Static_state() { return Static_state; };
	
	double Distance(TFParticle &pa){Vector d=loc-pa.loc;return d.TwoNorm();}

	TFluid* GetMaterial(){return mat;}
	TNonNewtonian* GetNonnewtonMaterial() { return Nonnewtonmat; }//////////////////////Nikooei-Non-Newtonian

	double GetP(){return p;}
	double GetTotalP() { return TotalP; }
	void SetTotalP(double totP) { TotalP = totP; }
	double GetU(){return u;}
	double GetDU(){return du;}

	double GetVDivergence(){return vdivergence;}

	double GetRHS2() { return RHS2; }
	double SetRHS2(double RightHand2) { RHS2=RightHand2; }

	friend std::ostream &operator<<(std::ostream &stream, TFParticle &par);
	
};
	std::ostream &operator<<(std::ostream &stream, TFParticle &par);
#endif
