//****************************************************************************************
//
// SPH Code
//
// created	85/12 By Fatehi
// Last change	86/09 By Fatehi
// Last change	87/07 By Fatehi
// Last change	87/08 By Fatehi
// Last change	88/01 By Fatehi
// Last change	88/02 By Fatehi
// changed	88/07 By Fatehi
// changed	88/11 By Fatehi
// changed	89/03 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
// World Header file
//	Difines the computational domain
//****************************************************************************************


#ifndef _WORLD_H_
#define _WORLD_H_

#include "options.h"
#include "particle.h"
#include "interaction.h"
#include "extremum.h"
#include "boundary_conditions.h"
#include "tensors3.h"
#include "iostream"


class World {
private:

	double time;			//Simulation time (s)
	double deltat;			// Time-step size (s)
	double OMEGA;
	Vector VELOCITY;
	Vector ACCEL;
	double E0,L0;


public:

	std::vector<TFParticle> particles;
	/* All Fluid particles are stored here. particles of different materials are percieved by type num of
	the "mat" member of TFParticle class (type is a public member of "mat"). "mat" is a pointer to "TFluid" class.
	a class which defines a fluid  by its properties such as vicosity, density and so, TFluid is a child of
	"TRMaterial" class (real material, i.e. has physical properties). grand father of all material type class
	is "TMaterial" class which has only a public variable named "type". this integer is used to identify particles
	of each material. so if we have a 3 phase problem there are particles with 3 type numbers. particle coordinates
	and other properties like "h" are read from the mesh file. and material properties of each material type are read
	from a seperate ".mat" file. */
	std::vector<Interaction> interactions;
	/* Interaction is a class used for storing respective properties of neighbor particles. infomation about their distance
	relative velocity, stress, value of kernel function and so. each object interaction has to intrges "iindex" and "jindex"
	defining interaction matrix, this means particles "iindex" and "jindex" have intraction and the information of interaction
	is stored in object.*/

	TOptions options; /* in TOptions class parameters such as integration method, read/write method (binary or asci),
					  search method and some other reside here. we read them from a file called "options.opt".*/

	TMat mat;
	/* This is a class for storing physical properties of all materials present in problem. "TMat" class has a member
	variable which is a vector of "pointers" to class TMaterial. we like to store all material properties of our prolem
	in a single vector. but TFluid and TSolid are different types and so it is not possible to define a vector of them.
	but here by using a pointer to base class type of "TFluid" and "TSolid" we can use runtime overloading and "TMaterial"
	pointers of this vector can point to "TFluid" or "TSolid".*/

	Extremum minimum;	//Minimum properties
	Extremum maximum;	//Maximum properties
						/* if the vaues of pressure, velocity, tempreture, ... exceed of values defined in "Extremum" SePeHr will warn
						you. good for cheching if the solution is sane.*/

//	std::vector<TPeriodicBCBase*> periodicBC;		//Periodic Boundary Condition
													/*This is a vector for handling of periodic BC's. each periodic BC should have two boundries with
													exactly the same properties. matching boundries are defined in TperiodicBC type(refer to).
													for example if we have a square boundry.suppose 2 vertical & horizantal boundries are periodic
													thed this vector (i.e. periodicBC) will contain 2 elements of type TperiodicBC. in first one the
													vertical & in second one horizantal periodic BC's are defined by a series of imaginery points of
													type TIparticle (refer to).*/
	std::vector<Vector> GradSxx;
	std::vector<Vector> GradSxy;
	std::vector<Vector> GradSyy;

	int STEP;

	int Solver_itr;



	static double si0;

	/// By S.J
	void SetSTEP(int stp) { STEP = stp; }
	int GetSTEP() { return STEP; }

	void SetITR(int ITR) { Solver_itr = ITR; }
	int GetITR() { return Solver_itr; }
	

	///////////

	void SetOMEGA(double ww) { OMEGA = ww; }
	double GetOMEGA() { return OMEGA; }

	void SetVELOCITY(Vector ww) { VELOCITY = ww; }
	Vector GetVELOCITY() { return VELOCITY; }

	void SetACCEL(Vector ww) { ACCEL = ww; }
	Vector GetACCEL() { return ACCEL; }

	std::vector<TFlowBC*> flowBC;		//flow Boundary Condition
										/* This is a vector which contains information about our boundry particles. "TPointBC" is a general class we use to
										define different BC's. it contains only a series of points which show boundry and a number called "index"(in TPointBC) shows
										which boundyr these set of particles belong to. For example "TPeriodicBC" class mentioned above is a child of "TPointBC".
										"TFlowBC" is also a child of this class we use to define boundry condition such as inflow, outflow, pressure inlet/outlet.
										in problem each of these BC's use one elemnt of this vector called "flowBC" to identify themselves. as you see "flowBC"
										is a vector of "TFlowBC" pointers. the reason is that we want to include different bc's like inflow & outflow in "flowBC"
										vector so we define a vector of pointers to parent class. after reading different bc's from the input file we can access
										memberfunctions of "TInflow" & "TOutflow" by using the "->" operator by runtime overloading.*/

										//*********************************************************** BY MAFB & MSS
	TwallBCVector wallBCVector;		//Wall Boundary Conditions
									//***********************************************************END BY MAFB & MSS

	double c2;		//Square of speed of sound
	double dx0;		//(m) Initial particle spacing.
	double Pref;		//Reference pressure

	World()				//Constructor
	{
		time = 0.0;
		Pref = 0.0;
	}

	~World() {}			//Destructor


	void ReadDataFile();
	/* This is a function reads particles, their properties such as their type (each fluid or solid has a number)
	, velocity, mass and other. reding of the file can be done in two manners ASCI or Binary.*/
	void Initialize();
	/* Here things like initializing of time, finding interactions for the first time step, adjusting flow boundry
	conditions and ...*/


	void SaveDataFile();
	/* again like reading of data file we can use Binary or ASCI method.*/
	void Report(int stepc, double t);
	
	void UpdateBoundary();


	void FindInteractions();
	/* this functions searchs foe each particle neighbors. currently 2 search methods of "directfind" and "linkedlist"
	are available. in directfind distance for each particle from all other particles is computed and if less than kernel
	distance renge the interacion between particles is set. this is a method for small number of particles because of
	huge computational time needs. */
	void UpdateInteractions();
	/* In each timestep first interactions are set by "Findinteractions" mentioned above. now by this new interaction matrix
	we should update interacion inforamtion such as distance of particles, value of kernel function, relative velocity and so
	on. this is what "UpdateInteractions" do.*/

	void EvalParticleNumDensity();
	void EvalInitialParticleNumDensity();
	void EvalMomentum();
	/* after finding interactions in new timestep and setting interacions parameters by "FindInteractions" & "UpdateInteractions"
	now we should compute forces acting on particles. stress forces, gravity, pressure, surface tension (if present), electro
	magnetic forces (if present) and any other forces you want to be added is computed and set here. so if you want a multi-physics
	SePeHr here is your place.*/
	
	void EvalVStar(); /////////////////////Nikooei


	void DoTimeIntegration();
	/* after computing all forces on particles we should march to next timestep. we can use different Integration methods like
	simple eulerian integration, implicit integration methods, predictor-corrector methods and ..., computing of new location and
	velocity of particles is done here. currently euler integration, two predictor-corrector schemes, and a projection method ( for
	incompressible fluids) is available. Enjoy!*/

	void March() { time += deltat; }

	void SetTime(double t) { time = t; }
	void SetDT(double dt) { deltat = dt; }

	double GetTime() { return time; }
	double GetDT() { return deltat; }

	void Set_E0 (double Init_E) { E0 = Init_E; }
	double Get_E0() { return E0; }


	void  Set_L_0(double L_0) { L0 = L_0; }
	double Get_L_0() { return L0; }

	/* for solving continuity equation we need to compute divergence of velocity at each timestep.*/
	//void EvalDivergence2();
	void EvalDivergence3();		//[Bonet 2004] [Randles, Libersky 1996]

	/* compute divergence of pressure gradient in a form that avoids the oscillations.*/
	void EvalTimeStep();

	void EvalMinimum();
	void EvalMaximum();
	void CheckLimits();
	/* after solving a timestep here some processing is done. we shold find new max & min of our neew domain (important for
	free surface problems). setting of the new max and mins for the domain is done by functions "EvalMinimum" "EvalMaximum".
	after finding new extremums of the domain "CheckLimits" function look if any parameter has exceeded its allowable limits
	which the user has set in options file. as you should remember these values are stored in "options" variable member of
	the world.*/


	double Error();
	void UpdateCorrectors();
	//void UpdateCorrectors2();
	void UpdateCorrectors(short, short);
	//Some temporary vectors
	Tensorvec muDV;			// mu grad V

	Tensorvec DelV;


	Vectorvec muD2V;			// mu grad grad V
	Vectorvec DPtorho;		//grad P / rho
	STensorvec PI;			//Surface Stress Tensor
	Vectorvec DPItoRho;


	STensorvec gradSum_r;
	Vectorvec gradSum;
	Vectorvec gradSum_Erosion;

	Vectorvec wallNormal;

	doublevec Wallneighbor;           //  Nikooei
	doublevec flowneighbor;            /////////////////////////For interface of flow and bed (entrainment)
	doublevec bedneighbor;         /////////////////////////For interface of flow and bed (entrainment)
	doublevec Walldistance;

	std::vector<TFParticle*> interfaceParticles;
	std::vector<Interaction> interfaceInteractions;

};

//void ReadBinaryFile(World & mw, const char* szFileName, unsigned long int loadbyte);
bool ReadASCIIFile(World & mw, const char *DFile, unsigned long int loadbyte);
//void ReadPoint(World & mw, int ptype, Vector r, Vector n);
/* these are a set of functions for reading the particles file of SePeHr refer to their declarations for further info*/

//void WriteBinaryFile(World & mw, const char* szFileName, const char* szTitle);
void WriteASCIIFile(World & mw, const char* DFile, const char* szTitle);
void ClearASCIIFile(const char* DFile);
void SavePosition(std::ofstream &ofs);
void NormalizeDensity(World & mw);
void UpdateWallDensity(World &mw);
void LinkedList(World &mw);
void AddXSPH(World &mw);			//[Monaghan 2005]
void ShiftParticles3(World &mw);			//[Farzin]                       //Nikooei  
void BodyForce(World &mw);
void PressureForceSymmetric(World &mw);
/////////////////////////////////////////////////////////////////////////////////Nikooei-Non-Newtonian

void Granularmaterial1(World &mw);  //Pressure-dependent Rheology           constant viscosity model (ETHA=cst)            Ionescu 2015   

void Projection3(World &mw);
void SetFreeSurface(std::vector<TFParticle> &particles, double  siInit, double FREESURFACE);
void SetFreeSurface(std::vector<TFParticle> &particles, doublevec& SI, double  siInit, double FREESURFACE);

#endif
