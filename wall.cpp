//****************************************************************************************
//
// SPH Code
//
// created	87/07 By US
// Last change	87/07 By US
// Last change	87/11 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
// Wall boundary file
//	Handles wall boundary conditions.
//****************************************************************************************

#include "boundary_conditions.h"
#include "kernel.h"



void TWallBC::EvalNormals()
{
	int pNo = points.size();
	for (register int i = 1; i < pNo - 1; i++)
	{

		points[i].normal = Vector(-(points[i + 1].GetLoc().GetY() - points[i - 1].GetLoc().GetY()) / (points[i + 1].GetLoc() - points[i - 1].GetLoc()).TwoNorm(),
			(points[i + 1].GetLoc().GetX() - points[i - 1].GetLoc().GetX()) / (points[i + 1].GetLoc() - points[i - 1].GetLoc()).TwoNorm());

	}
	points[0].normal = points[1].normal;
	points[pNo - 1].normal = points[pNo - 2].normal;
}



void TwallBCVector::Update(std::vector<TFParticle> &particles, bool SLIP)
{
	imagesIter.clear();
	int pNo=particles.size();

	int newpNo=particles.size();
	//for(register int i=pNo; i<newpNo; i++)
		//if (particles[i].GetBCindex()==wallimage)
			//imagesIter.push_back(particles.begin()+i);
}


void TwallBCVector::clearImages(std::vector<TFParticle> &particles)
{
	for (int i = imagesIter.size() - 1; i >= 0; i--)
		particles.erase(imagesIter[i]);
}



void TwallBCVector::EvalNormals()
{
	for (unsigned int i=0; i<this->size(); i++)
		(this->begin()+i)->EvalNormals();
}

