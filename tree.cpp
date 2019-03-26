//****************************************************************************************
//
// SPH Code
//
// created	87/04 By MAFB
// Last change	87/04 By MAFB
// Last change	87/08 By MAFB
// Last change	87/11 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
//  tree file
//	Finds interaction through a tree search.
//****************************************************************************************
#include "tree.h"
#include "kernel.h"



bool TtreeCell::inCell(const Vector &point)
{
return (point.GetX() >= lbloc.GetX()) &&
	   (point.GetX() < rtloc.GetX()) &&
	   (point.GetY() >= lbloc.GetY()) &&
	   (point.GetY() < rtloc.GetY());
}


void TtreeCell::divide()
{
branches[0]=new TtreeCell(mw, this, lbloc,
						  Vector( (lbloc.GetX()+rtloc.GetX())/2,(lbloc.GetY()+rtloc.GetY())/2 ));

branches[1]=new TtreeCell(mw, this, Vector( (lbloc.GetX()+rtloc.GetX())/2,lbloc.GetY() ),
						  Vector( rtloc.GetX(),(lbloc.GetY()+rtloc.GetY())/2 ));

branches[2]=new TtreeCell(mw, this, Vector( (lbloc.GetX()+rtloc.GetX())/2,(lbloc.GetY()+rtloc.GetY())/2 ),
						  rtloc);

branches[3]=new TtreeCell(mw, this, Vector( lbloc.GetX(),(lbloc.GetY()+rtloc.GetY())/2 ),
						  Vector( (lbloc.GetX()+rtloc.GetX())/2,rtloc.GetY() ));
}


TtreeCell *TtreeCell::findCell(TFParticle &particle)
{
	double h=particle.GetH();

	if (rtloc.GetX() - particle.GetLoc().GetX() <= h ||
		rtloc.GetY() - particle.GetLoc().GetY() <= h ||
		particle.GetLoc().GetX() - lbloc.GetX() <= h ||
		particle.GetLoc().GetY() - lbloc.GetY() <= h) return ( (parent) ? parent : this );
	else if ( branches[0]->inCell(particle.GetLoc()) ) return branches[0]->findCell(particle);
	else if ( branches[1]->inCell(particle.GetLoc()) ) return branches[1]->findCell(particle);
	else if ( branches[2]->inCell(particle.GetLoc()) ) return branches[2]->findCell(particle);
	else if ( branches[3]->inCell(particle.GetLoc()) ) return branches[3]->findCell(particle);
	else return NULL;
}

void TtreeCell::findNeigh(TFParticle &particle, intvec &neiparticles)
{
	double h=particle.GetH()*1.01;

	Vector centersdiff=0.5*(lbloc+rtloc)-particle.GetLoc();
	double halfsidex=0.5*(rtloc.GetX()-lbloc.GetX());
	double halfsidey=0.5*(rtloc.GetY()-lbloc.GetY());
	if ( centersdiff.TwoNorm()<=Vector(halfsidex+h,halfsidey+h).TwoNorm())

		  if (iparticles.size()>1)
			  for (int i=0; i<4; i++)
				  branches[i]->findNeigh(particle, neiparticles);

		  else if (iparticles.size()>0)
			  neiparticles.push_back(iparticles[0]);
}


void TreeSearch(World &mw)
{
	double hm;
	Vector dx;
	double dr;
	int pNo=mw.particles.size();

	//TtreeCell *cell;
	TtreeCell treeCell(&mw, NULL,
		mw.minimum.loc-Vector(mw.maximum.h,mw.maximum.h),
		mw.maximum.loc+Vector(mw.maximum.h,mw.maximum.h));

	mw.interactions.clear();

	for (int i=0; i<pNo; i++)
	{
		//cell=treeCell.findCell(mw.particles[i]);

		intvec neiparticles;
		treeCell.findNeigh(mw.particles[i], neiparticles);
		int nNo=neiparticles.size();
		for (int j=0; j<nNo; j++)
		{
			if (i>=neiparticles[j]) continue;

			hm=0.5*(mw.particles[neiparticles[j]].GetH()+ mw.particles[i].GetH());	//Average smoothing length

			dx=mw.particles[neiparticles[j]].GetLoc()-mw.particles[i].GetLoc();		//Distance Vector from j to i: ri-rj
			dr=dx.TwoNorm();

			if (dr<=hm)
				mw.interactions.push_back(Interaction(neiparticles[j], i, hm, dr, dx, Kernel(dr, hm, mw.options.KERNELTYPE), DKernel(dr, hm, mw.options.KERNELTYPE) ) );
		}
	}
}


