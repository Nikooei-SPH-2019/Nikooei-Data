//****************************************************************************************
//
// SPH Code
//
// created	87/04 By MAFB
// Last change	97/11 By Nikooei
// Last change	97/11 By Nikooei
//
//****************************************************************************************
//  tree file
//	Finds interaction through a tree search.
//****************************************************************************************


#ifndef _TREE_H_
#define _TREE_H_

#include "world.h"
#include "vectors.h"
#include "particle.h"



class TtreeCell {
public:

	World *mw;
	TtreeCell *parent;
	TtreeCell *branches[4];
	Vector lbloc,rtloc;
	std::vector<int> iparticles;

	TtreeCell(World *w, TtreeCell *parentCell, Vector lb, Vector rt) {
		mw=w;
		parent=parentCell;
		branches[0]=branches[1]=branches[2]=branches[3]=NULL;
		lbloc=lb; rtloc=rt;

		if (parentCell==NULL) {
			int pNo=mw->particles.size();
			for (int i=0;i<pNo;i++) iparticles.push_back(i);
		}
		else {
			int pNo=parentCell->iparticles.size();
			for (int i=0;i<pNo;i++)
				if ( inCell(mw->particles[parentCell->iparticles[i]].GetLoc()) )
					iparticles.push_back(parentCell->iparticles[i]);
		}

		if (iparticles.size()>1) divide();
	}

	~TtreeCell(){
		iparticles.clear();
		if (branches[0]) delete branches[0];
		if (branches[1]) delete branches[1];
		if (branches[2]) delete branches[2];
		if (branches[3]) delete branches[3];
	}

	bool inCell(const Vector &point);
	void divide();
	TtreeCell *findCell(TFParticle &particle);
	void findNeigh(TFParticle &particle, intvec &neiparticles);
};


#endif //_TREE_H_

