//****************************************************************************************
//
// SPH Code
//
// created	86/05 By Fatehi
// Last change	86/08 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
//  linkedlist file
//	Finds interaction through a grid-based search.
//****************************************************************************************


#include "world.h"
#include "linkedlist.h"

#include "kernel.h"


TGrid::TGrid(Vector p1, Vector p2, double hmax)
{
	lbloc = p1;
	rtloc = p2;

	std::vector<TCell> collomn;
	xNo = int((rtloc.GetX() - lbloc.GetX()) / hmax - 0.1);
	yNo = max(int((rtloc.GetY() - lbloc.GetY()) / hmax - 0.1), 1);
	delta = max((rtloc.GetX() - lbloc.GetX()) / xNo,
		(rtloc.GetY() - lbloc.GetY()) / yNo);
	collomn.assign(yNo, TCell());
	cells.assign(xNo, collomn);

}



void TGrid::Initialize()
{

	for (unsigned int i = 0; i<cells.size(); i++)
		for (unsigned int j = 0; j<cells[i].size(); j++)
		{
			cells[i][j].lbloc = lbloc + Vector(i*delta, j*delta);		//Assign the location of left-bottom point of each cell
			cells[i][j].iparticles.clear();
		}

}


void TGrid::AssignParticles(World &mw)
{
	int pNo = mw.particles.size();

	for (int i = 0; i<pNo; i++)
	{
		const Vector relposition = mw.particles[i].GetLoc() - lbloc;
		const int xindex = min(max(int(relposition.GetX() / delta), 0), xNo - 1);
		const int yindex = min(max(int(relposition.GetY() / delta), 0), yNo - 1);

		cells[xindex][yindex].iparticles.push_back(i);

	}

}


void TGrid::FindNeighbourCells()
{
	int iNo = cells.size();

	for (register int i = 0; i<iNo; i++)
	{
		int jNo = cells[i].size();
		for (register int j = 0; j<jNo; j++)
		{
			cells[i][j].neighbours.clear();

			if (j<jNo - 1)
				cells[i][j].neighbours.push_back(&cells[i][j + 1]);
			if (i<iNo - 1)									//Find and save neighbours of each cell.
			{
				cells[i][j].neighbours.push_back(&cells[i + 1][j]);
				if (j>0)
					cells[i][j].neighbours.push_back(&cells[i + 1][j - 1]);
				if (j<jNo - 1)
					cells[i][j].neighbours.push_back(&cells[i + 1][j + 1]);
			}
		}
	}

}




void LinkedList(World &mw)
{
	double hm;
	Vector dx;
	double dr;
	TGrid grid(mw.minimum.loc, mw.maximum.loc, mw.maximum.h);
	grid.Initialize();
	grid.AssignParticles(mw);
	grid.FindNeighbourCells();
	mw.interactions.clear();
	int pNo = mw.particles.size();
	for (int i = 0; i < pNo; i++)
		mw.particles[i].AllInter.clear();

	for (unsigned int i = 0; i<grid.cells.size(); i++)
		for (unsigned int j = 0; j<grid.cells[i].size(); j++)			//For each Cell in the grid
			for (unsigned int k = 0; k<grid.cells[i][j].iparticles.size(); k++)
			{
				const int kk = grid.cells[i][j].iparticles[k];

				for (unsigned int m = k + 1; m<grid.cells[i][j].iparticles.size(); m++)		//search the other perticles in the same cell
				{
					const int mm = grid.cells[i][j].iparticles[m];

					hm = 0.5*(mw.particles[kk].GetH() + mw.particles[mm].GetH());	//Average smoothing length

					dx = mw.particles[kk].GetLoc() - mw.particles[mm].GetLoc();		//Distance Vector from j to i: ri-rj
					dr = dx.TwoNorm();

					if (dr <= hm)
					{
						mw.interactions.push_back(Interaction(kk, mm, hm, dr, dx, Kernel(dr, hm, mw.options.KERNELTYPE), DKernel(dr, hm, mw.options.KERNELTYPE)));
						mw.particles[kk].AllInter.push_back(mm); 
						mw.particles[mm].AllInter.push_back(kk); 
					}



				}

				for (unsigned int n = 0; n<grid.cells[i][j].neighbours.size(); n++)			// For each neighbour of the cell
					for (unsigned int p = 0; p<grid.cells[i][j].neighbours[n]->iparticles.size(); p++)	//search the other perticles in the neighbouring cells
					{
						const int pp = grid.cells[i][j].neighbours[n]->iparticles[p];

						hm = 0.5*(mw.particles[kk].GetH() + mw.particles[pp].GetH());	//Average smoothing length

						dx = mw.particles[kk].GetLoc() - mw.particles[pp].GetLoc();		//Distance Vector from j to i: ri-rj
						dr = dx.TwoNorm();

						if (dr <= hm)
						{
							mw.interactions.push_back(Interaction(kk, pp, hm, dr, dx, Kernel(dr, hm, mw.options.KERNELTYPE), DKernel(dr, hm, mw.options.KERNELTYPE)));

							mw.particles[kk].AllInter.push_back(pp); 
							mw.particles[pp].AllInter.push_back(kk); 
						}

					}
			}

}

