/*
 * clustering.h
 *
 *  Created on: Sep 21, 2014
 *      Author: brunt
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_

#include "base.h"
#include<vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/clustering.h>

namespace brndan022 {
template <typename M> class ClusteringDecimator: public Simplifier<M> {
private:
	int CellNum;
	float CellSize;
	bool DupFace;

public:
	ClusteringDecimator(): CellNum(10000),CellSize(0),DupFace(false){};

	virtual void setParameters(int argc, char ** argv) override
	{
		  int i=4;
		  while(i<argc)
		  {
		    if(argv[i][0]!='-')
		    { printf("Error unable to parse option '%s'\n",argv[i]); exit(0); }
		    switch(argv[i][1])
		    {
		    case 'k' :	CellNum=atoi(argv[i+1]); ++i; printf("Using %i clustering cells\n",CellNum); break;
		    case 's' :	CellSize=atof(argv[i+1]); ++i; printf("Using %5f as clustering cell size\n",CellSize); break;
		    case 'd' :	DupFace=true; printf("Enabling the duplication of faces for double surfaces\n"); break;

		    default : {printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
		    }
		    ++i;
		  }
	};

	virtual void simplify(M &m) override
	{
		vcg::tri::UpdateBounding<M>::Box(m);
		vcg::tri::UpdateNormal<M>::PerFace(m);
		vcg::tri::Clustering<M, vcg::tri::AverageColorCell<M> > Grid;
		Grid.DuplicateFaceParam=DupFace;
		Grid.Init(m.bbox,CellNum,CellSize);

		printf("Clustering to %i cells\n",Grid.Grid.siz[0]*Grid.Grid.siz[1]*Grid.Grid.siz[2] );
		printf("Grid of %i x %i x %i cells\n",Grid.Grid.siz[0],Grid.Grid.siz[1],Grid.Grid.siz[2]);
		printf("with cells size of %.2f x %.2f x %.2f units\n",Grid.Grid.voxel[0],Grid.Grid.voxel[1],Grid.Grid.voxel[2]);

		Grid.AddMesh(m);
		Grid.ExtractMesh(m);
		printf("Output mesh vn:%i fn:%i\n",m.VN(),m.FN());
	};

	virtual	~ClusteringDecimator(){};
};

} /* namespace brndan022 */

#endif /* CLUSTERING_H_ */
