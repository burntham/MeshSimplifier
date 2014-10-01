/*
 * main.cpp
 *
 *  Created on: Sep 14, 2014
 *      Author: brunt
 */

#include <iostream>

#include "simplifier/quadricdecimator.h"
#include "simplifier/clustering.h"
#include "simplifier/mesh_properties.h"
#include "file_handler.h"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

using namespace std;
using namespace brndan022;
using namespace vcg;

int main(int argc, char**argv)
{
    if(argc<4)
    {
        printf("Invalid arguments have been entered\n"
                "Usage: simplify [method] filein.ply fileout.ply [opt/s]\n"
                "[method] can be 0 or 1\n "
                    "\t-0 for Quadric edge colapse \n"
                    "\t-1 for Clustering\n "
               );
        exit(-1);
    }


	MyMesh Mesh;
    Mesh.bn=0;



	Simplifier<MyMesh> *s;

    if(atoi(argv[1])==0)
        s=new QuadricDecimator<MyMesh>();
    else
        s = new ClusteringDecimator<MyMesh>();

    readFile<MyMesh>(Mesh,argv[2]);

    printf("xdim:%f: %f\n"
           "ydim:%f: %f\n"
           "zdim:%f: %f\n",
           Mesh.bbox.min.X(),Mesh.bbox.max.X(),
           Mesh.bbox.min.Y(),Mesh.bbox.max.Y(),
           Mesh.bbox.min.Z(),Mesh.bbox.max.Z());

    Point3f min = Mesh.bbox.min;
    Point3f max = Mesh.bbox.max;
    min.X()+=100;
    min.Y()+=100;
    min.Z()+=100;
    max.X()-=100;
    max.Y()-=100;
    max.Z()-=100;

    Mesh.workingBBox = Box3f(min,max);

    s->setParameters(argc, argv);
    s->simplify(Mesh);

    writeMesh(Mesh,argv[3]);

	printf("Successful run!");
	return 0;
}
