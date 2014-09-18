/*
 * main.cpp
 *
 *  Created on: Sep 14, 2014
 *      Author: brunt
 */

#include <iostream>

#include "simplifier/quadricdecimator.h"
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
	if(argc<3)
	{
		std::cout<<"Invalid arguments have been entered"<<std::endl;
		exit(-1);
	}

	//Test case!
	MyMesh m;

	readFile<MyMesh>(m,argv[1]);

	Simplifier<MyMesh> *s = new QuadricDecimator<MyMesh> ();
	s->setParameters(argc, argv);
	s->simplify(m);

	writeMesh(m,argv[2]);

	cout<<"Successful run!"<<endl;
	return 0;
}


