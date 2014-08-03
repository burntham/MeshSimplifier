#include <iostream>

//io
#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

using namespace vcg;
using namespace tri;

//Class prototypes for Mesh definition
class MyVertex;
class MyEdge;
class MyFace;

/*NB templated struct?*/
struct MyUsedTypes: public UsedTypes<Use<MyVertex>::AsVertexType, Use<MyEdge>::AsEdgeType,
		Use<MyFace>::AsFaceType>{};

//class MyVertex: public Vertex<MyUsedTypes, vertex::VFAdj, vertex>
//{}

int main()
{
	using namespace std;
	std::cout<<"Hello World!"<<"\n";
	return 5;
}

