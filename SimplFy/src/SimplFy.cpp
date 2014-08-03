//Initial methods copied from original Quadric error metric mesh simplification.

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
//Need to look up the structure of this struct - it's extended?
struct MyUsedTypes: public UsedTypes<Use<MyVertex>::AsVertexType, Use<MyEdge>::AsEdgeType,
		Use<MyFace>::AsFaceType>{

};

//Define vertex class, included the Quadriq
class MyVertex: public Vertex<MyUsedTypes, vertex::VFAdj,vertex::Coord3f,
vertex::Normal3f, vertex::Mark, vertex::BitFlags>{
public:
	math::Quadric<double> &Qd(){
		return q;
	}
private:
	math::Quadric<double> q;
};

class MyEdge: public Edge<MyUsedTypes>{
};

typedef BasicVertexPair<MyVertex> VertexPair;

class MyFace: public Face<MyUsedTypes, face::VFAdj, face::VertexRef, face::BitFlags>{
};

//Define the mesh
class MyMesh: public tri::TriMesh<std::vector<MyVertex>,std::vector<MyFace> >{
};

class MyTriEdgeCollapse: public tri::TriEdgeCollapseQuadric<MyMesh, VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex> >
{
public:
	typedef tri::TriEdgeCollapseQuadric<MyMesh, VertexPair, MyTriEdgeCollapse,QInfoStandard<MyVertex> > TECQ;
	typedef MyMesh::VertexType::EdgeType EdgeType;
	MyTriEdgeCollapse(const VertexPair &p, int i, BaseParameterClass *pp) :
			TECQ(p,i,pp){
	}
};



MyMesh mesh;

int main()
{
	using namespace std;
	std::cout<<"Hello World!"<<"\n";
	return 5;
}

