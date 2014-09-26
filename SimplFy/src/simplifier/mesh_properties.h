#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>
#include <wrap/io_trimesh/import.h>

#ifndef MESH_PROPERTIES_H_
#define MESH_PROPERTIES_H_


namespace brndan022
{
		class MyVertex;
		class MyEdge;
		class MyFace;

		class MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>::AsVertexType, vcg::Use<MyEdge>::AsEdgeType, vcg::Use<MyFace>::AsFaceType>{};

		class MyVertex : public vcg::Vertex < MyUsedTypes ,
			vcg::vertex::VFAdj,
			vcg::vertex::Coord3f,
            vcg::vertex::Normal3f,
			vcg::vertex::Mark,
			vcg::vertex::BitFlags,
            vcg::vertex::Color4b >
		{
			public:
				vcg::math::Quadric<double> & Qd() {return q;}
			private:
				vcg::math::Quadric<double> q;
		};

		class MyEdge : public vcg::Edge< MyUsedTypes> {};

		typedef vcg::tri::BasicVertexPair<MyVertex> VertexPair;

		class MyFace : public vcg::Face<MyUsedTypes,
			vcg::face::VFAdj,
			vcg::face::VertexRef,
            vcg::face::Normal3f,
            vcg::face::BitFlags >{};

		//my mesh class:
		class MyMesh : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> > {};
}
#endif
