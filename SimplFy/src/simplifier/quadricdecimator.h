/*
 * quadricdecimator.h
 *
 *  Created on: Sep 14, 2014
 *      Author: brunt
 */

#include "base.h"
#include "mesh_properties.h"
#include <iostream>
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

#ifndef QUADRICDECIMATOR_H_
#define QUADRICDECIMATOR_H_

namespace brndan022 {
using namespace vcg;
	template <typename M> class QuadricDecimator : public Simplifier<M> {
	private:
		TriEdgeCollapseQuadricParameter qparams;
		float TargetError=std::numeric_limits<float>::max();
		bool CleaningFlag = false;
		int FinalSize;

		class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric<MyMesh,
						VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex> > {
				public:
					typedef vcg::tri::TriEdgeCollapseQuadric<MyMesh, VertexPair,
							MyTriEdgeCollapse, QInfoStandard<MyVertex> > TECQ;
					typedef MyMesh::VertexType::EdgeType EdgeType;
					inline MyTriEdgeCollapse( const VertexPair &p, int i, BaseParameterClass *q) :TECQ(p, i, q){};
				};


	public:
		QuadricDecimator(){

		};

		void simplify(M &mesh) override
		{
			if (CleaningFlag) {
					int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
					int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
					printf("Removed %i duplicate and %i unreferenced vertices from mesh \n",
							dup, unref);
				}

				printf("reducing it to %i\n", FinalSize);

				vcg::tri::UpdateBounding<MyMesh>::Box(mesh);

				// decimator initialization
				vcg::LocalOptimization<MyMesh> DeciSession(mesh, &qparams);

				int t1 = clock();
				DeciSession.Init<MyTriEdgeCollapse>();
				int t2 = clock();
				printf("Initial Heap Size %i\n", int(DeciSession.h.size()));

				DeciSession.SetTargetSimplices(FinalSize);
				DeciSession.SetTimeBudget(0.5f);
				if (TargetError < std::numeric_limits<float>::max())
					DeciSession.SetTargetMetric(TargetError);

				while (DeciSession.DoOptimization() && mesh.fn > FinalSize
						&& DeciSession.currMetric < TargetError)
					printf("Current Mesh size %7i heap sz %9i err %9g \r", mesh.fn,
							int(DeciSession.h.size()), DeciSession.currMetric);

				int t3 = clock();
				printf("mesh  %d %d Error %g \n", mesh.vn, mesh.fn, DeciSession.currMetric);
				printf("\nCompleted in (%i+%i) msec\n", t2 - t1, t3 - t2);

		};

		void setParameters(int argc, char ** argv) override
		{
			FinalSize=atoi(argv[3]);

			std::cout<<"Quadric parameters called"<<std::endl;
				qparams.QualityThr = .3;
				float TargetError = std::numeric_limits<float>::max();
				bool CleaningFlag = false;
				// parse command line.
				for (int i = 4; i < argc;) {
					if (argv[i][0] == '-')
						switch (argv[i][1]) {
						case 'H':
							qparams.SafeHeapUpdate = true;
							printf("Using Safe heap option\n");
							break;
						case 'Q':
							if (argv[i][2] == 'y') {
								qparams.QualityCheck = true;
								printf("Using Quality Checking\n");
							} else {
								qparams.QualityCheck = false;
								printf("NOT Using Quality Checking\n");
							}
							break;
						case 'N':
							if (argv[i][2] == 'y') {
								qparams.NormalCheck = true;
								printf("Using Normal Deviation Checking\n");
							} else {
								qparams.NormalCheck = false;
								printf("NOT Using Normal Deviation Checking\n");
							}
							break;
						case 'O':
							if (argv[i][2] == 'y') {
								qparams.OptimalPlacement = true;
								printf("Using OptimalPlacement\n");
							} else {
								qparams.OptimalPlacement = false;
								printf("NOT Using OptimalPlacement\n");
							}
							break;
						case 'S':
							if (argv[i][2] == 'y') {
								qparams.ScaleIndependent = true;
								printf("Using ScaleIndependent\n");
							} else {
								qparams.ScaleIndependent = false;
								printf("NOT Using ScaleIndependent\n");
							}
							break;
						case 'B':
							if (argv[i][2] == 'y') {
								qparams.PreserveBoundary = true;
								printf("Preserving Boundary\n");
							} else {
								qparams.PreserveBoundary = false;
								printf("NOT Preserving Boundary\n");
							}
							break;
						case 'T':
							if (argv[i][2] == 'y') {
								qparams.PreserveTopology = true;
								printf("Preserving Topology\n");
							} else {
								qparams.PreserveTopology = false;
								printf("NOT Preserving Topology\n");
							}
							break;
						case 'q':
							qparams.QualityThr = atof(argv[i] + 2);
							printf("Setting Quality Thr to %f\n", atof(argv[i] + 2));
							break;
						case 'n':
							qparams.NormalThrRad = math::ToRad(atof(argv[i] + 2));
							printf("Setting Normal Thr to %f deg\n", atof(argv[i] + 2));
							break;
						case 'b':
							qparams.BoundaryWeight = atof(argv[i] + 2);
							printf("Setting Boundary Weight to %f\n", atof(argv[i] + 2));
							break;
						case 'e':
							TargetError = float(atof(argv[i] + 2));
							printf("Setting TargetError to %g\n", atof(argv[i] + 2));
							break;
						case 'P':
							CleaningFlag = true;
							printf("Cleaning mesh before simplification\n");
							break;

						default:
							printf("Unknown option '%s'\n", argv[i]);
							exit(0);
						}
					i++;
				}
		};

		//Nothing here just yet :)
		virtual ~QuadricDecimator(){};
	};

} /* namespace brndan022 */

#endif /* QUADRICDECIMATOR_H_ */
