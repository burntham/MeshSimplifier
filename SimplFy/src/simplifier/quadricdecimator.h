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
using namespace tri;

float TargetError =std::numeric_limits<float>::max();

template <typename M> class QuadricDecimator : public Simplifier<M> {
	private:
		TriEdgeCollapseQuadricParameter qparams;
		bool CleaningFlag;
		int FinalSize;
        Box3<float> t;

		/*
		 * Define class for simplifier object
		 */
		class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric<MyMesh,
						VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex> > {
				public:

            class ModFlag: public UpdateFlags<MyMesh>
            {
            public:
                typedef MyMesh MeshType;

                static void PEW(MeshType &m)
                {

                  RequirePerFaceFlags(m);
                  RequireVFAdjacency(m);
                  FaceClearB(m);

                  int visitedBit=VertexType::NewBitFlag();

                  const int BORDERFLAG[3]={FaceType::BORDER0, FaceType::BORDER1, FaceType::BORDER2};
                  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
                  {
                    if((!(*vi).IsD()) && (*vi).P().X()>0)
                    {


                    }
                      VertexType::DeleteBitFlag(visitedBit);
                }
}
            };

					typedef vcg::tri::TriEdgeCollapseQuadric<MyMesh, VertexPair,
							MyTriEdgeCollapse, QInfoStandard<MyVertex> > TECQ;
					typedef MyMesh::VertexType::EdgeType EdgeType;
                    typedef MyMesh TriMeshType;
                    typedef MyTriEdgeCollapse MYTYPE;
                    typedef typename TriEdgeCollapse< TriMeshType, VertexPair, MYTYPE>::HeapType HeapType;
                    inline MyTriEdgeCollapse( const VertexPair &p, int i, BaseParameterClass *q) :TECQ(p, i, q){}

                    static void Init(TriMeshType &m, HeapType &h_ret, BaseParameterClass *_pp)
                    {
                        printf("here we go again!");
                      QParameter *pp=(QParameter *)_pp;

                    typename 	TriMeshType::VertexIterator  vi;
                    typename 	TriMeshType::FaceIterator  pf;

                    pp->CosineThr=cos(pp->NormalThrRad);

                    vcg::tri::UpdateTopology<TriMeshType>::VertexFace(m);
                    //ModFlag::PEW(m);
                      vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromVF(m);

                    if(pp->FastPreserveBoundary)
                      {
                        for(pf=m.face.begin();pf!=m.face.end();++pf)
                        if( !(*pf).IsD() && (*pf).IsW() )
                          for(int j=0;j<3;++j)
                            if((*pf).IsB(j))
                            {
                              (*pf).V(j)->ClearW();
                              (*pf).V1(j)->ClearW();
                            }
                        }

                    if(pp->PreserveBoundary)
                      {
                        WV().clear();
                        for(pf=m.face.begin();pf!=m.face.end();++pf)
                        if( !(*pf).IsD() && (*pf).IsW() )
                          for(int j=0;j<3;++j)
                            if((*pf).IsB(j))
                            {
                              if((*pf).V(j)->IsW())  {(*pf).V(j)->ClearW(); WV().push_back((*pf).V(j));}
                              if((*pf).V1(j)->IsW()) {(*pf).V1(j)->ClearW();WV().push_back((*pf).V1(j));}
                            }
                        }

                      InitQuadric(m,pp);

                    // Initialize the heap with all the possible collapses
                      if(IsSymmetric(pp))
                      { // if the collapse is symmetric (e.g. u->v == v->u)
                        for(vi=m.vert.begin();vi!=m.vert.end();++vi)
                          if(!(*vi).IsD() && (*vi).IsRW())
                              {
                                  vcg::face::VFIterator<FaceType> x;
                                  for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++ x){
                                    x.V1()->ClearV();
                                    x.V2()->ClearV();
                                  }
                                  for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++x )
                                  {
                                    assert(x.F()->V(x.I())==&(*vi));
                                    if((x.V0()<x.V1()) && x.V1()->IsRW() && !x.V1()->IsV()){
                                          x.V1()->SetV();
                                          h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(),x.V1()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp )));
                                          }
                                    if((x.V0()<x.V2()) && x.V2()->IsRW()&& !x.V2()->IsV()){
                                          x.V2()->SetV();
                                          h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(),x.V2()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp )));
                                        }
                                  }
                              }
                      }
                          else
                      { // if the collapse is A-symmetric (e.g. u->v != v->u)
                              for(vi=m.vert.begin();vi!=m.vert.end();++vi)
                                  if(!(*vi).IsD() && (*vi).IsRW())
                                      {
                                          vcg::face::VFIterator<FaceType> x;
                                          UnMarkAll(m);
                                          for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++ x)
                                          {
                                              assert(x.F()->V(x.I())==&(*vi));
                                              if(x.V()->IsRW() && x.V1()->IsRW() && !IsMarked(m,x.F()->V1(x.I()))){
                                      h_ret.push_back( HeapElem( new MYTYPE( VertexPair (x.V(),x.V1()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp)));
                                                          }
                                              if(x.V()->IsRW() && x.V2()->IsRW() && !IsMarked(m,x.F()->V2(x.I()))){
                                      h_ret.push_back( HeapElem( new MYTYPE( VertexPair (x.V(),x.V2()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp)));
                                                      }
                                          }
                                      }
                        }
                  }
                };
	public:
		QuadricDecimator(){
			CleaningFlag = false;
			FinalSize=0;
        }


        /*
         * Configure the parameters of the simplifier
         */
        virtual void setParameters(int argc, char ** argv) override
        {
            FinalSize=atoi(argv[4]);

            std::cout<<"Quadric configuring"<<std::endl;
                qparams.QualityThr = .3;
                TargetError = std::numeric_limits<float>::max();
                CleaningFlag = false;
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

		/*
		 * Execute simplification
		 */
		virtual void simplify(M &mesh) override
		{
			if (CleaningFlag) {
					int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
					int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
					printf("Removed %i duplicate and %i unreferenced vertices from mesh \n",
							dup, unref);
				}

				printf("reducing it to %i\n", FinalSize);

                vcg::tri::UpdateBounding<MyMesh>::Box(mesh);

				// decimator initialization (Simplifeier object)
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

                if (CleaningFlag) {
                        int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
                        int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
                        printf("Removed %i duplicate and %i unreferenced vertices from mesh \n",
                                dup, unref);
                    }

				int t3 = clock();
				printf("mesh  %d %d Error %g \n", mesh.vn, mesh.fn, DeciSession.currMetric);
				printf("\nCompleted in (%i+%i) msec\n", t2 - t1, t3 - t2);

		};



		//Nothing here just yet :)
		virtual ~QuadricDecimator(){};
	};

} /* namespace brndan022 */

#endif /* QUADRICDECIMATOR_H_ */
