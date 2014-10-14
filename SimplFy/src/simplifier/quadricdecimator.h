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
#include<vcg/simplex/face/pos.h>

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
    bool defaultBoundaries;
    int FinalSize;
    Box3<float> t;

    /*
         * Define class for simplifier object
         */
    class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric<MyMesh,
            VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex> >
    {

    public:
        /**
         * @brief The ModFlag class
         * This class is to modify the default behavior of the update border flag
         */
        class ModFlag: public UpdateFlags<MyMesh>
        {
        public:
            typedef MyMesh MeshType;

            static void UpdateCustomBoundaryTriangles(MeshType &m)
            {

                RequirePerFaceFlags(m);
                RequireVFAdjacency(m);

                MeshType::FaceIterator pf;
                for (pf=m.face.begin();pf!=m.face.end();++pf)
                {   bool isBoundry = false;
                    for (int j=0; j<3; ++j)
                    {
                        if(! m.workingBBox.IsIn((*pf).V(j)->P()))
                        {
                            isBoundry=true;
                        }
                    }
                    if(isBoundry)
                    {
                        (*pf).V(0)->SetB();
                        (*pf).V(1)->SetB();
                        (*pf).V(2)->SetB();
                        (*pf).SetB(0);
                        (*pf).SetB(1);
                        (*pf).SetB(2);
                    }
                }

            }
            static void UpdateBoundaryCount(MeshType &m)
            {
                RequirePerFaceFlags(m);
                RequireVFAdjacency(m);

                MeshType::FaceIterator pf;
                for (pf=m.face.begin();pf!=m.face.end();++pf)
                {
                    if ((*pf).IsB(0) || (*pf).IsB(1) || (*pf).IsB(2) )
                    { ++m.bn;}
                }
            }
        };

        typedef vcg::tri::TriEdgeCollapseQuadric<MyMesh, VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex> > TECQ;
        typedef MyMesh::VertexType::EdgeType EdgeType;
        typedef MyMesh TriMeshType;
        typedef MyTriEdgeCollapse MYTYPE;
        typedef typename TECQ::QParameter QParameter;
        typedef typename TriMeshType::FaceType FaceType;
        typedef typename TriEdgeCollapse< TriMeshType, VertexPair, MYTYPE>::HeapType HeapType;
        typedef typename TriEdgeCollapse<TriMeshType, VertexPair, MYTYPE>::HeapElem HeapElem;


        inline MyTriEdgeCollapse( const VertexPair &p, int i, BaseParameterClass *q) :TECQ(p, i, q){ }

        static void Init(TriMeshType &m, HeapType &h_ret, BaseParameterClass *_pp)
        {

            QParameter *pp=(QParameter *)_pp;

            typename 	TriMeshType::VertexIterator  vi;
            typename 	TriMeshType::FaceIterator  pf;

            pp->CosineThr=cos(pp->NormalThrRad);

            vcg::tri::UpdateTopology<TriMeshType>::VertexFace(m);
            if (pp->PreserveBoundary&&pp->FastPreserveBoundary)
            {
                vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromVF(m);
                ModFlag::UpdateCustomBoundaryTriangles(m);
            } else if(pp->PreserveBoundary)
            {
                vcg::tri::UpdateFlags<TriMeshType>::FaceClearB(m);
                ModFlag::UpdateCustomBoundaryTriangles(m);
            }
            else
            {
                vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromVF(m);
            }

            ModFlag::UpdateBoundaryCount(m);


            if(pp->PreserveBoundary)
            {
                WV().clear();
                for(pf=m.face.begin();pf!=m.face.end();++pf)
                    if( !(*pf).IsD() && (*pf).IsW() )
                        for(int j=0;j<3;++j)
                            if((*pf).IsB(j))                            {
                                if((*pf).V(j)->IsW())  {(*pf).V(j)->ClearW(); WV().push_back((*pf).V(j));}
                                if((*pf).V1(j)->IsW()) {(*pf).V1(j)->ClearW();WV().push_back((*pf).V1(j));}
                            }
            }

            InitQuadric(m,pp);

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

        // Final Clean up after the end of the simplification process
        static void Finalize(TriMeshType &m, HeapType& /*h_ret*/, BaseParameterClass *_pp)
        {
          QParameter *pp=(QParameter *)_pp;

          if(pp->PreserveBoundary)
          {
            typename 	std::vector<typename TriMeshType::VertexPointer>::iterator wvi;
            for(wvi=WV().begin();wvi!=WV().end();++wvi)
              if(!(*wvi)->IsD()) (*wvi)->SetW();
          }
        }
    };

public:
    static float borderCount;
    Box3f WorkingBox;
    QuadricDecimator(){
        CleaningFlag = false;
        FinalSize=0;
    }

    //     static float HeapSimplexRatio(BaseParameterClass *_pp) {return 2.3f;}

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
                        if (argv[i][3]=='d')
                        {
                            qparams.FastPreserveBoundary=true;
                            printf("Default boundary preservation also enabled\n");
                        }
                        qparams.PreserveBoundary = true;
                        float minX, maxX, minY,maxY, minZ, maxZ;
                        if (i+6<argc){
                            try{
                                minX=std::stof(argv[++i]),maxX=stof(argv[++i]),minY=stof(argv[++i]),maxY=stof(argv[++i]),minZ=stof(argv[++i]),maxZ=stof(argv[++i]);
                            }
                            catch(...)
                            {
                                cout<<"invalid arguments for boundary box\n -By <minX> <maxX> <minY> <maxY> <minZ> <maxZ>\n";
                                break;
                            }
                        }else{
                            cout<<"Custom Boundary Preservation disabled\n";
                            qparams.PreserveBoundary = false;
                            break;
                        }
                        Point3f min(minX, minY, minZ);
                        Point3f max(maxX,maxY,maxZ);
                        WorkingBox = Box3f(min,max);
                        printf("Preserving Boundary defined as all faces with atleast 1 vertex outside of:\n");
                        printf("\tX from %f to %f \n",minX,maxX);
                        printf("\tY from %f to %f \n",minY,maxY);
                        printf("\tZ from %f to %f \n",minZ,maxZ);

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
                    TargetError = float(stof(argv[++i])); // Fixed the Target error option
                    printf("Setting TargetError to %g\n", TargetError);
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

    }

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

        vcg::tri::UpdateBounding<MyMesh>::Box(mesh);
        mesh.workingBBox = Box3f(WorkingBox);
        // decimator initialization (Simplifeier object)
        vcg::LocalOptimization<MyMesh> DeciSession(mesh, &qparams);

        int t1 = clock();
        //MyTriEdgeCollapse::defBound = true;
        DeciSession.Init<MyTriEdgeCollapse>();
        int t2 = clock();
        printf("Initial Heap Size %i\n", int(DeciSession.h.size()));

        if(qparams.PreserveBoundary)
        {
            int origTarget=FinalSize;
            if(FinalSize<=mesh.bn)
                FinalSize = FinalSize+mesh.bn;

            printf("Target Faces:%d\n"
                   "Initial Faces:%d\n"
                   "Boundary Triangles:%d\n"
                   "Adjusted target:%d\n",
                   origTarget,mesh.fn,mesh.bn,FinalSize);
        }

        printf("reducing it to %i\n", FinalSize);

        DeciSession.SetTargetSimplices(FinalSize);
        DeciSession.SetTimeBudget(0.5f);
        if (TargetError < std::numeric_limits<float>::max())
            DeciSession.SetTargetMetric(TargetError);

        printf("Initial MeshSize:%d ,Boundries:%d , ExpectedSize:%d",mesh.fn,mesh.bn,(FinalSize));

        printf("Mesh Size:%d, Borders:%d, Simplifiable faces:%d",mesh.fn, mesh.bn, (mesh.fn-mesh.bn));
        while (DeciSession.DoOptimization() && mesh.fn>FinalSize
               && DeciSession.currMetric < TargetError)
            printf("Current Mesh size %7i test %7i heap sz %9i err %9g \r", mesh.fn, (mesh.fn-mesh.bn) ,
                   int(DeciSession.h.size()), DeciSession.currMetric);

        if (CleaningFlag) {
            int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
            int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
            printf("Removed %i duplicate and %i unreferenced vertices from mesh \n",
                   dup, unref);
        }

        int t3 = clock();
        printf("mesh  vertices:%d faces:%d NonBoundryFaces:%d Error %g \n", mesh.vn, mesh.fn,(mesh.fn-mesh.bn), DeciSession.currMetric);
        printf("\nCompleted in (%i+%i) msec\n", t2 - t1, t3 - t2);

    }

    //Nothing here just yet :)
    virtual ~QuadricDecimator(){}
};

} /* namespace brndan022 */

#endif /* QUADRICDECIMATOR_H_ */
