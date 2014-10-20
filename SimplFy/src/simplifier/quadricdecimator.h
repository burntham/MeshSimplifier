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

    /**
     * @brief The MyTriEdgeCollapse class
     * This class is where most of the customizations are implemented
     */
    class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric<MyMesh,
            VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex> >
    {

    public:
        /**
         * @brief The ModFlag class
         * This class is to modify the default behavior of the update border flag
         * It also allows borders to be indentified
         */
        class ModFlag: public UpdateFlags<MyMesh>
        {
        public:
            typedef MyMesh MeshType;

            /**
             * @brief UpdateCustomBoundaryTriangles Identifies boundary triangles outside of a specified bounding box.
             * (Used to prevent seams from splitting,simplifying and stitching).
             * @param m
             */
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

            /**
             * @brief UpdateBoundaryCount Counts the number of identified boundary triangles, updating the boundary count in the mesh.
             * @param m - a reference to the mesh
             */
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

        /*
         * Custom edge collapse Constructor
         * */
        inline MyTriEdgeCollapse( const VertexPair &p, int i, BaseParameterClass *q) :TECQ(p, i, q){ }

        /**
         * @brief HeapSimplexRatio
         * Overides the heapSimplexRatio function in the vcglib implementation
         * This allows for a reduced memory cap (increased computation time though) as additional elements in the heap are removed.
         * If the heap is > heapRatio*meshsize then remove deleted elements from heap.
         * @return
         */
        static float HeapSimplexRatio(BaseParameterClass *_pp) {return 2.8f;}


        /**
         * @brief Init is method which overrides the default in QuadricDecimator -> note, vcglib requires this methods' signiature to remain the same.
         * @param m
         * @param h_ret
         * @param _pp
         */
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

            //Update the count of the boundar triangles (count is stored in the mesh object).
            ModFlag::UpdateBoundaryCount(m);

            if(pp->PreserveBoundary||pp->FastPreserveBoundary)
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

        /**
         * @brief Finalize is another overidden method, This has been overridden to prevent
         * operations happening when the fast-preserve boundary is flag is true. Since this flag has been repurposed(to retain default boundary preservation).
         *
         * @param m
         * @param _pp
         */
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
    Box3f WorkingBox;
    QuadricDecimator(): FinalSize(0),CleaningFlag(false){}

    /**
     * @brief setParameters This method is called to configure the parameters of the simplifier.
     * @param argc
     * @param argv
     */
    virtual void setParameters(int argc, char ** argv) override
    {
        FinalSize=atoi(argv[4]);
        qparams.QualityThr = .3;
        TargetError = std::numeric_limits<float>::max();
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
                {
                    int bPos = i;
                    if (argv[i][2] == 'y') {

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
                        printf("Classifying all faces as boundary faces if atleast 1 vertex outside of:\n");
                        printf("\tX from %f to %f \n",minX,maxX);
                        printf("\tY from %f to %f \n",minY,maxY);
                        printf("\tZ from %f to %f \n",minZ,maxZ);

                    } else {
                        qparams.PreserveBoundary = false;
                        printf("NOT Preserving Custom Boundary\n");
                    }

                    if (argv[bPos][3] == 'd')
                    {
                        qparams.FastPreserveBoundary=true;
                        printf("Default boundary preservation enabled\n");
                    }

                    break;
                }

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
                    printf("Mesh Cleaning enabled\n");
                    break;

                default:
                    printf("Unknown option '%s'\n", argv[i]);
                    exit(0);
                }
            i++;
        }
        printf("\n");
    }

    /**
     * @brief simplify Executes the simplification algorithm on the mesh
     * @param mesh
     */
    virtual void simplify(M &mesh) override
    {
        /*
         * Removed duplicate and unreferenced faces before the mesh is processed.
         * */
        if (CleaningFlag) {
            int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
            int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
            printf("Performing pre simplification clean\n"
                   "Removed %i duplicate and %i unreferenced vertices from mesh \n\n",
                   dup, unref);
        }

        vcg::tri::UpdateBounding<MyMesh>::Box(mesh);


        /*
         * Set a Custom bounding box used to identify boundry triangles
         * */
        mesh.workingBBox = Box3f(WorkingBox);

        vcg::LocalOptimization<MyMesh> DeciSession(mesh, &qparams);


        int t1 = clock();
        DeciSession.Init<MyTriEdgeCollapse>();
        int t2 = clock();

        /*
         *Modify the target face# if preserve boundary is enabled and taget is below boundary face#
         * */
        int origTarget=FinalSize;
        if(qparams.PreserveBoundary || qparams.FastPreserveBoundary)
        {

            if(FinalSize<=mesh.bn)
            {
                FinalSize = FinalSize+mesh.bn;
                printf("Boundary face# > target face#.\n"
                       "Target face# adjusted as boundary face# + target face#\n");
            }
        }

        printf("\tInitial Face#: %d\n"
               "\tBoundary Face#: %d\n"
               "\tOriginal Target#:%d\n"
               "\tAdjusted Target#: %d\n",
               mesh.fn, mesh.bn,origTarget,FinalSize);

        printf("\n");

        printf("Initial Heap Size %i\n", int(DeciSession.h.size()));

        DeciSession.SetTargetSimplices(FinalSize);
        DeciSession.SetTimeBudget(0.5f);

        if (TargetError < std::numeric_limits<float>::max())
            DeciSession.SetTargetMetric(TargetError);

        while (DeciSession.DoOptimization() && mesh.fn>FinalSize
               && DeciSession.currMetric < TargetError)
            printf("Current Mesh size %7i heap sz %9i err %9g \r", mesh.fn,int(DeciSession.h.size()), DeciSession.currMetric);

        printf("\n");

        if (CleaningFlag) {
            int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
            int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
            printf("Performing post simplification clean\n"
                   "Removed %i duplicate and %i unreferenced vertices from mesh \n\n",
                   dup, unref);
        }

        int t3 = clock();
        printf("Final stats:\n"
               "\tvertices:%d\n"
               "\tfaces:%d\n"
               "\tBoundary faces:%d\n"
               "\tNon-Boundry Faces:%d\n"
               "\tMesh has %d more faces then initial target:%d\n"
               "\tError:%g \n", mesh.vn, mesh.fn,mesh.bn,(mesh.fn-mesh.bn), \
               (mesh.fn-origTarget),origTarget, DeciSession.currMetric);

        printf("\nCompleted in (%i+%i) msec\n\n", t2 - t1, t3 - t2);

    }


    virtual ~QuadricDecimator(){}
};

} /* namespace brndan022 */

#endif /* QUADRICDECIMATOR_H_ */
