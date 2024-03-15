/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2024                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/*! \file OpenvdbRemeshing.cpp
\ingroup apps

\brief Remeshing a mesh using the OpenVDB library.

*/
#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>
#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/openvdb/OpenVDBAdapter.h>

using namespace vcg;
using namespace std;

class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<
	Use<MyVertex>   ::AsVertexType,
	Use<MyFace>     ::AsFaceType
>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Color4b, vertex::Normal3f, vertex::BitFlags>{};
class MyFace    : public Face<MyUsedTypes, face::VertexRef, vertex::BitFlags> {};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace>> {};

void testWindingAngle(MyMesh &m)
{
    tri::WindingNumber<MyMesh> windingNumber;
    windingNumber.init(m);
    printf("Testing Winding Angle\n");
    int t0=clock();
    int sampleNum=10000;
    MyMesh montecarloMesh;
    float offset= m.bbox.Diag() * 0.01f;
    typedef tri::MeshSampler<MyMesh> BaseSampler;
    tri::MeshSampler<MyMesh> mcSampler(montecarloMesh);
    tri::SurfaceSampling<MyMesh,BaseSampler>::Montecarlo(m, mcSampler, sampleNum);
    
    for(int i=0;i<montecarloMesh.VN();++i)
    {
        montecarloMesh.vert[i].P()[0]+=offset;
        Point3f p = {montecarloMesh.vert[i].P()[0],montecarloMesh.vert[i].P()[1],montecarloMesh.vert[i].P()[2]};

        auto wa = windingNumber.computeWindingNumber(p,2.0f);
        if (wa<0.5) montecarloMesh.vert[i].C()=Color4b(Color4b::Red);
        else montecarloMesh.vert[i].C()=Color4b(Color4b::Green);            
    }
    int t1=clock();
    printf("Evaluated %i samples in %f secs\n",sampleNum,(t1-t0)/float(CLOCKS_PER_SEC));
    
    // Save the mesh with montecarlo samples as ply with color
    tri::io::ExporterPLY<MyMesh>::Save(montecarloMesh,"montecarlo.ply",tri::io::Mask::IOM_VERTCOLOR);
    
}
int main( int argc, char **argv )
{
	MyMesh original,toremesh;
	if(argc<2)
	{
		printf("Usage: openvdb-remesher <filename> [voxelSize] [isovalue] [adaptivity]\n\tvoxelSize is a percentage of the bounding box diagonal\n");
		exit(0);
	}

	int lmask = 0;
	if(tri::io::ImporterOBJ<MyMesh>::Open(original, argv[1], lmask)!=0)
	{
		printf("Error reading file  %s\n",argv[1]);
		exit(0);
	}
    // Mesh cleaning
    tri::Clean<MyMesh>::RemoveUnreferencedVertex(original);
    vcg::tri::Allocator<MyMesh>::CompactEveryVector(original);
    tri::UpdateBounding<MyMesh>::Box(original);
    printf(" Input mesh %8i v %8i f\n",original.VN(),original.FN());
    testWindingAngle(original);
	assert(original.VN()>0);

	// OpenVDB remeshing parameters
	double targetLenPerc = 1.2;
	double isovalue = 0.0;
	double adaptivity = 0.0;
	bool useLevelSet = false;	// Use level set instead of volume
	if(argc>=3) targetLenPerc = atof(argv[2]);
	if(argc>=4) isovalue = atof(argv[3]);
	if(argc>=5) adaptivity = atof(argv[4]);
    if(argc>=6) useLevelSet = true;
	
	double voxelSize = targetLenPerc * (original.bbox.Diag() / 100.0);
    printf("Voxel Size %f\n",voxelSize);
    printf("Box size %.3f %.3f %.3f  - %i x %i x %i\n",original.bbox.DimX(),original.bbox.DimY(),original.bbox.DimZ(),
           int(original.bbox.DimX()/voxelSize),int(original.bbox.DimY()/voxelSize),int(original.bbox.DimZ()/voxelSize));
    

	// OpenVDB mesh to volume
	vcg::tri::OpenVDBAdapter<MyMesh> adapter;
	
	adapter.setMesh(&original);
	adapter.setVoxelSize(voxelSize);
    adapter.setIsovalue(isovalue);
    
    if(useLevelSet){
        printf("Building LevelSet\n");
		adapter.meshToLevelSet();
	}else{
        printf("Building Volume using winding number\n");
		adapter.meshToVolume();
	}

	// OpenVDB volume to mesh
	adapter.setIsovalue(isovalue);
	adapter.setAdaptivity(adaptivity);
	adapter.volumeToMesh(toremesh);

	tri::Clean<MyMesh>::RemoveUnreferencedVertex(toremesh);
	vcg::tri::Allocator<MyMesh>::CompactEveryVector(toremesh);

	printf("Output mesh %8i v %8i f\n",toremesh.VN(),toremesh.FN());
	vcg::tri::io::ExporterOBJ<MyMesh>::Save(toremesh, "remesh.obj", 0);

	return 0;
}
