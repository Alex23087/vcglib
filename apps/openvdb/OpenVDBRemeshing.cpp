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
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/openvdb/OpenVDBAdapter.cpp>

using namespace vcg;
using namespace std;

class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<
	Use<MyVertex>   ::AsVertexType,
	Use<MyFace>     ::AsFaceType
>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::BitFlags>{};
class MyFace    : public Face<MyUsedTypes, face::VertexRef, vertex::BitFlags> {};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace>> {};

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

	assert(original.VN()>0);

	// OpenVDB remeshing parameters
	double targetLenPerc = 0.2;
	double isovalue = 0.0;
	double adaptivity = 0.0;
	bool useLevelSet = false;	// Use level set instead of volume
	if(argc>=3) targetLenPerc = atof(argv[2]);
	if(argc>=4) isovalue = atof(argv[3]);
	if(argc>=5) adaptivity = atof(argv[4]);
	if(argc>=6) useLevelSet = atoi(argv[5]);
	
	double voxelSize = targetLenPerc * (original.bbox.Diag() / 100.0);

	// Mesh cleaning
	tri::Clean<MyMesh>::RemoveUnreferencedVertex(original);
	vcg::tri::Allocator<MyMesh>::CompactEveryVector(original);
	tri::UpdateBounding<MyMesh>::Box(original);

	// OpenVDB mesh to volume
	vcg::tri::OpenVDBAdapter<MyMesh> adapter;

	printf(" Input mesh %8i v %8i f\n",original.VN(),original.FN());

	adapter.setMesh(&original);
	adapter.setVoxelSize(voxelSize);
	if(useLevelSet){
		adapter.meshToLevelSet();
	}else{
		adapter.meshToVolume();
	}

	// OpenVDB volume to mesh
	adapter.setIsovalue(isovalue);
	adapter.setAdaptivity(adaptivity);
	adapter.volumeToMesh(toremesh);

	assert(toremesh.VN()>0);

	tri::Clean<MyMesh>::RemoveUnreferencedVertex(toremesh);
	vcg::tri::Allocator<MyMesh>::CompactEveryVector(toremesh);

	printf("Output mesh %8i v %8i f\n",toremesh.VN(),toremesh.FN());
	vcg::tri::io::ExporterOBJ<MyMesh>::Save(toremesh, "remesh.obj", 0);

	return 0;
}
