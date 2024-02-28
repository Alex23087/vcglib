#ifndef OPENVDB_VCGLIB_ADAPTER_H
#define OPENVDB_VCGLIB_ADAPTER_H

#include <vcg/complex/complex.h>
#include <vcg/complex/allocate.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/util/Util.h>

double isovalue = 0.01;
double adaptivity = 1;
double voxelSize = 0.005;

namespace vcg{
namespace tri{
    
template<class TRI_MESH_TYPE>
class OpenVDBAdapter {
    typedef TRI_MESH_TYPE MeshType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename FaceType::VertexType VertexType;
    typedef	typename VertexType::ScalarType ScalarType;

    private:
    std::vector<openvdb::v11_0::Vec3s> mVtx;
    std::vector<openvdb::v11_0::Vec3I> mTri;
    std::vector<openvdb::v11_0::Vec4I> mQuad;
    openvdb::FloatGrid::Ptr grid;

    //MARK: Parameters
    double isovalue = 0;
    double adaptivity = 0;
    double voxelSize = -1;

    public:
    OpenVDBAdapter(){
        openvdb::initialize();
    }

    void setIsovalue(double isovalue){
        this->isovalue = isovalue;
    }

    void setAdaptivity(double adaptivity){
        this->adaptivity = adaptivity;
    }

    void setVoxelSize(double voxelSize){
        this->voxelSize = voxelSize;
    }

    void loadMesh(MeshType& m){
        clearVectors();

        for (size_t i = 0; i < m.vert.size(); i++){
            mVtx.emplace_back(m.vert[i].P()[0], m.vert[i].P()[1], m.vert[i].P()[2]);
        }
        for (size_t i = 0; i < m.face.size(); i++){
            mTri.emplace_back(Index(m, m.face[i].V(0)), Index(m, m.face[i].V(1)), Index(m, m.face[i].V(2)));
        }
    }

    void meshToVolume(){
        assert(mVtx.size() > 0);
        assert(mTri.size() > 0);
        assert(voxelSize > 0);

        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(voxelSize);
        grid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*xform, mVtx, mTri, 200.0f);
    }

    void volumeToMesh(MeshType& m){
        assert(grid);
        
        openvdb::tools::volumeToMesh(*grid, mVtx, mTri, mQuad, isovalue, adaptivity);

        m.Clear();
        auto vi = Allocator<MeshType>::AddVertices(m, mVtx.size());
        auto fi = Allocator<MeshType>::AddFaces(m, mTri.size() + mQuad.size() * 2);
        for (size_t i = 0; i < mVtx.size(); vi++, i++){
            vi->P()[0] = (ScalarType) mVtx[i][0];
            vi->P()[1] = (ScalarType) mVtx[i][1];
            vi->P()[2] = (ScalarType) mVtx[i][2];
        }
        for (size_t i = 0; i < mTri.size(); fi++, i++){
            fi->V(0) = &m.vert[mTri[i][2]];
            fi->V(1) = &m.vert[mTri[i][1]];
            fi->V(2) = &m.vert[mTri[i][0]];
        }
        for (size_t i = 0; i < mQuad.size(); i++){
            fi->V(0) = &m.vert[mQuad[i][2]];
            fi->V(1) = &m.vert[mQuad[i][1]];
            fi->V(2) = &m.vert[mQuad[i][0]];
            fi++;
            fi->V(0) = &m.vert[mQuad[i][3]];
            fi->V(1) = &m.vert[mQuad[i][2]];
            fi->V(2) = &m.vert[mQuad[i][0]];
            fi++;
        }
    }

    private:
    void clearVectors(){
        mVtx.clear();
        mTri.clear();
        mQuad.clear();
    }
};
}
}

#endif