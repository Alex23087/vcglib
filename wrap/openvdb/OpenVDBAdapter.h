#pragma once

#ifndef OPENVDB_VCGLIB_ADAPTER_H
#define OPENVDB_VCGLIB_ADAPTER_H

#include <vcg/complex/complex.h>
#include <vcg/complex/allocate.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/util/Util.h>
#include <wrap/WindingNumber/WindingNumber.h>

namespace vcg::tri{
    
template<class TRI_MESH_TYPE>
class OpenVDBAdapter {
    private:
    typedef TRI_MESH_TYPE MeshType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename FaceType::VertexType VertexType;
    typedef	typename VertexType::ScalarType ScalarType;
    using GridType = openvdb::Grid<openvdb::v11_0::tree::Tree<openvdb::v11_0::tree::RootNode<openvdb::v11_0::tree::InternalNode<openvdb::v11_0::tree::InternalNode<openvdb::v11_0::tree::LeafNode<ScalarType, 3U>, 4U>, 5U>>>>;
    class MeshTypeDataAdapter;


    private:
    std::vector<openvdb::v11_0::math::Vec3<ScalarType>> mVtx;
    std::vector<openvdb::v11_0::Vec3I> mTri;
    std::vector<openvdb::v11_0::Vec4I> mQuad;
    MeshType *m = nullptr;
    typename GridType::Ptr grid;
    WindingNumber<MeshType> windingNumber;
    MeshTypeDataAdapter meshDataAdapter;

    //MARK: Parameters
    double isovalue = 0;
    double adaptivity = 0;
    double voxelSize = -1;

    public:
    OpenVDBAdapter(){
        openvdb::initialize();
        meshDataAdapter = MeshTypeDataAdapter();
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

    void setMesh(MeshType* m){
        this->m = m;
    }

    void meshToVolume(){
        assert(voxelSize > 0);

        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(voxelSize);

        meshDataAdapter.setMesh(m);
        meshDataAdapter.setTransform(xform);

        windingNumber.init(*m);
        auto interiorTest = [xform, &windingNumber = windingNumber, voxelSize = voxelSize](const openvdb::Coord &coord) -> bool
        {
            auto worldCoord = xform->indexToWorld(coord);
            // Using float instead of ScalarType because the winding number lib uses float internally.
            auto coordV = std::vector<float>{
                static_cast<float>(worldCoord.x()),
                static_cast<float>(worldCoord.y()),
                static_cast<float>(worldCoord.z())
            };
            auto wn = windingNumber.computeWindingNumber(coordV, voxelSize);
            return fabs(wn) >= 0.5f ? true : false;
        };

        float outerBand = (isovalue > 0 ? isovalue/voxelSize : 0.0f) + 0.5f;
        float innerBand = (isovalue < 0 ? (-isovalue)/voxelSize : 0.0f) + 0.5f;

        // NullInterrupter needed because openvdb does not correctly pass parameters otherwise.
        openvdb::util::NullInterrupter nullInterrupter;
        grid = openvdb::tools::meshToVolume<GridType>(nullInterrupter, meshDataAdapter, *xform, outerBand, innerBand, 0, nullptr, interiorTest, openvdb::tools::EVAL_EVERY_TILE);
    }

    void meshToLevelSet(){
        copyMeshVectors();

        assert(mVtx.size() > 0);
        assert(mTri.size() > 0);
        assert(voxelSize > 0);

        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(voxelSize);

        // The band half-width is dynamically computed as abs(isovalue/voxelSize) + 1.0f, which should be the smallest value to avoid artifacts.
        grid = openvdb::tools::meshToLevelSet<GridType>(*xform, mVtx, mTri, std::fabs(isovalue/voxelSize) + 1.0f);
    }

    void volumeToMesh(MeshType& outMesh){
        assert(grid);
        
        openvdb::tools::volumeToMesh(*grid, mVtx, mTri, mQuad, isovalue, adaptivity);

        outMesh.Clear();
        auto vi = Allocator<MeshType>::AddVertices(outMesh, mVtx.size());
        auto fi = Allocator<MeshType>::AddFaces(outMesh, mTri.size() + mQuad.size() * 2);
        for (size_t i = 0; i < mVtx.size(); vi++, i++){
            vi->P()[0] = (ScalarType) mVtx[i][0];
            vi->P()[1] = (ScalarType) mVtx[i][1];
            vi->P()[2] = (ScalarType) mVtx[i][2];
        }
        for (size_t i = 0; i < mTri.size(); fi++, i++){
            fi->V(0) = &outMesh.vert[mTri[i][2]];
            fi->V(1) = &outMesh.vert[mTri[i][1]];
            fi->V(2) = &outMesh.vert[mTri[i][0]];
        }
        for (size_t i = 0; i < mQuad.size(); i++){
            fi->V(0) = &outMesh.vert[mQuad[i][2]];
            fi->V(1) = &outMesh.vert[mQuad[i][1]];
            fi->V(2) = &outMesh.vert[mQuad[i][0]];
            fi++;
            fi->V(0) = &outMesh.vert[mQuad[i][3]];
            fi->V(1) = &outMesh.vert[mQuad[i][2]];
            fi->V(2) = &outMesh.vert[mQuad[i][0]];
            fi++;
        }
    }

    private:
    void copyMeshVectors(){
        clearVectors();

        for (size_t i = 0; i < m->vert.size(); i++){
            mVtx.emplace_back(m->vert[i].P()[0], m->vert[i].P()[1], m->vert[i].P()[2]);
        }
        for (size_t i = 0; i < m->face.size(); i++){
            mTri.emplace_back(Index(*m, m->face[i].V(0)), Index(*m, m->face[i].V(1)), Index(*m, m->face[i].V(2)));
        }
    }

    void clearVectors(){
        mVtx.clear();
        mTri.clear();
        mQuad.clear();
    }
};



template <class TRI_MESH_TYPE>
class OpenVDBAdapter<TRI_MESH_TYPE>::MeshTypeDataAdapter{
        private:
        MeshType *m = nullptr;
        openvdb::math::Transform::Ptr xform;

        public:
        MeshTypeDataAdapter(){}

        MeshType* getMesh(){
            return m;
        }

        void setMesh(MeshType* mesh){
            m = mesh;
        }

        void setTransform(openvdb::math::Transform::Ptr xform){
            this->xform = xform;
        }

        // Total number of polygons
        size_t polygonCount() const {
            return m->face.size();
        }

        // Total number of points
        size_t pointCount() const {
            return m->vert.size();
        }
        
        // Vertex count for polygon n
        size_t vertexCount(size_t n) const {
            return 3;
        }

        // Return position pos in local grid index space for polygon n and vertex v
        void getIndexSpacePoint(size_t n, size_t v, openvdb::Vec3d& pos) const {
            pos = openvdb::Vec3d(m->face[n].V(v)->P()[0], m->face[n].V(v)->P()[1], m->face[n].V(v)->P()[2]);
            pos = xform->worldToIndex(pos);
        }
    };
}

#endif