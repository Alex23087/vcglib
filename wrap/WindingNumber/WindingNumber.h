#pragma once

#ifndef WINDING_NUMBER_H
#define WINDING_NUMBER_H

#include <wrap/WindingNumber/UT_SolidAngle.h>
#include <vcg/complex/complex.h>

namespace vcg::tri{
template<class TRI_MESH_TYPE>
class WindingNumber {
    typedef TRI_MESH_TYPE MeshType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename FaceType::VertexType VertexType;
    typedef	typename VertexType::ScalarType ScalarType;

	private:
	HDK_Sample::UT_SolidAngle<float,float> solid_angle; // Only works with float

	public:
	void init(MeshType &m, int order = 2)
	{
		// Initialize vector of vertex positions
		std::vector<HDK_Sample::UT_Vector3T<float>> U(m.vert.size());
		for(int i = 0; i < m.vert.size(); i++)
		{
			for(int j = 0;j<3;j++)
			{
				U[i][j] = m.vert[i].P()[j];
			}
		}

		// Initialize vector of triangle indices
		std::vector<int> mTriIndex(m.face.size() * 3);
		for(int i = 0; i < m.face.size(); i++)
		{
			for(int j = 0;j<3;j++)
			{
				mTriIndex[i*3+j] = vcg::tri::Index(m, m.face[i].V(j));
			}
		}

		solid_angle.init(m.face.size(), mTriIndex.data(), m.vert.size(), &U[0], order);
	}

	ScalarType computeWindingNumber(std::vector<ScalarType> coordV, double accuracy_scale = 2.0){
		auto pt = HDK_Sample::UT_Vector3T<float>(coordV.data());
		return solid_angle.computeSolidAngle(pt, accuracy_scale) / (4.0 * M_PI);
	}
};
}

#endif