#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Sparse>

namespace AssistFunctions
{

	class AssistFunc
	{
	public:
		AssistFunc();
		~AssistFunc();

		template <typename T>
		static void DeleteRepeatElements(const std::vector<T> &valueVec)
		{
			std::sort(valueVec.begin(), valueVec.end());
			valueVec.erase(std::unique(valueVec.begin(), valueVec.end()), valueVec.end());
		}

		static void ComputeSVDDecomposition(double *inputMat, double *U, double &s1, double &s2, double *V);	//Calculate the SVD decomposition. Use the colume first express, and inputMat = U * sigmaMat * V^T
		static void ComputeMinVPos(double x0, double y0, double x1, double y1, double theta, double &res0, double &res1);
		static double ComputeAbsCross(const Eigen::Vector2d &vec0, const Eigen::Vector2d &vec1);
	};

}
