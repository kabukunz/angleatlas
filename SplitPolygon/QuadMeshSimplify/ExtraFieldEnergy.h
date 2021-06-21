//This function is very complicated. So this might be bad and useless

#pragma once
#include <ANN/ANN.h>
#include <Eigen/Sparse>

namespace DataOperation
{
	struct OptimizeInfo;

	typedef struct FaceVerts
	{
		uint32_t v[4];
	}FaceVerts;

	class ExtraFieldEnergy
	{
	public:
		ExtraFieldEnergy(OptimizeInfo *oi, ANNkd_tree *kdTree, std::vector<Eigen::Vector3d> &fields, double weight0, double weight1);
		~ExtraFieldEnergy();

		void ExtractQuadVerts();
		void ExtractQuadAreas();
		void ExtractFields();
		double E(Eigen::VectorXd &currVariables);
		void ExtraNablaHessionE(Eigen::VectorXd &nabla, std::vector<Eigen::Triplet<double>> &triList, Eigen::VectorXd &currVariables, Eigen::SparseMatrix<double> & hession);

	private:
		OptimizeInfo *oi_;
		ANNkd_tree *kdTree_;
		std::vector<Eigen::Vector3d> &fields_;
		std::vector<FaceVerts> quadVerts_;
		std::vector<Eigen::Vector2d> quadMids_;
		std::vector<Eigen::Vector3d> quadFields_;
		std::vector<double> quadAreas_;

		double weight0_ = 1, weight1_ = 1;
	};
}

