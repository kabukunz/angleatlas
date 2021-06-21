#pragma once
#include <ANN/ANN.h>
#include <Eigen/Sparse>
#include "QuadMesh.h"

namespace DataOperation
{
	struct OptimizeInfo;

	class ExtraFieldEnergyAlone
	{
	public:
		ExtraFieldEnergyAlone(BaseDataStructure::QuadMesh &qm, OptimizeInfo *oi, ANNkd_tree *kdTree, std::vector<Eigen::Vector3d> &fields, double weight0, double weight1, double targetEdgeLength);
		~ExtraFieldEnergyAlone();

		void ExtractQuadVerts();
		void ExtractQuadAreas();
		void ExtractFields();
		double E();
		void NablaE();
		void Optimize();
		void iterate(double nablaDotPk, double ek);	

		void BuildOiAndLocal();

		void Initialize();

		double QuadraticRoot1(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20);
		double QuadraticRoot2(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20);
		double LinearRoot(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20);
		double ComputeNoInverse(Eigen::VectorXd &pk);

		void SetValueFront(OptimizeInfo *oi);
		void SetValueBack(OptimizeInfo *oi);

	public:
		BaseDataStructure::QuadMesh &qm_;
		OptimizeInfo *oi_;
		ANNkd_tree *kdTree_;
		std::vector<Eigen::Vector3d> &fields_;
		std::vector<Eigen::Vector4i> quadVerts_;
		std::vector<Eigen::Vector2d> quadMids_;
		std::vector<Eigen::Vector3d> quadFields_;
		std::vector<double> quadAreas_;

		double weight0_ = 1, weight1_ = 1, targetEdgeLength_ = 0.3;
		
		std::vector<uint32_t> oiToLocal_;
		std::vector<uint32_t> localToOi_;

		Eigen::VectorXd currVariables_, nablaE_, backupVars_, pk_;
	};
}