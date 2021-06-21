#include "ExtraFieldEnergyAlone.h"
#include "AssistFunc.h"
#include "Pipeline.h"
using namespace BaseDataStructure;
const int OUT_LOOP_COUNT = 15;
const int INNER_LOOP_COUNT = 30;

namespace DataOperation{

	ExtraFieldEnergyAlone::ExtraFieldEnergyAlone(QuadMesh &qm, OptimizeInfo *oi, ANNkd_tree *kdTree, std::vector<Eigen::Vector3d> &fields, double weight0, double weight1, double targetEdgeLength) : qm_(qm), oi_(oi), kdTree_(kdTree), fields_(fields), weight0_(weight0), weight1_(weight1), targetEdgeLength_(targetEdgeLength)
	{
		BuildOiAndLocal();
		currVariables_.resize(localToOi_.size() * 2);
		currVariables_.setZero();
		nablaE_.resize(localToOi_.size() * 2);
		nablaE_.setZero();
	}


	ExtraFieldEnergyAlone::~ExtraFieldEnergyAlone()
	{
	}

	void ExtraFieldEnergyAlone::Initialize()
	{
		ExtractQuadVerts();
		ExtractQuadAreas();
		ExtractFields();
	}

	double ExtraFieldEnergyAlone::E()
	{
		std::function<double(double, double, double, double)> EachFieldEnergy = [&](double xx, double yy, double u00, double u10) -> double
		{
			double v0 = -4 * std::pow(u00, 3)*u10 + 4 * u00 * std::pow(u10, 3) + 4 * std::pow(xx, 3)*yy - 4 * xx*std::pow(yy, 3), v1 = -std::pow(u00, 4) + 6 * std::pow(u00, 2) * std::pow(u10, 2) - std::pow(u10, 4) + std::pow(xx, 4) - 6 * std::pow(xx, 2) * yy*yy + std::pow(yy, 4);
			return v0 * v0 + v1 * v1;
		};

		std::function<double(double, double, double, double)> RegularTerm = [&](double x0, double y0, double x1, double y1)->double
		{
			return x0 * x0*x1*x1 + y0 * y0*y1*y1 + 2 * x0*x1*y0*y1;
		};

		double sumV = 0;
		double x[4], y[4], x0, x1, y0, y1;
		for (int i = 0; i < quadVerts_.size(); ++i)
		{
			double u0 = targetEdgeLength_ * quadFields_[i][0], u1 = targetEdgeLength_ * quadFields_[i][1];
			for (int j = 0; j < 4; ++j)
			{
				if (oiToLocal_[quadVerts_[i][j]] == (uint32_t)-1)
				{
					x0 = oi_->localVPos[quadVerts_[i][j]][0];
					y0 = oi_->localVPos[quadVerts_[i][j]][1];
				}
				else
				{
					uint32_t realId = oiToLocal_[quadVerts_[i][j]];
					x0 = currVariables_[realId * 2];
					y0 = currVariables_[realId * 2 + 1];
				}

				if (oiToLocal_[quadVerts_[i][(j + 1) % 4]] == (uint32_t)-1)
				{
					x1 = oi_->localVPos[quadVerts_[i][(j + 1) % 4]][0];
					y1 = oi_->localVPos[quadVerts_[i][(j + 1) % 4]][1];
				}
				else
				{
					uint32_t realId = oiToLocal_[quadVerts_[i][(j+1)%4]];
					x1 = currVariables_[realId * 2];
					y1 = currVariables_[realId * 2 + 1];
				}

				x[j] = x1 - x0; y[j] = y1 - y0;
			}

			for (int j = 0; j < 4; ++j)
			{
				sumV += weight0_ * quadAreas_[i] *  EachFieldEnergy(x[j], y[j], u0, u1) + weight1_ * quadAreas_[i] * RegularTerm(x[j], y[j], x[(j+1)%4], y[(j+1)%4]) ;
			}
		}
		return sumV;
	}

	void ExtraFieldEnergyAlone::NablaE()
	{
		nablaE_.setZero();

		double x[4], y[4], x0, x1, y0, y1;
		for (int i = 0; i < quadVerts_.size(); ++i)
		{
			double u0 = targetEdgeLength_ * quadFields_[i][0], u1 = targetEdgeLength_ * quadFields_[i][1];
			for (int j = 0; j < 4; ++j)
			{
				if (oiToLocal_[quadVerts_[i][j]] == (uint32_t)-1)
				{
					x0 = oi_->localVPos[quadVerts_[i][j]][0];
					y0 = oi_->localVPos[quadVerts_[i][j]][1];
				}
				else
				{
					uint32_t realId = oiToLocal_[quadVerts_[i][j]];
					x0 = currVariables_[realId * 2];
					y0 = currVariables_[realId * 2 + 1];
				}

				if (oiToLocal_[quadVerts_[i][(j + 1) % 4]] == (uint32_t)-1)
				{
					x1 = oi_->localVPos[quadVerts_[i][(j + 1) % 4]][0];
					y1 = oi_->localVPos[quadVerts_[i][(j + 1) % 4]][1];
				}
				else
				{
					uint32_t realId = oiToLocal_[quadVerts_[i][(j + 1) % 4]];
					x1 = currVariables_[realId * 2];
					y1 = currVariables_[realId * 2 + 1];
				}

				x[j] = x1 - x0; y[j] = y1 - y0;
			}

			for (int j = 0; j < 4; ++j)
			{
				//uint32_t x0 = quadVerts_[i][(j + 1) % 4] * 2, x1 = quadVerts_[i][j] * 2, y0 = quadVerts_[i][(j + 1) % 4] * 2 + 1, y1 = quadVerts_[i][j] * 2 + 1;
				double temp0 = (-4 * std::pow(u0, 3)*u1 + 4 * u0 * std::pow(u1, 3) + 4 * std::pow(x[j], 3)*y[j] - 4 * x[j] * std::pow(y[j], 3)) * 2;
				double temp1 = (-std::pow(u0, 4) + 6 * std::pow(u0, 2) * std::pow(u1, 2) - std::pow(u1, 4) + std::pow(x[j], 4) - 6 * std::pow(x[j], 2) * y[j] * y[j] + std::pow(y[j], 4)) * 2;
				double temp2 = temp0 * (12 * x[j] * x[j] * y[j] - 4 * std::pow(y[j], 3)) + temp1 * (4 * std::pow(x[j], 3) - 12 * x[j] * y[j] * y[j]);
				double temp3 = temp0 * (4 * std::pow(x[j], 3) - 12 * x[j] * y[j] * y[j]) + temp1 * (-12 * x[j] * x[j] * y[j] + 4 * std::pow(y[j], 3));

				if (oiToLocal_[quadVerts_[i][j]] != (uint32_t)-1)
				{
					uint32_t realId = oiToLocal_[quadVerts_[i][j]];
					nablaE_[realId * 2] -= weight0_ * quadAreas_[i] * temp2;
					nablaE_[realId * 2 + 1] -= weight0_ * quadAreas_[i] * temp3;
				}
				if (oiToLocal_[quadVerts_[i][(j+1)%4]] != (uint32_t)-1)
				{
					uint32_t realId = oiToLocal_[quadVerts_[i][(j + 1) % 4]];
					nablaE_[realId * 2] += weight0_ * quadAreas_[i] * temp2;
					nablaE_[realId * 2 + 1] += weight0_ * quadAreas_[i] * temp3;
				}

			}

			for (int j = 0; j < 4; ++j)
			{
				uint32_t v0Id = quadVerts_[i][j], v1Id = quadVerts_[i][(j + 1) % 4], v2Id = quadVerts_[i][(j + 2) % 4];
				double commonV = (x[j] * x[(j + 1) % 4] + y[j] * y[(j + 1) % 4]) * 2;
				if (oiToLocal_[v0Id] != (uint32_t)-1)
				{
					uint32_t realId = oiToLocal_[v0Id];
					nablaE_[realId * 2] += weight1_ * quadAreas_[i]* commonV * (-x[(j+1)%4])/*(-2 * x[j] * x[(j + 1) % 4] * x[(j + 1) % 4] - 2 * x[(j+1)%4] * y[j] * y[(j + 1) % 4])*/;
					nablaE_[realId * 2 + 1] += weight1_ * quadAreas_[i]* commonV * (-y[(j+1)%4])/*( -2 * y[j] * y[(j + 1) % 4] * y[(j + 1) % 4] - 2 * x[j] * x[(j + 1) % 4] * y[(j + 1) % 4])*/;
				}
				if (oiToLocal_[v1Id] != (uint32_t)-1)
				{
					uint32_t realId = oiToLocal_[v1Id];
					nablaE_[realId * 2] += weight1_ * quadAreas_[i] * commonV * x[(j+1)%4]/*(2 * x[j] * x[(j + 1) % 4] * x[(j + 1) % 4] + 2 * x[(j + 1) % 4] * y[j] * y[(j + 1) % 4])*/;
					nablaE_[realId * 2 + 1] += weight1_ * quadAreas_[i] * commonV * y[(j+1)%4]/*(2 * y[j] * y[(j + 1) % 4] * y[(j + 1) % 4] + 2 * x[j] * x[(j + 1) % 4] * y[(j + 1) % 4])*/;

					nablaE_[realId * 2] += weight1_ * quadAreas_[i] * commonV * (-x[j])/*(-2 * x[j] * x[j] * x[(j + 1) % 4] - 2 * x[j] * y[j] * y[(j + 1) % 4])*/;
					nablaE_[realId * 2 + 1] += weight1_ * quadAreas_[i] * commonV * (-y[j])/*(-2 * y[j] * y[j] * y[(j + 1) % 4] - 2 * x[j] * x[(j + 1) % 4] * y[j])*/;
				}
				if (oiToLocal_[v2Id] != (uint32_t)-1)
				{
					uint32_t realId = oiToLocal_[v2Id];
					nablaE_[realId * 2] += weight1_ * quadAreas_[i] * commonV * x[j]/*(2 * x[j] * x[j] * x[(j + 1) % 4] + 2 * x[j] * y[j] * y[(j + 1) % 4])*/;
					nablaE_[realId * 2 + 1] += weight1_ * quadAreas_[i] * commonV * y[j]/*(2 * y[j] * y[j] * y[(j + 1) % 4] + 2 * x[j] * x[(j + 1) % 4] * y[j])*/;
				}
			}
		}
	}

	double ExtraFieldEnergyAlone::ComputeNoInverse(Eigen::VectorXd &pk)
	{
		std::vector<double> allAlpha;
		allAlpha.reserve(oi_->localTriIds.size());
		for (int i = 0; i < oi_->localTriIds.size(); ++i)
		{
			double tempAlpha = 0;
			int v0Id = oi_->localTriIds[i][0];
			int v1Id = oi_->localTriIds[i][1];
			int v2Id = oi_->localTriIds[i][2];

			double x0, y0, x1, y1, x2, y2, px0, py0, px1, py1, px2, py2;
			if (oiToLocal_[v0Id] != (uint32_t)-1)
			{
				x0 = currVariables_[oiToLocal_[v0Id] * 2];
				y0 = currVariables_[oiToLocal_[v0Id] * 2 + 1];
				px0 = pk[oiToLocal_[v0Id] * 2];
				py0 = pk[oiToLocal_[v0Id] * 2 + 1];
			}
			else
			{
				x0 = oi_->localVPos[v0Id][0];
				y0 = oi_->localVPos[v0Id][1];
				px0 = 0; py0 = 0;
			}
			if (oiToLocal_[v1Id] != (uint32_t)-1)
			{
				x1 = currVariables_[oiToLocal_[v1Id] * 2];
				y1 = currVariables_[oiToLocal_[v1Id] * 2 + 1];
				px1 = pk[oiToLocal_[v1Id] * 2];
				py1 = pk[oiToLocal_[v1Id] * 2 + 1];
			}
			else
			{
				x1 = oi_->localVPos[v1Id][0];
				y1 = oi_->localVPos[v1Id][1];
				px1 = 0; py1 = 0;
			}
			if (oiToLocal_[v2Id] != (uint32_t)-1)
			{
				x2 = currVariables_[oiToLocal_[v2Id] * 2];
				y2 = currVariables_[oiToLocal_[v2Id] * 2 + 1];
				px2 = pk[oiToLocal_[v2Id] * 2];
				py2 = pk[oiToLocal_[v2Id] * 2 + 1];
			}
			else
			{
				x2 = oi_->localVPos[v2Id][0];
				y2 = oi_->localVPos[v2Id][1];
				px2 = 0; py2 = 0;
			}
			double x10 = x1 - x0;
			double y10 = y1 - y0;
			double x20 = x2 - x0;
			double y20 = y2 - y0;
			double px10 = px1 - px0;
			double py10 = py1 - py0;
			double px20 = px2 - px0;
			double py20 = py2 - py0;

			double root0, root1;
			if (std::abs(px10*py20 - px20 * py10) < 1E-6)
			{
				root0 = LinearRoot(x10, y10, x20, y20, px10, py10, px20, py20);
				root1 = root0;
			}
			else
			{
				root0 = QuadraticRoot1(x10, y10, x20, y20, px10, py10, px20, py20);
				root1 = QuadraticRoot2(x10, y10, x20, y20, px10, py10, px20, py20);
			}
			if (isnan(root0) || isinf(root0))
				root0 = 1000;
			if (isnan(root1) || isinf(root1))
				root1 = 1000;
			if (root0 > root1)
				std::swap(root0, root1);	//make sure root0 <= root1.

			if (root0 > 0)
				tempAlpha = root0;
			else if (root1 > 0)
				tempAlpha = root1;
			else
				tempAlpha = 1;
			allAlpha.emplace_back(tempAlpha);
		}
		std::sort(allAlpha.begin(), allAlpha.end());
		return (allAlpha[0] > 0) ? allAlpha[0] * 0.8 : 1;
	}

	void ExtraFieldEnergyAlone::iterate(double nablaDotPk, double ek)
	{
		backupVars_ = currVariables_;
		pk_ = -nablaE_;
		double alpha = ComputeNoInverse(pk_);

		int loopCount = 0;
		currVariables_ = backupVars_ + alpha * pk_;
		double currE = E();

		while (currE > ek + 0.2 * alpha * nablaDotPk && loopCount < INNER_LOOP_COUNT)
		{
			alpha *= 0.8;
			currVariables_ = backupVars_ + alpha * pk_;
			currE = E();
			++loopCount;
		}
		if (loopCount == INNER_LOOP_COUNT && currE > ek)
		{
			std::cout << "AESolver::Iterate::INNER_LOOP_COUNT " << std::endl;
			currVariables_ = backupVars_;
		}
	}

	void ExtraFieldEnergyAlone::Optimize()
	{
		int outLoopCount = 0;
		NablaE();
		while (nablaE_.norm() > 1E-2 && outLoopCount < OUT_LOOP_COUNT)
		{
			++outLoopCount;

			double ek = E();
			double nablaENorm = nablaE_.norm();
//			std::cout << "In the " << outLoopCount << "th iterate, E = " << ek << ",	Nabla_E = " << nablaENorm << std::endl;
			
			double nablaDotPk = -nablaE_.dot(nablaE_);
//			std::cout << ", nablaDotPk is " << nablaDotPk << std::endl;

			iterate(nablaDotPk, ek);
			if (outLoopCount != OUT_LOOP_COUNT)
			{
				NablaE();
			}
		}
	}

	double ExtraFieldEnergyAlone::QuadraticRoot1(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20)
	{
		double a = px10 * py20 - px20 * py10;
		double b = px10 * y20 + py20 * x10 - px20 * y10 - py10 * x20;
		double c = x10 * y20 - x20 * y10;

		if (std::abs(a) < 1E-13 || b * b - 4 * a*c < 0)
			return 1;
		else
			return (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
	}

	double ExtraFieldEnergyAlone::QuadraticRoot2(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20)
	{
		double a = px10 * py20 - px20 * py10;
		double b = px10 * y20 + py20 * x10 - px20 * y10 - py10 * x20;
		double c = x10 * y20 - x20 * y10;

		if (std::abs(a) < 1E-13 || b * b - 4 * a*c<0)
			return 1;
		else
			return (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
	}

	double ExtraFieldEnergyAlone::LinearRoot(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20)
	{
		double b = px10 * y20 + py20 * x10 - px20 * y10 - py10 * x20;
		double c = x10 * y20 - x20 * y10;
		if (std::abs(b) < 1E-15)
			return 1;
		else
			return -c / b;
	}

	void ExtraFieldEnergyAlone::ExtractQuadVerts()
	{
		quadVerts_.clear();
		quadVerts_.reserve(oi_->localTriIds.size() / 4 + 1);
		for (int i = 0; i < oi_->localTriIds.size(); i += 4)
		{
			quadVerts_.emplace_back(Eigen::Vector4i());
			quadVerts_[i / 4][0] = oi_->localTriIds[i][0];
			quadVerts_[i / 4][1] = oi_->localTriIds[i][1];
			quadVerts_[i / 4][2] = oi_->localTriIds[i][2];
			quadVerts_[i / 4][3] = oi_->localTriIds[i+1][2];
		}
	}

	void ExtraFieldEnergyAlone::ExtractQuadAreas()
	{
		quadAreas_.clear();
		quadAreas_.reserve(quadVerts_.size());
		for (int i = 0; i < quadVerts_.size(); ++i)
		{
			quadAreas_.emplace_back(AssistFunctions::AssistFunc::ComputeAbsCross(oi_->localVPos[quadVerts_[i][1]] - oi_->localVPos[quadVerts_[i][0]], oi_->localVPos[quadVerts_[i][3]] - oi_->localVPos[quadVerts_[i][0]])*0.5 +
				AssistFunctions::AssistFunc::ComputeAbsCross(oi_->localVPos[quadVerts_[i][1]] - oi_->localVPos[quadVerts_[i][2]], oi_->localVPos[quadVerts_[i][3]] - oi_->localVPos[quadVerts_[i][2]])*0.5);
		}
	}

	void ExtraFieldEnergyAlone::ExtractFields()
	{
		quadMids_.clear();
		quadMids_.reserve(quadVerts_.size());
		for (int i = 0; i < quadVerts_.size(); ++i)
		{
			Eigen::Vector2d avePos(0, 0);
			for (int j = 0; j < 4; ++j)
			{
				avePos += oi_->localVPos[quadVerts_[i][j]];
			}
			quadMids_.emplace_back(avePos/4.0);
		}

		quadFields_.clear();
		quadFields_.reserve(quadVerts_.size());
		ANNpoint ap = annAllocPt(3);
		ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
		for (int i = 0; i < quadMids_.size(); ++i)
		{
			ap[0] = quadMids_[i][0]; ap[1] = quadMids_[i][1]; ap[2] = 0;
			kdTree_->annkSearch(ap, 1, nnIdx, dists);
			quadFields_.emplace_back(fields_[nnIdx[0]].normalized());
		}

		annDeallocPt(ap);
		delete[] nnIdx; delete[] dists;
	}

	void ExtraFieldEnergyAlone::BuildOiAndLocal()
	{
		bool *oToL = new bool[oi_->localToGlobalVVec.size()];
		std::memset(oToL, 0, oi_->localToGlobalVVec.size() * sizeof(bool));
		for (int i = 0; i < oi_->localToGlobalVVec.size(); ++i)
		{
			if (qm_.Vs_[oi_->localToGlobalVVec[i]].boundary)
				oToL[i] = true;
		}
		for (int i = 0; i < oi_->localVsNotMove.size(); ++i)
		{
			oToL[i] = true;
		}

		oiToLocal_.clear();
		localToOi_.clear();
		oiToLocal_.resize(oi_->localToGlobalVVec.size());
		localToOi_.reserve(oi_->localToGlobalVVec.size());

		uint32_t currC = 0;
		for (int i = 0; i < oi_->localToGlobalVVec.size(); ++i)
		{
			if (oToL[i])
			{
				oiToLocal_[i] = (uint32_t)-1;
			}
			else
			{
				oiToLocal_[i] = currC++;
				localToOi_.emplace_back(i);
			}
		}

		delete[] oToL;
	}

	void ExtraFieldEnergyAlone::SetValueFront(OptimizeInfo *oi)
	{
		for (int i = 0; i < oi->localVPos.size(); ++i)
		{
			uint32_t realId = oiToLocal_[i];
			if (realId != (uint32_t)-1)
			{
				currVariables_[2 * realId] = oi->localVPos[i][0];
				currVariables_[2 * realId + 1] = oi->localVPos[i][1];
			}
		}
	}

	void ExtraFieldEnergyAlone::SetValueBack(OptimizeInfo *oi)
	{
		for (int i = 0; i < oi->localVPos.size(); ++i)
		{
			uint32_t realId = oiToLocal_[i];
			if (realId != (uint32_t)-1)
			{
				oi->localVPos[i][0] = currVariables_[2 * realId];
				oi->localVPos[i][1] = currVariables_[2 * realId + 1];
			}
		}
	}

}