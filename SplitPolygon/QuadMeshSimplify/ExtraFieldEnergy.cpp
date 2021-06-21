#include "ExtraFieldEnergy.h"
#include "Pipeline.h"
#include "AssistFunc.h"

namespace DataOperation
{

	ExtraFieldEnergy::ExtraFieldEnergy(OptimizeInfo *oi, ANNkd_tree *kdTree, std::vector<Eigen::Vector3d> &fields, double weight0, double weight1) : oi_(oi), kdTree_(kdTree), fields_(fields), weight0_(weight0), weight1_(weight1)
	{
		ExtractQuadVerts();
		ExtractQuadAreas();
	}


	ExtraFieldEnergy::~ExtraFieldEnergy()
	{
	}

	void ExtraFieldEnergy::ExtractFields()
	{
		quadMids_.clear();
		quadMids_.reserve(quadVerts_.size());
		for (int i = 0; i < quadVerts_.size(); ++i)
		{
			Eigen::Vector2d avePos(0, 0);
			for (int j = 0; j < 4; ++j)
			{
				avePos += oi_->localVPos[quadVerts_[i].v[j]];
			}
			quadMids_.emplace_back(avePos / 4.0);
		}

		quadFields_.clear();
		quadFields_.reserve(quadVerts_.size());
		ANNpoint ap = annAllocPt(3);
		ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
		for (int i = 0; i < quadVerts_.size(); ++i)
		{
			ap[0] = quadMids_[i][0]; ap[1] = quadMids_[i][1]; ap[2] = 0;
			kdTree_->annkSearch(ap, 1, nnIdx, dists);
			quadFields_.emplace_back(fields_[nnIdx[0]]);
		}

		annDeallocPt(ap);
		delete[] nnIdx; delete[] dists;
	}

	void ExtraFieldEnergy::ExtractQuadVerts()
	{
		quadVerts_.clear();
		quadVerts_.reserve(oi_->localTriIds.size() / 4+1);
		for (int i = 0; i < oi_->localTriIds.size(); i += 4)
		{
			quadVerts_.emplace_back(FaceVerts());
			quadVerts_[i/4].v[0] = oi_->localTriIds[i][0]; 
			quadVerts_[i/4].v[1] = oi_->localTriIds[i][1];
			quadVerts_[i/4].v[2] = oi_->localTriIds[i][2];
			quadVerts_[i/4].v[3] = oi_->localTriIds[i + 1][2];
		}
	}

	void ExtraFieldEnergy::ExtractQuadAreas()
	{
		quadAreas_.clear();
		quadAreas_.reserve(quadVerts_.size());
		for (int i = 0; i < quadVerts_.size(); ++i)
		{
			quadAreas_.emplace_back(AssistFunctions::AssistFunc::ComputeAbsCross(oi_->localVPos[quadVerts_[i].v[1]] - oi_->localVPos[quadVerts_[i].v[0]], oi_->localVPos[quadVerts_[i].v[3]] - oi_->localVPos[quadVerts_[i].v[0]])*0.5 +
				AssistFunctions::AssistFunc::ComputeAbsCross(oi_->localVPos[quadVerts_[i].v[1]] - oi_->localVPos[quadVerts_[i].v[2]], oi_->localVPos[quadVerts_[i].v[3]] - oi_->localVPos[quadVerts_[i].v[2]])*0.5);
		}
	}

	double ExtraFieldEnergy::E(Eigen::VectorXd &currVariables)
	{
		std::function<double(double, double, double, double)> EachFieldEnergy = [&](double xx, double yy, double u00, double u10) -> double
		{
			double v0 = -4 * std::pow(u00, 3)*u10 + 4 * u00 * std::pow(u10, 3) + 4 * std::pow(xx, 3)*yy - 4 * xx*std::pow(yy, 3), v1 = -std::pow(u00,4) + 6 * std::pow(u00,2) * std::pow(u10,2) - std::pow(u10,4) + std::pow(xx,4) - 6 * std::pow(xx,2) * yy*yy + std::pow(yy,4);
			return v0 * v0 + v1 * v1;
		};


		double sumV = 0;
		for (int i = 0; i < quadVerts_.size(); ++i)
		{
			double x[4], y[4];
			//double  
			for (int j = 0; j < 4; ++j)
			{
				x[j] = currVariables[quadVerts_[i].v[(j+1)%4] * 2] - currVariables[quadVerts_[i].v[j] * 2];
				y[j] = currVariables[quadVerts_[i].v[(j+1)%4] * 2 + 1] - currVariables[quadVerts_[i].v[j] * 2 + 1];
			}

			double u0 = quadFields_[i][0], u1 = quadFields_[i][1];
			for (int j = 0; j < 4; ++j)
			{
				sumV += weight0_ * quadAreas_[i] * EachFieldEnergy(x[j], y[j], u0, u1);
			}
		}
		
		return sumV;
	}

	void ExtraFieldEnergy::ExtraNablaHessionE(Eigen::VectorXd &nabla, std::vector<Eigen::Triplet<double>> &triList, Eigen::VectorXd &currVariables, Eigen::SparseMatrix<double> &hession)
	{
		//std::vector<Eigen::Triplet<double>> currTriList;
		//currTriList = triList;

		for (int i = 0; i < quadVerts_.size(); ++i)
		{
			double x[4], y[4];
			//double  
			for (int j = 0; j < 4; ++j)
			{
				x[j] = currVariables[quadVerts_[i].v[(j + 1) % 4] * 2] - currVariables[quadVerts_[i].v[j] * 2];
				y[j] = currVariables[quadVerts_[i].v[(j + 1) % 4] * 2 + 1] - currVariables[quadVerts_[i].v[j] * 2 + 1];
			}

			double u0 = quadFields_[i][0], u1 = quadFields_[i][1];
			for (int j = 0; j < 4; ++j)
			{
				uint32_t x0 = quadVerts_[i].v[(j + 1) % 4] * 2, x1 = quadVerts_[i].v[j] * 2, y0 = quadVerts_[i].v[(j + 1) % 4] * 2 + 1, y1 = quadVerts_[i].v[j] * 2 + 1;
				double temp0 = (-4 * std::pow(u0, 3)*u1 + 4 * u0 * std::pow(u1, 3) + 4 * std::pow(x[j], 3)*y[j] - 4 * x[j]*std::pow(y[j], 3)) * 2;
				double temp1 = (-std::pow(u0, 4) + 6 * std::pow(u0, 2) * std::pow(u1, 2) - std::pow(u1, 4) + std::pow(x[j], 4) - 6 * std::pow(x[j], 2) * y[j]*y[j] + std::pow(y[j], 4)) * 2;
				double temp2 = temp0 * (12 * x[j] * x[j] * y[j] - 4 * std::pow(y[j], 3)) + temp1 * (4 * std::pow(x[j], 3) - 12 * x[j] * y[j] * y[j]);
				nabla[x0] += weight0_ * quadAreas_[i] * temp2;
				nabla[x1] -= weight0_ * quadAreas_[i] * temp2;

				double temp3 = temp0 * (4 * std::pow(x[j], 3) - 12 * x[j] * y[j] * y[j]) + temp1 * (-12 * x[j] * x[j] * y[j] + 4 * std::pow(y[j], 3));
				nabla[y0] += weight0_ * quadAreas_[i] * temp3;
				nabla[y1] -= weight0_ * quadAreas_[i] * temp3;

				/*double hX0X0 = temp0 * (24 * x[j] * y[j])
					+ 2 * std::pow(12 * x[j] * x[j] * y[j] - 4 * std::pow(y[j], 3), 2)
					+ temp1 * (12 * std::pow(x[j], 2) - 12 * std::pow(y[j], 2))
					+ 2 * std::pow(4 * std::pow(x[j], 3) - 12 * x[j] * y[j] * y[j], 2);

				double hX0X1 = -hX0X0;
				double hX0Y0 = temp0 * (12 * std::pow(x[j], 2) - 12 * std::pow(y[j], 2))
					+ 2 * (4 * pow(x[j], 3) - 12 * x[j] * std::pow(y[j], 2)) * (12 * x[j] * x[j] * y[j] - 4 * std::pow(y[j], 3))
					+ temp1 * (-24 * x[j] * y[j])
					+ 2 * (-12 * std::pow(x[j], 2) * y[j] + 4 * std::pow(y[j], 3)) * (4 * std::pow(x[j], 3) - 12 * x[j] * y[j] * y[j]);
				double hX0Y1 = -hX0Y0, hX1X1 = hX0X0, hX1Y0 = -hX0Y0, hX1Y1 = hX0Y0;
				double hY0Y0 = temp0 * (-24 * x[j] * y[j])
					+ 2 * std::pow(4 * std::pow(x[j], 3) - 12 * x[j] * std::pow(y[j], 2), 2)
					+ temp1 * (-12 * std::pow(x[j], 2) + 12 * std::pow(y[j], 3))
					+ 2 * std::pow(-12 * std::pow(x[j], 2) * y[j] + 4 * std::pow(y[j], 3), 2);
				double hY0Y1 = -hY0Y0, hY1Y1 = hY0Y0;

				triList.emplace_back(Eigen::Triplet<double>(x0, x0, weight0_ * quadAreas_[i] * hX0X0));
				triList.emplace_back(Eigen::Triplet<double>(x1, x1, weight0_ * quadAreas_[i] * hX1X1));
				triList.emplace_back(Eigen::Triplet<double>(y0, y0, weight0_ * quadAreas_[i] * hY0Y0));
				triList.emplace_back(Eigen::Triplet<double>(y1, y1, weight0_ * quadAreas_[i] * hY1Y1));

				triList.emplace_back(Eigen::Triplet<double>(x0, x1, weight0_ * quadAreas_[i] * hX0X1));
				triList.emplace_back(Eigen::Triplet<double>(x1, x0, weight0_ * quadAreas_[i] * hX0X1));
				triList.emplace_back(Eigen::Triplet<double>(x0, y0, weight0_ * quadAreas_[i] * hX0Y0));
				triList.emplace_back(Eigen::Triplet<double>(y0, x0, weight0_ * quadAreas_[i] * hX0Y0));
				triList.emplace_back(Eigen::Triplet<double>(x0, y1, weight0_ * quadAreas_[i] * hX0Y1));
				triList.emplace_back(Eigen::Triplet<double>(y1, x0, weight0_ * quadAreas_[i] * hX0Y1));

				triList.emplace_back(Eigen::Triplet<double>(x1, y0, weight0_ * quadAreas_[i] * hX1Y0));
				triList.emplace_back(Eigen::Triplet<double>(y0, x1, weight0_ * quadAreas_[i] * hX1Y0));
				triList.emplace_back(Eigen::Triplet<double>(x1, y1, weight0_ * quadAreas_[i] * hX1Y1));
				triList.emplace_back(Eigen::Triplet<double>(y1, x1, weight0_ * quadAreas_[i] * hX1Y1));

				triList.emplace_back(Eigen::Triplet<double>(y0, y1, weight0_ * quadAreas_[i] * hY0Y1));
				triList.emplace_back(Eigen::Triplet<double>(y1, y0, weight0_ * quadAreas_[i] * hY0Y1));*/
			}
		}

		//hession.setFromTriplets(triList.begin(), triList.end());
		//double tempV = nabla.transpose() * hession * nabla;
		//if (1/*tempV < 0*/)
		//{
		//	for (int j = 0; j < oi_->localVPos.size() * 2; ++j)
		//	{
		//		currTriList.emplace_back(Eigen::Triplet<double>(j, j, 1));
		//		hession.setFromTriplets(currTriList.begin(), currTriList.end());
		//	}
		//}
		hession.setFromTriplets(triList.begin(), triList.end());
		hession.makeCompressed();
	}
}
