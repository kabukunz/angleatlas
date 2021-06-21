#include "AESolver.h"
#include "Pipeline.h"
#include "AssistFunc.h"
#include <Eigen/SVD>

#define USE_NON_FIXED_TARGET_ENERGY 1
#define USE_FIELD_ENERGY 0
#define USE_ALPHA_VEC 0

namespace DataOperation
{
	const int OUT_LOOP_COUNT = 1;
	const int INNER_LOOP_COUNT = 30;
	double LAMBDA_D = 1;
	double LAMBDA_B = 1E6;
	double LAMBDA_T = 1E4;
	const double DEGEN_TRI_THRESHOLD = 1E-16;
	AESolver::AESolver(OptimizeInfo *oi, ANNkd_tree *kdTree, std::vector<Eigen::Vector3d> &fields) : oi_(oi)//, efe_(oi, kdTree, fields, 10000, 1)
	{
		eleNum_ = oi_->localVPos.size() * 2 + oi_->localFc.ids_L.size();
		vVarNum_ = oi_->localVPos.size() * 2;
		currentVariables_.resize(eleNum_);
		currentVariables_.setZero();
		for (int i = 0; i < oi_->localVPos.size(); ++i)
		{
			currentVariables_[2 * i] = oi_->localVPos[i][0];
			currentVariables_[2 * i + 1] = oi_->localVPos[i][1];
		}

		standardInverseMatrix_ = new double[oi_->localTriIds.size()][4];
		//sim1Num_ = oi_->localTriIds.size();
		BuildStandardInverseMatrix();

		triAreas_ = new double[oi_->localTriIds.size()];
		BuildTriAreas();

		nablaE_.resize(eleNum_);
		hessionE_.resize(eleNum_, eleNum_);
		std::fill(leftMat[0], leftMat[0] + 6 * 4, 0);
#if USE_ALPHA_VEC
		alphaPtr_ = new double[eleNum_];
#endif
	}


	AESolver::~AESolver()
	{
		//SafeDeletePtr(triAreas_);
		if (triAreas_ != NULL)
		{
			delete[] triAreas_;
			triAreas_ = NULL;
		}
		//SafeDeletePtr(standardInverseMatrix_);
		if (standardInverseMatrix_ != NULL)
		{
			delete[] standardInverseMatrix_;
			standardInverseMatrix_ = NULL;
		}
	}

	double AESolver::E()
	{
		et_ = LAMBDA_T * E_T();
		eb_ = LAMBDA_B * E_B();
		ed_ = LAMBDA_D * E_D();
#if USE_FIELD_ENERGY
		ef_ = efe_.E(currentVariables_);
		return ed_ + eb_ + et_ + ef_;
#else
		return ed_ + eb_ + et_;
#endif
	}

	double AESolver::E_B()
	{
		double sum = 0;
		for (int i = 0; i < oi_->localFc.ids_C.size(); ++i)
		{
			uint32_t currC = oi_->localFc.ids_C[i];
			sum += (currentVariables_[2 * currC] - oi_->localFc.C[i][0])*(currentVariables_[2 * currC] - oi_->localFc.C[i][0])
				+ (currentVariables_[2 * currC + 1] - oi_->localFc.C[i][1])*(currentVariables_[2 * currC + 1] - oi_->localFc.C[i][1]);
		}
		for (int i = 0; i < oi_->localFc.ids_L.size(); ++i)
		{
			int currVId = oi_->localFc.ids_L[i];
			double x0 = (currentVariables_[2 * currVId] - oi_->localFc.origin_L[i][0] - currentVariables_[vVarNum_ + i] * oi_->localFc.axa_L[i][0]);
			double y0 = (currentVariables_[2 * currVId + 1] - oi_->localFc.origin_L[i][1] - currentVariables_[vVarNum_ + i] * oi_->localFc.axa_L[i][1]);
			sum += x0 * x0 + y0 * y0;
		}
		return sum;
	}

	double AESolver::E_D()
	{
		double sum = 0;
		for (int i = 0; i < oi_->localTriIds.size(); ++i)
		{
			double currSyDir = ComputeSyDir(i);
			sum += triAreas_[i] * currSyDir;
		}
		return sum;
	}

	double AESolver::E_T()
	{
		double sum = 0;
#if !USE_NON_FIXED_TARGET_ENERGY
		for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			for (int j = 0; j < currVsGroup.size(); ++j)
			{
				sum += (currentVariables_[2 * currVsGroup[j]] - oi_->localVsCoords[i][0])*(currentVariables_[2 * currVsGroup[j]] - oi_->localVsCoords[i][0])
					+ (currentVariables_[2 * currVsGroup[j] + 1] - oi_->localVsCoords[i][1])*(currentVariables_[2 * currVsGroup[j] + 1] - oi_->localVsCoords[i][1]);
			}
		}

#else
#pragma region Old_Version_Without_Boundary_Judge
		/*for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			for (int k = 0; k < currVsGroup.size(); ++k)
			{
				for (int j = 0; j < currVsGroup.size(); ++j)
				{

					if (j != k)
					{
						sum += (currentVariables_[2 * currVsGroup[j]] - currentVariables_[2 * currVsGroup[k]]) * (currentVariables_[2 * currVsGroup[j]] - currentVariables_[2 * currVsGroup[k]])
							+ (currentVariables_[2 * currVsGroup[j] + 1] - currentVariables_[2 * currVsGroup[k] + 1]) * (currentVariables_[2 * currVsGroup[j] + 1] - currentVariables_[2 * currVsGroup[k] + 1]);
					}
				}
			}
		}*/
#pragma endregion
		for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			if (false/*oi_->vsGroupHasBoundary[i]*/)
			{
				for (int j = 0; j < currVsGroup.size(); ++j)
				{
					sum += (currentVariables_[2 * currVsGroup[j]] - oi_->localVsCoords[i][0])*(currentVariables_[2 * currVsGroup[j]] - oi_->localVsCoords[i][0])
						+ (currentVariables_[2 * currVsGroup[j] + 1] - oi_->localVsCoords[i][1])*(currentVariables_[2 * currVsGroup[j] + 1] - oi_->localVsCoords[i][1]);
				}
			}
			else
			{
				/*for (int j = 1; j < currVsGroup.size(); ++j)
				{
					sum += (currentVariables_[2 * currVsGroup[j]] - currentVariables_[2 * currVsGroup[0]])*(currentVariables_[2 * currVsGroup[j]] - currentVariables_[2 * currVsGroup[0]])
						+ (currentVariables_[2 * currVsGroup[j] + 1] - currentVariables_[2 * currVsGroup[0] + 1])*(currentVariables_[2 * currVsGroup[j] + 1] - currentVariables_[2 * currVsGroup[0] + 1]);
				}*/

				for (int k = 0; k < currVsGroup.size(); ++k)
				{
					for (int j = 0; j < currVsGroup.size(); ++j)
					{

						if (j != k)
						{
							sum += (currentVariables_[2 * currVsGroup[j]] - currentVariables_[2 * currVsGroup[k]]) * (currentVariables_[2 * currVsGroup[j]] - currentVariables_[2 * currVsGroup[k]])
								+ (currentVariables_[2 * currVsGroup[j] + 1] - currentVariables_[2 * currVsGroup[k] + 1]) * (currentVariables_[2 * currVsGroup[j] + 1] - currentVariables_[2 * currVsGroup[k] + 1]);
						}
					}
				}
			}
		}
#endif
		for (int i = 0; i < oi_->localVsNotMove.size(); ++i)
		{
			sum += (currentVariables_[2 * oi_->localVsNotMove[i]] - oi_->localVsNotMovePos[i][0]) * (currentVariables_[2 * oi_->localVsNotMove[i]] - oi_->localVsNotMovePos[i][0])
				+ (currentVariables_[2 * oi_->localVsNotMove[i] + 1] - oi_->localVsNotMovePos[i][1]) *(currentVariables_[2 * oi_->localVsNotMove[i] + 1] - oi_->localVsNotMovePos[i][1]);
		}
		return sum;
	}

	void AESolver::Nabla_Hession_E()
	{
		nablaE_.setZero();
		hessionE_.setZero();
		std::vector<Eigen::Triplet<double>> triList;
		uint32_t groupVNum = 0;
		for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			groupVNum += oi_->localVsGroups[i].size();
		}
		triList.reserve(2 * groupVNum + 2 * oi_->localVsNotMove.size() + 2 * oi_->localFc.ids_C.size() + 7 * oi_->localFc.ids_L.size() + 36 * oi_->localTriIds.size());

		//nablaE_
		//E_T
#if !USE_NON_FIXED_TARGET_ENERGY
		for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			for (int j = 0; j < currVsGroup.size(); ++j)
			{
				double x0 = currentVariables_[2 * currVsGroup[j]], y0 = currentVariables_[2 * currVsGroup[j] + 1];
				double x1 = oi_->localVsCoords[i][0], y1 = oi_->localVsCoords[i][1];
				nablaE_[2 * currVsGroup[j]] += LAMBDA_T * 2 * (x0 - x1);
				nablaE_[2 * currVsGroup[j] + 1] += LAMBDA_T * 2 * (y0 - y1);
			}
		}
#else
#pragma region Old_Version_Without_Boundary_Judge
		/*for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			for (int k = 0; k < currVsGroup.size(); ++k)
			{
				double x0 = currentVariables_[2 * currVsGroup[k]], y0 = currentVariables_[2 * currVsGroup[k] + 1];
				for (int j = 0; j < currVsGroup.size(); ++j)
				{
					if (j != k)
					{
						double x1 = currentVariables_[2 * currVsGroup[j]], y1 = currentVariables_[2 * currVsGroup[j] + 1];
						nablaE_[2 * currVsGroup[j]] += LAMBDA_T * 2 * (x1 - x0);
						nablaE_[2 * currVsGroup[j] + 1] += LAMBDA_T * 2 * (y1 - y0);
						nablaE_[2 * currVsGroup[k]] += LAMBDA_T * 2 * (x0 - x1);
						nablaE_[2 * currVsGroup[k] + 1] += LAMBDA_T * 2 * (y0 - y1);
					}
				}
			}
		}*/
#pragma endregion
		for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			if (false/*oi_->vsGroupHasBoundary[i]*/)
			{
				for (int j = 0; j < currVsGroup.size(); ++j)
				{
					double x0 = currentVariables_[2 * currVsGroup[j]], y0 = currentVariables_[2 * currVsGroup[j] + 1];
					double x1 = oi_->localVsCoords[i][0], y1 = oi_->localVsCoords[i][1];
					nablaE_[2 * currVsGroup[j]] += LAMBDA_T * 2 * (x0 - x1);
					nablaE_[2 * currVsGroup[j] + 1] += LAMBDA_T * 2 * (y0 - y1);
				}
			}
			else
			{
				/*double x1 = currentVariables_[2 * currVsGroup[0]], y1 = currentVariables_[2 * currVsGroup[0] + 1];
				for (int j = 1; j < currVsGroup.size(); ++j)
				{
					double x0 = currentVariables_[2 * currVsGroup[j]], y0 = currentVariables_[2 * currVsGroup[j] + 1];
					nablaE_[2 * currVsGroup[0]] += LAMBDA_T * 2 * (x1 - x0);
					nablaE_[2 * currVsGroup[0] + 1] += LAMBDA_T * 2 * (y1 - y0);
					nablaE_[2 * currVsGroup[j]] += LAMBDA_T * 2 * (x0 - x1);
					nablaE_[2 * currVsGroup[j] + 1] += LAMBDA_T * 2 * (y0 - y1);
				}*/
				for (int k = 0; k < currVsGroup.size(); ++k)
				{
					double x0 = currentVariables_[2 * currVsGroup[k]], y0 = currentVariables_[2 * currVsGroup[k] + 1];
					for (int j = 0; j < currVsGroup.size(); ++j)
					{
						if (j != k)
						{
							double x1 = currentVariables_[2 * currVsGroup[j]], y1 = currentVariables_[2 * currVsGroup[j] + 1];
							nablaE_[2 * currVsGroup[j]] += LAMBDA_T * 2 * (x1 - x0);
							nablaE_[2 * currVsGroup[j] + 1] += LAMBDA_T * 2 * (y1 - y0);
							nablaE_[2 * currVsGroup[k]] += LAMBDA_T * 2 * (x0 - x1);
							nablaE_[2 * currVsGroup[k] + 1] += LAMBDA_T * 2 * (y0 - y1);
						}
					}
				}
			}
		}
#endif
		for (int i = 0; i < oi_->localVsNotMove.size(); ++i)
		{
			uint32_t currV = oi_->localVsNotMove[i];
			double x0 = currentVariables_[2 * currV], y0 = currentVariables_[2 * currV + 1];
			double x1 = oi_->localVsNotMovePos[i][0], y1 = oi_->localVsNotMovePos[i][1];
			nablaE_[2 * currV] += LAMBDA_T * 2 * (x0 - x1);
			nablaE_[2 * currV + 1] += LAMBDA_T * 2 * (y0 - y1);
		}

		//E_B
		for (int i = 0; i < oi_->localFc.ids_C.size(); ++i)
		{
			uint32_t currV = oi_->localFc.ids_C[i];
			nablaE_[2 * currV] += 2 * LAMBDA_B * (currentVariables_[2 * currV] - oi_->localFc.C[i][0]);
			nablaE_[2 * currV + 1] += 2 * LAMBDA_B * (currentVariables_[2 * currV + 1] - oi_->localFc.C[i][1]);
		}

		for (int i = 0; i < oi_->localFc.ids_L.size(); ++i)
		{
			uint32_t currV = oi_->localFc.ids_L[i];
			double x0 = currentVariables_[2 * currV], y0 = currentVariables_[2 * currV + 1];
			double x1 = oi_->localFc.origin_L[i][0], y1 = oi_->localFc.origin_L[i][1];
			double x2 = oi_->localFc.axa_L[i][0], y2 = oi_->localFc.axa_L[i][1];
			nablaE_[2 * currV] += 2 * LAMBDA_B * (x0 - x1 - currentVariables_[vVarNum_ + i] * x2);
			nablaE_[2 * currV + 1] += 2 * LAMBDA_B * (y0 - y1 - currentVariables_[vVarNum_ + i] * y2);
			nablaE_[vVarNum_ + i] += LAMBDA_B * 2 * (-x2) * (x0 - x1 - currentVariables_[vVarNum_ + i] * x2);
			nablaE_[vVarNum_ + i] += LAMBDA_B * 2 * (-y2) * (y0 - y1 - currentVariables_[vVarNum_ + i] * y2);
		}

		//hessionE_
		//E_T
#if !USE_NON_FIXED_TARGET_ENERGY
		for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			for (int j = 0; j < currVsGroup.size(); ++j)
			{
				triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j], 2 * currVsGroup[j], 2 * LAMBDA_T));
				triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j]+1, 2 * currVsGroup[j]+1, 2 * LAMBDA_T));
			}
		}
#else
#pragma region Old_Version_Without_Boundary_Judge
		/*for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			for (int k = 0; k < currVsGroup.size(); ++k)
			{
				for (int j = 0; j < currVsGroup.size(); ++j)
				{
					if (j != k)
					{
						triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j], 2 * currVsGroup[j], 2 * LAMBDA_T));
						triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j] + 1, 2 * currVsGroup[j] + 1, 2 * LAMBDA_T));
						triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[k], 2 * currVsGroup[k], 2 * LAMBDA_T));
						triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[k] + 1, 2 * currVsGroup[k] + 1, 2 * LAMBDA_T));
					}
				}
			}
		}*/
#pragma endregion
		for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			if (false/*oi_->vsGroupHasBoundary[i]*/)
			{
				for (int j = 0; j < currVsGroup.size(); ++j)
				{
					triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j], 2 * currVsGroup[j], 2 * LAMBDA_T));
					triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j] + 1, 2 * currVsGroup[j] + 1, 2 * LAMBDA_T));
				}
			}
			else
			{
				/*for (int j = 1; j < currVsGroup.size(); ++j)
				{
					triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j], 2 * currVsGroup[j], 2 * LAMBDA_T));
					triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j] + 1, 2 * currVsGroup[j] + 1, 2 * LAMBDA_T));
					triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[0], 2 * currVsGroup[0], 2 * LAMBDA_T));
					triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[0] + 1, 2 * currVsGroup[0] + 1, 2 * LAMBDA_T));
				}*/
				for (int k = 0; k < currVsGroup.size(); ++k)
				{
					for (int j = 0; j < currVsGroup.size(); ++j)
					{
						if (j != k)
						{
							triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j], 2 * currVsGroup[j], 2 * LAMBDA_T));
							triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j] + 1, 2 * currVsGroup[j] + 1, 2 * LAMBDA_T));
							triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[k], 2 * currVsGroup[k], 2 * LAMBDA_T));
							triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[k] + 1, 2 * currVsGroup[k] + 1, 2 * LAMBDA_T));
						}
					}
				}
			}
		}
#endif
		for (int i = 0; i < oi_->localVsNotMove.size(); ++i)
		{
			uint32_t currV = oi_->localVsNotMove[i];
			triList.emplace_back(Eigen::Triplet<double>(2 * currV, 2 * currV, 2 * LAMBDA_T));
			triList.emplace_back(Eigen::Triplet<double>(2 * currV + 1, 2 * currV + 1, 2 * LAMBDA_T));
		}

		//E_B
		for (int i = 0; i < oi_->localFc.ids_C.size(); ++i)
		{
			uint32_t currV = oi_->localFc.ids_C[i];
			triList.emplace_back(Eigen::Triplet<double>(2 * currV, 2 * currV, 2 * LAMBDA_B));
			triList.emplace_back(Eigen::Triplet<double>(2 * currV + 1, 2 * currV + 1, 2 * LAMBDA_B));
		}

		for (int i = 0; i < oi_->localFc.ids_L.size(); ++i)
		{
			uint32_t currV = oi_->localFc.ids_L[i];
			double x0 = currentVariables_[2 * currV], y0 = currentVariables_[2 * currV + 1];
			double x1 = oi_->localFc.origin_L[i][0], y1 = oi_->localFc.origin_L[i][1];
			double x2 = oi_->localFc.axa_L[i][0], y2 = oi_->localFc.axa_L[i][1];

			triList.emplace_back(Eigen::Triplet<double>(2*currV, 2*currV, 2*LAMBDA_B));
			triList.emplace_back(Eigen::Triplet<double>(2 * currV + 1, 2 * currV + 1, 2 * LAMBDA_B));
			triList.emplace_back(Eigen::Triplet<double>(2 * currV, vVarNum_ + i, -2 * LAMBDA_B * x2));
			triList.emplace_back(Eigen::Triplet<double>(2 * currV + 1, vVarNum_ + i, -2 * LAMBDA_B * y2));
			triList.emplace_back(Eigen::Triplet<double>(vVarNum_ + i, 2 * currV, -2 * LAMBDA_B * x2));
			triList.emplace_back(Eigen::Triplet<double>(vVarNum_ + i, 2 * currV + 1, -2 * LAMBDA_B * y2));
			triList.emplace_back(Eigen::Triplet<double>(vVarNum_ + i, vVarNum_ + i, 2 * LAMBDA_B*(x2*x2 + y2 * y2) ));
		}

		//E_D
		for (int i = 0; i < oi_->localTriIds.size(); ++i)
		{
			std::fill(sdfq, sdfq + 4, 0);
			double *standInverse = standardInverseMatrix_[i];
			double currArea = triAreas_[i];

			int v0Index = oi_->localTriIds[i][0];
			int v1Index = oi_->localTriIds[i][1];
			int v2Index = oi_->localTriIds[i][2];

			ComputeJacobi(i, ptrJ);
			/*currJacobi << ptrJ[0], ptrJ[2],
				ptrJ[1], ptrJ[3];

			Eigen::JacobiSVD<Eigen::MatrixXd> svd(currJacobi, Eigen::ComputeFullU | Eigen::ComputeFullV);
			const Eigen::Matrix2d &currentU = svd.matrixU();
			const Eigen::Matrix2d &currentV = svd.matrixV();*/

			/*double s0 = svd.singularValues()(0, 0);
			double s1 = svd.singularValues()(1, 0);*/
			double s0, s1;

			AssistFunctions::AssistFunc::ComputeSVDDecomposition(ptrJ, currU, s0, s1, currV);
			/*if (s0==0 || s1==0)
			{
				AssistFunctions::AssistFunc::ComputeSVDDecomposition(ptrJ, currU, s0, s1, currV);
				std::cout.precision(18);
				std::cout << s00 << " " << s0 << std::endl;
				std::cout << s11 << " " << s1 << std::endl;
				std::cout << (s00 - s0) / s00 << std::endl;
				std::cout << (s11 - s1) / s11 << std::endl;
				exit(0);
			}*/
			/*if ((s00 - s0)/s00 > 0.1 || (s11 - s1)/s11 > 0.1)
			{
				AssistFunctions::AssistFunc::ComputeSVDDecomposition(ptrJ, currU, s0, s1, currV);
				std::cout.precision(18);
				std::cout << s00 << " " << s0 << std::endl;
				std::cout << s11 << " " << s1 << std::endl;
				std::cout << (s00 - s0) / s00 << std::endl;
				std::cout << (s11 - s1) / s11 << std::endl;
				exit(0);
			}*/
			/*Eigen::Matrix2d currentU;
			Eigen::Matrix2d currentV;
			currentU << -currU[0], currU[2],
				-currU[1], currU[3];
			currentV << -currV[0], currV[2],
				-currV[1], currV[3];*/

			for (int j = 0; j < 2; ++j)
			{
				for (int k = 0; k < 2; ++k)
				{
					//currentF2[j * 2 + k] = 2 * currJacobi(k, j);
					currentF2[j * 2 + k] = 2 * ptrJ[j*2+k];
					//currentG[j * 2 + k] = s1 * currentU(k, 0)*currentV(j, 0) + s0 * currentU(k, 1) * currentV(j, 1);
					currentG[j * 2 + k] = s1 * currU[k]*currV[j] + s0 * currU[k+2] * currV[j+2];
				}
			}

			double i1 = s0 + s1;
			double i2 = s0 * s0 + s1 * s1;
			double i3 = s0 * s1;
			if (std::abs(i3) < DEGEN_TRI_THRESHOLD)
			{
				//std::cout << "Error!" << std::endl;
				isInverse = true;
				return;
			}

			double sd_i2 = SD_I2(i1, i2, i3);
			double sd_i3 = SD_I3(i1, i2, i3);
			double *i2_f = currentF2;
			double *i3_f = currentG;
			for (int j = 0; j < 4; ++j)
			{
				sdfq[j] += sd_i2 * i2_f[j];
				sdfq[j] += sd_i3 * i3_f[j];
			}

			for (int j = 0; j < 2; ++j)
			{
				nablaE_[2 * v0Index + j] += LAMBDA_D * currArea * (standInverse[1] * sdfq[j] + standInverse[3] * sdfq[j + 2]);
				nablaE_[2 * v1Index + j] += LAMBDA_D * currArea * ((-standInverse[0] - standInverse[1])*sdfq[j] + (-standInverse[2] - standInverse[3])*sdfq[j + 2]);
				nablaE_[2 * v2Index + j] += LAMBDA_D * currArea * (standInverse[0] * sdfq[j] + standInverse[2] * sdfq[j + 2]);
			}

			//hession
			double coff = 1.41421 / 2;
			for (int j = 0; j < 2; ++j)
			{
				for (int k = 0; k < 2; ++k)
				{
					/*d0[j * 2 + k] = currentU(k, 0) * currentV(j, 0);
					d1[j * 2 + k] = currentU(k, 1) * currentV(j, 1);
					l0[j * 2 + k] = coff * (currentU(k, 1) * currentV(j, 0) + currentU(k, 0)*currentV(j, 1));
					t0[j * 2 + k] = coff * (currentU(k, 1) * currentV(j, 0) - currentU(k, 0)*currentV(j, 1));*/
					d0[j * 2 + k] = currU[k] * currV[j];
					d1[j * 2 + k] = currU[k+2] * currV[j+2];
					l0[j * 2 + k] = coff * (currU[k+2] * currV[j] + currU[k]*currV[j+2]);
					t0[j * 2 + k] = coff * (currU[k+2] * currV[j] - currU[k]*currV[j+2]);
				}
			}
			if (std::abs(s0) < std::sqrt(DEGEN_TRI_THRESHOLD) || std::abs(s1) < std::sqrt(DEGEN_TRI_THRESHOLD) || std::abs(i3) < DEGEN_TRI_THRESHOLD)
			{
				std::cout << "AESolver::Nabla_Hession_E " << std::endl;
				isInverse = true;
				return;
			}

			lambdaVec[0] = 1.0 + 3 / std::pow(s0, 4);
			lambdaVec[1] = 1.0 + 3 / std::pow(s1, 4);
			lambdaVec[2] = 1.0 + 1.0 / (i3*i3) + i2 / (i3*i3*i3);
			lambdaVec[3] = 1.0 + 1.0 / (i3*i3) - i2 / (i3*i3*i3);

			for (int j = 0; j < 4; ++j)
			{
				if (lambdaVec[j] < 0.1)
					lambdaVec[j] = 0.1;
			}

			for (int j = 0; j < 4; ++j)
			{
				for (int k = 0; k <= j; ++k)
				{
					sd2_fq_fq[j][k] = 2 * (lambdaVec[0] * d0[j] * d0[k] + lambdaVec[1] * d1[j] * d1[k] + lambdaVec[2] * l0[j] * l0[k] + lambdaVec[3] * t0[j] * t0[k]);
				}
			}
			for (int j = 0; j < 4; ++j)
			{
				for (int k = j + 1; k < 4; ++k)
				{
					sd2_fq_fq[j][k] = sd2_fq_fq[k][j];
				}
			}

			//std::fill(leftMat[0], leftMat[0] + 6*4, 0);
			mapV[0] = 2 * v0Index; mapV[1] = 2 * v0Index + 1;
			mapV[2] = 2 * v1Index; mapV[3] = 2 * v1Index + 1;
			mapV[4] = 2 * v2Index; mapV[5] = 2 * v2Index + 1;

			leftMat[0][0] = standInverse[1]; leftMat[0][2] = standInverse[3];
			leftMat[1][1] = leftMat[0][0]; leftMat[1][3] = leftMat[0][2];
			leftMat[2][0] = -standInverse[0] - standInverse[1];
			leftMat[2][2] = -standInverse[2] - standInverse[3];
			leftMat[3][1] = leftMat[2][0]; leftMat[3][3] = leftMat[2][2];
			leftMat[4][0] = standInverse[0]; leftMat[4][2] = standInverse[2];
			leftMat[5][1] = leftMat[4][0]; leftMat[5][3] = leftMat[4][2];

			for (int j = 0; j < 6; ++j)
			{
				for (int k = 0; k < 6; ++k)
				{
					double tempValue = LAMBDA_D * (currArea * (leftMat[j][j % 2] * sd2_fq_fq[j % 2][k % 2] + leftMat[j][j % 2 + 2] * sd2_fq_fq[j % 2 + 2][k % 2] ) * leftMat[k][k % 2]
						+ currArea * (leftMat[j][j % 2] * sd2_fq_fq[j % 2][k % 2 + 2] + leftMat[j][j % 2 + 2] * sd2_fq_fq[j % 2 + 2][k % 2 + 2]) * leftMat[k][k % 2+ 2]);

					triList.emplace_back(Eigen::Triplet<double>(mapV[j], mapV[k], tempValue));
				}
			}
		}

#if USE_FIELD_ENERGY
		efe_.ExtraNablaHessionE(nablaE_, triList, currentVariables_, hessionE_);
#else
		hessionE_.setFromTriplets(triList.begin(), triList.end());
		hessionE_.makeCompressed();
#endif
	}

	double AESolver::SD_I2(double i1, double i2, double i3)
	{
		return 1.0 + 1.0 / (i3*i3);
	}

	double AESolver::SD_I3(double i1, double i2, double i3)
	{
		return -2.0 * i2 / (i3*i3*i3);
	}

	void AESolver::ComputeJacobi(int index, double *ptrJ)
	{
		double *currStanMat = standardInverseMatrix_[index];
		const Eigen::Vector3i &currVs = oi_->localTriIds[index];

		double x0 = currentVariables_[2 * currVs[0]];
		double y0 = currentVariables_[2 * currVs[0] + 1];
		double x1 = currentVariables_[2 * currVs[1]];
		double y1 = currentVariables_[2 * currVs[1] + 1];
		double x2 = currentVariables_[2 * currVs[2]];
		double y2 = currentVariables_[2 * currVs[2] + 1];

		tempMat[0] = x2 - x1; tempMat[2] = x0 - x1;
		tempMat[1] = y2 - y1; tempMat[3] = y0 - y1;

		ptrJ[0] = tempMat[0] * currStanMat[0] + tempMat[2] * currStanMat[1];
		ptrJ[2] = tempMat[0] * currStanMat[2] + tempMat[2] * currStanMat[3];
		ptrJ[1] = tempMat[1] * currStanMat[0] + tempMat[3] * currStanMat[1];
		ptrJ[3] = tempMat[1] * currStanMat[2] + tempMat[3] * currStanMat[3];
	}

	void AESolver::BuildStandardInverseMatrix()
	{
		for (int i = 0; i < oi_->localStanTriPos.size(); ++i)
		{
			std::vector<Eigen::Vector2d> &currVPos = oi_->localStanTriPos[i];
			const Eigen::Vector2d &stanVec0 = currVPos[2] - currVPos[1];
			const Eigen::Vector2d &stanVec1 = currVPos[0] - currVPos[1];

			double detM = stanVec0[0] * stanVec1[1] - stanVec0[1] * stanVec1[0];
			if (std::abs(detM) < DEGEN_TRI_THRESHOLD)
			{
				std::cout << "degenerate quad! " << std::endl;
				isInverse = true;
				standardInverseMatrix_[i][0] = 1;
				standardInverseMatrix_[i][2] = 0;
				standardInverseMatrix_[i][1] = 0;
				standardInverseMatrix_[i][3] = 1;
				return;
			}
			standardInverseMatrix_[i][0] = stanVec1[1] / detM;
			standardInverseMatrix_[i][2] = -stanVec1[0] / detM;
			standardInverseMatrix_[i][1] = -stanVec0[1] / detM;
			standardInverseMatrix_[i][3] = stanVec0[0] / detM;
		}
	}

	void AESolver::BuildTriAreas()
	{
		for (int i = 0; i < oi_->localTriIds.size(); ++i)
		{
			Eigen::Vector3i &currVIds = oi_->localTriIds[i];
			const Eigen::Vector2d &vPos0 = oi_->localVPos[currVIds[0]];
			const Eigen::Vector2d &vPos1 = oi_->localVPos[currVIds[1]];
			const Eigen::Vector2d &vPos2 = oi_->localVPos[currVIds[2]];
			const Eigen::Vector2d &vec0 = vPos1 - vPos0;
			const Eigen::Vector2d &vec1 = vPos2 - vPos0;
			triAreas_[i] = std::abs(vec0[0] * vec1[1] - vec0[1] * vec1[0])/2;
		}
	}

	void AESolver::MatrixInverse(double *input, double *inverse)
	{
		double detM = input[0]*input[3] - input[1]*input[2];
		if (std::abs(detM) < DEGEN_TRI_THRESHOLD)
		{
			std::cout << "Singular matrix! " << std::endl;
			isInverse = true;
			inverse[0] = 1;
			inverse[1] = 0;
			inverse[2] = 0;
			inverse[3] = 1;
			return;
		}
		inverse[0] = input[3] / detM;
		inverse[1] = -input[1] / detM;
		inverse[2] = -input[2] / detM;
		inverse[3] = input[0] / detM;
	}

	double AESolver::ComputeSyDir(int index)
	{
		double *currStanMat = standardInverseMatrix_[index];
		const Eigen::Vector3i &currVs = oi_->localTriIds[index];

		double x0 = currentVariables_[2 * currVs[0]];
		double y0 = currentVariables_[2 * currVs[0] + 1];
		double x1 = currentVariables_[2 * currVs[1]];
		double y1 = currentVariables_[2 * currVs[1] + 1];
		double x2 = currentVariables_[2 * currVs[2]];
		double y2 = currentVariables_[2 * currVs[2] + 1];

		tempMat[0] = x2 - x1; tempMat[2] = x0 - x1;
		tempMat[1] = y2 - y1; tempMat[3] = y0 - y1;

		tempMat2[0] = tempMat[0] * currStanMat[0] + tempMat[2] * currStanMat[1];
		tempMat2[2] = tempMat[0] * currStanMat[2] + tempMat[2] * currStanMat[3];
		tempMat2[1] = tempMat[1] * currStanMat[0] + tempMat[3] * currStanMat[1];
		tempMat2[3] = tempMat[1] * currStanMat[2] + tempMat[3] * currStanMat[3];

		MatrixInverse(tempMat2, tempMat3);
		
		double sum = 0;
		for (int i = 0; i < 4; ++i)
		{
			sum += tempMat2[i] * tempMat2[i] + tempMat3[i] * tempMat3[i];
		}
		return sum;
	}

	void AESolver::Optimize()
	{
		int outLoopCount = 0;
		double preE = -1000, currE = 0;
		bool firstIteFlag = true;
		int sameEValueCount = 0;

		Nabla_Hession_E();
		if (isInverse)
			return;
		while (nablaE_.norm() > 1E-2 && outLoopCount < OUT_LOOP_COUNT)
		{
			++outLoopCount;

			Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > solver;
			solver.compute(hessionE_);
			Eigen::VectorXd pk = solver.solve(-nablaE_);
			/*Eigen::VectorXd guass(eleNum_);
			guass.setZero();
			Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
			solver.setTolerance(1e-8);
			solver.compute(hessionE_);
			Eigen::VectorXd pk = solver.solveWithGuess(-nablaE_, guass);*/

			double ek = E();
			currE = ek;
			if (std::abs(currE - preE) < 1E-2)
				++sameEValueCount;
			else
				sameEValueCount = 0;

			if (sameEValueCount == 2)
				break;

			double nablaENorm = nablaE_.norm();
			std::cout << "In the " << outLoopCount << "th iterate, E = " << ek << ",	Nabla_E = " << nablaENorm << std::endl;
			std::cout << " E_T = " << et_ << ",	E_B = " << eb_ << ", E_D = " << ed_ << ", E_F = " << ef_ << std::endl;
			double nablaDotPk = nablaE_.dot(pk);
			std::cout << "The length of pk is " << pk.norm() << ", nablaDotPk is " << nablaDotPk << std::endl;
			if (nablaDotPk > 0)
			{
				std::cout << "It is not a falling direction. " << std::endl;
				//return;
			}

			Iterate(pk, nablaDotPk, ek);
			preE = currE;
			if (outLoopCount != OUT_LOOP_COUNT)
			{
				Nabla_Hession_E();
				if (isInverse)
					return;
			}
		}
	}

	double AESolver::ComputeNoInverseMinAlpha(Eigen::VectorXd &pk)
	{
		std::vector<double> allAlpha;
		allAlpha.reserve(oi_->localTriIds.size());
		for (int i = 0; i < oi_->localTriIds.size(); ++i)
		{
			double tempAlpha = 0;
			int v0Id = oi_->localTriIds[i][0];
			int v1Id = oi_->localTriIds[i][1];
			int v2Id = oi_->localTriIds[i][2];

			double x10 = currentVariables_[2 * v1Id] - currentVariables_[2 * v0Id];
			double y10 = currentVariables_[2 * v1Id + 1] - currentVariables_[2 * v0Id + 1];
			double x20 = currentVariables_[2 * v2Id] - currentVariables_[2 * v0Id];
			double y20 = currentVariables_[2 * v2Id + 1] - currentVariables_[2 * v0Id + 1];
			double px10 = pk[2 * v1Id] - pk[2 * v0Id];
			double py10 = pk[2 * v1Id + 1] - pk[2 * v0Id + 1];
			double px20 = pk[2 * v2Id] - pk[2 * v0Id];
			double py20 = pk[2 * v2Id + 1] - pk[2 * v0Id + 1];

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
		return (allAlpha[0] > 0) ? allAlpha[0]*0.8 : 1;
	}

	void AESolver::ComputeNoInverseAlpha(Eigen::VectorXd &pk, double *alphaPtr)
	{
		std::fill(alphaPtr, alphaPtr + pk.size(), 1);
		for (int i = 0; i < oi_->localTriIds.size(); ++i)
		{
			double tempAlpha = 0;
			int v0Id = oi_->localTriIds[i][0];
			int v1Id = oi_->localTriIds[i][1];
			int v2Id = oi_->localTriIds[i][2];

			double x10 = currentVariables_[2 * v1Id] - currentVariables_[2 * v0Id];
			double y10 = currentVariables_[2 * v1Id + 1] - currentVariables_[2 * v0Id + 1];
			double x20 = currentVariables_[2 * v2Id] - currentVariables_[2 * v0Id];
			double y20 = currentVariables_[2 * v2Id + 1] - currentVariables_[2 * v0Id + 1];
			double px10 = pk[2 * v1Id] - pk[2 * v0Id];
			double py10 = pk[2 * v1Id + 1] - pk[2 * v0Id + 1];
			double px20 = pk[2 * v2Id] - pk[2 * v0Id];
			double py20 = pk[2 * v2Id + 1] - pk[2 * v0Id + 1];

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

			if (tempAlpha < alphaPtr[v0Id * 2])
			{
				alphaPtr[v0Id * 2] = tempAlpha;
				alphaPtr[v0Id * 2 + 1] = tempAlpha;
			}
			if (tempAlpha < alphaPtr[v1Id * 2])
			{
				alphaPtr[v1Id * 2] = tempAlpha;
				alphaPtr[v1Id * 2 + 1] = tempAlpha;
			}
			if (tempAlpha < alphaPtr[v2Id * 2])
			{
				alphaPtr[v2Id * 2] = tempAlpha;
				alphaPtr[v2Id * 2 + 1] = tempAlpha;
			}
		}
	}

	void AESolver::ComputeNewCurrVariables(Eigen::VectorXd &backupVars, double *alphaPtr, Eigen::VectorXd &pk, Eigen::VectorXd &outVars)
	{
		for (int i = 0; i < backupVars.size(); ++i)
		{
			outVars[i] = backupVars[i] + alphaPtr[i] * pk[i];
		}
	}

	void AESolver::Iterate(Eigen::VectorXd &pk, double nablaDotPk, double ek) 
	{
		backupVars_ = currentVariables_;
#if USE_ALPHA_VEC
		double alpha = 1;
		ComputeNoInverseAlpha(pk, alphaPtr_);
#else
		double alpha = ComputeNoInverseMinAlpha(pk);
#endif
		int loopCount = 0;

#if USE_ALPHA_VEC
		ComputeNewCurrVariables(backupVars_, alphaPtr_, pk, currentVariables_);
#else
		currentVariables_ = backupVars_ + alpha * pk;
#endif
		double currE = E();

		while (currE > ek + 0.2*alpha*nablaDotPk && loopCount < INNER_LOOP_COUNT)
		{
			alpha *= 0.8;
#if USE_ALPHA_VEC
			for (int j = 0; j < pk.size(); ++j)
			{
				alphaPtr_[j] *= 0.8;
			}
			ComputeNewCurrVariables(backupVars_, alphaPtr_, pk, currentVariables_);
#else
			currentVariables_ = backupVars_ + alpha * pk;
#endif
			currE = E();
			++loopCount;
		}

		if (loopCount == INNER_LOOP_COUNT && currE>ek)
		{
			std::cout << "AESolver::Iterate::INNER_LOOP_COUNT " << std::endl;
			currentVariables_ = backupVars_;
		}

	}

	double AESolver::QuadraticRoot1(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20)
	{
		double a = px10 * py20 - px20 * py10;
		double b = px10 * y20 + py20 * x10 - px20 * y10 - py10 * x20;
		double c = x10 * y20 - x20 * y10;

		if (std::abs(a) < 1E-13 || b * b - 4 * a*c < 0)
			return 1;
		else
		return (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
	}
	double AESolver::QuadraticRoot2(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20)
	{
		double a = px10 * py20 - px20 * py10;
		double b = px10 * y20 + py20 * x10 - px20 * y10 - py10 * x20;
		double c = x10 * y20 - x20 * y10;

		if (std::abs(a) < 1E-13 || b * b - 4 * a*c < 0)
			return 1;
		else
		return (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
	}
	double AESolver::LinearRoot(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20)
	{
		double b = px10 * y20 + py20 * x10 - px20 * y10 - py10 * x20;
		double c = x10 * y20 - x20 * y10;
		if (std::abs(b) < 1E-15)
			return 1;
		else
			return -c / b;
	}

	void AESolver::SetValueBack(OptimizeInfo *oi)
	{
		for (int i = 0; i < oi->localVPos.size(); ++i)
		{
			oi->localVPos[i][0] = currentVariables_[2 * i];
			oi->localVPos[i][1] = currentVariables_[2 * i + 1];
		}
	}

	void AESolver::SetValueFront(OptimizeInfo *oi)
	{
		for (int i = 0; i < oi->localVPos.size(); ++i)
		{
			/*oi->localVPos[i][0] = currentVariables_[2 * i];
			oi->localVPos[i][1] = currentVariables_[2 * i + 1];*/
			currentVariables_[2 * i] = oi->localVPos[i][0];
			currentVariables_[2 * i+1] = oi->localVPos[i][1];
		}
	}
}
