#include "AESolverSquare.h"
#include "Pipeline.h"
#include "AssistFunc.h"
#include <Eigen/SVD>

#define USE_FIELD_ENERGY 0
using namespace BaseDataStructure;

namespace DataOperation
{
	const int OUT_LOOP_COUNT = 2;
	const int INNER_LOOP_COUNT = 100;
	double LAMBDA_D_S = 1;
	double LAMBDA_B_S = 1E9;
	double LAMBDA_T_S = 1E5;
	const double DEGEN_TRI_THRESHOLD = 1E-20;
	AESolverSquare::AESolverSquare(OptimizeInfo *oi, MeshFeatures &mf, ANNkd_tree *kdTree, std::vector<Eigen::Vector3d> &fields, bool proPreInv) : oi_(oi), mf_(mf), efe_(oi, kdTree, fields, 10, 1), proPreInv_(proPreInv)
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
	}


	AESolverSquare::~AESolverSquare()
	{
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

	double AESolverSquare::E()
	{
#if USE_FIELD_ENERGY
		double ef = efe_.E(currentVariables_);
		return LAMBDA_D_S * E_D() + LAMBDA_B_S * E_B() + LAMBDA_T_S * E_T() + ef;
#else
		return LAMBDA_D_S * E_D() + LAMBDA_B_S * E_B() + LAMBDA_T_S * E_T();
#endif
	}

	double AESolverSquare::E_B()
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

	double AESolverSquare::E_D()
	{
		double sum = 0;
		for (int i = 0; i < oi_->localTriIds.size(); ++i)
		{
			double currSyDir = ComputeSyDir(i);
			sum += triAreas_[i] * currSyDir * currSyDir;	//Change: 能量函数变为对称迪利克雷的平方。
		}
		return sum;
	}

	double AESolverSquare::E_T()
	{
		double sum = 0;
		for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			for (int j = 0; j < currVsGroup.size(); ++j)
			{
				sum += (currentVariables_[2 * currVsGroup[j]] - oi_->localVsCoords[i][0])*(currentVariables_[2 * currVsGroup[j]] - oi_->localVsCoords[i][0])
					+ (currentVariables_[2 * currVsGroup[j] + 1] - oi_->localVsCoords[i][1])*(currentVariables_[2 * currVsGroup[j] + 1] - oi_->localVsCoords[i][1]);
			}
		}

		for (int i = 0; i < oi_->localVsNotMove.size(); ++i)
		{
			sum += (currentVariables_[2 * oi_->localVsNotMove[i]] - oi_->localVsNotMovePos[i][0]) * (currentVariables_[2 * oi_->localVsNotMove[i]] - oi_->localVsNotMovePos[i][0])
				+ (currentVariables_[2 * oi_->localVsNotMove[i] + 1] - oi_->localVsNotMovePos[i][1]) *(currentVariables_[2 * oi_->localVsNotMove[i] + 1] - oi_->localVsNotMovePos[i][1]);
		}
		return sum;
	}

	void AESolverSquare::Nabla_Hession_E()
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
		for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			for (int j = 0; j < currVsGroup.size(); ++j)
			{
				double x0 = currentVariables_[2 * currVsGroup[j]], y0 = currentVariables_[2 * currVsGroup[j] + 1];
				double x1 = oi_->localVsCoords[i][0], y1 = oi_->localVsCoords[i][1];
				nablaE_[2 * currVsGroup[j]] += LAMBDA_T_S * 2 * (x0 - x1);
				nablaE_[2 * currVsGroup[j] + 1] += LAMBDA_T_S * 2 * (y0 - y1);
			}
		}
		for (int i = 0; i < oi_->localVsNotMove.size(); ++i)
		{
			uint32_t currV = oi_->localVsNotMove[i];
			double x0 = currentVariables_[2 * currV], y0 = currentVariables_[2 * currV + 1];
			double x1 = oi_->localVsNotMovePos[i][0], y1 = oi_->localVsNotMovePos[i][1];
			nablaE_[2 * currV] += LAMBDA_T_S * 2 * (x0 - x1);
			nablaE_[2 * currV + 1] += LAMBDA_T_S * 2 * (y0 - y1);
		}

		//E_B
		for (int i = 0; i < oi_->localFc.ids_C.size(); ++i)
		{
			uint32_t currV = oi_->localFc.ids_C[i];
			nablaE_[2 * currV] += 2 * LAMBDA_B_S * (currentVariables_[2 * currV] - oi_->localFc.C[i][0]);
			nablaE_[2 * currV + 1] += 2 * LAMBDA_B_S * (currentVariables_[2 * currV + 1] - oi_->localFc.C[i][1]);
		}

		for (int i = 0; i < oi_->localFc.ids_L.size(); ++i)
		{
			uint32_t currV = oi_->localFc.ids_L[i];
			double x0 = currentVariables_[2 * currV], y0 = currentVariables_[2 * currV + 1];
			double x1 = oi_->localFc.origin_L[i][0], y1 = oi_->localFc.origin_L[i][1];
			double x2 = oi_->localFc.axa_L[i][0], y2 = oi_->localFc.axa_L[i][1];
			nablaE_[2 * currV] += 2 * LAMBDA_B_S * (x0 - x1 - currentVariables_[vVarNum_ + i] * x2);
			nablaE_[2 * currV + 1] += 2 * LAMBDA_B_S * (y0 - y1 - currentVariables_[vVarNum_ + i] * y2);
			nablaE_[vVarNum_ + i] += LAMBDA_B_S * 2 * (-x2) * (x0 - x1 - currentVariables_[vVarNum_ + i] * x2);
			nablaE_[vVarNum_ + i] += LAMBDA_B_S * 2 * (-y2) * (y0 - y1 - currentVariables_[vVarNum_ + i] * y2);
		}

		//hessionE_
		//E_T
		for (int i = 0; i < oi_->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oi_->localVsGroups[i];
			for (int j = 0; j < currVsGroup.size(); ++j)
			{
				triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j], 2 * currVsGroup[j], 2 * LAMBDA_T_S));
				triList.emplace_back(Eigen::Triplet<double>(2 * currVsGroup[j] + 1, 2 * currVsGroup[j] + 1, 2 * LAMBDA_T_S));
			}
		}
		for (int i = 0; i < oi_->localVsNotMove.size(); ++i)
		{
			uint32_t currV = oi_->localVsNotMove[i];
			triList.emplace_back(Eigen::Triplet<double>(2 * currV, 2 * currV, 2 * LAMBDA_T_S));
			triList.emplace_back(Eigen::Triplet<double>(2 * currV + 1, 2 * currV + 1, 2 * LAMBDA_T_S));
		}

		//E_B
		for (int i = 0; i < oi_->localFc.ids_C.size(); ++i)
		{
			uint32_t currV = oi_->localFc.ids_C[i];
			triList.emplace_back(Eigen::Triplet<double>(2 * currV, 2 * currV, 2 * LAMBDA_B_S));
			triList.emplace_back(Eigen::Triplet<double>(2 * currV + 1, 2 * currV + 1, 2 * LAMBDA_B_S));
		}

		for (int i = 0; i < oi_->localFc.ids_L.size(); ++i)
		{
			uint32_t currV = oi_->localFc.ids_L[i];
			double x0 = currentVariables_[2 * currV], y0 = currentVariables_[2 * currV + 1];
			double x1 = oi_->localFc.origin_L[i][0], y1 = oi_->localFc.origin_L[i][1];
			double x2 = oi_->localFc.axa_L[i][0], y2 = oi_->localFc.axa_L[i][1];

			triList.emplace_back(Eigen::Triplet<double>(2 * currV, 2 * currV, 2 * LAMBDA_B_S));
			triList.emplace_back(Eigen::Triplet<double>(2 * currV + 1, 2 * currV + 1, 2 * LAMBDA_B_S));
			triList.emplace_back(Eigen::Triplet<double>(2 * currV, vVarNum_ + i, -2 * LAMBDA_B_S * x2));
			triList.emplace_back(Eigen::Triplet<double>(2 * currV + 1, vVarNum_ + i, -2 * LAMBDA_B_S * y2));
			triList.emplace_back(Eigen::Triplet<double>(vVarNum_ + i, 2 * currV, -2 * LAMBDA_B_S * x2));
			triList.emplace_back(Eigen::Triplet<double>(vVarNum_ + i, 2 * currV + 1, -2 * LAMBDA_B_S * y2));
			triList.emplace_back(Eigen::Triplet<double>(vVarNum_ + i, vVarNum_ + i, 2 * LAMBDA_B_S*(x2*x2 + y2 * y2)));
		}

		//E_D
		for (int i = 0; i < oi_->localTriIds.size(); ++i)
		{
			std::fill(sdfq, sdfq + 4, 0);
			double currSyDir = ComputeSyDir(i);
			double *standInverse = standardInverseMatrix_[i];
			double currArea = triAreas_[i];

			int v0Index = oi_->localTriIds[i][0];
			int v1Index = oi_->localTriIds[i][1];
			int v2Index = oi_->localTriIds[i][2];

			ComputeJacobi(i, ptrJ);
			currJacobi << ptrJ[0], ptrJ[2],
				ptrJ[1], ptrJ[3];

			Eigen::JacobiSVD<Eigen::MatrixXd> svd(currJacobi, Eigen::ComputeFullU | Eigen::ComputeFullV);
			const Eigen::Matrix2d &currentU = svd.matrixU();
			const Eigen::Matrix2d &currentV = svd.matrixV();

			double s0 = svd.singularValues()(0, 0);
			double s1 = svd.singularValues()(1, 0);

			for (int j = 0; j < 2; ++j)
			{
				for (int k = 0; k < 2; ++k)
				{
					currentF2[j * 2 + k] = 2 * currJacobi(k, j);
					currentG[j * 2 + k] = s1 * currentU(k, 0)*currentV(j, 0) + s0 * currentU(k, 1) * currentV(j, 1);
				}
			}

			double i1 = s0 + s1;
			double i2 = s0 * s0 + s1 * s1;
			double i3 = s0 * s1;
			if (std::abs(i3) < std::sqrt(DEGEN_TRI_THRESHOLD))
			{
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
				currSourceFDX[j] = (standInverse[1] * sdfq[j] + standInverse[3] * sdfq[j + 2]);
				nablaE_[2 * v0Index + j] += LAMBDA_D_S * currArea * currSourceFDX[j] * 2 * currSyDir;
				currSourceFDX[2+j] = ((-standInverse[0] - standInverse[1])*sdfq[j] + (-standInverse[2] - standInverse[3])*sdfq[j + 2]);
				nablaE_[2 * v1Index + j] += LAMBDA_D_S * currArea * currSourceFDX[2+j] * 2 * currSyDir;
				currSourceFDX[4+j] = (standInverse[0] * sdfq[j] + standInverse[2] * sdfq[j + 2]);
				nablaE_[2 * v2Index + j] += LAMBDA_D_S * currArea * currSourceFDX[4+j] * 2 * currSyDir;
			}

			//hession
			double coff = 1.41421 / 2;
			for (int j = 0; j < 2; ++j)
			{
				for (int k = 0; k < 2; ++k)
				{
					d0[j * 2 + k] = currentU(k, 0) * currentV(j, 0);
					d1[j * 2 + k] = currentU(k, 1) * currentV(j, 1);
					l0[j * 2 + k] = coff * (currentU(k, 1) * currentV(j, 0) + currentU(k, 0)*currentV(j, 1));
					t0[j * 2 + k] = coff * (currentU(k, 1) * currentV(j, 0) - currentU(k, 0)*currentV(j, 1));
				}
			}
			if (std::abs(s0) < std::sqrt(DEGEN_TRI_THRESHOLD) || std::abs(s1) < std::sqrt(DEGEN_TRI_THRESHOLD) || std::abs(i3) < std::sqrt(DEGEN_TRI_THRESHOLD))
			{
				std::cout << "AESolverSquare::Nabla_Hession_E " << std::endl;
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

			//std::fill(leftMat[0], leftMat[0] + 6 * 4, 0);
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
					double tempValue = LAMBDA_D_S * (currArea * (leftMat[j][j % 2] * sd2_fq_fq[j % 2][k % 2] + leftMat[j][j % 2 + 2] * sd2_fq_fq[j % 2 + 2][k % 2]) * leftMat[k][k % 2]
						+ currArea * (leftMat[j][j % 2] * sd2_fq_fq[j % 2][k % 2 + 2] + leftMat[j][j % 2 + 2] * sd2_fq_fq[j % 2 + 2][k % 2 + 2]) * leftMat[k][k % 2 + 2]);

					double finalTempValue = 2 * currSourceFDX[j] * currSourceFDX[k] * LAMBDA_D_S * currArea + 2 * currSyDir * tempValue;
					triList.emplace_back(Eigen::Triplet<double>(mapV[j], mapV[k], finalTempValue));
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

	double AESolverSquare::SD_I2(double i1, double i2, double i3)
	{
		return 1.0 + 1.0 / (i3*i3);
	}

	double AESolverSquare::SD_I3(double i1, double i2, double i3)
	{
		return -2.0 * i2 / (i3*i3*i3);
	}

	void AESolverSquare::ComputeJacobi(int index, double *ptrJ)
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

	void AESolverSquare::BuildStandardInverseMatrix()
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

	void AESolverSquare::BuildTriAreas()
	{
		for (int i = 0; i < oi_->localTriIds.size(); ++i)
		{
			Eigen::Vector3i &currVIds = oi_->localTriIds[i];
			const Eigen::Vector2d &vPos0 = oi_->localVPos[currVIds[0]];
			const Eigen::Vector2d &vPos1 = oi_->localVPos[currVIds[1]];
			const Eigen::Vector2d &vPos2 = oi_->localVPos[currVIds[2]];
			const Eigen::Vector2d &vec0 = vPos1 - vPos0;
			const Eigen::Vector2d &vec1 = vPos2 - vPos0;
			triAreas_[i] = std::abs(vec0[0] * vec1[1] - vec0[1] * vec1[0]) / 2;
		}
	}

	void AESolverSquare::MatrixInverse(double *input, double *inverse)
	{
		double detM = input[0] * input[3] - input[1] * input[2];
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

	double AESolverSquare::ComputeSyDir(int index)
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

	void AESolverSquare::Optimize()
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
			double nablaDotPk = nablaE_.dot(pk);
			std::cout << "The length of pk is " << pk.norm() << ", nablaDotPk is " << nablaDotPk << std::endl;
			if (nablaDotPk > 0)
				std::cout << "It is not a falling direction. " << std::endl;

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

	double AESolverSquare::ComputeNoInverseMinAlphaForWholeQM(Eigen::VectorXd &pk)
	{
		std::vector<double> allAlpha;
		allAlpha.reserve(qm_->Fs_.size() * 4);

		double x0, x1, x2, y0, y1, y2, px0, px1, px2, py0, py1, py2;
		for (int i = 0; i < qm_->Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm_->Fs_[i].vs;
			for (int j = 0; j < 4; ++j)
			{
				double tempAlpha = 0;
				int gV0Id = fvs[quadTriTable[j][0]];
				int gV1Id = fvs[quadTriTable[j][1]];
				int gV2Id = fvs[quadTriTable[j][2]];
				uint32_t v0Id = oi_->globalToLocalVVec[gV0Id];
				uint32_t v1Id = oi_->globalToLocalVVec[gV1Id];
				uint32_t v2Id = oi_->globalToLocalVVec[gV2Id];
				if (v0Id == (uint32_t)-1 && v1Id == (uint32_t)-1 && v2Id == (uint32_t)-1)
				{
					allAlpha.emplace_back(1);
					continue;
				}

				if (v0Id == (uint32_t)-1)
				{
					x0 = qm_->V_[gV0Id][0];
					y0 = qm_->V_[gV0Id][1];
					px0 = 0; py0 = 0;
				}
				else
				{
					x0 = currentVariables_[2 * v0Id];
					y0 = currentVariables_[2 * v0Id + 1];
					px0 = pk[2 * v0Id];
					py0 = pk[2 * v0Id + 1];
				}
				if (v1Id == (uint32_t)-1)
				{
					x1 = qm_->V_[gV1Id][0];
					y1 = qm_->V_[gV1Id][1];
					px1 = 0; py1 = 0;
				}
				else
				{
					x1 = currentVariables_[2 * v1Id];
					y1 = currentVariables_[2 * v1Id + 1];
					px1 = pk[2 * v1Id];
					py1 = pk[2 * v1Id + 1];
				}
				if (v2Id == (uint32_t)-1)
				{
					x2 = qm_->V_[gV2Id][0];
					y2 = qm_->V_[gV2Id][1];
					px2 = 0; py2 = 0;
				}
				else
				{
					x2 = currentVariables_[2 * v2Id];
					y2 = currentVariables_[2 * v2Id + 1];
					px2 = pk[2 * v2Id];
					py2 = pk[2 * v2Id + 1];
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
		}
		std::sort(allAlpha.begin(), allAlpha.end());
		return (allAlpha[0] > 0) ? allAlpha[0] * 0.8 : 1;
	}

	double AESolverSquare::ComputeNoInverseMinAlpha(Eigen::VectorXd &pk)
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

	void AESolverSquare::ComputeNoInverseAlpha(Eigen::VectorXd &pk, double *alphaPtr)
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
			
			if (tempAlpha < alphaPtr[v0Id])
				alphaPtr[v0Id] = tempAlpha;
			if (tempAlpha < alphaPtr[v1Id])
				alphaPtr[v1Id] = tempAlpha;
			if (tempAlpha < alphaPtr[v2Id])
				alphaPtr[v2Id] = tempAlpha;
		}
	}

	bool AESolverSquare::JudgeJacobiInOptimization(MeshFeatures &mf)
	{
		Eigen::VectorXd copyVariable_ = currentVariables_;
		
		ProjectInOptimization(oi_, mf, copyVariable_);

		//judge if inverse 
		double minJacobi = 10000;
		Eigen::Vector2d vPos0(copyVariable_[boundaryTris_[0][0] * 2], copyVariable_[boundaryTris_[0][0] * 2 + 1]),
			vPos1(copyVariable_[boundaryTris_[0][1] * 2], copyVariable_[boundaryTris_[0][1] * 2 + 1]),
			vPos2(copyVariable_[boundaryTris_[0][2] * 2], copyVariable_[boundaryTris_[0][2] * 2 + 1]);
		Eigen::Vector2d vec0 = (vPos2 - vPos1).normalized();
		Eigen::Vector2d vec1 = (vPos0 - vPos1).normalized();
		double startSign = vec0[0] * vec1[1] - vec0[1] * vec1[0];
		for (int i = 0; i < boundaryTris_.size(); ++i)
		{
			vPos0[0] = copyVariable_[boundaryTris_[i][0] * 2]; vPos0[1] = copyVariable_[boundaryTris_[i][0] * 2 + 1];
			vPos1[0] = copyVariable_[boundaryTris_[i][1] * 2]; vPos1[1] = copyVariable_[boundaryTris_[i][1] * 2 + 1];
			vPos2[0] = copyVariable_[boundaryTris_[i][2] * 2]; vPos2[1] = copyVariable_[boundaryTris_[i][2] * 2 + 1];
			vec0 = (vPos2 - vPos1).normalized();
			vec1 = (vPos0 - vPos1).normalized();
			double currSign = vec0[0] * vec1[1] - vec0[1] * vec1[0];
			/*if (std::abs(currSign) < minJacobi_ && JudgeSign(startSign, currSign))
				minJacobi_ = std::abs(currSign);*/
			if (!BaseDataStructure::QuadMesh::JudgeSign(startSign, currSign))
			{
				return false;
			}
			else if (std::abs(currSign) < minJacobi)
				minJacobi = std::abs(currSign);
		}
		if (minJacobi < 0.1)
		{
			return false;
		}
		return true;
	}

	void AESolverSquare::ProjectInOptimization(OptimizeInfo *oi, MeshFeatures &mf, Eigen::VectorXd &copyVariables)
	{
		std::vector<uint32_t> &currIdsC = oi->localFc.ids_C;
		for (int i = 0; i < currIdsC.size(); ++i)
		{
			//oi->localVPos[currIdsC[i]] = oi->localFc.C[i];
			copyVariables[currIdsC[i] * 2] = oi->localFc.C[i][0];
			copyVariables[currIdsC[i] * 2 + 1] = oi->localFc.C[i][1];
		}
		std::vector<uint32_t> &currIdsL = oi->localFc.ids_L;
		for (int i = 0; i < currIdsL.size(); ++i)
		{
			const Eigen::Vector2d currP(copyVariables[currIdsL[i] * 2], copyVariables[currIdsL[i] * 2+1]);
			Point segP(currP[0], currP[1], 0);
			uint32_t onWhichL = oi->localFc.on_which_L[i];
			
			Point closeP = mf.aabbTrees[onWhichL]->closest_point(segP);
			copyVariables[currIdsL[i] * 2] = closeP[0]; copyVariables[currIdsL[i] * 2 + 1] = closeP[1];
		}
	}

	void AESolverSquare::ExtractBoundaryTriangles(BaseDataStructure::QuadMesh *qm)
	{
		//std::vector<Eigen::Vector3i>().swap(boundaryTris_);
		boundaryTris_.clear();
		boundaryTris_.reserve(oi_->localTriIds.size());
		for (int i = 0; i < oi_->localTriIds.size(); ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				if (qm->Vs_[oi_->localToGlobalVVec[oi_->localTriIds[i][j]]].boundary)
				{
					boundaryTris_.emplace_back(oi_->localTriIds[i]);
					break;
				}
			}
		}
	}

	void AESolverSquare::Iterate(Eigen::VectorXd &pk, double nablaDotPk, double ek)
	{
		backupVars_ = currentVariables_;
		double alpha;
		if (qm_ == NULL)
			alpha = ComputeNoInverseMinAlpha(pk);
		else
			alpha = ComputeNoInverseMinAlphaForWholeQM(pk);
		
		int loopCount = 0;
		currentVariables_ = backupVars_ + alpha * pk;
		double currE = E();

		while ((currE > ek + 0.2*alpha*nablaDotPk && loopCount < INNER_LOOP_COUNT) || (proPreInv_ && !JudgeJacobiInOptimization(mf_) && loopCount < INNER_LOOP_COUNT))
		{
			alpha *= 0.8;
			currentVariables_ = backupVars_ + alpha * pk;
			currE = E();
			++loopCount;
		}

		if (loopCount == INNER_LOOP_COUNT)
		{
			std::cout << "AESolverSquare::Iterate::INNER_LOOP_COUNT " << std::endl;
		}

	}

	double AESolverSquare::QuadraticRoot1(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20)
	{
		double a = px10 * py20 - px20 * py10;
		double b = px10 * y20 + py20 * x10 - px20 * y10 - py10 * x20;
		double c = x10 * y20 - x20 * y10;

		return (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
	}
	double AESolverSquare::QuadraticRoot2(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20)
	{
		double a = px10 * py20 - px20 * py10;
		double b = px10 * y20 + py20 * x10 - px20 * y10 - py10 * x20;
		double c = x10 * y20 - x20 * y10;

		return (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
	}
	double AESolverSquare::LinearRoot(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20)
	{
		double b = px10 * y20 + py20 * x10 - px20 * y10 - py10 * x20;
		double c = x10 * y20 - x20 * y10;
		if (std::abs(b) < 1E-15)
			return 1;
		else
			return -c / b;
	}

	void AESolverSquare::SetValueBack(OptimizeInfo *oi)
	{
		for (int i = 0; i < oi->localVPos.size(); ++i)
		{
			oi->localVPos[i][0] = currentVariables_[2 * i];
			oi->localVPos[i][1] = currentVariables_[2 * i + 1];
		}
	}

	void AESolverSquare::SetValueFront(OptimizeInfo *oi)
	{
		for (int i = 0; i < oi->localVPos.size(); ++i)
		{
			/*oi->localVPos[i][0] = currentVariables_[2 * i];
			oi->localVPos[i][1] = currentVariables_[2 * i + 1];*/
			currentVariables_[2 * i] = oi->localVPos[i][0];
			currentVariables_[2 * i + 1] = oi->localVPos[i][1];
		}
	}
}
