#pragma once
#include <Eigen/Sparse>
#include <vector>
#include <ANN/ANN.h>
#include "ExtraFieldEnergy.h"

namespace DataOperation
{
	struct OptimizeInfo;
	struct FeatureConstraints;
	struct MeshFeatures;
	class AESolver
	{
	public:
		AESolver(OptimizeInfo *oi, ANNkd_tree *kdTree, std::vector<Eigen::Vector3d> &fields);
		~AESolver();

		template<typename T> 
		void SafeDeletePtr(T *ptrT)
		{
			if (ptrT != NULL)
			{
				delete ptrT;
				ptrT = NULL;
			}
		}
		template<typename T>
		void SafeDeletePtr(T **ptrT, int ptr1Num)
		{
			if (ptrT != NULL)
			{
				for (int i = 0; i < ptr1Num; ++i)
				{
					if (ptrT[i] != NULL)
					{
						delete[] ptrT[i];
						ptrT[i] = NULL;
					}
				}
				delete[] ptrT;
				ptrT = NULL;
			}
		}

		void Optimize();
		void Iterate(Eigen::VectorXd &pk, double nablaDotPk, double ek);
		double ComputeNoInverseMinAlpha(Eigen::VectorXd &pk);
		void ComputeNoInverseAlpha(Eigen::VectorXd &pk, double *alphaPtr);
		void ComputeNewCurrVariables(Eigen::VectorXd &backupVars, double *alphaPtr, Eigen::VectorXd &pk, Eigen::VectorXd &outVars);

		double E();
		double E_D();
		double E_T();
		double E_B();

		void Nabla_Hession_E();

		void BuildStandardInverseMatrix();
		void BuildTriAreas();
		void MatrixInverse(double *input, double *inverse);
		double ComputeSyDir(int index);
		void ComputeJacobi(int index, double *ptrJ);

		//SD_I1没必要有，因为后面的计算用不到
		double SD_I2(double i1, double i2, double i3);
		double SD_I3(double i1, double i2, double i3);

		double QuadraticRoot1(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20);
		double QuadraticRoot2(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20);
		double LinearRoot(double x10, double y10, double x20, double y20, double px10, double py10, double px20, double py20);

		void SetValueBack(OptimizeInfo *oi);
		void SetValueFront(OptimizeInfo *oi);

	public:
		OptimizeInfo *oi_ = NULL;
		Eigen::VectorXd currentVariables_;
		Eigen::VectorXd backupVars_;
		double (*standardInverseMatrix_)[4] = NULL;	//每一个standardInverseMatrix_中的元素都是一个2*2的矩阵，行主元保存的。
		//int sim1Num_ = 0;
		double *triAreas_ = NULL;

		int vVarNum_ = 0;
		int eleNum_ = 0;

		Eigen::VectorXd nablaE_;
		Eigen::SparseMatrix<double> hessionE_;

		//temp values,为了防止频繁分配内存而预先分配好的。
	public:
		double tempMat[4];
		double tempMat2[4];
		double tempMat3[4];

		double ptrJ[4];
		double currU[4], currV[4];
		Eigen::Matrix2d currJacobi;

		double currentF2[4], currentG[4];
		double d0[4], d1[4], l0[4], t0[4];
		double lambdaVec[4];
		double sdfq[4];
		double sd2_fq_fq[4][4];
		double leftMat[6][4];
		int mapV[6];

		bool isInverse = false;
		double et_, eb_, ed_, ef_;
		double *alphaPtr_;

		//ExtraFieldEnergy efe_;
	};

	extern double LAMBDA_D;
	extern double LAMBDA_B;
	extern double LAMBDA_T;
}