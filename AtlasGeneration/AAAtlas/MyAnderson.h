#pragma once

#include <Eigen/Core>
#include"Eigen/Dense"
#include "Eigen/SVD"
#include "Eigen/QR"

using namespace Eigen;

class Anderson
{
public:
	Anderson();
	~Anderson();
public:
	void replace(VectorXd & u)
	{
		current_u_ = u;
	}

	const VectorXd& compute(VectorXd & G)
	{
		assert(iter_ >= 0);

		current_F_ = G - current_u_;

		if (iter_ == 0)
		{
			prev_dF_.col(0) = -current_F_;
			prev_dG_.col(0) = -G;
			current_u_ = G;
		}
		else
		{
			prev_dF_.col(col_idx_) += current_F_;
			prev_dG_.col(col_idx_) += G;

			double eps = 1e-14;
			double scale = std::max(eps, prev_dF_.col(col_idx_).norm());
			dF_scale_(col_idx_) = scale;
			prev_dF_.col(col_idx_) /= scale;

			int m_k = std::min(m_, iter_);

			if (m_k == 1)
			{
				theta_(0) = 0;
				double dF_sqrnorm = prev_dF_.col(col_idx_).squaredNorm();
				M_(0, 0) = dF_sqrnorm;
				double dF_norm = std::sqrt(dF_sqrnorm);

				if (dF_norm > eps){
					// compute theta = (dF * F) / (dF * dF)
					theta_(0) = (prev_dF_.col(col_idx_) / dF_norm).dot(current_F_ / dF_norm);
				}
			}
			else
			{
				// Update the normal equation matrix, for the column and row corresponding to the new dF column
				VectorXd new_inner_prod = (prev_dF_.leftCols(m_k).transpose())* prev_dF_.col(col_idx_);
				M_.block(0, col_idx_, m_k, 1) = new_inner_prod;
				M_.block(col_idx_, 0, 1, m_k) = new_inner_prod.transpose();
				// Solve normal equation
				cod_.compute(M_.block(0, 0, m_k, m_k));
				theta_.head(m_k) = cod_.solve(prev_dF_.leftCols(m_k).transpose() * current_F_);
			}

			// Use rescaled theata to compute new u
			for (int k = 0; k < m_k; ++k)
			{
				theta_[k] = theta_[k] / dF_scale_[k];
			}
			current_u_ = G - prev_dG_.leftCols(m_k) * theta_.head(m_k);

			col_idx_ = (col_idx_ + 1) % m_;
			prev_dF_.col(col_idx_) = -current_F_;
			prev_dG_.col(col_idx_) = -G;
		}

		iter_++;
		return current_u_;
	}

	// m: number of previous iterations used
	// d: dimension of variables
	// u0: initial variable values
	void init(int m, int d, VectorXd & u0)
	{
		assert(m > 0);
		m_ = m;
		dim_ = d;
		v_n = u0.size();
		current_u_ = u0;
		prev_dG_.resize(v_n, m);
		prev_dF_.resize(v_n, m);
		M_.resize(m, m);
		theta_.resize(m);
		dF_scale_.resize(m);
		iter_ = 0;
		col_idx_ = 0;
	}

private:
	VectorXd current_u_;
	VectorXd current_F_;
	MatrixXd prev_dG_;
	MatrixXd prev_dF_;
	MatrixXd M_;		// Normal equations matrix for the computing theta
	VectorXd theta_;	// theta value computed from normal equations
	VectorXd dF_scale_;		// The scaling factor for each column of prev_dF
	Eigen::CompleteOrthogonalDecomposition<MatrixXd> cod_;

	int m_;		// Number of previous iterates used for Andreson Acceleration
	int dim_;	// Dimension of variables
	int v_n;	// Number of vertices * dim_
	int iter_;	// Iteration count since initialization
	int col_idx_;	// Index for history matrix column to store the next value
};

