#pragma once
#include "ChartDeformation.h"
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

class ChartUntangle
{
public:
	ChartUntangle(ChartDeformation& _parent);
	~ChartUntangle();

	bool calc();

private:
	ChartDeformation& parent;
	const std::vector<double>& ref_uv;
	std::vector<double>& uv;

	struct solution
	{
		Eigen::VectorXd u, v;

		solution() = default;
		solution(int n) { resize(n); };
		void resize(int n) { u.resize(n); v.resize(n); };
		void setZero() { u.setZero(); v.setZero(); };
	};

	void get_uv(solution& sol);
	void set_uv(const solution& sol);
	std::vector<int> sol2kkt;

	void LSCM_uv(solution& sol);

	void get_fv(int f, int& v0, int& v1, int& v2);

	double sum_det;
	std::vector<double> ref_det;//signed
	std::vector<Eigen::Matrix2d> inv_S;
	bool prepare_mapping();
	bool check_inversion(const solution& sol);

	//using Solver = Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>;
	using Solver = Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>;
	Solver hessian_solver;
	solution LSCM_b;
	bool assemble_hessian();
	void solve_step(const solution& grad, solution& step);

	double untangling_eps;

	bool calc_eps(const solution& sol);

	template <bool GRADIENT>
	double calc_energy(const solution& sol, solution& grad = solution());

	double calc_next_sol(const solution& sol0, solution& sol1, double& step_size, int& count);
};

