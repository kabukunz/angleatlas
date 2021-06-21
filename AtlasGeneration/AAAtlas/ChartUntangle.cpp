#include <iostream>
#include "ChartUntangle.h"

ChartUntangle::ChartUntangle(ChartDeformation& _parent)
	:parent(_parent), ref_uv(*(_parent.uv_x_vec.end() - 2)), uv(_parent.uv_x)
{
}


ChartUntangle::~ChartUntangle()
{
}

bool ChartUntangle::calc()
{
	int nv = ref_uv.size() / 2;

	int max_iter = 2000;
	double step_size = 1.0;
	double theta = (1 - 0.4) / (1 + 0.4);
	bool no_inversion = false;
	
	std::cout << "--------------------Untangling--------------------" << std::endl;
	if (!prepare_mapping())
	{
		std::cout << "Input Error!" << std::endl;
		return false;
	}
	if (!assemble_hessian())
	{
		std::cout << "Hessian Error!" << std::endl;
		return false;
	}
	
	solution sol;
//	get_uv(sol);
	LSCM_uv(sol);
	std::cout << "LSCM Done";
	if (calc_eps(sol))
	{
		set_uv(sol);
		return true;
	}
	
	std::cout << std::endl << "EPS: " << untangling_eps << std::endl;
	solution prev_sol = sol;
	solution curr_sol(nv);

	int iter_count = 0;
	for (; iter_count < max_iter; iter_count++)
	{
		if (iter_count >= 1)
		{
			sol.u = (1 + theta) * curr_sol.u - theta * prev_sol.u;
			sol.v = (1 + theta) * curr_sol.v - theta * prev_sol.v;
			prev_sol = curr_sol;
		}

		int count_line = 0;
		double e1 = calc_next_sol(sol, curr_sol, step_size, count_line);

		if ((iter_count < 10) || (iter_count < 100 && (iter_count % 10 == 0)) || (iter_count % 100 == 0))
		{
			printf("[%4d] SS:%5.4e, E:%9.8e, %d\n", iter_count, step_size, e1, count_line);
		}

		if (check_inversion(curr_sol))
		{
			no_inversion = true;
			break;
		}
	}

	std::cout << "Iters: " << iter_count;

	set_uv(curr_sol);

	return no_inversion;
}

void ChartUntangle::get_uv(solution& sol)
{
	int nv = ref_uv.size() / 2;
	sol.resize(nv);
	for (int i = 0; i < nv; i++)
	{
		sol.u[i] = uv[2 * i + 0];
		sol.v[i] = uv[2 * i + 1];
	}
}

void ChartUntangle::set_uv(const solution& sol)
{
	int nv = ref_uv.size() / 2;
	for (int i = 0; i < nv; i++)
	{
		uv[2 * i + 0] = sol.u[i];
		uv[2 * i + 1] = sol.v[i];
	}
}

void ChartUntangle::LSCM_uv(solution& sol)
{
	int nv = ref_uv.size() / 2;

	solution LSCM_kkt;
	LSCM_kkt.u = hessian_solver.solve(LSCM_b.u);
	LSCM_kkt.v = hessian_solver.solve(LSCM_b.v);

	sol.resize(nv);
	for (int i = 0; i < nv; i++)
	{
		if (sol2kkt[i] < 0)
		{
			sol.u[i] = uv[2 * i + 0];
			sol.v[i] = uv[2 * i + 1];
		}
		else
		{
			sol.u[i] = LSCM_kkt.u[sol2kkt[i]];
			sol.v[i] = LSCM_kkt.v[sol2kkt[i]];
		}
	}
}

void ChartUntangle::get_fv(int f, int & v0, int & v1, int & v2)
{
	const auto& finfo = parent.parent.mesh_face_info[parent.mesh_faces[f]];
	v0 = parent.parent.v_chart[finfo.uv0].second;
	v1 = parent.parent.v_chart[finfo.uv1].second;
	v2 = parent.parent.v_chart[finfo.uv2].second;
}

bool ChartUntangle::prepare_mapping()
{
	double factor = 1.0 / parent.goal_length;

	sum_det = 0.0;
	inv_S.reserve(parent.mesh_faces.size());
	ref_det.reserve(parent.mesh_faces.size());
	for (int i = 0; i < parent.mesh_faces.size(); i++)
	{
		int v0, v1, v2;
		get_fv(i, v0, v1, v2);

		inv_S.emplace_back();
		Eigen::Matrix2d& inv = inv_S.back();
		inv(0, 0) = ref_uv[2 * v2 + 1] - ref_uv[2 * v0 + 1];// S(1, 1)
		inv(1, 0) = ref_uv[2 * v0 + 1] - ref_uv[2 * v1 + 1];//-S(1, 0)
		inv(0, 1) = ref_uv[2 * v0 + 0] - ref_uv[2 * v2 + 0];//-S(0, 1)
		inv(1, 1) = ref_uv[2 * v1 + 0] - ref_uv[2 * v0 + 0];// S(0, 0)

		inv *= factor;
		double det = inv.determinant();
		if ((parent.z_pos && det < 0.0) || (!parent.z_pos && det > 0.0)) return false;

		sum_det += std::abs(det);
		ref_det.push_back(det);
		inv /= det;
	}

	std::vector<bool> fixed_var(ref_uv.size() / 2, false);
	for (int i = 0; i < parent.n_boundary_edges; i++)
	{
		int v0 = parent.boundary_h_vert[i].first;
		fixed_var[v0] = true;
	}
	int n_kkt = 0;
	sol2kkt.resize(ref_uv.size() / 2);
	for (int i = 0; i < ref_uv.size() / 2; i++)
	{
		sol2kkt[i] = fixed_var[i] ? -1 : (n_kkt++);
	}
	return true;
}

bool ChartUntangle::check_inversion(const solution& sol)
{
	for (int i = 0; i < parent.mesh_faces.size(); i++)
	{
		int v0, v1, v2;
		get_fv(i, v0, v1, v2);

		Eigen::Matrix2d S;
		S(0, 0) = sol.u(v1) - sol.u(v0);
		S(1, 0) = sol.v(v1) - sol.v(v0);
		S(0, 1) = sol.u(v2) - sol.u(v0);
		S(1, 1) = sol.v(v2) - sol.v(v0);

		if (std::signbit(S.determinant()) != std::signbit(ref_det[i])) return false;
	}
	return true;
}

bool ChartUntangle::assemble_hessian()
{
	int nv = ref_uv.size() / 2;
	int nf = parent.mesh_faces.size();
	int nb = parent.n_boundary_edges;

	LSCM_b.resize(nv - nb);
	LSCM_b.setZero();
	std::vector<Eigen::Triplet<double>> hessian_tlist;
	auto set_triplet = [&](int i, int j, double val)
	{
		int kkt_i = sol2kkt[i];
		int kkt_j = sol2kkt[j];
		if (kkt_i >= kkt_j && kkt_j >= 0)
		{
			hessian_tlist.emplace_back(kkt_i, kkt_j, val);
		}
		else if (kkt_i >= 0 && kkt_j < 0)
		{
			LSCM_b.u[kkt_i] -= val * uv[2 * j + 0];
			LSCM_b.v[kkt_i] -= val * uv[2 * j + 1];
		}
	};

	hessian_tlist.reserve(9 * nf);

	for (int i = 0; i < nf; ++i)
	{
		int v0, v1, v2;
		get_fv(i, v0, v1, v2);

		const auto& T = inv_S[i];
		double det = std::abs(ref_det[i]);

		double h00 = det * ((T(0, 0) + T(1, 0))*(T(0, 0) + T(1, 0)) + (T(0, 1) + T(1, 1))*(T(0, 1) + T(1, 1)));
		double h01 = det * (-T(0, 0) * T(0, 0) - T(0, 0) * T(1, 0) - T(0, 1) * T(0, 1) - T(0, 1) * T(1, 1));
		double h02 = det * (-T(1, 0) * T(1, 0) - T(0, 0) * T(1, 0) - T(1, 1) * T(1, 1) - T(0, 1) * T(1, 1));
		double h11 = det * (T(0, 0)*T(0, 0) + T(0, 1) * T(0, 1));
		double h12 = det * (T(0, 0)*T(1, 0) + T(0, 1) * T(1, 1));
		double h22 = det * (T(1, 0)*T(1, 0) + T(1, 1) * T(1, 1));

		set_triplet(v0, v0, h00);
		set_triplet(v0, v1, h01);
		set_triplet(v0, v2, h02);

		set_triplet(v1, v0, h01);
		set_triplet(v1, v1, h11);
		set_triplet(v1, v2, h12);

		set_triplet(v2, v0, h02);
		set_triplet(v2, v1, h12);
		set_triplet(v2, v2, h22);
	}

	Eigen::SparseMatrix<double> hessian;
	hessian.resize(nv - nb, nv - nb);
	hessian.setFromTriplets(hessian_tlist.begin(), hessian_tlist.end());
	hessian.makeCompressed();
	
	hessian_solver.compute(hessian);
	return hessian_solver.info() == Eigen::Success;
}

void ChartUntangle::solve_step(const solution& grad, solution& step)
{
	int nv = ref_uv.size() / 2;
	int n_var = nv - parent.n_boundary_edges;
	solution grad_kkt(n_var), step_kkt(n_var);

	for (int i = 0; i < nv; i++)
	{
		if (sol2kkt[i] < 0) continue;
		grad_kkt.u[sol2kkt[i]] = -grad.u[i];
		grad_kkt.v[sol2kkt[i]] = -grad.v[i];
	}

	step_kkt.u = hessian_solver.solve(grad_kkt.u);
	step_kkt.v = hessian_solver.solve(grad_kkt.v);

// 	std::cout << "U diff " << (hessian * step_kkt.u - grad_kkt.u).norm() / grad_kkt.u.norm() << std::endl;
// 	std::cout << "V diff " << (hessian * step_kkt.v - grad_kkt.v).norm() / grad_kkt.v.norm() << std::endl;

	step.setZero();
	for (int i = 0; i < nv; i++)
	{
		if (sol2kkt[i] < 0) continue;
		step.u[i] = step_kkt.u[sol2kkt[i]];
		step.v[i] = step_kkt.v[sol2kkt[i]];
	}
}

bool ChartUntangle::calc_eps(const solution& sol)
{
	double min_det = std::numeric_limits<double>::max();
	for (int i = 0; i < parent.mesh_faces.size(); ++i)
	{
		int v0, v1, v2;
		get_fv(i, v0, v1, v2);

		Eigen::Matrix2d S;
		S(0, 0) = sol.u(v1) - sol.u(v0);
		S(1, 0) = sol.v(v1) - sol.v(v0);
		S(0, 1) = sol.u(v2) - sol.u(v0);
		S(1, 1) = sol.v(v2) - sol.v(v0);

		S *= inv_S[i];
		double det = S.determinant();

		min_det = std::min(det, min_det);
	}

	double eps = (min_det < 0) ? std::max(1.0e-10 * min_det * min_det, -2.0e-8 * min_det) : 1.0e-10;

	untangling_eps = eps;
	return (min_det > 0.0);
}

template <bool GRADIENT>
double ChartUntangle::calc_energy(const solution& sol, solution& grad)
{
	double e = 0.0;

	if constexpr (GRADIENT) grad.setZero();

	for (int i = 0; i < parent.mesh_faces.size(); ++i)
	{
		int v0, v1, v2;
		get_fv(i, v0, v1, v2);
		double det_factor = std::abs(ref_det[i]) / sum_det;

		Eigen::Matrix2d S;
		S(0, 0) = sol.u(v1) - sol.u(v0);
		S(1, 0) = sol.v(v1) - sol.v(v0);
		S(0, 1) = sol.u(v2) - sol.u(v0);
		S(1, 1) = sol.v(v2) - sol.v(v0);

		const auto& inv = inv_S[i];
		Eigen::Matrix2d J = S * inv;
		double Frob_J = J.squaredNorm();
		double det_J = J.determinant();

		double untangling_det_J = 0.5 * (det_J + std::sqrt(det_J * det_J + untangling_eps));
		double inv_untangling_det_J = 1.0 / std::max(1.0e-30, untangling_det_J);

		e += 0.5 * det_factor * Frob_J * inv_untangling_det_J;

		if constexpr (GRADIENT)
		{
			double d_J_00_x0 = -inv(0, 0) - inv(1, 0); double d_J_00_x1 = inv(0, 0); double d_J_00_x2 = inv(1, 0);
			double d_J_01_x0 = -inv(0, 1) - inv(1, 1); double d_J_01_x1 = inv(0, 1); double d_J_01_x2 = inv(1, 1);
			double d_J_10_y0 = -inv(0, 0) - inv(1, 0); double d_J_10_y1 = inv(0, 0); double d_J_10_y2 = inv(1, 0);
			double d_J_11_y0 = -inv(0, 1) - inv(1, 1); double d_J_11_y1 = inv(0, 1); double d_J_11_y2 = inv(1, 1);

			double d_JF_x0 = 2 * J(0, 0)*d_J_00_x0 + 2 * J(0, 1)*d_J_01_x0;
			double d_JF_x1 = 2 * J(0, 0)*d_J_00_x1 + 2 * J(0, 1)*d_J_01_x1;
			double d_JF_x2 = 2 * J(0, 0)*d_J_00_x2 + 2 * J(0, 1)*d_J_01_x2;
			double d_JF_y0 = 2 * J(1, 0)*d_J_10_y0 + 2 * J(1, 1)*d_J_11_y0;
			double d_JF_y1 = 2 * J(1, 0)*d_J_10_y1 + 2 * J(1, 1)*d_J_11_y1;
			double d_JF_y2 = 2 * J(1, 0)*d_J_10_y2 + 2 * J(1, 1)*d_J_11_y2;

			double d_det_J_x0 = d_J_00_x0 * J(1, 1) - d_J_01_x0 * J(1, 0);
			double d_det_J_x1 = d_J_00_x1 * J(1, 1) - d_J_01_x1 * J(1, 0);
			double d_det_J_x2 = d_J_00_x2 * J(1, 1) - d_J_01_x2 * J(1, 0);
			double d_det_J_y0 = J(0, 0)*d_J_11_y0 - J(0, 1)*d_J_10_y0;
			double d_det_J_y1 = J(0, 0)*d_J_11_y1 - J(0, 1)*d_J_10_y1;
			double d_det_J_y2 = J(0, 0)*d_J_11_y2 - J(0, 1)*d_J_10_y2;

			double temp_det = 1.0 / std::sqrt(det_J*det_J + untangling_eps);
			double d_U_det_J_x0 = 0.5*(d_det_J_x0 + det_J * d_det_J_x0 * temp_det);
			double d_U_det_J_x1 = 0.5*(d_det_J_x1 + det_J * d_det_J_x1 * temp_det);
			double d_U_det_J_x2 = 0.5*(d_det_J_x2 + det_J * d_det_J_x2 * temp_det);
			double d_U_det_J_y0 = 0.5*(d_det_J_y0 + det_J * d_det_J_y0 * temp_det);
			double d_U_det_J_y1 = 0.5*(d_det_J_y1 + det_J * d_det_J_y1 * temp_det);
			double d_U_det_J_y2 = 0.5*(d_det_J_y2 + det_J * d_det_J_y2 * temp_det);

			double d_e_x0 = 0.5 * det_factor * ((d_JF_x0 - Frob_J * d_U_det_J_x0 * inv_untangling_det_J) * inv_untangling_det_J);
			double d_e_x1 = 0.5 * det_factor * ((d_JF_x1 - Frob_J * d_U_det_J_x1 * inv_untangling_det_J) * inv_untangling_det_J);
			double d_e_x2 = 0.5 * det_factor * ((d_JF_x2 - Frob_J * d_U_det_J_x2 * inv_untangling_det_J) * inv_untangling_det_J);
			double d_e_y0 = 0.5 * det_factor * ((d_JF_y0 - Frob_J * d_U_det_J_y0 * inv_untangling_det_J) * inv_untangling_det_J);
			double d_e_y1 = 0.5 * det_factor * ((d_JF_y1 - Frob_J * d_U_det_J_y1 * inv_untangling_det_J) * inv_untangling_det_J);
			double d_e_y2 = 0.5 * det_factor * ((d_JF_y2 - Frob_J * d_U_det_J_y2 * inv_untangling_det_J) * inv_untangling_det_J);

			grad.u(v0) += d_e_x0;
			grad.v(v0) += d_e_y0;
			grad.u(v1) += d_e_x1;
			grad.v(v1) += d_e_y1;
			grad.u(v2) += d_e_x2;
			grad.v(v2) += d_e_y2;
		}
	}

	return e;
}

double ChartUntangle::calc_next_sol(const solution& sol0, solution& sol1, double& step_size, int& count)
{
	double min_step_size = 1e-16;
	int nv = ref_uv.size() / 2;
	solution grad(nv);
	solution step(nv);

	double e0 = calc_energy<true>(sol0, grad);
	solve_step(grad, step);

	//line search
	sol1.u = sol0.u + step_size * step.u;
	sol1.v = sol0.v + step_size * step.v;

	double e1 = calc_energy<false>(sol1);
	if (e1 <= e0)
	{
		step_size *= 2.0;
		return e1;
	}

	count = 0;
	while (e1 > e0)
	{
		if (step_size < min_step_size)
		{
			sol1 = sol0;
			return e0;
		}
		step_size *= 0.5;
		sol1.u = sol0.u + step_size * step.u;
		sol1.v = sol0.v + step_size * step.v;
		e1 = calc_energy<false>(sol1);

		count++;
	}

	return e1;
}