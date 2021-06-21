#include <iostream>
#include <memory>
#include "PolySquareDeformation.h"

#include "Optimization/Numeric/fmath.hpp"
#include "Optimization/HLBFGS/HLBFGS.h"
#include "Eigen/Eigen"

#define EXP_ENERGY

void energy_func(const size_t N, const std::vector<double>& x, double& f, std::vector<double>& g, void* user_supply)
{
	auto ptr_param = static_cast<std::pair<size_t, PolySquareDeformation*>*>(user_supply);
	(ptr_param->second)->deformation_energy_evalfunc(ptr_param->first, x, f, g);
};

void energy_iter(const size_t niter, const size_t call_iter, const size_t n_vars, const std::vector<double>& variables, const double& func_value, const std::vector<double>& gradient, const double& gnorm, void* user_pointer)
{
	if ((niter < 10) || (niter < 100 && (niter % 10 == 0)) || (niter % 100 == 0))
		std::cout << niter << " " << call_iter << " " << func_value << " " << gnorm << std::endl;
}

void PolySquareDeformation::calc_deformation()
{
	std::vector<std::vector<double>> chart_x;
	chart_x.reserve(atlas.size());
	for (int i = 0; i < atlas.size(); i++)
	{
		chart_x.emplace_back(atlas[i].n_vertices * 2);
	}
	for (int i = 0; i < para.n_vertices(); i++)
	{
		int cid = v_chart[i].first;
		int vid = v_chart[i].second;

		chart_x[cid][2 * vid + 0] = vertex_uv[i][0];
		chart_x[cid][2 * vid + 1] = vertex_uv[i][1];
	}
	
	for (int i = 0; i < atlas.size(); i++)
	{
		int n_vars = chart_x[i].size();

		HLBFGS energy_solver;
		energy_solver.set_number_of_variables(n_vars);
		energy_solver.set_verbose(true);
		energy_solver.set_func_callback(energy_func, 0, 0, /*energy_iter*/0, 0);

		energy_lambda = 1.0;
		std::pair<size_t, PolySquareDeformation*> rotation_param = { i, this };

		energy_solver.optimize_without_constraints(chart_x[i].data(), 2000, &rotation_param);
	}

	for (int i = 0; i < para.n_vertices(); i++)
	{
		int cid = v_chart[i].first;
		int vid = v_chart[i].second;

		vertex_uv[i][0] = chart_x[cid][2 * vid + 0];
		vertex_uv[i][1] = chart_x[cid][2 * vid + 1];
	}
}

void PolySquareDeformation::deformation_energy_evalfunc(int chart_id, const std::vector<double>& x, double & f, std::vector<double>& g)
{
	double align_energy = 0.0;
	const double align_factor = energy_lambda / atlas[chart_id].total_boundary_length;
	for (int i = 0; i < boundary_h[chart_id].size(); i++)
	{
		double length0 = boundary_h[chart_id][i].length;
		
		int v0 = v_chart[boundary_h[chart_id][i].v0].second;
		int v1 = v_chart[boundary_h[chart_id][i].v1].second;
		OpenMesh::Vec2d h_vec(x[2 * v1 + 0] - x[2 * v0 + 0], x[2 * v1 + 1] - x[2 * v0 + 1]);
		OpenMesh::Vec2d h_tag = tag_direction[boundary_h_tag[chart_id][i]];

		double len = h_vec.norm();
		h_vec /= len;

		align_energy += 2.0 * align_factor * length0 * (1 - OpenMesh::dot(h_vec, h_tag));

		double g_delta = -2.0 * align_factor * length0 * (h_tag[0] * h_vec[1] - h_tag[1] * h_vec[0]) / len;

		g[2 * v0 + 0] -= h_vec[1] * g_delta;
		g[2 * v1 + 0] += h_vec[1] * g_delta;
		g[2 * v0 + 1] += h_vec[0] * g_delta;
		g[2 * v1 + 1] -= h_vec[0] * g_delta;
	}

 	double mips_energy = 0.0;
	const double mips_factor = 1.0 / atlas[chart_id].total_area;
	const double exp_factor = energy_exp_factor;
	for (auto fid : chart_faces[chart_id])
	{
		double alpha = 0.5;
		double T[2][2];
		
		const auto& finfo = mesh_face_info[fid];
		int v0 = v_chart[finfo.v0].second;
		int v1 = v_chart[finfo.v1].second;
		int v2 = v_chart[finfo.v2].second;

		T[0][0] = x[2 * v1 + 0] - x[2 * v0 + 0];
		T[0][1] = x[2 * v2 + 0] - x[2 * v0 + 0];
		T[1][0] = x[2 * v1 + 1] - x[2 * v0 + 1];
		T[1][1] = x[2 * v2 + 1] - x[2 * v0 + 1];

		double l2_p1 = finfo.l2_p1;
		double l2_p2 = finfo.l2_p2;
		double dot_p = finfo.dot_p;
		double det_p = finfo.det_p;

		double det_u = T[0][0] * T[1][1] - T[0][1] * T[1][0];
		det_u = (finfo.normal_towards == 1) ? det_u : -det_u;
		double det = det_u / det_p;

		if (det < 1.0e-12 || abs(det_u) < 1.0e-20)
		{
			mips_energy += 1.0e+100;
			continue;
		}

		double xHx = l2_p2 * (T[0][0] * T[0][0] + T[1][0] * T[1][0])
				   + l2_p1 * (T[0][1] * T[0][1] + T[1][1] * T[1][1])
			 - 2.0 * dot_p * (T[0][0] * T[0][1] + T[1][0] * T[1][1]);
		double delta_det = 0.5 * (det + 1.0 / det);
		double delta_conf = 0.5 * (xHx / det_p / det_u);
		double delta = (1 - alpha) * delta_det + alpha * delta_conf;

		double area = det_p / 2.0;

		double delta_exp;
#ifdef EXP_ENERGY
		delta_exp = fmath::expd(exp_factor * delta);
		mips_energy += mips_factor * std::max(delta_exp - fmath::expd(exp_factor), 0.0) * area;
		delta_exp *= exp_factor;
#else
		delta_exp = 1.0;
		mips_energy += mips_factor * std::max(delta - 1.0, 0.0) * area;
#endif
		delta_exp = mips_factor * delta_exp * area;

		double d_xHx_u1 = -2.0 * (l2_p2 * T[0][0] - dot_p * T[0][1]);
		double d_xHx_u2 = -2.0 * (l2_p1 * T[0][1] - dot_p * T[0][0]);
		double d_xHx_v1 = -2.0 * (l2_p2 * T[1][0] - dot_p * T[1][1]);
		double d_xHx_v2 = -2.0 * (l2_p1 * T[1][1] - dot_p * T[1][0]);

		double factor_xHx = 0.5 * alpha / det_u / det_p * delta_exp;
		double factor_det = 0.5 * ((1 - alpha) * (1.0 / det_p - det_p / det_u / det_u) - alpha * xHx / det_p / det_u / det_u) * delta_exp;

		double g_u1 = factor_xHx * d_xHx_u1 - factor_det * T[1][1];
		double g_u2 = factor_xHx * d_xHx_u2 + factor_det * T[1][0];
		double g_v1 = factor_xHx * d_xHx_v1 + factor_det * T[0][1];
		double g_v2 = factor_xHx * d_xHx_v2 - factor_det * T[0][0];

		g[2 * v0 + 0] += g_u1 + g_u2;
		g[2 * v0 + 1] += g_v1 + g_v2;
		g[2 * v1 + 0] -= g_u1;
		g[2 * v1 + 1] -= g_v1;
		g[2 * v2 + 0] -= g_u2;
		g[2 * v2 + 1] -= g_v2;
	}

//	std::cout << "Ea " << align_energy << ", Em " << mips_energy << std::endl;
	f += align_energy + mips_energy;
}