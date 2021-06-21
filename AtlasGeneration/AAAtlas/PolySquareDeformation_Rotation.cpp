#include <iostream>
#include "PolySquareDeformation.h"

#include "Optimization/HLBFGS/HLBFGS.h"

void angle_func(const size_t N, const std::vector<double>& x, double& f, std::vector<double>& g, void* user_supply)
{
	auto ptr_param = static_cast<std::pair<size_t, PolySquareDeformation*>*>(user_supply);
	(ptr_param->second)->rotating_angle_evalfunc(ptr_param->first, x, f, g);
};

void PolySquareDeformation::calc_global_rotation()
{
	HLBFGS angle_solver;
	angle_solver.set_number_of_variables(1);
	angle_solver.set_verbose(false);
	angle_solver.set_func_callback(angle_func, 0, 0, 0, 0);

	std::vector<double> global_theta(boundary_h.size(), 0.0);
	for (int i = 0; i < boundary_h.size(); i++)
	{
		std::pair<size_t, PolySquareDeformation*> rotation_param = { i, this };

		double theta[1] = { global_theta[i] };
		angle_solver.optimize_without_constraints(theta, 100, &rotation_param);
		global_theta[i] = theta[0];
	}
	
	std::vector<OpenMesh::Vec2d> chart_center(boundary_h.size(), OpenMesh::Vec2d(0.0, 0.0));
	for (int i = 0; i < para.n_vertices(); i++)
	{
		chart_center[v_chart[i].first] += vertex_uv[i];
	}

	std::vector<double> cos_theta(boundary_h.size()), sin_theta(boundary_h.size());
	for (int i = 0; i < boundary_h.size(); i++)
	{
		chart_center[i] /= atlas[i].n_vertices;
		cos_theta[i] = std::cos(global_theta[i]);
		sin_theta[i] = std::sin(global_theta[i]);
	}
	for (int i = 0; i < para.n_vertices(); i++)
	{
		int chart_id = v_chart[i].first;
		vertex_uv[i] -= chart_center[chart_id];
		rotate_uv(vertex_uv[i], cos_theta[chart_id], sin_theta[chart_id]);
		vertex_uv[i] += chart_center[chart_id];
	}

	vertex_uv_vec.push_back(vertex_uv);
}

void PolySquareDeformation::rotating_angle_evalfunc(int chart_id, const std::vector<double>& x, double & f, std::vector<double>& g)
{
	double cos_theta_t = cos(x[0]);
	double sin_theta_t = sin(x[0]);

	f = 0.0;
	g[0] = 0.0;

	for (int i = 0; i < boundary_h[chart_id].size(); i++)
	{
		OpenMesh::Vec2d h_vec = boundary_h_svec[chart_id][i];
		double length = boundary_h[chart_id][i].length;

		rotate_uv(h_vec, cos_theta_t, sin_theta_t);

		f += length * h_vec[0]* h_vec[0] * h_vec[1] * h_vec[1];
		g[0] += 2.0 * length * h_vec[0] * h_vec[1] * (h_vec[1] * h_vec[1] - h_vec[0] * h_vec[0]);
	}
}