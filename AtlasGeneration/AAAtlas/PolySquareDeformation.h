#pragma once

#include "MeshViewer/MeshDefinition.h"

#include <set>

class ChartUntangle;
class ChartDeformation;
class PolySquareDeformation
{
public:
	PolySquareDeformation(const Mesh& _mesh, Mesh& _para, std::map<int, int>& _vk);
	virtual ~PolySquareDeformation();

	bool calc(double kernel, double goal, double exp, std::string str_);
	void set_phase(int i) { set_vertices(vertex_uv_vec[i]); };
	void set_path(const char* str) { str_path = std::string(str); }

	double get_goal_length() { return goal_length; };
	void gather_vk();

	int get_n_phase() { return vertex_uv_vec.size(); };

private:
	const Mesh& mesh;
	Mesh& para;
	std::map<int, int>& viewer_vk;

	friend class ChartDeformation;
	friend class ChartUntangle;

	std::string str_path;
	bool amips_exp = true;
	double kernel_width = 4.0;
//	double goal_factor = 0.1;
	double energy_exp_factor = 1.0;

	double goal_length;

	OpenMesh::Vec2d tag_direction[4] = { { 1.0, 0.0 },{ 0.0, 1.0 },{ -1.0, 0.0 },{ 0.0, -1.0 } };

	std::vector<std::vector<OpenMesh::Vec2d>> vertex_uv_vec;

	struct face_info
	{
		int normal_towards;
		int uv0, uv1, uv2;
		double l2_p1, l2_p2, dot_p, det_p;
	};
	inline double face_area(int f) { return std::abs(mesh_face_info[f].det_p) / 2.0; };
	std::vector<face_info> mesh_face_info;

	bool load_face_info();
	void calc_face_vectors(int f, OpenMesh::Vec2d& p1, OpenMesh::Vec2d& p2, const std::vector<OpenMesh::Vec2d>& points);

	void build_boundary_cycles();
	void build_charts_info();

	std::map<int, std::pair<int, int>> h_boundary;

	std::vector<ChartDeformation> atlas;
	std::vector<std::pair<int, int>> v_chart;
	std::vector<int> f_chart;

	inline void rotate_uv(OpenMesh::Vec2d& u, const double& cos_theta, const double& sin_theta) const
	{
		double t = cos_theta * u[0] + sin_theta * u[1];
		u[1] = -sin_theta * u[0] + cos_theta * u[1];
		u[0] = t;
	}

	void set_uv();
	void get_uv();
	void check_flip(const std::vector<OpenMesh::Vec2d>& vertex_uv);

	void set_vertices(const std::vector<OpenMesh::Vec2d>& vertex_uv);
	void color_edges();
};