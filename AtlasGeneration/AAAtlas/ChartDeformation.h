#pragma once

#include <set>
#include "PolySquareDeformation.h"

#include "Eigen\Eigen"
#include "Eigen\Sparse"

#include "Scaffold/TriangleInterface.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <iostream>

class ChartUntangle;
class ChartDeformation
{
public:
	ChartDeformation(PolySquareDeformation& _parent);
	~ChartDeformation();

	void build_chart();
	void build_vertices();
	void build_faces();

	void deformation(std::string str1_);

	int n_vertices, n_edges, n_faces;
	int n_boundary_edges;
	double total_area, total_boundary_length;

	double kernel_width, energy_lambda;
	double align_alpha, amips_alpha;
	double energy_align, energy_amips, max_angle_align;
	double goal_length;

	bool z_pos;

	std::vector<int> boundary_h_meshid;
	std::vector<int> boundary_h_tag;
	std::vector<std::pair<int, int>> boundary_h_vert;
	std::vector<double> boundary_h_len0;
	std::vector<OpenMesh::Vec2d> boundary_h_svec;

	std::vector<int> boundary_v_k;
	std::vector<double> boundary_v_sangle;
	std::vector<double> boundary_v_angle0;

	std::vector<int> mesh_faces;
	std::vector<PolySquareDeformation::face_info> chart_face_info;
	std::set<int> fixed_uv_x;

	std::vector<double> uv_x;
	std::vector<std::vector<double>> uv_x_vec;

	std::vector<std::vector<double>> lbfgs_vec;

	void rotating_angle_evalfunc(const std::vector<double>& x, double& f, std::vector<double>& g = std::vector<double>());

	template <bool GRADIENT>
	double align_energy_evalfunc(const std::vector<double>& x, double& f, std::vector<double>& g = std::vector<double>());
	template <bool EXP_ENERGY, bool GRADIENT, bool FIX_POINTS>
	double amips_energy_evalfunc(const std::vector<double>& x, double& f, std::vector<double>& g = std::vector<double>());
	template <bool GRADIENT, bool HESSIAN_P, bool FIX_POINTS>
	double dirichlet_energy_evalfunc(const std::vector<double>& x, double& f, std::vector<double>& g = std::vector<double>(), Eigen::SparseMatrix<double>& h = Eigen::SparseMatrix<double>());

	//My 
	std::vector<double> uv_x_plus;
	std::vector<PolySquareDeformation::face_info> chart_face_info_plus;
	int n_vertices_plus;
	double total_area_plus;
	std::vector<OpenMesh::Vec3d> build_frame(Mesh &mesh_);
	void build_boundary(Mesh &mesh_, std::vector<OpenMesh::Vec3d> &boundary_, std::vector<int> &boundary_size_, std::vector<int> &boundary_id_);
	void build_scaffold(Mesh &mesh_);
	void uvplus2uv();
	void uv2uvplus();
	template <bool EXP_ENERGY, bool GRADIENT, bool FIX_POINTS>
	double amips_energy_evalfunc1(const std::vector<double>& x, double& f, std::vector<double>& g = std::vector<double>());
	template <bool GRADIENT, bool HESSIAN_P, bool FIX_POINTS>
	double dirichlet_energy_evalfunc1(const std::vector<double>& x, double& f, std::vector<double>& g = std::vector<double>(), Eigen::SparseMatrix<double>& h = Eigen::SparseMatrix<double>());
	Mesh mesh_temp;
	void get_align_direction();
	void write_txt(const std::string & outstr);
	void write_obj(const std::string & outstr);
	bool detect_intersect(OpenMesh::Vec3d start_point1, OpenMesh::Vec3d end_point1, OpenMesh::Vec3d start_point2, OpenMesh::Vec3d end_point2);
	std::vector<int> corner_vid_list;
	std::vector<int> boundary_cc_list;
	std::vector<int> boundary_cce_list;
	std::vector<std::vector<int>> segment_list;
	std::vector<int> origin_corner_list;
	void construct_segment();
	double polar_angle(OpenMesh::Vec2d vec2);
	double get_error(OpenMesh::Vec2d e1, OpenMesh::Vec2d e2);
	OpenMesh::Vec2d get_normal(int edge_);
	OpenMesh::Vec2d get_proxy(std::vector<int> edge_part);
	int get_seed(std::vector<int> edge_part, OpenMesh::Vec2d proxy_);
	void merge_hole(std::vector<std::vector<int>> &segment_part_, std::vector<int> hole_);
	void merge_segment_part(std::vector<std::vector<int>> &segment_part_, int merge_id_);
	std::vector<int> get_new_segmentpart(std::vector<int> hole_);
	std::vector<double> lloyd_iteration();
	std::vector<double> align_direction;
	double pre_angle(int edge_id_);
	double next_angle(int edge_id_);
	template <bool GRADIENT>
	double align_energy_evalfunc1(const std::vector<double>& x, double& f, std::vector<double>& g = std::vector<double>());
	void get_align_direction_corner();
	std::vector<double> h_tags_pre;
	double energy_amips_threshold;

	template <bool GRADIENT>
	double line_energy_evalfunc(const std::vector<double>& x, double & f, std::vector<double>& g = std::vector<double>());
	template <bool EXP_ENERGY>
	void calc_line_energy();
	template <bool EXP_ENERGY>
	void calc_line_deformation(int max_iter_times);
	void construct_line();
	std::vector<std::vector<int>> line_list;
	double angle_lambda;
	double line_lambda;
	

private:
	PolySquareDeformation& parent;

	inline int boundary(int bid);

	template <bool EXP_ENERGY>
	void calc_align_energy();

	void calc_boundary_directions();

	void calc_global_rotation();
	template <bool EXP_ENERGY>
	void calc_align_deformation(int max_iter_times);
	template <bool EXP_ENERGY>
	void calc_inner_deformation(int max_iter_times);
	template <bool EXP_ENERGY>
	void calc_final_deformation(int max_iter_times);

	bool CM_inner_deformation(int max_iter_times);
	//My
	bool CM_inner_deformation1(int max_iter_times);
	double calc_max_step1(const std::vector<double>& uv, const std::vector<double>& p);
	double backtracking_line_search1(const std::vector<double>& uv, const std::vector<double>& g, const std::vector<double>& p, double& alpha);
	template <bool EXP_ENERGY>
	void calc_align_energy1();
	template <bool EXP_ENERGY>
	void calc_align_deformation1(int max_iter_times);
	//
	double backtracking_line_search(const std::vector<double>& uv, const std::vector<double>& g, const std::vector<double>& p, double& alpha);
	double calc_max_step(const std::vector<double>& uv, const std::vector<double>& p);

	void find_corners();
	void move_corners();
	void segment_flattening();
	bool mosek_flattening(std::vector<double>& sol_x);

	void polysquare_post_deformation();

	std::vector<int> uv2kkt;
	Eigen::Matrix<double, 6, 6> H_det;
	void hessian_preparation();
	void update_uv2kkt();

	void update_face_info();

	void scale_finfo(double scale);
	void check_vk();

	struct segment_info
	{
		int begin, end;
		int tag;
		int size;
		double coordinate;
		double length0;
	};
	std::vector<segment_info> segments;

	double interior_angle(int next_boundary_h, const std::vector<double>& uv);

	inline OpenMesh::Vec2d get_vec(int v0, int v1, const std::vector<double>& uv) { return OpenMesh::Vec2d(uv[2 * v1 + 0] - uv[2 * v0 + 0], uv[2 * v1 + 1] - uv[2 * v0 + 1]); };
	inline double tag_dot(const OpenMesh::Vec2d& p, int t) { return (t & 2) ? -p[t & 1] : p[t & 1]; };

	std::vector<std::pair<int, int>> aligned_uv;

private:
	struct corner_info
	{
		int bvid;
		int status;//bit01: tag, bit2: deleted
		double seg_length;

		corner_info() : bvid(-1), status(0), seg_length(0.0) {};
	};
	std::list<corner_info> boundary_corners;
	using corner = std::list<corner_info>::iterator;
	corner next_corner(corner it);
	corner prev_corner(corner it);
	void calc_seglength(corner it, const std::vector<double>& bh_len);
	void update_segs_len(int uv_id = -1);

 	bool modify_short_segments(double thres_factor = 1.0);
 	bool conflicting_segments();
	void tag_from_segments();
	void built_segments();

	friend class ChartUntangle;

	double max_sigma = 20.0;
	double min_sigma = 0.0;
};