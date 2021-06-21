#ifndef PARAQUADCUTTING_H
#define PARAQUADCUTTING_H

#include "MeshViewer/MeshDefinition.h"
#include "ParaQuadChartDecomposition.h"
#include "Viewer2D_Parameterization.h"
#include <Eigen/Eigen>
#include "Eigen\Dense"

#include "Eigen\Sparse"
#include <memory>

class ParaQuadCutting
{
public:
	struct retangle_type
	{
		retangle_type() :min(std::numeric_limits<int>::max(), std::numeric_limits<int>::max()), max(std::numeric_limits<int>::min(), std::numeric_limits<int>::min()) {};
		retangle_type(const OpenMesh::Vec2i& _min, const OpenMesh::Vec2i& _max) :min(_min), max(_max) {};

		OpenMesh::Vec2i min, max;

		inline int width() const { return max[0] - min[0]; };
		inline int height() const { return max[1] - min[1]; };
		inline int area() const { return width() * height(); };
	};
	ParaQuadCutting(Mesh& _mesh, const char* path, const char* file);
	~ParaQuadCutting() = default;

	void set_path(const char* str) { str_path = std::string(str); }
	void set_file(const char* str) { str_file = std::string(str); }
	const std::string& get_path() { return str_path; }

	void cutting(double thres);
	void save_para(const char* filename);

	double get_split_thres() { return split_thres; };

	void set_seleted_vertices(const std::vector<int>& mesh_v);
	void set_seleted_edges(const std::vector<int>& mesh_e);
	void set_seleted_faces(const std::vector<int>& mesh_f);
	void set_seleted(const std::vector<int>& mesh_v, const std::vector<int>& mesh_e, const std::vector<int>& mesh_f);

	void set_viewer_title(const QString& title) { viewer->setWindowTitle(title); };
	void toggle_viewer_window() { viewer->setHidden(!viewer->isHidden()); };

	void get_textured_mesh(Mesh& tar);
	void trans_textured(Mesh& tar);
	void update_textured_mesh(Mesh& tar, bool use_new = false);
	void split_edges(const std::vector<int>& mesh_e);

	Mesh& get_origin_para() { return origin_para; };
	void update_para(double factor = 1.0);
	void flip_neg_charts();

	const OpenMesh::Vec3d& get_BB_Max() { return BB_Max; };
	const OpenMesh::Vec3d& get_BB_Min() { return BB_Min; };

	std::unique_ptr<Viewer2D_Parameterization> viewer;

	const std::set<int>& get_flipped_faces() { return flipped_faces; };

	std::map<int, OpenMesh::Vec3d> old_point;

	double calc_distortion(bool silence);
	double get_distortion() { return para_distortion; };

	void get_scaf_info(Eigen::MatrixXd& v_pos, Eigen::MatrixXd& uv_v_pos, Eigen::MatrixXi& fv_id, Eigen::MatrixXi& uv_fv_id, Eigen::MatrixXd& f_n);
	void load_from_scaf(const Eigen::MatrixXd& uv_v_pos);

	//my
	void split_poly();
	void read_cut();
	void read_cut1();
	void read_cut2();
	bool is_point_on_edge(OpenMesh::Vec3d p_, OpenMesh::EdgeHandle edge_, Mesh &mesh_, double average_length);
	bool is_point_in_face(OpenMesh::Vec3d p_, OpenMesh::FaceHandle face_, Mesh &mesh_);
	OpenMesh::HalfedgeHandle find_split_uv(OpenMesh::Vec3d p_, OpenMesh::FaceHandle face_, Mesh &mesh_, double &u, double &v);
	double find_split_t(OpenMesh::Vec3d p_, OpenMesh::HalfedgeHandle he_);
	void cut_segment(std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> &my_cut_segments_);
	void correct_segment(OpenMesh::Vec3d &segment1_, OpenMesh::Vec3d &segment2_, std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> &my_cut_segments_);
	bool is_intersect(OpenMesh::Vec3d p1_, OpenMesh::Vec3d p2_, OpenMesh::Vec3d p3_, OpenMesh::Vec3d p4_);
	OpenMesh::Vec3d compute_intersect_point(OpenMesh::Vec3d p1_, OpenMesh::Vec3d p2_, OpenMesh::Vec3d p3_, OpenMesh::Vec3d p4_);
	OpenMesh::HalfedgeHandle find_next_point(OpenMesh::VertexHandle vertex_, OpenMesh::Vec3d segment1_, OpenMesh::Vec3d segment2_, double &t_);
	void process_point_in_face(OpenMesh::Vec3d p_, std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> &my_cut_segments_);
	void update_cut_segments(OpenMesh::Vec3d p_old_, OpenMesh::Vec3d p_new_ , std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> &my_cut_segments_);
	void mark_e_cut(OpenMesh::VertexHandle v_pre_, OpenMesh::VertexHandle v_now_);
	void my_decomposition();
	void my_cutting(double thres);
	void my_write_obj(Mesh mesh_, std::string str_);
	void set_boundary(int ii_, std::vector<OpenMesh::VertexHandle> &cc_boundary_, std::vector<OpenMesh::Vec3d> &new_bpoint_);
	void uniform_tutte();
	void set_boundary1(int ii_, std::vector<OpenMesh::VertexHandle> &cc_boundary_, std::vector<OpenMesh::Vec3d> &new_bpoint_);
	void uniform_tutte1(bool is_reconstruct_);
	double my_retangle_packing(std::vector<retangle_type> my_ret_candidate_);
	void my_packing();
	void find_BB(Mesh mesh_, OpenMesh::Vec3d &b_max_, OpenMesh::Vec3d &b_min_);
	void my_get_scaf_info(Eigen::MatrixXd& v_pos, Eigen::MatrixXd& uv_v_pos, Eigen::MatrixXi& fv_id, Eigen::MatrixXi& uv_fv_id, Mesh para_, Mesh mesh_);
	void find_new_cut(bool my_switch);
	void build_rotation_message(Mesh mesh_);
	void my_calc_global_rotation();
	void my_rotating_angle_evalfunc(const std::vector<double>& x, double& f, std::vector<double>& g);
	inline void my_rotate_uv(OpenMesh::Vec2d& u, const double& cos_theta, const double& sin_theta) const
	{
		double t = cos_theta * u[0] + sin_theta * u[1];
		u[1] = -sin_theta * u[0] + cos_theta * u[1];
		u[0] = t;
	}

private:
	Mesh& origin_mesh;
	Mesh origin_para;
	Mesh new_mesh, new_para;
	Mesh mesh_unpack;
	std::vector<Mesh> para_atlas;

	OpenMesh::Vec3d BB_Max, BB_Min;

	double para_distortion = 0.0;

	double cut_length = 0.0;
	double split_thres = 0.02;

	std::set<int> flipped_faces;

	std::vector<std::vector<int>> v_para2layer;
	std::vector<int> face_layer;

	bool get_para_mesh();

	OpenMesh::EPropHandleT<int> e_segment;
	OpenMesh::HPropHandleT<std::pair<int, int>> h_mesh2para, h_para2mesh;

	std::vector<int> origin_h_mesh2para;

	OpenMesh::HalfedgeHandle get_mesh2para(OpenMesh::HalfedgeHandle h_mesh)
	{
		const auto& h_pair = new_mesh.property(h_mesh2para, h_mesh);
		return new_para.find_halfedge(new_para.vertex_handle(h_pair.first), new_para.vertex_handle(h_pair.second));
	};
	OpenMesh::HalfedgeHandle get_para2mesh(OpenMesh::HalfedgeHandle h_para)
	{
		const auto& h_pair = new_para.property(h_para2mesh, h_para);
		return new_mesh.find_halfedge(new_mesh.vertex_handle(h_pair.first), new_mesh.vertex_handle(h_pair.second));
	};

	void cut_on_boundary();
	void cut_segments();
	void decomposition();
	void get_packing_result();

	std::vector<int> get_seed_triangles();

	void init();
	std::string str_path, str_file;
	std::unique_ptr<ParaQuadChartDecomposition> charts_decomposition;
	void save_tri_charts();
	void save_quad_charts();

	std::vector<int> v_quad2para;

	void find_cut_path(OpenMesh::HalfedgeHandle& h_prev, int k, int dst_quad_v, int seg_id);

	OpenMesh::VertexHandle split_edge(OpenMesh::HalfedgeHandle h_para, double t);
	OpenMesh::VertexHandle split_face(OpenMesh::FaceHandle f_para, double u, double v);

	double find_split_t(OpenMesh::HalfedgeHandle h, int val_u, int val_v);
	OpenMesh::HalfedgeHandle find_split_uv(OpenMesh::FaceHandle f, int val_u, int val_v, double& u, double& v);

	void find_replace(OpenMesh::HalfedgeHandle& h0_p, OpenMesh::VertexHandle v_hf_p, OpenMesh::VertexHandle v_ht_p, OpenMesh::HalfedgeHandle& h0_m, OpenMesh::VertexHandle v_hf_m, OpenMesh::VertexHandle v_ht_m);

	inline uint get_tag(const OpenMesh::Vec3d& vec);

	//my
	OpenMesh::EPropHandleT<int> e_cut;
	std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> my_cut_segments;
	std::vector<Mesh> my_para_atlas;
	std::vector<Mesh> my_mesh_atlas;
	std::vector<Mesh> my_origin_atlas;
	std::vector<std::vector<int>> f_atlas2mesh;
	std::vector<std::vector<int>> f_atlas2para;
	std::vector<std::vector<int>> v_atlas2mesh;

	std::vector<retangle_type> my_ret_candidate;
	std::vector<int> my_cf_candidate;
	std::vector<OpenMesh::Vec2d> my_ct_candidate;
	double my_scale;
	std::vector<int> vertex_on_cut;
	double my_total_area;
	std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> my_cut_segments_new;
	
    //rotation
	int n_vertices;
	int n_boundary_edges;
	std::vector<OpenMesh::Vec2d> my_boundary_h_svec;
	std::vector<double> my_boundary_h_len0;
	std::vector<double> my_uv_x;

public:
	double calc_signal_error_3D();
};

#endif // PARAQUADCUTTING_H