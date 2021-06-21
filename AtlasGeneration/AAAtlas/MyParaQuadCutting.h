#pragma once
#include "MeshViewer/MeshDefinition.h"
#include "Viewer2D_Parameterization.h"

#include <Eigen/Eigen>
#include <memory>
#include <io.h>
#include <iostream>
#include <string>
#include <fstream>

class MyParaQuadCutting
{
public:
	MyParaQuadCutting(Mesh& _mesh, const char* path, const char* file);
	~MyParaQuadCutting() = default;

	void set_path(const char* str) { str_path = std::string(str); }
	void set_file(const char* str) { str_file = std::string(str); }
	const std::string& get_path() { return str_path; }

	void read_cut();
	bool is_point_on_edge(OpenMesh::Vec3d p_, OpenMesh::EdgeHandle edge_, Mesh &mesh_);
	bool is_point_in_face(OpenMesh::Vec3d p_, OpenMesh::FaceHandle face_, Mesh &mesh_);
	OpenMesh::HalfedgeHandle find_split_uv(OpenMesh::Vec3d p_, OpenMesh::FaceHandle face_, Mesh &mesh_, double &u, double &v);
	double find_split_t(OpenMesh::Vec3d p_, OpenMesh::HalfedgeHandle he_);
	void find_replace(OpenMesh::HalfedgeHandle& h0_p, OpenMesh::VertexHandle v_hf_p, OpenMesh::VertexHandle v_ht_p, OpenMesh::HalfedgeHandle& h0_m, OpenMesh::VertexHandle v_hf_m, OpenMesh::VertexHandle v_ht_m);
	OpenMesh::VertexHandle split_edge(OpenMesh::HalfedgeHandle h_para, double t);
	OpenMesh::VertexHandle split_face(OpenMesh::FaceHandle f_para, double u, double v);
	void trans_textured(Mesh& tar);
	void cut_segment();
	bool is_intersect(OpenMesh::Vec3d p1_, OpenMesh::Vec3d p2_, OpenMesh::Vec3d p3_, OpenMesh::Vec3d p4_);
	OpenMesh::Vec3d compute_intersect_point(OpenMesh::Vec3d p1_, OpenMesh::Vec3d p2_, OpenMesh::Vec3d p3_, OpenMesh::Vec3d p4_);
	OpenMesh::HalfedgeHandle find_next_point(OpenMesh::VertexHandle vertex_, OpenMesh::Vec3d segment1_, OpenMesh::Vec3d segment2_, double &t_);
	void process_point_in_face(OpenMesh::Vec3d p_);
	void update_cut_segments(OpenMesh::Vec3d p_old_, OpenMesh::Vec3d p_new_);
	void mark_e_cut(OpenMesh::VertexHandle v_pre_, OpenMesh::VertexHandle v_now_);
	void para_decomposition();

private:
	Mesh& origin_mesh;
	Mesh origin_para;
	Mesh new_mesh, new_para;
	Mesh mesh_unpack;
	std::vector<Mesh> para_atlas;
	std::string str_path, str_file;

	OpenMesh::Vec3d BB_Max, BB_Min;

	double cut_length = 0.0;
	double split_thres = 0.02;

	bool get_para_mesh();

	OpenMesh::EPropHandleT<int> e_cut;
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

	std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> cut_segments;
	std::vector<Mesh> my_para_atlas;
	std::vector<std::vector<int>> f_atlas2mesh;

	void init();
};