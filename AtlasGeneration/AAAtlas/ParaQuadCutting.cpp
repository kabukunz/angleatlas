#include "ParaQuadCutting.h"
#include "MeshViewer/Mesh_doubleIO.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <queue>

#include <Eigen/Eigen>
#include <Eigen/CholmodSupport>

#include <QTime>

#include "finders_interface.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "Scaffold/StateManager.h"

#include "Common\CommonFunctions.h"
#include "Optimization\Numeric\FastMath.h"
#include "Mosek\ConvexQuadOptimization.h"
#include "Optimization\HLBFGS\HLBFGS.h"

#include "MyParafun.h"

#ifdef MY_DEBUG
#define MY_DOUBT(cond, msg) if (cond) std::cout << msg << std::endl;
#else
#define MY_DOUBT(cond, msg)
#endif

OpenMesh::Vec3uc boundary_edge_color = { 229, 156,  59 };
OpenMesh::Vec3uc interior_edge_color = { 100, 100, 100 };

ParaQuadCutting::ParaQuadCutting(Mesh& _mesh, const char* path, const char* file)
	: origin_mesh(_mesh)
{
	str_path = std::string(path);
	str_file = std::string(file);
	init();
}

void ParaQuadCutting::init()
{
	Mesh_doubleIO::copy_mesh(origin_mesh, new_mesh);

	if (!get_para_mesh()) return;

	viewer = std::make_unique<Viewer2D_Parameterization>(origin_para);
	viewer->show();
	viewer->adjust_range();
	
	new_para.add_property(e_segment);
	for (auto e_h : new_para.edges())
	{
		new_para.property(e_segment, e_h) = -1;
	}
}

uint ParaQuadCutting::get_tag(const OpenMesh::Vec3d& vec)
{
	uint b = (abs(vec[0]) <= abs(vec[1]));
	uint a = (vec[b] <= 0.0);

	return (a << 1) + b;
}

double ParaQuadCutting::calc_signal_error_3D()
{
	constexpr int split_n = 8;
	OpenMesh::Vec3d mesh_BB_max, mesh_BB_min;
	mesh_BB_max = mesh_BB_min = origin_mesh.point(origin_mesh.vertex_handle(0));
	for (auto v_h : origin_mesh.vertices())
	{
		const auto& p = origin_mesh.point(v_h);
		mesh_BB_max = mesh_BB_max.maximize(p);
		mesh_BB_min = mesh_BB_min.minimize(p);
	}

	OpenMesh::Vec3d grid_length = (mesh_BB_max - mesh_BB_min) / 100.0;

	auto get_color = [&](const OpenMesh::Vec3d& pt)
	{
		auto pt1 = (pt - mesh_BB_min) / grid_length;
		bool block = (((std::lround(std::floor(pt1[0])) + std::lround(std::floor(pt1[1])) + std::lround(std::floor(pt1[2]))) & 1) == 0);

		return block ? OpenMesh::Vec3d(1.0, 1.0, 1.0) : OpenMesh::Vec3d(0.0, 0.0, 0.0);
	};

	auto get_interior_point = [](const std::vector<OpenMesh::Vec3d>& v_list, double u, double v)
	{
		return u * v_list[0] + v * v_list[1] + (1.0 - u - v) * v_list[2];
	};

	double delta_u = (BB_Max - BB_Min).max() / 4096.0;
	double delta_v = (BB_Max - BB_Min).max() / 4096.0;
	auto signal_error = [&](const std::vector<OpenMesh::Vec3d>& v_list, const std::vector<OpenMesh::Vec3d>& uv_list, int i0, int j0, int i1, int j1, int i2, int j2)
	{
		std::vector<OpenMesh::Vec3d> pt(3), pt_color(3);
		double divisor = 1.0 / split_n;

		pt[0] = get_interior_point(uv_list, i0 * divisor, j0 * divisor);
		pt[1] = get_interior_point(uv_list, i1 * divisor, j1 * divisor);
		pt[2] = get_interior_point(uv_list, i2 * divisor, j2 * divisor);

		pt_color[0] = get_color(get_interior_point(v_list, i0 * divisor, j0 * divisor));
		pt_color[1] = get_color(get_interior_point(v_list, i1 * divisor, j1 * divisor));
		pt_color[2] = get_color(get_interior_point(v_list, i2 * divisor, j2 * divisor));

		Eigen::Matrix3d para_mat;
		Eigen::Matrix3d signal_mat;

		for (int i = 0; i < 3; i++)
		{
			para_mat(0, i) = pt[i][0] / delta_u;
			para_mat(1, i) = pt[i][1] / delta_v;
			para_mat(2, i) = 1.0;

			signal_mat(0, i) = pt_color[i][0];
			signal_mat(1, i) = pt_color[i][1];
			signal_mat(2, i) = pt_color[i][2];
		}

		Eigen::Matrix3d J_mat = signal_mat * para_mat.inverse();
		return J_mat.col(0).squaredNorm() + J_mat.col(1).squaredNorm();
	};

	double error = 0.0;
	double total_area = 0.0;
	for (auto f_h : new_mesh.faces())
	{
		double s_area = new_mesh.calc_sector_area(*new_mesh.fh_begin(f_h));
		total_area += s_area;
		s_area /= split_n * split_n;

		std::map<int, int> v_mesh2uv;

		std::vector<OpenMesh::Vec3d> fv_p, fuv_p;
		fv_p.reserve(3);
		fuv_p.reserve(3);
		for (auto fh_h : new_mesh.fh_range(f_h))
		{
			fv_p.push_back(new_mesh.point(new_mesh.to_vertex_handle(fh_h)));
			fuv_p.push_back(new_para.point(new_para.to_vertex_handle(get_mesh2para(fh_h))));
		}

		double f_error = 0.0;
		for (int j = 0; j <= split_n; j++)
		{
			for (int k = 0; k <= split_n - j; k++)
			{
				if (j + k < split_n) f_error += s_area * signal_error(fv_p, fuv_p, j, k, j + 1, k, j, k + 1);

				if (j > 0 && k > 0) f_error += s_area * signal_error(fv_p, fuv_p, j, k, j - 1, k, j, k - 1);
			}
		}

		error += f_error;
	}

	error /= 12.0 * total_area;

	std::cout << "3D Checkerboard Signal Error " << error << std::endl;

	return error;
}

bool ParaQuadCutting::get_para_mesh()
{
	origin_para.clear();
	new_para.clear();

	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
	OpenMesh::HPropHandleT<int> hvt_index;

	if (!origin_mesh.get_property_handle(mvt_list, "mvt_list") || !origin_mesh.get_property_handle(hvt_index, "hvt_index"))
	{
		std::cout << "Texture data is invalid." << std::endl;
		return false;
	}

	for (int i = 0; i < origin_mesh.property(mvt_list).size(); i++)
	{
		origin_para.add_vertex(Mesh::Point(origin_mesh.property(mvt_list)[i][0], origin_mesh.property(mvt_list)[i][1], 0.0));
		new_para.add_vertex(Mesh::Point(origin_mesh.property(mvt_list)[i][0], origin_mesh.property(mvt_list)[i][1], 0.0));
	}

	for (int i = 0; i < origin_mesh.n_faces(); i++)
	{
		auto f_h = origin_mesh.face_handle(i);
		std::vector<OpenMesh::VertexHandle> para_f;
		para_f.reserve(3);

		for (auto fh_h = origin_mesh.cfh_begin(f_h); fh_h != origin_mesh.cfh_end(f_h); fh_h++)
		{
			para_f.push_back(origin_para.vertex_handle(origin_mesh.property(hvt_index, *fh_h)));
		}

		origin_para.add_face(para_f);
		new_para.add_face(para_f);
	}

	origin_h_mesh2para.resize(origin_mesh.n_halfedges(), -1);
	new_mesh.add_property(h_mesh2para);
	new_para.add_property(h_para2mesh);
	for (int i = 0; i < origin_mesh.n_faces(); i++)
	{
		auto f0 = origin_mesh.face_handle(i);
		auto f1 = new_mesh.face_handle(i);
		auto f2 = new_para.face_handle(i);
		
		std::map<int, int> v_para2mesh;
		for (auto f0_h : origin_mesh.fh_range(f0))
		{
			v_para2mesh[origin_mesh.property(hvt_index, f0_h)] = origin_mesh.to_vertex_handle(f0_h).idx();
		}

		std::map<int, int> mesh_v2h;
		for (auto f1_h : new_mesh.fh_range(f1))
		{
			mesh_v2h[new_mesh.to_vertex_handle(f1_h).idx()] = f1_h.idx();
		}

		for (auto f2_h : new_para.fh_range(f2))
		{
			auto f1_h = new_mesh.halfedge_handle(mesh_v2h[v_para2mesh[new_para.to_vertex_handle(f2_h).idx()]]);

			new_mesh.property(h_mesh2para, f1_h) = std::make_pair(new_para.from_vertex_handle(f2_h).idx(), new_para.to_vertex_handle(f2_h).idx());
			new_para.property(h_para2mesh, f2_h) = std::make_pair(new_mesh.from_vertex_handle(f1_h).idx(), new_mesh.to_vertex_handle(f1_h).idx());

			origin_h_mesh2para[f1_h.idx()] = f2_h.idx();
		}
	}

	cut_length = 0.0;
	OpenMesh::EPropHandleT<bool> e_oncut;
	origin_mesh.get_property_handle(e_oncut, "e_oncut");
	for (auto e_h : origin_mesh.edges())
	{
		if (origin_mesh.is_boundary(e_h))
		{
			cut_length += origin_mesh.calc_edge_length(e_h);
			continue;
		}
		
		int to0 = origin_mesh.property(hvt_index, origin_mesh.halfedge_handle(e_h, 0));
		int to1 = origin_mesh.property(hvt_index, origin_mesh.halfedge_handle(e_h, 1));
		int from0 = origin_mesh.property(hvt_index, origin_mesh.prev_halfedge_handle(origin_mesh.halfedge_handle(e_h, 0)));
		int from1 = origin_mesh.property(hvt_index, origin_mesh.prev_halfedge_handle(origin_mesh.halfedge_handle(e_h, 1)));

		bool is_cut = ((to0 != from1) || (to1 != from0));
		origin_mesh.property(e_oncut, e_h) = is_cut;
		if (is_cut) cut_length += origin_mesh.calc_edge_length(e_h) * 2.0;
	}

	origin_mesh.request_vertex_normals();

	new_para.request_edge_status();
	new_para.request_face_status();

	new_mesh.request_edge_status();
	new_mesh.request_face_status();

	origin_para.request_face_normals();
	origin_para.update_face_normals();

	origin_para.request_edge_colors();
	for (auto e_h : origin_para.edges())
	{
		origin_para.set_color(e_h, origin_para.is_boundary(e_h) ? boundary_edge_color : interior_edge_color);
	}
//	std::cout << "Euler of Para : " << origin_para.n_vertices() + origin_para.n_faces() - origin_para.n_edges() << std::endl;

	BB_Max = BB_Min = origin_para.point(origin_para.vertex_handle(0));
	for (auto v_h : origin_para.vertices())
	{
		const auto& p = origin_para.point(v_h);
		BB_Max = BB_Max.maximize(p);
		BB_Min = BB_Min.minimize(p);
	}
//	std::cout << "BB of Para : X [" << BB_Min[0] << ", " << BB_Max[0] << "] Y [" << BB_Min[1] << ", " << BB_Max[1] << "]" << std::endl;
	return true;
}

double ParaQuadCutting::calc_distortion(bool silence)
{
	int nf = origin_mesh.n_faces();

	int n_pos = 0;
	int n_neg = 0;
	double total_area = 0.0;
	double total_uv_area = 0.0;

	OpenMesh::Vec3d mesh_BB_max, mesh_BB_min;
	std::vector<double> face_area(nf);
	mesh_BB_max = mesh_BB_min = origin_mesh.point(origin_mesh.vertex_handle(0));
	for (auto f_h : origin_mesh.faces())
	{
		OpenMesh::Vec3d mesh_p[3];
		OpenMesh::Vec3d para_p[3];

		auto fh_iter = origin_mesh.cfh_begin(f_h);
		for (int i = 0; i < 3; i++, fh_iter++)
		{
			mesh_p[i] = origin_mesh.point(origin_mesh.to_vertex_handle(*fh_iter));
			para_p[i] = origin_para.point(origin_para.to_vertex_handle(get_mesh2para(*fh_iter)));

			mesh_BB_max = mesh_BB_max.maximize(mesh_p[i]);
			mesh_BB_min = mesh_BB_min.minimize(mesh_p[i]);
		}

		mesh_p[1] -= mesh_p[0];
		mesh_p[2] -= mesh_p[0];
		para_p[1] -= para_p[0];
		para_p[2] -= para_p[0];

		face_area[f_h.idx()] = OpenMesh::cross(mesh_p[1], mesh_p[2]).norm() / 2.0;

		total_area += face_area[f_h.idx()];
		total_uv_area += OpenMesh::cross(para_p[1], para_p[2]).norm() / 2.0;
	}
	
	double factor = std::sqrt(total_uv_area / total_area);
	std::vector<double> distor(nf);
	for (auto f_h : origin_mesh.faces())
	{
		OpenMesh::Vec3d mesh_p[3];
		OpenMesh::Vec3d para_p[3];

		auto fh_iter = origin_mesh.cfh_begin(f_h);
		for (int i = 0; i < 3; i++, fh_iter++)
		{
			mesh_p[i] = origin_mesh.point(origin_mesh.to_vertex_handle(*fh_iter));
			para_p[i] = origin_para.point(origin_para.to_vertex_handle(get_mesh2para(*fh_iter)));
		}
		Eigen::Matrix2d mesh_M, para_M;

		mesh_p[1] -= mesh_p[0];
		mesh_p[2] -= mesh_p[0];
		para_p[1] -= para_p[0];
		para_p[2] -= para_p[0];

		mesh_p[1] *= factor;
		mesh_p[2] *= factor;

		OpenMesh::Vec3d e1 = mesh_p[1].normalized();
		OpenMesh::Vec3d e2 = OpenMesh::cross(origin_mesh.normal(f_h), e1);

		mesh_M(0, 0) = mesh_p[1].norm();
		mesh_M(1, 0) = 0.0;
		mesh_M(0, 1) = OpenMesh::dot(mesh_p[2], e1);
		mesh_M(1, 1) = OpenMesh::dot(mesh_p[2], e2);

		para_M(0, 0) = para_p[1][0];
		para_M(1, 0) = para_p[1][1];
		para_M(0, 1) = para_p[2][0];
		para_M(1, 1) = para_p[2][1];

		double det_p = para_M.determinant();

		(det_p > 0 ? n_pos : n_neg)++;

		if (det_p <= 0)
		{
			para_M.row(0) = -para_M.row(0);
			flipped_faces.insert(f_h.idx());
		}

		Eigen::Matrix2d J = para_M * mesh_M.inverse();
		Eigen::JacobiSVD<Eigen::Matrix2d> SVD_solver;

		SVD_solver.compute(J);
		Eigen::Vector2d singulars = SVD_solver.singularValues();

		double det_J = J.determinant();
		double s_max = singulars.maxCoeff();
		double s_min = singulars.minCoeff();

		distor[f_h.idx()] = (s_max * s_max + s_min * s_min + 1.0 / s_max / s_max + 1.0 / s_min / s_min) * 0.25;
	}
	//std::cout << distor[0] << std::endl;
	double x_avg_w = 0.0;
	for (int i = 0; i < distor.size(); i++) x_avg_w += distor[i] * face_area[i];
	x_avg_w /= total_area;

	
	para_distortion = x_avg_w;

	if (!silence)
	{
		std::cout << "PE " << total_uv_area / viewer->get_bb_area() << std::endl;
		std::cout << "BL " << cut_length / (mesh_BB_max - mesh_BB_min).norm() << std::endl;
		std::cout << "ED " << x_avg_w << std::endl;
	}
	return total_uv_area / viewer->get_bb_area();
}

void ParaQuadCutting::get_scaf_info(Eigen::MatrixXd& v_pos, Eigen::MatrixXd& uv_v_pos, Eigen::MatrixXi& fv_id, Eigen::MatrixXi& uv_fv_id, Eigen::MatrixXd& f_n)
{
	v_pos.resize(new_mesh.n_vertices(), 3);
	uv_v_pos.resize(new_para.n_vertices(), 2);
	fv_id.resize(new_mesh.n_faces(), 3);
	uv_fv_id.resize(new_para.n_faces(), 3);
	f_n.resize(new_mesh.n_faces(), 3);

	for (int i = 0; i < new_mesh.n_vertices(); i++)
	{
		const auto& pt = new_mesh.point(new_mesh.vertex_handle(i));
		v_pos(i, 0) = pt[0];
		v_pos(i, 1) = pt[1];
		v_pos(i, 2) = pt[2];
	}

	for (int i = 0; i < new_para.n_vertices(); i++)
	{
		const auto& pt = new_para.point(new_para.vertex_handle(i));
		uv_v_pos(i, 0) = pt[0];
		uv_v_pos(i, 1) = pt[1];
	}

	for (int i = 0; i < new_mesh.n_faces(); i++)
	{
		int j = 0;
		for (auto fh_it = new_mesh.cfh_begin(new_mesh.face_handle(i)); fh_it != new_mesh.cfh_end(new_mesh.face_handle(i)); fh_it++)
		{
			fv_id(i, j++) = new_mesh.to_vertex_handle(fh_it.handle()).idx();
		}
	}

	for (int i = 0; i < new_para.n_faces(); i++)
	{
		int j = 0;
		for (auto fh_it = new_mesh.cfh_begin(new_mesh.face_handle(i)); fh_it != new_mesh.cfh_end(new_mesh.face_handle(i)); fh_it++)
		{
			uv_fv_id(i, j++) = new_para.to_vertex_handle(get_mesh2para(fh_it.handle())).idx();
		}
	}

	for (int i = 0; i < origin_mesh.n_faces(); i++)
	{
		OpenMesh::Vec3d pt = origin_mesh.normal(origin_mesh.face_handle(i));
		f_n(i, 0) = pt[0];
		f_n(i, 1) = pt[1];
		f_n(i, 2) = pt[2];
	}
}

void ParaQuadCutting::load_from_scaf(const Eigen::MatrixXd& uv_v_pos)
{
	for (int i = 0; i < new_para.n_vertices(); i++)
	{
		auto& pt = new_para.point(new_para.vertex_handle(i));
		pt[0] = uv_v_pos(i, 0);
		pt[1] = uv_v_pos(i, 1);
	}
}

void ParaQuadCutting::flip_neg_charts()
{
	std::vector<bool> f_visited(origin_para.n_faces(), false);
	for (auto f_h : origin_para.faces())
	{
		if (f_visited[f_h.idx()] || origin_para.normal(f_h)[2] > 0) continue;

		std::set<int> comp_v;
		std::queue<int> queue_f;

		queue_f.push(f_h.idx());
		while (!queue_f.empty())
		{
			int cur_f = queue_f.front();
			queue_f.pop();

			if (f_visited[cur_f]) continue;
			auto cur_fh = origin_para.face_handle(cur_f);
			f_visited[cur_f] = true;

			for (auto fv : origin_para.fv_range(cur_fh))
			{
				comp_v.insert(fv.idx());
			}

			for (auto fh : origin_para.fh_range(cur_fh))
			{
				auto next_fh = origin_para.opposite_face_handle(fh);
				if (next_fh.is_valid() && !f_visited[next_fh.idx()])
				{
					queue_f.push(next_fh.idx());
				}
			}
		}

		double x_max = -std::numeric_limits<double>::max();
		double x_min = std::numeric_limits<double>::max();
		for (int vid : comp_v)
		{
			double v_x = origin_para.point(origin_para.vertex_handle(vid))[0];
			x_max = std::max(v_x, x_max);
			x_min = std::min(v_x, x_min);
		}

		for (int vid : comp_v)
		{
			double& v_x = new_para.point(new_para.vertex_handle(vid))[0];
			v_x = x_max + x_min - v_x;
		}
	}
}

void ParaQuadCutting::set_seleted_vertices(const std::vector<int>& mesh_v)
{
	viewer->selected_v.clear();
	for (int v : mesh_v)
	{
		auto v_h = origin_mesh.vertex_handle(v);
		for (auto vih : origin_mesh.vih_range(v_h))
		{
			if (origin_mesh.is_boundary(vih)) continue;

			int para_v = origin_para.to_vertex_handle(origin_para.halfedge_handle(origin_h_mesh2para[vih.idx()])).idx();
			viewer->selected_v.insert(para_v);
		}
	}

	viewer->updateGL();
}

void ParaQuadCutting::set_seleted_edges(const std::vector<int>& mesh_e)
{
	viewer->selected_e.clear();
	for (int e : mesh_e)
	{
		viewer->selected_e.insert(origin_h_mesh2para[2 * e + 0] / 2);
		viewer->selected_e.insert(origin_h_mesh2para[2 * e + 1] / 2);
	}

	viewer->updateGL();
}

void ParaQuadCutting::set_seleted_faces(const std::vector<int>& mesh_f)
{
	viewer->selected_f.clear();
	viewer->selected_f.insert(mesh_f.begin(), mesh_f.end());

	viewer->updateGL();
}

void ParaQuadCutting::set_seleted(const std::vector<int>& mesh_v, const std::vector<int>& mesh_e, const std::vector<int>& mesh_f)
{
	set_seleted_vertices(mesh_v);
	set_seleted_edges(mesh_e);
	set_seleted_faces(mesh_f);
}

void ParaQuadCutting::get_textured_mesh(Mesh& tar)
{
	Mesh_doubleIO::copy_mesh(new_mesh, tar);

	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
	OpenMesh::HPropHandleT<int> hvt_index;

	if (!tar.get_property_handle(mvt_list, "mvt_list") || !tar.get_property_handle(hvt_index, "hvt_index"))
	{
		std::cout << "Texture data is invalid." << std::endl;
		return;
	}

	tar.property(mvt_list).resize(new_para.n_vertices());
	for (auto v_h : new_para.vertices())
	{
		auto p = new_para.point(v_h);
		tar.property(mvt_list)[v_h.idx()][0] = p[0];
		tar.property(mvt_list)[v_h.idx()][1] = p[1];
	}

	for (auto f_h : new_mesh.faces())
	{
		std::map<int, int> v_mesh2para;
		for (auto fh_h : new_mesh.fh_range(f_h))
		{
			int v0 = new_mesh.to_vertex_handle(fh_h).idx();
			int v1 = new_mesh.property(h_mesh2para, fh_h).second;

			v_mesh2para[v0] = v1;
		}

		for (auto fh_h : tar.fh_range(tar.face_handle(f_h.idx())))
		{
			tar.property(hvt_index, fh_h) = v_mesh2para.at(tar.to_vertex_handle(fh_h).idx());
		}
	}
}

void ParaQuadCutting::trans_textured(Mesh& tar)
{
	Mesh_doubleIO::copy_mesh(new_mesh, tar);

	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list0, mvt_list1;
	OpenMesh::HPropHandleT<int> hvt_index0, hvt_index1;

	tar.get_property_handle(mvt_list0, "mvt_list");
	tar.get_property_handle(hvt_index0, "hvt_index");
	new_mesh.get_property_handle(mvt_list1, "mvt_list");
	new_mesh.get_property_handle(hvt_index1, "hvt_index");

	tar.property(mvt_list0) = std::move(new_mesh.property(mvt_list1));
	for (auto h1 : new_mesh.halfedges())
	{
		auto h0 = tar.find_halfedge(new_mesh.from_vertex_handle(h1), new_mesh.to_vertex_handle(h1));
		tar.property(hvt_index0, h0) = new_mesh.property(hvt_index1, h1);
	}
}

void ParaQuadCutting::update_textured_mesh(Mesh& tar, bool use_new /*= false*/)
{
	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
	OpenMesh::HPropHandleT<int> hvt_index;

	if (!tar.get_property_handle(mvt_list, "mvt_list") || !tar.get_property_handle(hvt_index, "hvt_index"))
	{
		std::cout << "Texture data is invalid." << std::endl;
		return;
	}

	const Mesh& src = use_new ? new_para : origin_para;

	if (tar.property(mvt_list).size() != src.n_vertices())
	{
		std::cout << "Para vertices do not match." << std::endl;
	}

	for (auto v_h : src.vertices())
	{
		auto p = src.point(v_h);
		tar.property(mvt_list)[v_h.idx()][0] = p[0];
		tar.property(mvt_list)[v_h.idx()][1] = p[1];
	}
}

void ParaQuadCutting::split_edges(const std::vector<int>& mesh_e)
{
	for (int i : mesh_e)
	{
		split_edge(get_mesh2para(new_mesh.halfedge_handle(2 * i)), 0.5);
	}

	new_mesh.garbage_collection();
	new_para.garbage_collection();
}

void ParaQuadCutting::update_para(double factor /*= 1.0*/)
{
	viewer->update();
	viewer->para_scale(factor);
//	viewer->adjust_range();

	if (new_para.n_vertices() != origin_para.n_vertices())
	{
		std::cout << "Para vertices do not match." << std::endl;
	}

	for (auto v_h : new_para.vertices())
	{
		new_para.point(v_h) = origin_para.point(origin_para.vertex_handle(v_h.idx()));
	}
}

void ParaQuadCutting::cutting(double thres)
{
	QTime timer;
	split_thres = thres;

	charts_decomposition = std::make_unique<ParaQuadChartDecomposition>(new_para);

	timer.start();
	charts_decomposition->decomposition();
	
	cut_segments();
	decomposition();
	std::cout << "Time: " << timer.elapsed() / 1000.0 << "s" << std::endl;

	get_packing_result();

// 	save_quad_charts();
// 	save_tri_charts();
}

OpenMesh::VertexHandle ParaQuadCutting::split_edge(OpenMesh::HalfedgeHandle h_para, double t)
{
//	MY_DOUBT(t < split_thres || t > 1.0 - split_thres, "Split Edge Error");

	if (new_para.status(h_para).deleted())
	{
		return OpenMesh::VertexHandle();
	}

	auto h0_mesh = get_para2mesh(h_para);

	bool is_bh_mesh[2];
	OpenMesh::VertexHandle v_new[2];
	OpenMesh::HalfedgeHandle h_split[2];
	OpenMesh::HalfedgeHandle h_split_mesh[2] = { h0_mesh, new_mesh.opposite_halfedge_handle(h0_mesh) };

	//boundary may change after a face is splited, so boundary info need to save
	is_bh_mesh[0] = new_mesh.is_boundary(h_split_mesh[0]);
	is_bh_mesh[1] = new_mesh.is_boundary(h_split_mesh[1]);

	auto insert_p_mesh = (1 - t) * new_mesh.point(new_mesh.to_vertex_handle(h_split_mesh[1])) + t * new_mesh.point(new_mesh.to_vertex_handle(h_split_mesh[0]));
	auto v_new_mesh = new_mesh.add_vertex(insert_p_mesh);

	bool is_e_cut = 0;
	OpenMesh::EdgeHandle my_edge = new_para.edge_handle(h_para);
	int fromv_id, tov_id;
	if (new_para.get_property_handle(e_cut, "e_cut"))
	{
		if (new_para.property(e_cut, my_edge) == 1)
		{
			is_e_cut = 1;
			fromv_id = new_para.from_vertex_handle(h_para).idx();
			tov_id = new_para.to_vertex_handle(h_para).idx();
		}
	}
	

	for (int i = 0; i < 2; i++)
	{
		if (is_bh_mesh[i]) continue;

		h_split[i] = get_mesh2para(h_split_mesh[i]);

		OpenMesh::HalfedgeHandle h_face[2], h_new[2];
		OpenMesh::VertexHandle v_hf[2], v_ht[2];

		OpenMesh::HalfedgeHandle h_face_mesh[2], h_new_mesh[2];
		OpenMesh::VertexHandle v_hf_mesh[2], v_ht_mesh[2];

		if (i == 0)
		{
			auto v_hsf = new_para.from_vertex_handle(h_split[i]);
			auto v_hst = new_para.to_vertex_handle(h_split[i]);

			auto insert_p = (1.0 - t) * new_para.point(v_hsf) + t * new_para.point(v_hst);
			v_new[0] = new_para.add_vertex(insert_p);
		}
		else if (i == 1)
		{
			if (new_para.edge_handle(h_split[0]) != new_para.edge_handle(h_split[1]))
			{
				auto v_hsf = new_para.from_vertex_handle(h_split[i]);
				auto v_hst = new_para.to_vertex_handle(h_split[i]);

				auto insert_p = t * new_para.point(v_hsf) + (1.0 - t) * new_para.point(v_hst);
				v_new[1] = new_para.add_vertex(insert_p);
			}
			else
			{
				v_new[1] = v_new[0];
			}
		}

		h_face[0] = new_para.next_halfedge_handle(h_split[i]);
		h_face[1] = new_para.prev_halfedge_handle(h_split[i]);
		auto f_h = new_para.face_handle(h_split[i]);

		h_face_mesh[0] = new_mesh.next_halfedge_handle(h_split_mesh[i]);
		h_face_mesh[1] = new_mesh.prev_halfedge_handle(h_split_mesh[i]);
		auto f_h_mesh = new_mesh.face_handle(h_split_mesh[i]);

		for (int j = 0; j < 2; j++)
		{
			v_hf[j] = new_para.from_vertex_handle(h_face[j]);
			v_ht[j] = new_para.to_vertex_handle(h_face[j]);

			v_hf_mesh[j] = new_mesh.from_vertex_handle(h_face_mesh[j]);
			v_ht_mesh[j] = new_mesh.to_vertex_handle(h_face_mesh[j]);
		}

		new_para.delete_face(f_h, false);
		new_mesh.delete_face(f_h_mesh, false);

		new_para.add_face(v_hf[0], v_ht[0], v_new[i]);
		new_mesh.add_face(v_hf_mesh[0], v_ht_mesh[0], v_new_mesh);
		new_para.add_face(v_hf[1], v_ht[1], v_new[i]);
		new_mesh.add_face(v_hf_mesh[1], v_ht_mesh[1], v_new_mesh);

		for (int j = 0; j < 2; j++)
		{
			find_replace(h_face[j], v_hf[j], v_ht[j], h_face_mesh[j], v_hf_mesh[j], v_ht_mesh[j]);
		}

		new_para.property(h_para2mesh, new_para.next_halfedge_handle(h_face[0])) = std::make_pair(v_ht_mesh[0].idx(), v_new_mesh.idx());
		new_para.property(h_para2mesh, new_para.prev_halfedge_handle(h_face[1])) = std::make_pair(v_new_mesh.idx(), v_hf_mesh[1].idx());

		new_mesh.property(h_mesh2para, new_mesh.next_halfedge_handle(h_face_mesh[0])) = std::make_pair(v_ht[0].idx(), v_new[i].idx());
		new_mesh.property(h_mesh2para, new_mesh.prev_halfedge_handle(h_face_mesh[1])) = std::make_pair(v_new[i].idx(), v_hf[1].idx());

		new_para.property(h_para2mesh, new_para.prev_halfedge_handle(h_face[0])) = std::make_pair(v_new_mesh.idx(), v_hf_mesh[0].idx());
		new_para.property(h_para2mesh, new_para.next_halfedge_handle(h_face[1])) = std::make_pair(v_ht_mesh[1].idx(), v_new_mesh.idx());

		new_mesh.property(h_mesh2para, new_mesh.prev_halfedge_handle(h_face_mesh[0])) = std::make_pair(v_new[i].idx(), v_hf[0].idx());
		new_mesh.property(h_mesh2para, new_mesh.next_halfedge_handle(h_face_mesh[1])) = std::make_pair(v_ht[1].idx(), v_new[i].idx());
	
		if (new_para.get_property_handle(e_cut, "e_cut"))
		{
			if (is_e_cut)
			{
				OpenMesh::HalfedgeHandle he1 = new_para.find_halfedge(v_new[i], new_para.vertex_handle(fromv_id));
				OpenMesh::HalfedgeHandle he2 = new_para.find_halfedge(v_new[i], new_para.vertex_handle(tov_id));
				OpenMesh::EdgeHandle e1 = new_para.edge_handle(he1.idx() / 2);
				OpenMesh::EdgeHandle e2 = new_para.edge_handle(he2.idx() / 2);
				new_para.property(e_cut, e1) = 1;
				new_para.property(e_cut, e2) = 1;
			}
		}
	}

	return v_new[0];
}

OpenMesh::VertexHandle ParaQuadCutting::split_face(OpenMesh::FaceHandle f_para, double u, double v)
{
	MY_DOUBT(u < split_thres || v < split_thres || u + v > 1.0 - split_thres, "Split Face Error");

	std::vector<OpenMesh::HalfedgeHandle> h_face, h_face_mesh;
	std::vector<OpenMesh::VertexHandle> v_ht, v_hf, v_ht_mesh, v_hf_mesh;
	double w = 1.0 - u - v;

	v_ht.reserve(3);
	v_hf.reserve(3);
	h_face.reserve(3);

	v_ht_mesh.reserve(3);
	v_hf_mesh.reserve(3);
	h_face_mesh.reserve(3);

	bool skip_f = new_para.status(f_para).deleted();
	for (auto fh_h : new_para.fh_range(f_para))
	{
		if (skip_f || new_para.status(fh_h).deleted())
		{
			skip_f = true;
			break;
		}
		h_face.push_back(fh_h);
		v_ht.push_back(new_para.to_vertex_handle(fh_h));
		v_hf.push_back(new_para.from_vertex_handle(fh_h));

		const auto& fh_h_pair = new_para.property(h_para2mesh, fh_h);
		auto vf_mesh = new_mesh.vertex_handle(fh_h_pair.first);
		auto vt_mesh = new_mesh.vertex_handle(fh_h_pair.second);

		h_face_mesh.push_back(new_mesh.find_halfedge(vf_mesh, vt_mesh));
		v_ht_mesh.push_back(vt_mesh);
		v_hf_mesh.push_back(vf_mesh);
	}

	if (skip_f)
	{
		return OpenMesh::VertexHandle();
	}

	auto insert_p = u * new_para.point(v_ht[0]) + v * new_para.point(v_ht[1]) + w * new_para.point(v_ht[2]);
	auto insert_p_mesh = u * new_mesh.point(v_ht_mesh[0]) + v * new_mesh.point(v_ht_mesh[1]) + w * new_mesh.point(v_ht_mesh[2]);

	auto v_new = new_para.add_vertex(insert_p);
	auto v_new_mesh = new_mesh.add_vertex(insert_p_mesh);

	new_para.delete_face(f_para, false);
	new_mesh.delete_face(new_mesh.face_handle(f_para.idx()), false);
	for (int i = 0; i < 3; i++)
	{
		new_para.add_face(v_hf[i], v_ht[i], v_new);
		new_mesh.add_face(v_hf_mesh[i], v_ht_mesh[i], v_new_mesh);

		find_replace(h_face[i], v_hf[i], v_ht[i], h_face_mesh[i], v_hf_mesh[i], v_ht_mesh[i]);
	}

	for (auto ve : new_para.ve_range(v_new))
	{
		new_para.property(e_segment, ve) = -1;
	}

	for (int i = 0; i < 3; i++)
	{
		auto h0_para = new_para.find_halfedge(v_new, v_ht[i]);
		auto h1_para = new_para.opposite_halfedge_handle(h0_para);

		new_para.property(h_para2mesh, h0_para) = std::make_pair(v_new_mesh.idx(), v_ht_mesh[i].idx());
		new_para.property(h_para2mesh, h1_para) = std::make_pair(v_ht_mesh[i].idx(), v_new_mesh.idx());

		auto h0_mesh = new_mesh.find_halfedge(v_new_mesh, v_ht_mesh[i]);
		auto h1_mesh = new_mesh.opposite_halfedge_handle(h0_mesh);

		new_mesh.property(h_mesh2para, h0_mesh) = std::make_pair(v_new.idx(), v_ht[i].idx());
		new_mesh.property(h_mesh2para, h1_mesh) = std::make_pair(v_ht[i].idx(), v_new.idx());
	}

	return v_new;
}

void ParaQuadCutting::find_replace(OpenMesh::HalfedgeHandle& h0_p, OpenMesh::VertexHandle v_hf_p, OpenMesh::VertexHandle v_ht_p, OpenMesh::HalfedgeHandle& h0_m, OpenMesh::VertexHandle v_hf_m, OpenMesh::VertexHandle v_ht_m)
{
	if (!new_para.status(h0_p).deleted()) return;

	auto new_h0_p = new_para.find_halfedge(v_hf_p, v_ht_p);
	auto new_h0_m = new_mesh.find_halfedge(v_hf_m, v_ht_m);

	new_para.property(h_para2mesh, new_h0_p) = new_para.property(h_para2mesh, h0_p);
	new_mesh.property(h_mesh2para, new_h0_m) = new_mesh.property(h_mesh2para, h0_m);

	h0_p = new_h0_p;
	h0_m = new_h0_m;
}

double ParaQuadCutting::find_split_t(OpenMesh::HalfedgeHandle h, int val_u, int val_v)
{
	const auto& uv0 = new_para.point(new_para.from_vertex_handle(h));
	const auto& uv1 = new_para.point(new_para.to_vertex_handle(h));
	OpenMesh::Vec3d p(val_u, val_v, 0.0);

	auto v_e = uv1 - uv0;
	auto v_p = p - uv0;

	return OpenMesh::dot(v_e, v_p) / v_e.sqrnorm();
}

OpenMesh::HalfedgeHandle ParaQuadCutting::find_split_uv(OpenMesh::FaceHandle f, int val_u, int val_v, double& u, double& v)
{
	OpenMesh::Vec3d tar = { (double)val_u, (double)val_v, 0.0 };
	std::vector<OpenMesh::HalfedgeHandle> fh_h;
	std::vector<OpenMesh::Vec3d> v_ht;
	fh_h.reserve(3);
	v_ht.reserve(3);

	for (auto fh : new_para.fh_range(f))
	{
		fh_h.push_back(fh);
		v_ht.push_back(new_para.point(new_para.to_vertex_handle(fh)));
	}

	tar -= v_ht[2];
	v_ht[0] -= v_ht[2];
	v_ht[1] -= v_ht[2];

	u =  tar[0] * v_ht[1][1] - tar[1] * v_ht[1][0];
	v = -tar[0] * v_ht[0][1] + tar[1] * v_ht[0][0];

	double det_f = v_ht[0][0] * v_ht[1][1] - v_ht[0][1] * v_ht[1][0];
	u /= det_f;
	v /= det_f;

	if (u <= split_thres) return fh_h[2];
	else if (v <= split_thres) return fh_h[0];
	else if (u + v >= 1.0 - split_thres) return fh_h[1];
	else return OpenMesh::HalfedgeHandle();
}

void ParaQuadCutting::find_cut_path(OpenMesh::HalfedgeHandle& h_prev, int k, int dst_quad_v, int seg_id)
{
	const auto& quad_mesh = charts_decomposition->get_quad_mesh();
	MY_DOUBT(k < 1 || k > 4, "K Error");

	const int tag_step = 1;
	int tag_k = (get_tag(-new_para.calc_edge_vector(h_prev)) - k * tag_step) & 3;
	int gtag = (tag_k - tag_step) & 3;
	
	auto& uv0 = new_para.point(new_para.to_vertex_handle(h_prev));
	double tracing_val = std::round(uv0[gtag & 1]);
	double dst_val = std::round(quad_mesh.point(quad_mesh.vertex_handle(dst_quad_v))[tag_k & 1]);

	struct path_node
	{
		path_node(OpenMesh::HalfedgeHandle _h, double _t, double _u, double _v)
			:h(_h), t(_t), tar_p(_u, _v, 0.0) {};

		double t;
		OpenMesh::HalfedgeHandle h;
		OpenMesh::Vec3d tar_p;
	};
	bool spilt_end = false;
	std::vector<path_node> pass_by;
	OpenMesh::VertexHandle dst_v;

	std::set<int> dst_tris; 
	if (v_quad2para[dst_quad_v] != -1)
	{
		for (auto f_h : new_para.vf_range(new_para.vertex_handle(v_quad2para[dst_quad_v])))
		{
			dst_tris.insert(f_h.idx());
		}

		for (auto f_h : new_para.vf_range(new_para.to_vertex_handle(h_prev)))
		{
			if (dst_tris.count(f_h.idx()) > 0)
			{
				dst_v = new_para.vertex_handle(v_quad2para[dst_quad_v]);
				break;
			}
		}
	}

	auto h_iter = h_prev;

	while (true)
	{
		if (dst_v.is_valid()) break;

		if (k == 4 && h_iter == h_prev)
		{
			h_iter = new_para.opposite_halfedge_handle(new_para.next_halfedge_handle(h_iter));
			continue;
		}

		auto v1 = new_para.point(new_para.from_vertex_handle(h_iter)) - uv0;
		auto v2 = new_para.point(new_para.to_vertex_handle(new_para.next_halfedge_handle(h_iter))) - uv0;

		if ((!(gtag & 2) && (v1[gtag & 1] <= 0) && (v2[gtag & 1] >= 0)) || ((gtag & 2) && (v1[gtag & 1] >= 0) && (v2[gtag & 1] <= 0)))
		{
			//tracing
			auto tracing_h = new_para.prev_halfedge_handle(h_iter);

			MY_DOUBT((new_para.point(new_para.from_vertex_handle(tracing_h))[gtag & 1] - tracing_val) * (new_para.point(new_para.to_vertex_handle(tracing_h))[gtag & 1] - tracing_val) > 0,
				"Tracing Error 1, v " << new_para.to_vertex_handle(h_prev).idx() << ", k " << k);

			while (true)
			{
				tracing_h = new_para.opposite_halfedge_handle(tracing_h);
				auto tracing_f = new_para.face_handle(tracing_h);

				const auto& p1 = new_para.point(new_para.from_vertex_handle(tracing_h));
				const auto& p2 = new_para.point(new_para.to_vertex_handle(tracing_h));
				double c1 = p1[gtag & 1];
				double c2 = p2[gtag & 1];

				double split_t = (tracing_val - c1) / (c2 - c1);
				OpenMesh::Vec3d split_p;
				split_p[tag_k & 1] = p1[tag_k & 1] + split_t * (p2[tag_k & 1] - p1[tag_k & 1]);
				split_p[gtag & 1] = tracing_val;
				split_p[2] = 0.0;

				if (v_quad2para[dst_quad_v] == -1 && (((tag_k & 2) && split_p[tag_k & 1] <= dst_val) || (!(tag_k & 2) && split_p[tag_k & 1] >= dst_val)))
				{
					auto split_f = new_para.opposite_face_handle(tracing_h);
					OpenMesh::Vec3d dst_i;
					dst_i[tag_k & 1] = dst_val;
					dst_i[gtag & 1] = tracing_val;
					dst_i[2] = 0.0;

					double split_u, split_v;
					auto t_h = find_split_uv(split_f, dst_i[0], dst_i[1], split_u, split_v);

					if (!t_h.is_valid())
					{
						dst_v = split_face(split_f, split_u, split_v);
						v_quad2para[dst_quad_v] = dst_v.idx();
						new_para.point(dst_v) = dst_i;
					}
					else 
					{
						if (pass_by.empty() || t_h != pass_by.back().h)
						{
							double t_t = find_split_t(t_h, dst_i[0], dst_i[1]);
							pass_by.emplace_back(t_h, t_t, dst_i[0], dst_i[1]);
						}
						spilt_end = true;
					}

					break;
				}

				pass_by.emplace_back(tracing_h, split_t, split_p[0], split_p[1]);
				if (dst_tris.count(tracing_f.idx()) > 0)
				{
					dst_v = new_para.vertex_handle(v_quad2para[dst_quad_v]);
					break;
				}

				double c0 = new_para.point(new_para.to_vertex_handle(new_para.next_halfedge_handle(tracing_h)))[gtag & 1];

				MY_DOUBT((c1 - tracing_val) * (c2 - tracing_val) > 0, "Tracing Error 2, v " << new_para.to_vertex_handle(h_prev).idx() << ", k " << k);

				auto t1t0 = new_para.prev_halfedge_handle(tracing_h);
				auto t2t0 = new_para.next_halfedge_handle(tracing_h);
				if (c0 != tracing_val && c1 != tracing_val && c2 != tracing_val)
				{
					tracing_h = ((c0 < tracing_val) == (c1 < tracing_val)) ? t2t0 : t1t0;
				}
				else if (c0 != tracing_val && c1 == tracing_val && c2 != tracing_val)
				{
					tracing_h = ((c0 < tracing_val) == (c2 < tracing_val)) ? t1t0 : t2t0;
				}
				else if (c0 != tracing_val && c1 != tracing_val && c2 == tracing_val)
				{
					tracing_h = ((c0 < tracing_val) == (c1 < tracing_val)) ? t2t0 : t1t0;
				}
				else if (c0 == tracing_val && c1 == tracing_val && c2 != tracing_val)
				{
					tracing_h = t2t0;
				}
				else if (c0 == tracing_val && c1 != tracing_val && c2 == tracing_val)
				{
					tracing_h = t1t0;
				}
				else //FTT or TFF
				{
					double d_t1 = std::abs(new_para.point(new_para.from_vertex_handle(tracing_h))[(gtag & 1) ^ 1] - uv0[(gtag & 1) ^ 1]);
					double d_t2 = std::abs(new_para.point(new_para.to_vertex_handle(tracing_h))[(gtag & 1) ^ 1] - uv0[(gtag & 1) ^ 1]);

					tracing_h = (d_t1 <= d_t2) ? t2t0 : t1t0;
				}
			}

			break;
		}
 	
		h_iter = new_para.opposite_halfedge_handle(new_para.next_halfedge_handle(h_iter));
	}

	std::vector<OpenMesh::VertexHandle> path_v = { new_para.to_vertex_handle(h_prev) };

	std::set<int> pass_by_v;
	std::vector<int> split_status(pass_by.size());
	for (int i = 0; i < pass_by.size(); i++)
	{
		const auto& split_info = pass_by[i];
		split_status[i] = -1;

		if (split_info.t <= split_thres)
		{
			auto already_v = new_para.from_vertex_handle(split_info.h);
			if (!new_para.is_boundary(already_v))
			{
				split_status[i] = already_v.idx();
				pass_by_v.insert(split_status[i]);
			}
		}
		else if (split_info.t >= 1.0 - split_thres)
		{
			auto already_v = new_para.to_vertex_handle(split_info.h);
			if (!new_para.is_boundary(already_v))
			{
				split_status[i] = already_v.idx();
				pass_by_v.insert(split_status[i]);
			}
		}
	}

	for (int i = 0; i < pass_by.size(); i++)
	{
		if (split_status[i] == -1)
		{
			int v_from = new_para.from_vertex_handle(pass_by[i].h).idx();
			int v_to = new_para.to_vertex_handle(pass_by[i].h).idx();

			if (!(i == pass_by.size() - 1 && spilt_end) && (pass_by_v.count(v_from) == 1 || pass_by_v.count(v_to) == 1))
			{
				split_status[i] = -2;
			}
		}
	}

	for (int kkk = 0; kkk < 0; kkk++)
	{
		for (int i = 1; i < pass_by.size(); i++)
		{
			if (split_status[i] < 0) continue;
			auto v1 = new_para.vertex_handle(split_status[i]);

			int iter_i;

			iter_i = i - 1;
			while (iter_i >= 0 && (split_status[iter_i] == -2 || split_status[iter_i] == v1.idx())) iter_i--;
			if (iter_i < 0 || split_status[iter_i] < 0) continue;
			auto v0 = new_para.vertex_handle(split_status[iter_i]);

			iter_i = i + 1;
			while (iter_i < pass_by.size() && (split_status[iter_i] == -2 || split_status[iter_i] == v1.idx())) iter_i++;
			if (iter_i < pass_by.size() && split_status[iter_i] < 0) continue;
			if (iter_i == pass_by.size() && !dst_v.is_valid()) continue;
			auto v2 = (iter_i < pass_by.size()) ? new_para.vertex_handle(split_status[iter_i]) : dst_v;

			std::set<int> f_v012;
			auto h0 = new_para.find_halfedge(v0, v1);
			int aaa = new_para.face_handle(h0).idx();
			f_v012.insert(new_para.face_handle(h0).idx());
			aaa = new_para.opposite_face_handle(h0).idx();
			f_v012.insert(new_para.opposite_face_handle(h0).idx());

			bool flat_triangle = false;
			auto h1 = new_para.find_halfedge(v1, v2);
			if (f_v012.count(new_para.face_handle(h1).idx()) == 1) flat_triangle = true;
			if (f_v012.count(new_para.opposite_face_handle(h1).idx()) == 1) flat_triangle = true;

			if (!flat_triangle) continue;

			for (int j = i; j < iter_i; j++)
			{
				split_status[j] = -2;
				pass_by[j].h = new_para.find_halfedge(v0, v2);
				pass_by[j].t = 0.0;
			}
		}
	}

	for (int i = 0; i < pass_by.size(); i++)
	{
		const auto& split_info = pass_by[i];

		OpenMesh::VertexHandle new_v;
		if (split_status[i] >= 0)
		{
			new_v = new_para.vertex_handle(split_status[i]);
		}
		else if (split_status[i] == -1)
		{
			new_v = split_edge(split_info.h, split_info.t);
		}

		if (new_v.is_valid() && new_v != path_v.back())
		{
			old_point[new_v.idx()] = new_para.point(new_v);
			for (int k : {0, 1})
			{
				double val = split_info.tar_p[k];
				if (val - std::round(val) > 1e-4) continue;

				new_para.point(new_v)[k] = val;
			}
			
			path_v.push_back(new_v);
		}
	}

	if (spilt_end)
	{
		new_para.point(path_v.back()) = quad_mesh.point(quad_mesh.vertex_handle(dst_quad_v));
		v_quad2para[dst_quad_v] = path_v.back().idx();
	}

	if (dst_v.is_valid() && dst_v != path_v.back())
	{
		path_v.push_back(dst_v);
	}

	OpenMesh::HalfedgeHandle path_h;
	for (int i = 0; i < path_v.size() - 1; i++)
	{
		path_h = new_para.find_halfedge(path_v[i], path_v[i + 1]);
		new_para.property(e_segment, new_para.edge_handle(path_h)) = seg_id;
	}

	h_prev = path_h;
}

void ParaQuadCutting::save_quad_charts()
{
	if (!str_path.empty())
	{
		Mesh_doubleIO::save_mesh(charts_decomposition->get_quad_mesh(), (str_path + "/quad_comp/quad_comp.obj").c_str(), true);
		Mesh_doubleIO::save_mesh(charts_decomposition->get_packing_mesh(), (str_path + "/quad_comp/packing.obj").c_str());

		for (int i = 0; i < charts_decomposition->get_atlas().size(); i++)
		{
			std::stringstream comp_filename;
			comp_filename << str_path << "/quad_comp/quad_comp" << i << ".obj";
			Mesh_doubleIO::save_mesh(charts_decomposition->get_atlas()[i], comp_filename.str().c_str());
		}
	}
}

void ParaQuadCutting::save_tri_charts()
{
// 	Mesh_doubleIO::save_mesh(new_mesh, (str_path + "/tri_comp/mesh.obj").c_str());
// 	Mesh_doubleIO::save_mesh(new_para, (str_path + "/tri_comp/para.obj").c_str());

	if (!str_path.empty())
	{
		Mesh_doubleIO::save_mesh(new_para, (str_path + "/tri_comp/tri_comp.obj").c_str(), true);
		Mesh_doubleIO::save_mesh(new_mesh, (str_path + "/mesh_result.obj").c_str(), true);
		Mesh_doubleIO::save_mesh(mesh_unpack, (str_path + "/tri_comp/unpack.obj").c_str(), true);
//		Mesh_doubleIO::save_uv_mesh(new_mesh, (str_path + "/para_result.obj").c_str());

// 		for (int i = 0; i < para_atlas.size(); i++)
// 		{
// 			std::stringstream comp_filename;
// 			comp_filename << str_path << "/tri_comp/tri_comp" << i << ".obj";
// 			Mesh_doubleIO::save_mesh(para_atlas[i], comp_filename.str().c_str());
// 		}
	}

// 	std::ofstream result_log((str_path + "/result_log.txt").c_str());
// 	result_log << "#Charts : " << para_atlas.size() << std::endl;
// 	result_log << "#Charts : " << para_atlas.size() << std::endl;
// 	result_log.close();
}

void ParaQuadCutting::save_para(const char* filename)
{
	if (!str_path.empty())
	{
		Mesh_doubleIO::save_mesh(new_para, filename);

		std::ofstream para_log(str_path + "/para_log.txt");
		para_log << "Euler of Para : " << origin_para.n_vertices() + origin_para.n_faces() - origin_para.n_edges() << std::endl;
		para_log.close();
	}
}

void ParaQuadCutting::cut_on_boundary()
{
	const auto& quad_mesh = charts_decomposition->get_quad_mesh();

	const auto& boundary_h = charts_decomposition->get_boundary_h();
	const auto& bh_seg = charts_decomposition->get_bh_segment();
	for (int i = 0; i < boundary_h.size(); i++)
	{
		auto e_h = new_para.edge_handle(boundary_h[i] / 2);
		new_para.property(e_segment, e_h) = bh_seg[i];
	}

	v_quad2para.assign(quad_mesh.n_vertices(), -1); 
	const auto& vert_corner = charts_decomposition->get_vert_corner();
	for (int i = 0; i < vert_corner.size(); i++)
	{
		v_quad2para[i] = vert_corner[i];
	}

	const auto& boundary_cutting_nodes = charts_decomposition->get_boundary_cutting_nodes();
	for (const auto& bc_cutting : boundary_cutting_nodes)
	{
		if (v_quad2para[bc_cutting.first] != -1) continue;

		int seg_id = bc_cutting.second;
		int vert0 = charts_decomposition->get_segments()[seg_id].vert0;
		int vert1 = charts_decomposition->get_segments()[seg_id].vert1;
		int tag = charts_decomposition->get_segments()[seg_id].tag;

		double x_tar = std::round(quad_mesh.point(quad_mesh.vertex_handle(bc_cutting.first))[tag & 1]);

		OpenMesh::HalfedgeHandle h_iter;
		for (auto voh : new_para.voh_range(new_para.vertex_handle(vert0)))
		{
			if (new_para.property(e_segment, new_para.edge_handle(voh)) == seg_id)
			{
				h_iter = voh;
				break;
			}
		}

		MY_DOUBT(!new_para.is_boundary(h_iter), "Boundary Error");

		double split_t;
		while (true)
		{
			double x0 = new_para.point(new_para.from_vertex_handle(h_iter))[tag & 1];
			double x1 = new_para.point(new_para.to_vertex_handle(h_iter))[tag & 1];

			if ((!(tag & 2) && x0 <= x_tar && x1 > x_tar) || ((tag & 2) && x1 <= x_tar && x0 > x_tar))
			{
				split_t = (x_tar - x0) / (x1 - x0);
				break;
			}

			if (new_para.to_vertex_handle(h_iter).idx() == vert1)
			{
				std::cout << "Cant find h in seg." << std::endl;
				break;
			}
			h_iter = new_para.next_halfedge_handle(h_iter);
		}

		OpenMesh::VertexHandle v_new;
		if (split_t <= split_thres)
		{
			v_new = new_para.from_vertex_handle(h_iter);
		}
		else if (split_t >= 1.0 - split_thres)
		{
			v_new = new_para.to_vertex_handle(h_iter);
		}
		else
		{
			auto h_inside = new_para.opposite_halfedge_handle(h_iter);
			auto vf = new_para.from_vertex_handle(h_inside);
			auto vt = new_para.to_vertex_handle(h_inside);

			v_new = split_edge(h_inside, 1.0 - split_t);

			auto h0 = new_para.find_halfedge(vf, v_new);
			auto h0_mesh = get_para2mesh(h0);
			auto h1_mesh = new_mesh.opposite_halfedge_handle(h0_mesh);
			if (!new_mesh.is_boundary(h1_mesh))
			{
				auto h1 = get_mesh2para(h1_mesh);
				auto v_new2 = new_para.from_vertex_handle(h1);

				int seg_id = new_para.property(e_segment, new_para.edge_handle(h1));
				const auto& seg2 = charts_decomposition->get_segments()[seg_id];
				const auto& seg2_nodes = charts_decomposition->get_seg_nodes()[seg_id];

				const auto& quad_vertices = charts_decomposition->get_quad_vertices();

				MY_DOUBT(std::lround(new_para.point(v_new2)[(seg2.tag & 1) ^ 1]) != seg2.coord, "Seg2 Error1");
			}
		}

		v_quad2para[bc_cutting.first] = v_new.idx();
		new_para.point(v_new)[tag & 1] = x_tar;
		new_para.point(v_new)[(tag & 1) ^ 1] = charts_decomposition->get_segments()[seg_id].coord;
		new_para.point(v_new)[2] = 0.0;
	}
}

void ParaQuadCutting::cut_segments()
{
	cut_on_boundary();
	
	const auto& quad_mesh = charts_decomposition->get_quad_mesh();
	const auto& cut_path_h = charts_decomposition->get_cut_path_h();
	const auto& cut_path_v = charts_decomposition->get_cut_path_v();
	const auto& cut_path_c = charts_decomposition->get_cut_path_seg();
	int n_cut_path = cut_path_v.size();

	for (int i = 0; i < n_cut_path; i++)
	{
		int h_quad_prev = cut_path_h[i][0];
		int v_quad = cut_path_v[i][0].first;
	
		MY_DOUBT(v_quad2para[v_quad] == -1, "Cut Path Error");

		OpenMesh::HalfedgeHandle h_prev;
		OpenMesh::VertexHandle v0 = new_para.vertex_handle(v_quad2para[v_quad]);
		for (auto vih : new_para.vih_range(v0))
		{
			int vih_seg = new_para.property(e_segment, new_para.edge_handle(vih));
			auto h_quad_seg = charts_decomposition->get_cutting_edge_seg().find(h_quad_prev / 2);

			if (h_quad_seg != charts_decomposition->get_cutting_edge_seg().end() && vih_seg == h_quad_seg->second 
				&& get_tag(new_para.calc_edge_vector(vih)) == get_tag(quad_mesh.calc_edge_vector(quad_mesh.halfedge_handle(h_quad_prev))))
			{
				h_prev = vih;
				break;
			}
		}
		
		for (int j = 0; j < cut_path_v[i].size() - 1; j++)
		{
			find_cut_path(h_prev, cut_path_v[i][j].second, cut_path_v[i][j + 1].first, cut_path_c[i][j]);
		}
	}

	new_mesh.garbage_collection();
	new_para.garbage_collection();
}

void ParaQuadCutting::decomposition()
{
	para_atlas.clear();

	face_layer.assign(new_para.n_faces(), -1);

	std::vector<int> seed_triangles = get_seed_triangles();

	for (int i = 0; i < seed_triangles.size(); i++)
	{
		auto f_h = new_para.face_handle(seed_triangles[i]);
		int seed_test = face_layer[f_h.idx()];
		bool aaaa = face_layer[f_h.idx()] >= 0;
		MY_DOUBT(face_layer[f_h.idx()] >= 0, "SEED TRIANGLE ERROR");

		para_atlas.emplace_back();
		Mesh& comp = para_atlas.back();
		std::vector<int> vertex_index(new_para.n_vertices(), -1);

		std::queue<int> face_queue;
		face_queue.push(f_h.idx());

		while (!face_queue.empty())
		{
			int fid_cur = face_queue.front();
			face_queue.pop();
			if (face_layer[fid_cur] >= 0) continue;

			face_layer[fid_cur] = para_atlas.size() - 1;
			auto f_cur = new_para.face_handle(fid_cur);

			std::vector<OpenMesh::VertexHandle> v_fcur;
			for (auto fv : new_para.fv_range(f_cur))
			{
				if (vertex_index[fv.idx()] == -1)
				{
					vertex_index[fv.idx()] = comp.n_vertices();
					comp.add_vertex(new_para.point(fv));
				}
				v_fcur.emplace_back(vertex_index[fv.idx()]);
			}
			comp.add_face(v_fcur);

			for (auto fh : new_para.fh_range(f_cur))
			{
				if (new_para.property(e_segment, new_para.edge_handle(fh)) == -1 && face_layer[new_para.opposite_face_handle(fh).idx()] == -1)
				{
					face_queue.push(new_para.opposite_face_handle(fh).idx());
				}
			}
		}

		v_para2layer.emplace_back(std::move(vertex_index));
	}
 
 	std::cout << "Components " << para_atlas.size() << std::endl;
}

void ParaQuadCutting::get_packing_result()
{
	std::vector<int> nv_offset(para_atlas.size());

	nv_offset[0] = 0;
	for (int i = 0; i < para_atlas.size() - 1; i++)
	{
		nv_offset[i + 1] = nv_offset[i] + para_atlas[i].n_vertices();
	}

	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
	OpenMesh::HPropHandleT<int> hvt_index;

	new_mesh.add_property(mvt_list, "mvt_list");
	new_mesh.add_property(hvt_index, "hvt_index");

	new_mesh.property(mvt_list).resize(nv_offset.back() + para_atlas.back().n_vertices());

	const auto& chart_translation = charts_decomposition->get_chart_translation();
	const auto& chart_flipped = charts_decomposition->get_chart_flipped();
	for (int i = 0; i < para_atlas.size(); i++)
	{
		double offset_u = chart_translation[i][0];
		double offset_v = chart_translation[i][1];
		for (int j = 0; j < para_atlas[i].n_vertices(); j++)
		{
			const auto& point_v = para_atlas[i].point(para_atlas[i].vertex_handle(j));
			if (chart_flipped[i])
			{
				new_mesh.property(mvt_list)[nv_offset[i] + j][0] = -point_v[1] + offset_u;
				new_mesh.property(mvt_list)[nv_offset[i] + j][1] = point_v[0] + offset_v;
			}
			else
			{
				new_mesh.property(mvt_list)[nv_offset[i] + j][0] = point_v[0] + offset_u;
				new_mesh.property(mvt_list)[nv_offset[i] + j][1] = point_v[1] + offset_v;
			}
		}
	}

	for (int i = 0; i < new_mesh.n_faces(); i++)
	{
		auto f_h = new_mesh.face_handle(i);
		int layer_id = face_layer[i];

		for (auto fh_h = new_mesh.cfh_begin(f_h); fh_h != new_mesh.cfh_end(f_h); fh_h++)
		{
			auto h_para = get_mesh2para(*fh_h);
			if (!h_para.is_valid())
			{
				std::cout << new_mesh.from_vertex_handle(*fh_h).idx() << " " << new_mesh.to_vertex_handle(*fh_h).idx() << " " << i << std::endl;
				continue;
			}

			new_mesh.property(hvt_index, *fh_h) = v_para2layer[layer_id][new_para.to_vertex_handle(h_para).idx()] + nv_offset[layer_id];
		}
	}

	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list_para;
	OpenMesh::HPropHandleT<int> hvt_index_para;

	new_para.add_property(mvt_list_para, "mvt_list");
	new_para.add_property(hvt_index_para, "hvt_index");

	new_para.property(mvt_list_para) = new_mesh.property(mvt_list);
	for (int i = 0; i < new_para.n_faces(); i++)
	{
		auto f_h = new_para.face_handle(i);
		int layer_id = face_layer[i];

		for (auto fh_h : new_para.fh_range(f_h))
		{
			new_para.property(hvt_index_para, fh_h) = v_para2layer[layer_id][new_para.to_vertex_handle(fh_h).idx()] + nv_offset[layer_id];
		}
	}

	Mesh_doubleIO::copy_mesh(new_mesh, mesh_unpack);
	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list_unpack;
	OpenMesh::HPropHandleT<int> hvt_index_unpack;

	mesh_unpack.add_property(mvt_list_unpack, "mvt_list");
	mesh_unpack.add_property(hvt_index_unpack, "hvt_index");

	mesh_unpack.property(mvt_list_unpack).resize(nv_offset.back() + para_atlas.back().n_vertices());

	for (int i = 0; i < para_atlas.size(); i++)
	{
		for (int j = 0; j < para_atlas[i].n_vertices(); j++)
		{
			const auto& point_v = para_atlas[i].point(para_atlas[i].vertex_handle(j));
			mesh_unpack.property(mvt_list_unpack)[nv_offset[i] + j][0] = point_v[0];
			mesh_unpack.property(mvt_list_unpack)[nv_offset[i] + j][1] = point_v[1];
		}
	}

	for (auto h_h : mesh_unpack.halfedges())
	{
		auto v0 = mesh_unpack.from_vertex_handle(h_h);
		auto v1 = mesh_unpack.to_vertex_handle(h_h);

		auto h_mesh = new_mesh.find_halfedge(v0, v1);
		mesh_unpack.property(hvt_index_unpack, h_h) = new_mesh.property(hvt_index, h_mesh);
	}
}

std::vector<int> ParaQuadCutting::get_seed_triangles()
{
	const auto& quad_mesh = charts_decomposition->get_quad_mesh();
	const auto& quad_atlas = charts_decomposition->get_atlas();
	const auto& v_chart2quad = charts_decomposition->get_v_chart2quad();

	std::vector<int> seed_triangles(quad_atlas.size(), -1);
	for (int i = 0; i < quad_atlas.size(); i++)
	{
		const auto& quad_chart = quad_atlas[i];

		OpenMesh::VertexHandle corner;
		for (auto v_h : quad_chart.vertices())
		{
			if (quad_chart.is_boundary(v_h) && quad_chart.valence(v_h) != 3)
			{
				corner = v_h;
				break;
			}
		}

		OpenMesh::HalfedgeHandle inner_out_h;
		for (auto vih : quad_chart.vih_range(corner))
		{
			if (quad_chart.is_boundary(vih))
			{
				inner_out_h = quad_chart.opposite_halfedge_handle(vih);
				break;
			}
		}

		int tag0 = get_tag(quad_chart.calc_edge_vector(inner_out_h));
		int tag1 = get_tag(-quad_chart.calc_edge_vector(quad_chart.prev_halfedge_handle(inner_out_h)));

		auto vf_quad = quad_mesh.vertex_handle(v_chart2quad[i][quad_chart.from_vertex_handle(inner_out_h).idx()]);
		auto vt_quad = quad_mesh.vertex_handle(v_chart2quad[i][quad_chart.to_vertex_handle(inner_out_h).idx()]);
		auto h_quad = quad_mesh.find_halfedge(vf_quad, vt_quad);
		int seg_id = charts_decomposition->get_cutting_edge_seg().at(h_quad.idx() / 2);

		OpenMesh::HalfedgeHandle h_tri;
		OpenMesh::VertexHandle v_tri = new_para.vertex_handle(v_quad2para[v_chart2quad[i][corner.idx()]]);
		MY_DOUBT(!v_tri.is_valid(), "Cut Path Error 2");

		for (auto voh : new_para.voh_range(v_tri))
		{
			if (get_tag(new_para.calc_edge_vector(voh)) == tag0 && new_para.property(e_segment, new_para.edge_handle(voh)) == seg_id)
			{
				h_tri = voh;
				break;
			}
		}

		OpenMesh::Vec3d h_vec;

		if (!new_para.is_boundary(h_tri))
		{
			h_vec = -new_para.calc_edge_vector(new_para.prev_halfedge_handle(h_tri));
			h_vec[tag0 & 1] = 0.0;
			if (get_tag(h_vec) == tag1) seed_triangles[i] = new_para.face_handle(h_tri).idx();
		}

		if (!new_para.is_boundary(new_para.opposite_halfedge_handle(h_tri)))
		{
			h_vec = new_para.calc_edge_vector(new_para.next_halfedge_handle(new_para.opposite_halfedge_handle(h_tri)));
			h_vec[tag0 & 1] = 0.0;
			if (get_tag(h_vec) == tag1) seed_triangles[i] = new_para.opposite_face_handle(h_tri).idx();
		}
	}

	return seed_triangles;
}

void ParaQuadCutting::split_poly()
{
	std::string outstr = str_file.substr(0, str_file.find_last_of('.')) + "_axis.obj";
	OpenMesh::IO::write_mesh(new_para, outstr);
	std::string outstr_temp = outstr;
	bool is_replace_ok = 1;
	while (is_replace_ok)
	{
		outstr_temp = outstr_temp.replace(outstr_temp.find("/"), 1, "\\");
		if (strstr(outstr_temp.c_str(), "/") == NULL)
		{
			is_replace_ok = 0;
		}
	}
	//std::string new_file = str_file.substr(0, str_file.find_last_of('.')) + "_result";
	//replace(new_file.begin(),new_file.end(), "/", "\\");
	/*std::cout <<"new_file: " <<new_file << std::endl;
	is_replace_ok = 1;
	while (is_replace_ok)
	{
		new_file = new_file.replace(new_file.find("/"), 1, "\\");
		if (strstr(new_file.c_str(), "/") == NULL)
		{
			is_replace_ok = 0;
		}
	}
	std::cout << "new_file: " << new_file << std::endl;
	std::string file_name = new_file.substr(new_file.find_last_of("\\")+1 , new_file.length());
	std::cout << "new_file: " << file_name << std::endl;
	std::string cmd = "mkdir " + new_file;
	system(cmd.data());
	std::cout << "npppppp" << std::endl;*/
	std::string str_path_temp = str_path;
	//replace(str_path_temp.begin(), str_path_temp.end(), "/", "\\");
	is_replace_ok = 1;
	while (is_replace_ok)
	{
		str_path_temp = str_path_temp.replace(str_path_temp.find("/"), 1, "\\");
		if (strstr(str_path_temp.c_str(), "/") == NULL)
		{
			is_replace_ok = 0;
		}
	}
	std::string tmp = str_path_temp + "\\SplitPoly.exe " + outstr_temp + " " + str_path_temp + "\\ 50";    /*"\\SplitPoly.exe D1_01164_input_result\\D1_01164_input_result.obj D1_01164_input_result\\ 30";*/
	system(tmp.data());
}

void ParaQuadCutting::read_cut()
{
	std::string file_name0 = str_file.substr(str_file.find_last_of("/") + 1);
	std::string file_name = file_name0.substr(0, file_name0.find_last_of('.'));
	std::string instr = str_file.substr(0, str_file.find_last_of('.')) + "_axis/" + file_name + "_axis.txt";
	std::vector<double> data_temp;
	std::ifstream data(instr);
	double d;
	while (data >> d)
	{
		data_temp.push_back(d);
	}
	data.close();
	for (int i = 0; i < data_temp.size() / 4; i++)
	{
		std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> cut_segments_temp;
		OpenMesh::Vec3d segment1_temp(0, 0, 0);
		OpenMesh::Vec3d segment2_temp(0, 0, 0);
		segment1_temp[0] = data_temp[4 * i];
		segment1_temp[1] = data_temp[4 * i + 1];
		segment2_temp[0] = data_temp[4 * i + 2];
		segment2_temp[1] = data_temp[4 * i + 3];
		cut_segments_temp.first = segment1_temp;
		cut_segments_temp.second = segment2_temp;
		my_cut_segments.push_back(cut_segments_temp);
	}
	new_para.add_property(e_cut);
	for (auto e_h : new_mesh.edges())
	{
		new_para.property(e_cut, e_h) = 0;
	}
}

void ParaQuadCutting::read_cut2()
{
	std::string file_name0 = str_file.substr(str_file.find_last_of("/") + 1);
	std::string file_name = file_name0.substr(0, file_name0.find_last_of('.'));
	std::string instr = str_file.substr(0, str_file.find_last_of('.')) + "_axis/" + file_name + "_axis.txt";

	std::vector<OpenMesh::Vec3d> segment_point;
	std::ifstream ifs(instr.c_str());

	int boundVNum;
	ifs >> boundVNum;
	for (int i = 0; i < boundVNum; ++i)
	{
		int v0, v1;
		double lambda;
		ifs >> v0 >> v1 >> lambda;

		OpenMesh::Vec3d point_temp1 = new_para.point(new_para.vertex_handle(v0));
		OpenMesh::Vec3d point_temp2 = new_para.point(new_para.vertex_handle(v1));
		OpenMesh::Vec3d point_temp = point_temp1 + lambda * (point_temp2 - point_temp1);
		segment_point.push_back(point_temp);
	}

	int innerVNum;
	ifs >> innerVNum;
	for (int i = 0; i < innerVNum; ++i)
	{
		int v0, v1;
		double lambda;
		ifs >> v0 >> v1 >> lambda;

		OpenMesh::Vec3d point_temp1 = segment_point[v0];
		OpenMesh::Vec3d point_temp2 = segment_point[v1];
		OpenMesh::Vec3d point_temp = point_temp1 + lambda * (point_temp2 - point_temp1);
		segment_point.push_back(point_temp);
	}

	int segmentNum;
	ifs >> segmentNum;
	for (int i = 0; i < segmentNum; ++i)
	{
		int v0, v1;
		ifs >> v0 >> v1;

		std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> cut_segments_temp;
		cut_segments_temp.first = segment_point[v0];
		cut_segments_temp.second = segment_point[v1];
		my_cut_segments.push_back(cut_segments_temp);
	}
	new_para.add_property(e_cut,"e_cut");
	for (auto e_h : new_mesh.edges())
	{
		new_para.property(e_cut, e_h) = 0;
	}
}

void ParaQuadCutting::read_cut1()
{
	std::string file_name0 = str_file.substr(str_file.find_last_of("/") + 1);
	std::string file_name = file_name0.substr(0, file_name0.find_last_of('.'));
	std::string instr = str_file.substr(0, str_file.find_last_of('.')) + "_axis/" + file_name + "_axis.txt";

	FILE *fid = fopen(instr.c_str(), "r");

	int LINE_MAX = 1024;
	char line[1024];

	std::vector<OpenMesh::Vec3d> segment_point;
	int type_temp = 0;
	while (fgets(line, LINE_MAX, fid) != NULL)
	{
		if (strlen(line) <= 2)
		{
			type_temp++;
			continue;
		}
		if (type_temp == 1)
		{
			int id1, id2;
			double lambda;
			int count = sscanf(line, "%ld %ld %lf\n", &id1, &id2, &lambda);
			if (count == 3)
			{
				OpenMesh::Vec3d point_temp1 = new_para.point(new_para.vertex_handle(id1));
				OpenMesh::Vec3d point_temp2 = new_para.point(new_para.vertex_handle(id2));
				OpenMesh::Vec3d point_temp = point_temp1 + lambda * (point_temp2 - point_temp1);
				segment_point.push_back(point_temp);
			}
			else
			{
				//std::cout << "read Wrong1" << std::endl;
				return;
			}
		}
		if (type_temp == 2)
		{
			int id1, id2;
			double lambda;
			int count = sscanf(line, "%ld %ld %lf\n", &id1, &id2, &lambda);
			if (count == 3)
			{
				OpenMesh::Vec3d point_temp1 = segment_point[id1];
				OpenMesh::Vec3d point_temp2 = segment_point[id2];
				OpenMesh::Vec3d point_temp = point_temp1 + lambda * (point_temp2 - point_temp1);
				segment_point.push_back(point_temp);
			}
			else
			{
				//std::cout << "read Wrong2" << std::endl;
				return;
			}
		}
		if (type_temp == 3)
		{
			int id1, id2;
			int count = sscanf(line, "%ld %ld\n", &id1, &id2);
			if (count == 2)
			{
				std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> cut_segments_temp;
				OpenMesh::Vec3d point_temp1 = segment_point[id1];
				OpenMesh::Vec3d point_temp2 = segment_point[id2];
				cut_segments_temp.first = point_temp1;
				cut_segments_temp.second = point_temp2;
				my_cut_segments.push_back(cut_segments_temp);
			}
			else
			{
				//std::cout << "read Wrong3" << std::endl;
				return;
			}
		}
	}
	fclose(fid);
	new_para.add_property(e_cut);
	for (auto e_h : new_mesh.edges())
	{
		new_para.property(e_cut, e_h) = 0;
	}
}

bool ParaQuadCutting::is_point_on_edge(OpenMesh::Vec3d p_, OpenMesh::EdgeHandle edge_, Mesh &mesh_,double average_length)
{
	OpenMesh::HalfedgeHandle he = mesh_.halfedge_handle(edge_, 0);
	OpenMesh::VertexHandle fromv = mesh_.from_vertex_handle(he);
	OpenMesh::VertexHandle tov = mesh_.to_vertex_handle(he);
	OpenMesh::Vec3d q1 = mesh_.point(fromv);
	OpenMesh::Vec3d q2 = mesh_.point(tov);
	OpenMesh::Vec3d vec1 = (p_ - q1).normalize();
	OpenMesh::Vec3d vec2 = (q2 - q1).normalize();
	if (std::fabs((vec1%vec2)[2]) < 1e-3*average_length
		&& std::min(q1[0], q2[0])- 1e-3*average_length <= p_[0] && p_[0] <= std::max(q1[0], q2[0])+ 1e-3*average_length
		&& std::min(q1[1], q2[1])- 1e-3*average_length <= p_[1] && p_[1] <= std::max(q1[1], q2[1])+ 1e-3*average_length)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool ParaQuadCutting::is_point_in_face(OpenMesh::Vec3d p_, OpenMesh::FaceHandle face_, Mesh &mesh_)
{
	OpenMesh::Vec3d tar = p_;
	std::vector<OpenMesh::HalfedgeHandle> fh_h;
	std::vector<OpenMesh::Vec3d> v_ht;
	fh_h.reserve(3);
	v_ht.reserve(3);

	for (auto fh : mesh_.fh_range(face_))
	{
		fh_h.push_back(fh);
		v_ht.push_back(mesh_.point(mesh_.to_vertex_handle(fh)));
	}

	tar -= v_ht[2];
	v_ht[0] -= v_ht[2];
	v_ht[1] -= v_ht[2];

	double u = tar[0] * v_ht[1][1] - tar[1] * v_ht[1][0];
	double v = -tar[0] * v_ht[0][1] + tar[1] * v_ht[0][0];

	double det_f = v_ht[0][0] * v_ht[1][1] - v_ht[0][1] * v_ht[1][0];
	u /= det_f;
	v /= det_f;

	if (u < -1e-3 || u > 1+1e-3)
	{
		return false;
	}

	if (v < -1e-3 || v > 1 + 1e-3)
	{
		return false;
	}

	return u + v < 1 + 1e-3;

	/*std::vector<OpenMesh::Vec3d> face_vertex;
	for (auto fv_h : mesh_.fv_range(face_))
	{
		face_vertex.push_back(mesh_.point(fv_h));
	}
	OpenMesh::Vec3d A = face_vertex[0];
	OpenMesh::Vec3d B = face_vertex[1];
	OpenMesh::Vec3d C = face_vertex[2];
	OpenMesh::Vec3d vec0 = C - A;
	OpenMesh::Vec3d vec1 = B - A;
	OpenMesh::Vec3d vec2 = p_ - A;

	double dot00 = vec0 | vec0;
	double dot01 = vec0 | vec1;
	double dot02 = vec0 | vec2;
	double dot11 = vec1 | vec1;
	double dot12 = vec1 | vec2;

	double inverDeno = 1 / (dot00 * dot11 - dot01 * dot01);

	double u = (dot11*dot02 - dot01 * dot12)*inverDeno;
	if (u <= 0 || u >= 1)
	{
		return false;
	}

	double v = (dot00 * dot12 - dot01 * dot02)*inverDeno;
	if (v <= 0 || v >= 1)
	{
		return false;
	}

	return u + v < 1;*/
}

OpenMesh::HalfedgeHandle ParaQuadCutting::find_split_uv(OpenMesh::Vec3d p_, OpenMesh::FaceHandle face_, Mesh &mesh_, double &u, double &v)
{
	OpenMesh::Vec3d tar = p_;
	std::vector<OpenMesh::HalfedgeHandle> fh_h;
	std::vector<OpenMesh::Vec3d> v_ht;
	fh_h.reserve(3);
	v_ht.reserve(3);

	for (auto fh : mesh_.fh_range(face_))
	{
		fh_h.push_back(fh);
		v_ht.push_back(mesh_.point(mesh_.to_vertex_handle(fh)));
	}

	tar -= v_ht[2];
	v_ht[0] -= v_ht[2];
	v_ht[1] -= v_ht[2];

	u = tar[0] * v_ht[1][1] - tar[1] * v_ht[1][0];
	v = -tar[0] * v_ht[0][1] + tar[1] * v_ht[0][0];

	double det_f = v_ht[0][0] * v_ht[1][1] - v_ht[0][1] * v_ht[1][0];
	u /= det_f;
	v /= det_f;

	if (u <= split_thres) return fh_h[2];
	else if (v <= split_thres) return fh_h[0];
	else if (u + v >= 1.0 - split_thres) return fh_h[1];
	else return OpenMesh::HalfedgeHandle();
}

double ParaQuadCutting::find_split_t(OpenMesh::Vec3d p_, OpenMesh::HalfedgeHandle he_)
{
	const auto& uv0 = new_para.point(new_para.from_vertex_handle(he_));
	const auto& uv1 = new_para.point(new_para.to_vertex_handle(he_));

	auto v_e = uv1 - uv0;
	auto v_p = p_ - uv0;

	return OpenMesh::dot(v_e, v_p) / v_e.sqrnorm();
}

void ParaQuadCutting::cut_segment(std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> &my_cut_segments_)
{
	double neighbor_length = 0;
	for (auto e_h : new_para.edges())
	{
		neighbor_length += new_para.calc_edge_length(e_h);
	}
	neighbor_length = neighbor_length / new_para.n_edges();
	for (int i = 0; i < my_cut_segments_.size(); i++)
	{
		OpenMesh::Vec3d segment1 = my_cut_segments_[i].first;
		OpenMesh::Vec3d segment2 = my_cut_segments_[i].second;
		process_point_in_face(segment1, my_cut_segments_);
		process_point_in_face(segment2, my_cut_segments_);
		
		segment1 = my_cut_segments_[i].first;
		segment2 = my_cut_segments_[i].second;
		OpenMesh::EdgeHandle edge_init;
		bool can_find_edge = 0;
		for (auto e_h : new_para.edges())
		{
			bool edge_init_flag = is_point_on_edge(segment1, e_h, new_para, neighbor_length);
			if (edge_init_flag)
			{
				can_find_edge = 1;
				edge_init = e_h;
				break;
			}
		}
		if (!can_find_edge)
		{
			//std::cout << "can not find edge" << std::endl;
		}
		//std::cout << "find edge ok" << std::endl;
		if ((segment1 - segment2).norm() < 1e-6)
		{
			continue;
		}

		OpenMesh::HalfedgeHandle he_temp = new_para.halfedge_handle(edge_init, 0);
		OpenMesh::Vec3d fpoint = new_para.point(new_para.from_vertex_handle(he_temp));
		OpenMesh::Vec3d tpoint = new_para.point(new_para.to_vertex_handle(he_temp));
		double t_temp = (segment1 - fpoint).norm() / (tpoint - fpoint).norm();
		OpenMesh::VertexHandle new_vertex_pre_temp;
		if (t_temp < split_thres)
		{
			new_vertex_pre_temp = new_para.from_vertex_handle(he_temp);
			update_cut_segments(segment1, fpoint, my_cut_segments_);
		}
		else if (t_temp > 1.0 - split_thres)
		{
			new_vertex_pre_temp = new_para.to_vertex_handle(he_temp);
			update_cut_segments(segment1, tpoint, my_cut_segments_);
		}
		else
		{
			new_vertex_pre_temp = split_edge(he_temp, t_temp);
		}
		for (auto vv_h : new_para.vv_range(new_vertex_pre_temp))
		{
			OpenMesh::HalfedgeHandle new_he_temp = new_para.find_halfedge(new_vertex_pre_temp, vv_h);
			OpenMesh::EdgeHandle  new_edge_temp = new_para.edge_handle(new_he_temp.idx() / 2);
			if (new_para.property(e_cut, new_edge_temp) != 1)
			{
				new_para.property(e_cut, new_edge_temp) = 0;
			}
		}
		OpenMesh::VertexHandle new_vertex_temp = new_vertex_pre_temp;
		he_temp = find_next_point(new_vertex_temp, segment1, segment2, t_temp);

		//std::cout << "init ok" << std::endl;

		OpenMesh::Vec3d point_temp = segment1;
		bool cut_over = 0;
		int count_temp = 0;
		if (!new_para.is_valid_handle(he_temp))
		{
			cut_over = 1;
		}
		while (!cut_over)
		{
			count_temp++;
			if (t_temp < split_thres)
			{
				new_vertex_temp = new_para.from_vertex_handle(he_temp);
				mark_e_cut(new_vertex_pre_temp, new_vertex_temp);
			}
			else if (t_temp > 1.0 - split_thres)
			{
				new_vertex_temp = new_para.to_vertex_handle(he_temp);
				mark_e_cut(new_vertex_pre_temp, new_vertex_temp);
			}
			else
			{
				new_vertex_temp = split_edge(he_temp, t_temp);
				mark_e_cut(new_vertex_pre_temp, new_vertex_temp);
			}
			OpenMesh::Vec3d new_point_temp = new_para.point(new_vertex_temp);
			he_temp = find_next_point(new_vertex_temp, segment1, segment2, t_temp);
			if (!new_para.is_valid_handle(he_temp))
			{
				//std::cout << "wrong over" << std::endl;
				break;
			}
			fpoint = new_para.point(new_para.from_vertex_handle(he_temp));
			tpoint = new_para.point(new_para.to_vertex_handle(he_temp));
			point_temp = (1 - t_temp)*fpoint + t_temp * tpoint;
			if (new_para.is_boundary(new_vertex_temp))
			{
				//std::cout << "over0" << std::endl;
				cut_over = 1;
				//if (count_temp < 5)
				//{
				//cut_over = 0;
				//}
			}
			if ((segment1 - point_temp).norm() <= (segment1 - new_point_temp).norm())
			{
				//std::cout << "over1" << std::endl;
				cut_over = 1;
				//if (count_temp < 5)
				//{
				//cut_over = 0;
			//	}
			}
			if ((segment1 - point_temp).norm() > 1.01*(segment1 - segment2).norm())
			{
				//std::cout << "over2" << std::endl;
				cut_over = 1;
			}
			if ((new_point_temp - segment2).norm() < 1e-3*neighbor_length)
			{
				//std::cout << "over3" << std::endl;
				cut_over = 1;
			}
			new_vertex_pre_temp = new_vertex_temp;
		}
	}
	for (auto e_h : new_para.edges())
	{
		if (new_para.is_boundary(e_h))
		{
			new_para.property(e_cut, e_h) = 1;
		}
	}

	vertex_on_cut.clear();
	vertex_on_cut.resize(new_para.n_vertices());
	for (auto e_h : new_para.edges())
	{
		if (new_para.property(e_cut, e_h) == 1)
		{
			OpenMesh::HalfedgeHandle he_temp = new_para.halfedge_handle(e_h, 0);
			OpenMesh::VertexHandle fromv = new_para.from_vertex_handle(he_temp);
			OpenMesh::VertexHandle tov = new_para.to_vertex_handle(he_temp);

			vertex_on_cut[fromv.idx()] = 1;
			vertex_on_cut[tov.idx()] = 1;
		}	
	}
	for (auto e_h : new_para.edges())
	{
		OpenMesh::HalfedgeHandle he_temp = new_para.halfedge_handle(e_h, 0);
		OpenMesh::VertexHandle fromv = new_para.from_vertex_handle(he_temp);
		OpenMesh::VertexHandle tov = new_para.to_vertex_handle(he_temp);
		if (new_para.property(e_cut, e_h) != 1 && vertex_on_cut[fromv.idx()] == 1 && vertex_on_cut[tov.idx()] == 1)
		{
			double t_temp = 0.5;
			OpenMesh::VertexHandle vertex_temp = split_edge(he_temp, t_temp);
			for (auto vv_h : new_para.vv_range(vertex_temp))
			{
				OpenMesh::HalfedgeHandle new_he_temp = new_para.find_halfedge(vertex_temp, vv_h);
				OpenMesh::EdgeHandle  new_edge_temp = new_para.edge_handle(new_he_temp.idx() / 2);
				if (new_para.property(e_cut, new_edge_temp) != 1)
				{
					new_para.property(e_cut, new_edge_temp) = 0;
				}
			}
		}
	}

	new_para.garbage_collection();
	new_mesh.garbage_collection();	
}

void ParaQuadCutting::correct_segment(OpenMesh::Vec3d &segment1_, OpenMesh::Vec3d &segment2_, std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> &my_cut_segments_)
{
	bool segment1_flag = 0;
	bool segment2_flag = 0;
	int segment1_num = 0;
	int segment2_num = 0;
	for (int i = 0; i < my_cut_segments_.size(); i++)
	{
		if (my_cut_segments_[i].first == segment1_)
		{
			segment1_num++;
		}
		if (my_cut_segments_[i].first == segment2_)
		{
			segment2_num++;
		}
		if (my_cut_segments_[i].second == segment1_)
		{
			segment1_num++;
		}
		if (my_cut_segments_[i].second == segment2_)
		{
			segment2_num++;
		}
	}
	if (segment1_num > 1)
	{
		segment1_flag = 1;
	}
	if (segment2_num > 1)
	{
		segment2_flag = 1;
	}

	double neighbor_length = 0;
	for (auto e_h : new_para.edges())
	{
		neighbor_length += new_para.calc_edge_length(e_h);
	}
	neighbor_length = neighbor_length / new_para.n_edges();
	for (auto v_h : new_para.vertices())
	{
		if (vertex_on_cut.size() > 0)
		{
			if (vertex_on_cut[v_h.idx()] == 1)
			{
				OpenMesh::Vec3d bpoint = new_para.point(v_h);
				if ((segment1_ - bpoint).norm() < 1e-3*neighbor_length)
				{
					segment1_flag = 1;
					segment1_ = bpoint;
				}
				if ((segment2_ - bpoint).norm() < 1e-3*neighbor_length)
				{
					segment2_flag = 1;
					segment2_ = bpoint;
				}
			}
		}
		else
		{
			if (new_para.is_boundary(v_h))
			{
				OpenMesh::Vec3d bpoint = new_para.point(v_h);
				if ((segment1_ - bpoint).norm() < 1e-3*neighbor_length)
				{
					segment1_flag = 1;
					segment1_ = bpoint;
				}
				if ((segment2_ - bpoint).norm() < 1e-3*neighbor_length)
				{
					segment2_flag = 1;
					segment2_ = bpoint;
				}
			}
		}
	}
	std::vector<OpenMesh::Vec3d> intersect_point;
	for (auto he_h : new_para.halfedges())
	{
		if (new_para.is_boundary(he_h))
		{
			OpenMesh::Vec3d fpoint = new_para.point(new_para.from_vertex_handle(he_h));
			OpenMesh::Vec3d tpoint = new_para.point(new_para.to_vertex_handle(he_h));
			OpenMesh::Vec3d vec0 = segment2_ - segment1_;
			OpenMesh::Vec3d vec1 = fpoint - segment1_;
			OpenMesh::Vec3d vec2 = tpoint - segment1_;
			bool is_intersect_temp = 0;
			if ((vec0%vec1)[2] * (vec0%vec2)[2] <= 0)
			{
				is_intersect_temp = 1;
			}
			if (is_intersect_temp)
			{
				OpenMesh::Vec3d segment1_temp = segment1_ /*+ 0.1*(segment1_ - segment2_)*/;
				OpenMesh::Vec3d segment2_temp = segment2_ /*+ 0.1*(segment2_ - segment1_)*/;
				double x1 = segment1_temp[0]; double y1 = segment1_temp[1];
				double x2 = segment2_temp[0]; double y2 = segment2_temp[1];
				double x3 = fpoint[0]; double y3 = fpoint[1];
				double x4 = tpoint[0]; double y4 = tpoint[1];
				double t = ((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)) / ((x2 - x1)*(y3 - y4) - (x3 - x4)*(y2 - y1));
				if (t > 1)
				{
					t = 1;
				}
				if (t < 0)
				{
					t = 0;
				}
				OpenMesh::Vec3d new_point_temp(x3 + t * (x4 - x3), y3 + t * (y4 - y3), 0);
				intersect_point.push_back(new_point_temp);
			}
		}
	}
	if (!segment1_flag)
	{
		std::vector<double> point_d;
		for (int i = 0; i < intersect_point.size(); i++)
		{
			OpenMesh::Vec3d segment1_temp = segment1_ /*+ 0.1*(segment2_ - segment1_)*/;
			double d_temp = (intersect_point[i] - segment1_temp).norm();
			point_d.push_back(d_temp);
		}
		auto smallest = std::min_element(point_d.begin(), point_d.end());
		int smallest_p = std::distance(point_d.begin(), smallest);
		segment1_ = intersect_point[smallest_p];
	}
	if (!segment2_flag)
	{
		std::vector<double> point_d;
		for (int i = 0; i < intersect_point.size(); i++)
		{
			OpenMesh::Vec3d segment2_temp = segment2_/* + 0.1*(segment1_ - segment2_)*/;
			double d_temp = (intersect_point[i] - segment2_temp).norm();
			point_d.push_back(d_temp);
		}
		auto smallest = std::min_element(point_d.begin(), point_d.end());
		int smallest_p = std::distance(point_d.begin(), smallest);
		segment2_ = intersect_point[smallest_p];
	}
}

bool ParaQuadCutting::is_intersect(OpenMesh::Vec3d p1_, OpenMesh::Vec3d p2_, OpenMesh::Vec3d p3_, OpenMesh::Vec3d p4_)
{
	p2_ = p2_ + (p2_ - p1_) / 100;
	bool flag1, flag2;
	OpenMesh::Vec3d vec12 = p2_ - p1_;
	OpenMesh::Vec3d vec13 = p3_ - p1_;
	OpenMesh::Vec3d vec14 = p4_ - p1_;
	if (((vec13%vec12)[2]*(vec14%vec12)[2]) <= 0)
	{
		flag1 = 1;
	}
	else
	{
		flag1 = 0;
	}

	OpenMesh::Vec3d vec34 = p4_ - p3_;
	OpenMesh::Vec3d vec31 = p1_ - p3_;
	OpenMesh::Vec3d vec32 = p2_ - p3_;
	if (((vec31%vec34)[2]*(vec32%vec34)[2]) <= 0)
	{
		flag2 = 1;
	}
	else
	{
		flag2 = 0;
	}

	if (flag1 && flag2)
	{
		return true;
	}
	else
	{
		return false;
	}
}

OpenMesh::Vec3d ParaQuadCutting::compute_intersect_point(OpenMesh::Vec3d p1_, OpenMesh::Vec3d p2_, OpenMesh::Vec3d p3_, OpenMesh::Vec3d p4_)
{
	OpenMesh::Vec3d intersect_point(0, 0, 0);
	double b1 = (p2_[1] - p1_[1])*p1_[0] + (p1_[0] - p2_[0])*p1_[1];
	double b2 = (p4_[1] - p3_[1])*p3_[0] + (p3_[0] - p4_[0])*p3_[1];
	double D = (p2_[0] - p1_[0])*(p4_[1] - p3_[1]) - (p4_[0] - p3_[0])*(p2_[1] - p1_[1]);
	double D1 = b2 * (p2_[0] - p1_[0]) - b1 * (p4_[0] - p3_[0]);
	double D2 = b2 * (p2_[1] - p1_[1]) - b1 * (p4_[1] - p3_[1]);
	intersect_point[0] = D1 / D;
	intersect_point[1] = D2 / D;
	return intersect_point;
}

OpenMesh::HalfedgeHandle ParaQuadCutting::find_next_point(OpenMesh::VertexHandle vertex_, OpenMesh::Vec3d segment1_, OpenMesh::Vec3d segment2_, double &t_)
{
	std::vector<OpenMesh::HalfedgeHandle> he_list;
	std::vector<double> t_list;
	std::vector<double> length_list;
	std::vector<OpenMesh::Vec3d> point_list;
	std::vector<OpenMesh::HalfedgeHandle> he_list1;
	for (auto voh_h : new_para.voh_range(vertex_))
	{
		OpenMesh::HalfedgeHandle he_temp = new_para.next_halfedge_handle(voh_h);
		OpenMesh::VertexHandle fromv = new_para.from_vertex_handle(he_temp);
		OpenMesh::VertexHandle tov = new_para.to_vertex_handle(he_temp);
		OpenMesh::Vec3d fpoint = new_para.point(fromv);
		OpenMesh::Vec3d tpoint = new_para.point(tov);
		bool flag_temp = is_intersect(segment1_, segment2_, fpoint, tpoint);
		if (flag_temp)
		{
			he_list.push_back(he_temp);
		}
	}
	for (int i = 0; i < he_list.size(); i++)
	{
		OpenMesh::HalfedgeHandle he_temp = he_list[i];
		OpenMesh::VertexHandle fromv = new_para.from_vertex_handle(he_temp);
		OpenMesh::VertexHandle tov = new_para.to_vertex_handle(he_temp);
		OpenMesh::Vec3d fpoint = new_para.point(fromv);
		OpenMesh::Vec3d tpoint = new_para.point(tov);
		OpenMesh::Vec3d new_point_temp = compute_intersect_point(segment1_, segment2_, fpoint, tpoint);
		double t_temp = (new_point_temp - fpoint).norm() / (tpoint - fpoint).norm();
		if (t_temp <= 1.1)
		{
			t_list.push_back(t_temp);
			point_list.push_back(new_point_temp);
			double length_temp = (segment1_ - new_point_temp).norm();
			length_list.push_back(length_temp);
			he_list1.push_back(he_temp);
		}
	}
	if (length_list.size() > 0)
	{
		auto biggest = std::max_element(length_list.begin(), length_list.end());
		int biggest_position = std::distance(length_list.begin(), biggest);
		t_ = t_list[biggest_position];
		return he_list1[biggest_position];
	}
	else
	{
		//cout << "compute wrong" << std::endl;
		return OpenMesh::HalfedgeHandle();
	}
}

void ParaQuadCutting::process_point_in_face(OpenMesh::Vec3d p_, std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> &my_cut_segments_)
{
	OpenMesh::FaceHandle face_temp;
	bool flag_temp = 0;
	for (auto f_h : new_para.faces())
	{
		flag_temp = is_point_in_face(p_, f_h, new_para);
		if (flag_temp)
		{
			face_temp = f_h;
			break;
		}
	}
	if (!flag_temp)
	{
		return;
	}
	else
	{
		double u_temp, v_temp;
		OpenMesh::HalfedgeHandle face_he_temp = find_split_uv(p_, face_temp, new_para, u_temp, v_temp);

		if (!face_he_temp.is_valid())
		{
			auto dst_v = split_face(face_temp, u_temp, v_temp);
		}
		else
		{
			double t_t = find_split_t(p_, face_he_temp);
			OpenMesh::Vec3d new_point_temp;
			if (t_t < split_thres)
			{
				new_point_temp = new_para.point(new_para.from_vertex_handle(face_he_temp));
				update_cut_segments(p_, new_point_temp, my_cut_segments_);
			}
			else if (t_t > 1.0 - split_thres)
			{
				new_point_temp = new_para.point(new_para.to_vertex_handle(face_he_temp));
				update_cut_segments(p_, new_point_temp, my_cut_segments_);
			}
			else
			{
				new_point_temp = (1 - t_t)*new_para.point(new_para.from_vertex_handle(face_he_temp)) + t_t * new_para.point(new_para.to_vertex_handle(face_he_temp));
				update_cut_segments(p_, new_point_temp, my_cut_segments_);
			}
		}
	}
}

void ParaQuadCutting::update_cut_segments(OpenMesh::Vec3d p_old_, OpenMesh::Vec3d p_new_, std::vector<std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>> &my_cut_segments_)
{
	for (int i = 0; i < my_cut_segments_.size(); i++)
	{
		OpenMesh::Vec3d segment1 = my_cut_segments_[i].first;
		OpenMesh::Vec3d segment2 = my_cut_segments_[i].second;
		if (segment1 == p_old_)
		{
			my_cut_segments_[i].first = p_new_;
		}
		if (segment2 == p_old_)
		{
			my_cut_segments_[i].second = p_new_;
		}
	}
}

void ParaQuadCutting::mark_e_cut(OpenMesh::VertexHandle v_pre_, OpenMesh::VertexHandle v_now_)
{
	for (auto vv_h : new_para.vv_range(v_now_))
	{
		OpenMesh::HalfedgeHandle new_he_temp = new_para.find_halfedge(v_now_, vv_h);
		OpenMesh::EdgeHandle new_edge_temp = new_para.edge_handle(new_he_temp.idx() / 2);
		if (vv_h == v_pre_)
		{
			new_para.property(e_cut, new_edge_temp) = 1;
		}
		else
		{
			if (new_para.property(e_cut, new_edge_temp) != 1)
			{
				new_para.property(e_cut, new_edge_temp) = 0;
			}
		}
	}
}

void ParaQuadCutting::my_decomposition()
{
	my_para_atlas.clear();
	my_mesh_atlas.clear();
	f_atlas2mesh.clear();
	v_atlas2mesh.clear();
	face_layer.clear();
	v_para2layer.clear();
	std::vector<int> f_para2mesh;
	f_para2mesh.resize(new_para.n_faces());
	for (auto h_h : new_para.halfedges())
	{
		if (!new_para.is_boundary(h_h))
		{
			OpenMesh::FaceHandle para_face_temp = new_para.face_handle(h_h);
			auto v_pair = new_para.property(h_para2mesh, h_h);
			OpenMesh::HalfedgeHandle mesh_he_temp = new_mesh.find_halfedge(new_mesh.vertex_handle(v_pair.first), new_mesh.vertex_handle(v_pair.second));
			OpenMesh::FaceHandle mesh_face_temp = new_mesh.face_handle(mesh_he_temp);
			f_para2mesh[para_face_temp.idx()] = mesh_face_temp.idx();
		}
	}
	std::vector<int> v_para2mesh;
	v_para2mesh.resize(new_para.n_vertices());
	for (auto h_h : new_para.halfedges())
	{
		OpenMesh::VertexHandle fvertex = new_para.from_vertex_handle(h_h);
		OpenMesh::VertexHandle tvertex = new_para.to_vertex_handle(h_h);
		auto v_pair = new_para.property(h_para2mesh, h_h);
		v_para2mesh[fvertex.idx()] = v_pair.first;
		v_para2mesh[tvertex.idx()] = v_pair.second;
	}

	face_layer.assign(new_para.n_faces(), -1);
	std::vector<int> is_check_face;
	is_check_face.resize(new_para.n_faces());

	std::vector<int> is_check_face0;
	is_check_face0.resize(new_para.n_faces());
	while (std::accumulate(is_check_face0.begin(), is_check_face0.end(), 0) != new_para.n_faces())
	{
		std::vector<OpenMesh::FaceHandle> atlas_face;
		std::queue<OpenMesh::FaceHandle> myqueue;
		OpenMesh::FaceHandle root;
		for (auto f_h : new_para.faces())
		{
			if (!is_check_face0[f_h.idx()])
			{
				root = f_h;
				break;
			}
		}
		myqueue.push(root);
		is_check_face0[root.idx()] = 1;
		atlas_face.push_back(root);

		while (!myqueue.empty())
		{
			OpenMesh::FaceHandle face_temp = myqueue.front();
			myqueue.pop();

			for (auto fe_h : new_para.fe_range(face_temp))
			{
				if (new_para.property(e_cut, fe_h) == 0)
				{
					OpenMesh::HalfedgeHandle he_temp = new_para.halfedge_handle(fe_h, 0);
					OpenMesh::FaceHandle neighbor_face_temp = new_para.face_handle(he_temp) == face_temp ? new_para.opposite_face_handle(he_temp) : new_para.face_handle(he_temp);
					if (is_check_face0[neighbor_face_temp.idx()] == 0)
					{
						myqueue.push(neighbor_face_temp);
						is_check_face0[neighbor_face_temp.idx()] = 1;
						atlas_face.push_back(neighbor_face_temp);
					}
				}
			}

		}
		if (atlas_face.size() < 50)
		{
			bool find_e_cut = 0;
			for (int i = 0; i < atlas_face.size(); i++)
			{
				OpenMesh::FaceHandle face_temp = atlas_face[i];
				for (auto fe : new_para.fe_range(face_temp))
				{
					if (new_para.property(e_cut, fe) == 1)
					{
						new_para.property(e_cut, fe) = 0;
						find_e_cut = 1;
						break;
					}
				}
				if (find_e_cut)
				{
					break;
				}
			}
		}
	}

	while (std::accumulate(is_check_face.begin(), is_check_face.end(), 0) != new_para.n_faces())
	{
		my_para_atlas.emplace_back();
		Mesh& atlas_temp = my_para_atlas.back();
		std::vector<int> vertex_index(new_para.n_vertices(), -1);

		std::vector<OpenMesh::FaceHandle> atlas_face;
		std::vector<int> f_atlas2mesh_temp;
		std::vector<int> v_atlas2mesh_temp;
		std::queue<OpenMesh::FaceHandle> myqueue;
		OpenMesh::FaceHandle root;
		for (auto f_h : new_para.faces())
		{
			if (!is_check_face[f_h.idx()])
			{
				root = f_h;
				break;
			}
		}
		myqueue.push(root);
		is_check_face[root.idx()] = 1;
		atlas_face.push_back(root);

		std::vector<int> old_faces;
		std::vector<int> new_faces;
		std::vector<int> old_vertex;
		std::vector<int> new_vertex;
		while (!myqueue.empty())
		{
			OpenMesh::FaceHandle face_temp = myqueue.front();
			int fid_cur = face_temp.idx();
			myqueue.pop();
			if (face_layer[fid_cur] >= 0) continue;

			face_layer[fid_cur] = my_para_atlas.size() - 1;

			std::vector<OpenMesh::VertexHandle> v_fcur;
			for (auto fv : new_para.fv_range(face_temp))
			{
				if (vertex_index[fv.idx()] == -1)
				{
					vertex_index[fv.idx()] = atlas_temp.n_vertices();
					OpenMesh::VertexHandle newv = atlas_temp.add_vertex(new_para.point(fv));
					new_vertex.push_back(newv.idx());
					old_vertex.push_back(fv.idx());
				}
				v_fcur.emplace_back(vertex_index[fv.idx()]);
			}
			OpenMesh::FaceHandle newf = atlas_temp.add_face(v_fcur);
			new_faces.push_back(newf.idx());
			old_faces.push_back(face_temp.idx());

			for (auto fe_h : new_para.fe_range(face_temp))
			{
				if (new_para.property(e_cut, fe_h) == 0)
				{
					OpenMesh::HalfedgeHandle he_temp = new_para.halfedge_handle(fe_h, 0);
					OpenMesh::FaceHandle neighbor_face_temp = new_para.face_handle(he_temp) == face_temp ? new_para.opposite_face_handle(he_temp) : new_para.face_handle(he_temp);
					if (is_check_face[neighbor_face_temp.idx()] == 0)
					{
						myqueue.push(neighbor_face_temp);
						is_check_face[neighbor_face_temp.idx()] = 1;
						atlas_face.push_back(neighbor_face_temp);
					}
				}
			}	
		}
		f_atlas2mesh_temp.resize(atlas_face.size());
		for (int i = 0; i < atlas_face.size(); i++)
		{
			f_atlas2mesh_temp[new_faces[i]] = /*f_para2mesh[*/old_faces[i]/*]*/;
		}
		v_atlas2mesh_temp.resize(new_vertex.size());
		for (int i = 0; i < new_vertex.size(); i++)
		{
			v_atlas2mesh_temp[new_vertex[i]] = /*v_para2mesh[*/old_vertex[i]/*]*/;
		}
		
		f_atlas2mesh.push_back(f_atlas2mesh_temp);
		v_atlas2mesh.push_back(v_atlas2mesh_temp);
		v_para2layer.emplace_back(std::move(vertex_index));
	}

	for (int i = 0; i < my_para_atlas.size(); i++)
	{
		Mesh mesh_temp;
		Mesh atlas_temp = my_para_atlas[i];
		std::vector<OpenMesh::VertexHandle> v2newv;
		v2newv.resize(new_mesh.n_vertices());
		for (auto v_it = atlas_temp.vertices_begin(); v_it != atlas_temp.vertices_end(); v_it++)
		{
			OpenMesh::VertexHandle mesh_v = new_mesh.vertex_handle(v_atlas2mesh[i][v_it.handle().idx()]);
			OpenMesh::Vec3d mesh_p = new_mesh.point(mesh_v);
			OpenMesh::VertexHandle newv = mesh_temp.add_vertex(mesh_p);
			v2newv[mesh_v.idx()] = newv;
		}
		for (auto f_it = atlas_temp.faces_begin(); f_it != atlas_temp.faces_end(); f_it++)
		{
			OpenMesh::FaceHandle mesh_f = new_mesh.face_handle(f_atlas2mesh[i][f_it.handle().idx()]);
			std::vector<int> face_vertex;
			for (auto cfv_it = new_mesh.cfv_begin(mesh_f); cfv_it != new_mesh.cfv_end(mesh_f); cfv_it++)
			{
				face_vertex.push_back(cfv_it.handle().idx());
			}
			mesh_temp.add_face(v2newv[face_vertex[0]], v2newv[face_vertex[1]], v2newv[face_vertex[2]]);
		}
		my_mesh_atlas.push_back(mesh_temp);
	}
}

void ParaQuadCutting::my_cutting(double thres)
{
	QTime timer;
	split_thres = thres;

	timer.start();

	split_poly();
	read_cut2();
	cut_segment(my_cut_segments);

	Mesh new_para0 = new_para;
	Mesh new_mesh0 = new_mesh;

	my_decomposition();
	bool is_reconstruct = 0;
	para_atlas = my_para_atlas;
	
	uniform_tutte1(is_reconstruct);

	double my_pe = my_retangle_packing(my_ret_candidate);
	std::vector<double> my_pe_list;
	int last_para_num = 0;
	bool my_switch = 0;
	my_pe_list.push_back(my_pe);
	while (my_pe < 0.85)
	{
		if (para_atlas.size() >= 20)
		{
			break;
		}

		find_new_cut(my_switch);
		//my_switch = !my_switch;
		if (my_cut_segments_new.size() < 1)
		{
			//std::cout << "segment_new ERROR!!" << std::endl;
			break;
		}
		my_ret_candidate.clear();
		cut_segment(my_cut_segments_new);
		my_decomposition();
		is_reconstruct = 0;
		para_atlas = my_para_atlas;
		uniform_tutte1(is_reconstruct);
		my_pe = my_retangle_packing(my_ret_candidate);
		my_pe_list.push_back(my_pe);
		my_cut_segments_new.clear();
		if (para_atlas.size() == last_para_num)
		{
			//std::cout << "Cut segment_new ERROR!" << std::endl;
			//my_switch = !my_switch;
			break;
		}
		last_para_num = para_atlas.size();
	}
	
	if (my_pe < 0.85)
	{
		new_para = new_para0;
		new_mesh = new_mesh0;
		auto biggest = std::max_element(my_pe_list.begin(), my_pe_list.end());
		int new_cut_num = std::distance(my_pe_list.begin(), biggest);
		if (new_cut_num == 0)
		{
			my_decomposition();
			my_ret_candidate.clear();
			para_atlas = my_para_atlas;
			uniform_tutte1(is_reconstruct);
			my_pe = my_retangle_packing(my_ret_candidate);
		}

		else
		{
			my_decomposition();
			my_ret_candidate.clear();
			para_atlas = my_para_atlas;
			uniform_tutte1(is_reconstruct);
			my_pe = my_retangle_packing(my_ret_candidate);
			for (int i = 0; i < new_cut_num; i++)
			{
				find_new_cut(my_switch);
				//my_switch = !my_switch;
				if (my_cut_segments_new.size() < 1)
				{
					//std::cout << "segment_new ERROR!!" << std::endl;
					break;
				}
				my_ret_candidate.clear();
				cut_segment(my_cut_segments_new);
				my_decomposition();
				is_reconstruct = 0;
				para_atlas = my_para_atlas;
				uniform_tutte1(is_reconstruct);
				my_pe = my_retangle_packing(my_ret_candidate);
				my_cut_segments_new.clear();
			}
		}
	}

	const auto& chart_translation0 = my_ct_candidate;
	const auto& chart_flipped0 = my_cf_candidate;
	for (int i = 0; i < para_atlas.size(); i++)
	{
		double offset_u = chart_translation0[i][0];
		double offset_v = chart_translation0[i][1];
		for (int j = 0; j < para_atlas[i].n_vertices(); j++)
		{
			const auto& point_v = para_atlas[i].point(para_atlas[i].vertex_handle(j));
			const auto& point_v1 = my_origin_atlas[i].point(my_origin_atlas[i].vertex_handle(j));
			if (chart_flipped0[i] == 1)
			{
				OpenMesh::Vec3d point_temp(-point_v[1] + offset_u, point_v[0] + offset_v, 0);
				para_atlas[i].set_point(para_atlas[i].vertex_handle(j), point_temp);
				OpenMesh::Vec3d point_temp1(-point_v1[1] + offset_u, point_v1[0] + offset_v, 0);
				my_origin_atlas[i].set_point(my_origin_atlas[i].vertex_handle(j), point_temp1);
			}
			else
			{
				OpenMesh::Vec3d point_temp(point_v[0] + offset_u, point_v[1] + offset_v, 0);
				para_atlas[i].set_point(para_atlas[i].vertex_handle(j), point_temp);
				OpenMesh::Vec3d point_temp1(point_v1[0] + offset_u, point_v1[1] + offset_v, 0);
				my_origin_atlas[i].set_point(my_origin_atlas[i].vertex_handle(j), point_temp1);
			}
		}
	}


	for (int i = 0; i < para_atlas.size(); i++)
	{
		bool my_is_flip = 0;
		Mesh &para_temp = para_atlas[i];
		for (auto f_h : para_temp.faces())
		{
			std::vector<OpenMesh::Vec3d> face_point;
			for (auto fv_h : para_temp.fv_range(f_h))
			{
				face_point.push_back(para_temp.point(fv_h));
			}
			OpenMesh::Vec3d vec1 = face_point[1] - face_point[0];
			OpenMesh::Vec3d vec2 = face_point[2] - face_point[0];
			double my_det_p = OpenMesh::cross(vec1, vec2).norm();
			if (my_det_p <= 0)
			{
				my_is_flip = 1;
				break;
			}
		}
		if (my_is_flip)
		{
			para_atlas[i] = my_origin_atlas[i];
		}
	}

	my_packing();
	std::cout << "Time: " << timer.elapsed() / 1000.0 << "s" << std::endl;
}

void ParaQuadCutting::my_write_obj(Mesh mesh_, std::string str_)
{
	std::vector<std::vector<int>> face_vertex_list;
	std::vector<OpenMesh::Vec3d> vertex_point_list;
	face_vertex_list.resize(mesh_.n_faces());
	vertex_point_list.resize(mesh_.n_vertices());
	for (auto f_h : mesh_.faces())
	{
		for (auto fv_h : mesh_.fv_range(f_h))
		{
			face_vertex_list[f_h.idx()].push_back(fv_h.idx());
		}
	}
	for (auto v_h : mesh_.vertices())
	{
		OpenMesh::Vec3d point_temp = mesh_.point(v_h);
		vertex_point_list[v_h.idx()] = point_temp;
	}

	std::ofstream of_obj(str_, std::ios::trunc);

	for (int i = 0; i < vertex_point_list.size(); i++)
	{
		of_obj << "v " << vertex_point_list[i][0] << " " << vertex_point_list[i][1] << " " << vertex_point_list[i][2] << std::endl;
	}

	for (int i = 0; i < face_vertex_list.size(); i++)
	{
		of_obj << "f " << face_vertex_list[i][0] + 1 << " " << face_vertex_list[i][1] + 1 << " " << face_vertex_list[i][2] + 1 << std::endl;
	}
	of_obj.close();
}


void ParaQuadCutting::set_boundary(int ii_, std::vector<OpenMesh::VertexHandle> &cc_boundary_, std::vector<OpenMesh::Vec3d> &new_bpoint_)
{
	//std::cout << "iter======================================================================" << ii_ << std::endl;
	my_scale = 3000.0 / (BB_Max - BB_Min).norm();
	Mesh atlas_temp = my_para_atlas[ii_];
	std::vector<int> f_atlas2mesh_temp = f_atlas2mesh[ii_];
	double atlas_area = 0;
	for (auto f_h : atlas_temp.faces())
	{
		OpenMesh::FaceHandle mesh_face_temp = new_mesh.face_handle(f_atlas2mesh_temp[f_h.idx()]);
		std::vector<OpenMesh::VertexHandle> face_vertex;
		for (auto cfv_it = new_mesh.cfv_begin(mesh_face_temp); cfv_it != new_mesh.cfv_end(mesh_face_temp); cfv_it++)
		{
			face_vertex.push_back(cfv_it.handle());
		}
		OpenMesh::Vec3d A = new_mesh.point(face_vertex[0]);
		OpenMesh::Vec3d B = new_mesh.point(face_vertex[1]);
		OpenMesh::Vec3d C = new_mesh.point(face_vertex[2]);
		//double area_temp = std::fabs((A[0] * B[1] + B[0] * C[1] + C[0] * A[1] - A[0] * C[1] - B[0] * A[1] - C[0] * B[1]) / 2);
		double area_temp = std::fabs(((B - A) % (C - A)).length() / 2.0);
		atlas_area += area_temp;
	}

	OpenMesh::HalfedgeHandle hedge_init;
	OpenMesh::HalfedgeHandle hedge_temp;
	OpenMesh::VertexHandle vertex_init;
	OpenMesh::VertexHandle vertex_temp;
	for (auto he_h : atlas_temp.halfedges())
	{
		if (atlas_temp.is_boundary(he_h))
		{
			hedge_init = he_h;
			break;
		}
	}
	vertex_init = atlas_temp.to_vertex_handle(hedge_init);
	hedge_temp = hedge_init;
	vertex_temp = vertex_init;
	do
	{
		cc_boundary_.push_back(vertex_temp);
		vertex_temp = atlas_temp.from_vertex_handle(hedge_temp);
		hedge_temp = atlas_temp.prev_halfedge_handle(hedge_temp);
	} while (vertex_temp != vertex_init);

	std::vector<int> quad_corner_candidate;
	std::vector<double> angle_list;
	for (int i = 0; i < cc_boundary_.size(); i++)
	{
		double angle_temp = 0;
		OpenMesh::Vec3d p0 = atlas_temp.point(cc_boundary_[i]);
		OpenMesh::HalfedgeHandle h_iter = atlas_temp.find_halfedge(cc_boundary_[(i + cc_boundary_.size() - 1) % cc_boundary_.size()], cc_boundary_[i]);
		while (!atlas_temp.is_boundary(h_iter))
		{
			OpenMesh::Vec3d p1 = atlas_temp.point(atlas_temp.from_vertex_handle(h_iter));
			h_iter = atlas_temp.next_halfedge_handle(h_iter);
			OpenMesh::Vec3d p2 = atlas_temp.point(atlas_temp.to_vertex_handle(h_iter));

			OpenMesh::Vec3d vec1 = p1 - p0;
			OpenMesh::Vec3d vec2 = p2 - p0;
			//angle_temp += std::atan2(vec1[0] * vec2[1] - vec1[1] * vec2[0], vec1[0] * vec2[0] + vec1[1] * vec2[1]);
			angle_temp += std::acos((vec1 | vec2) / (vec1.norm()*vec2.norm()));

			h_iter = atlas_temp.opposite_halfedge_handle(h_iter);
		}
		angle_list.push_back(std::fabs(angle_temp));
	}
	for (int i = 0; i < angle_list.size(); i++)
	{
		if (angle_list[i] < 8.0 / 9.0*M_PI)
		{
			quad_corner_candidate.push_back(i);
		}
	}
	bool is_can_delete = 1;
	while (is_can_delete)
	{
		is_can_delete = 0;
		for (int i = 0; i < quad_corner_candidate.size(); i++)
		{
			OpenMesh::Vec3d candidate_point0 = atlas_temp.point(cc_boundary_[quad_corner_candidate[i]]);
			OpenMesh::Vec3d candidate_point1 = atlas_temp.point(cc_boundary_[quad_corner_candidate[(i + quad_corner_candidate.size() - 1) % quad_corner_candidate.size()]]);
			OpenMesh::Vec3d candidate_point2 = atlas_temp.point(cc_boundary_[quad_corner_candidate[(i + 1) % quad_corner_candidate.size()]]);
			double x = ((candidate_point1 - candidate_point0) | (candidate_point2 - candidate_point0)) / ((candidate_point1 - candidate_point0).norm()*(candidate_point2 - candidate_point0).norm());
			if (std::fabs(x) >= 1)
			{
				x = x > 0 ? 1 : -1;
			}
			double angle_temp = std::acos(x);
			if (angle_temp > 8.0 / 9.0*M_PI)
			{
				quad_corner_candidate.erase(quad_corner_candidate.begin() + i);
				is_can_delete = 1;
				break;
			}
		}
	}

	std::vector<int> quad_corner;
	std::vector<double> quad_edge_length;
	double boundary_length = 0;
	for (int i = 0; i < cc_boundary_.size(); i++)
	{
		OpenMesh::Vec3d point_temp1 = atlas_temp.point(cc_boundary_[i]);
		OpenMesh::Vec3d point_temp2 = atlas_temp.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
		boundary_length += (point_temp1 - point_temp2).norm();
	}
	std::vector<double> quad_length_temp;
	int corner_temp;
	int corner_num = quad_corner_candidate.size();
	if (corner_num < 4)
	{
		Mesh atlas_temp1 = atlas_temp;
		build_rotation_message(atlas_temp);
		my_calc_global_rotation();
		for (auto v_h : atlas_temp1.vertices())
		{
			OpenMesh::Vec3d point_temp(my_uv_x[2 * v_h.idx()], my_uv_x[2 * v_h.idx() + 1], 0);
			atlas_temp1.set_point(v_h, point_temp);
		}
		OpenMesh::Vec3d b_max, b_min;
		find_BB(atlas_temp1, b_max, b_min);

		double b_width_temp = b_max[0] - b_min[0];
		double b_height_temp = b_max[1] - b_min[1];
		double b_ratio_hw = b_height_temp / b_width_temp;
		double width_temp = std::sqrt(atlas_area / b_ratio_hw);
		double height_temp = width_temp * b_ratio_hw;

		int corner1;
		int corner2;
		int corner3;
		int corner4;
		OpenMesh::Vec3d start_point(b_min[0], b_max[1], 0);
		std::vector<double> start_length;
		double boundary_length = 0;
		for (int i = 0; i < cc_boundary_.size(); i++)
		{
			OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i]);
			OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
			start_length.push_back((point_temp - start_point).norm());
			boundary_length += (point_temp1 - point_temp).norm();
		}
		auto smallest = std::min_element(start_length.begin(), start_length.end());
		corner1 = std::distance(start_length.begin(), smallest);

		double boundary_length1 = boundary_length * (height_temp / (2 * (height_temp + width_temp)));
		double boundary_length2 = boundary_length * (width_temp / (2 * (height_temp + width_temp)));
		double edge_length_temp = 0;
		for (int i = corner1; i < corner1 + cc_boundary_.size(); i++)
		{
			OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i%cc_boundary_.size()]);
			OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
			edge_length_temp += (point_temp1 - point_temp).norm();
			if (edge_length_temp > boundary_length1)
			{
				corner2 = (i + 1) % cc_boundary_.size();
				break;
			}
		}
		edge_length_temp = 0;
		for (int i = corner2; i < corner2 + cc_boundary_.size(); i++)
		{
			OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i%cc_boundary_.size()]);
			OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
			edge_length_temp += (point_temp1 - point_temp).norm();
			if (edge_length_temp > boundary_length2)
			{
				corner3 = (i + 1) % cc_boundary_.size();
				break;
			}
		}
		edge_length_temp = 0;
		for (int i = corner3; i < corner3 + cc_boundary_.size(); i++)
		{
			OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i%cc_boundary_.size()]);
			OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
			edge_length_temp += (point_temp1 - point_temp).norm();
			if (edge_length_temp > boundary_length1)
			{
				corner4 = (i + 1) % cc_boundary_.size();
				break;
			}
		}
		quad_corner.push_back(corner1);
		quad_corner.push_back(corner2);
		quad_corner.push_back(corner3);
		quad_corner.push_back(corner4);
		quad_edge_length.push_back(height_temp);
		quad_edge_length.push_back(width_temp);
		quad_edge_length.push_back(height_temp);
		quad_edge_length.push_back(width_temp);

		OpenMesh::Vec3d center_point = (b_max + b_min) / 2;
		b_max = b_max - center_point;
		b_min = b_min - center_point;
		int width_int = std::lround(width_temp * (my_scale / 2));
		int height_int = std::lround(height_temp * (my_scale / 2));
		double ratio1 = 2 * width_int / (b_max[0] - b_min[0]);
		double ratio2 = 2 * height_int / (b_max[1] - b_min[1]);
		double ratio = ratio1 > ratio2 ? ratio2 : ratio1;
		for (auto v_h : atlas_temp1.vertices())
		{
			OpenMesh::Vec3d new_point = atlas_temp1.point(v_h) - center_point;
			new_point = new_point * ratio;
			atlas_temp1.set_point(v_h, new_point);
		}
		my_origin_atlas.push_back(atlas_temp1);
	}
		
	if (corner_num == 4)
	{
		//std::cout << "case four corner" << std::endl;
		quad_corner.clear();
		if (my_para_atlas.size() <= 20)
		{
			quad_corner = quad_corner_candidate;
			//std::cout << "quad corner candidate" << quad_corner_candidate[0] << " " << quad_corner_candidate[1] << " " << quad_corner_candidate[2] << " " << quad_corner_candidate[3] << std::endl;
			double length1 = 0;
			double length2 = 0;
			double length3 = 0;
			double length4 = 0;
			for (int i = quad_corner_candidate[0]; i < quad_corner_candidate[1]; i++)
			{
				//OpenMesh::Vec3d point_temp1 = atlas_temp.point(cc_boundary_[i%cc_boundary_.size()]);
				//OpenMesh::Vec3d point_temp2 = atlas_temp.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				OpenMesh::Vec3d point_temp1 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[i%cc_boundary_.size()].idx()]));
				OpenMesh::Vec3d point_temp2 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[(i + 1) % cc_boundary_.size()].idx()]));
				length1 += (point_temp1 - point_temp2).norm();
			}
			for (int i = quad_corner_candidate[1]; i < quad_corner_candidate[2]; i++)
			{
				//OpenMesh::Vec3d point_temp1 = atlas_temp.point(cc_boundary_[i%cc_boundary_.size()]);
				//OpenMesh::Vec3d point_temp2 = atlas_temp.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				OpenMesh::Vec3d point_temp1 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[i%cc_boundary_.size()].idx()]));
				OpenMesh::Vec3d point_temp2 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[(i + 1) % cc_boundary_.size()].idx()]));
				length2 += (point_temp1 - point_temp2).norm();
			}
			for (int i = quad_corner_candidate[2]; i < quad_corner_candidate[3]; i++)
			{
				//OpenMesh::Vec3d point_temp1 = atlas_temp.point(cc_boundary_[i%cc_boundary_.size()]);
				//OpenMesh::Vec3d point_temp2 = atlas_temp.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				OpenMesh::Vec3d point_temp1 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[i%cc_boundary_.size()].idx()]));
				OpenMesh::Vec3d point_temp2 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[(i + 1) % cc_boundary_.size()].idx()]));
				length3 += (point_temp1 - point_temp2).norm();
			}
			for (int i = quad_corner_candidate[3]; i < quad_corner_candidate[0] + cc_boundary_.size(); i++)
			{
				//OpenMesh::Vec3d point_temp1 = atlas_temp.point(cc_boundary_[i%cc_boundary_.size()]);
				//OpenMesh::Vec3d point_temp2 = atlas_temp.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				OpenMesh::Vec3d point_temp1 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[i%cc_boundary_.size()].idx()]));
				OpenMesh::Vec3d point_temp2 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[(i + 1) % cc_boundary_.size()].idx()]));
				length4 += (point_temp1 - point_temp2).norm();
			}
			quad_edge_length.push_back(length1);
			quad_edge_length.push_back(length2);
			quad_edge_length.push_back(length3);
			quad_edge_length.push_back(length4);

			Mesh atlas_temp1 = atlas_temp;
			build_rotation_message(atlas_temp);
			my_calc_global_rotation();
			for (auto v_h : atlas_temp1.vertices())
			{
				OpenMesh::Vec3d point_temp(my_uv_x[2 * v_h.idx()], my_uv_x[2 * v_h.idx() + 1], 0);
				atlas_temp1.set_point(v_h, point_temp);
			}
			OpenMesh::Vec3d b_max, b_min;
			find_BB(atlas_temp1, b_max, b_min);

			OpenMesh::Vec3d center_point = (b_max + b_min) / 2;
			b_max = b_max - center_point;
			b_min = b_min - center_point;
			int width_int = std::lround(((length2 + length4) / 2) * (my_scale / 2));
			int height_int = std::lround(((length1 + length3) / 2) * (my_scale / 2));
			double ratio1 = 2 * width_int / (b_max[0] - b_min[0]);
			double ratio2 = 2 * height_int / (b_max[1] - b_min[1]);
			double ratio = ratio1 > ratio2 ? ratio2 : ratio1;
			for (auto v_h : atlas_temp1.vertices())
			{
				OpenMesh::Vec3d new_point = atlas_temp1.point(v_h) - center_point;
				new_point = new_point * ratio;
				atlas_temp1.set_point(v_h, new_point);
			}
			my_origin_atlas.push_back(atlas_temp1);
		}
		else
		{
			Mesh atlas_temp1 = atlas_temp;
			build_rotation_message(atlas_temp);
			my_calc_global_rotation();
			for (auto v_h : atlas_temp1.vertices())
			{
				OpenMesh::Vec3d point_temp(my_uv_x[2 * v_h.idx()], my_uv_x[2 * v_h.idx() + 1], 0);
				atlas_temp1.set_point(v_h, point_temp);
			}
			OpenMesh::Vec3d b_max, b_min;
			find_BB(atlas_temp1, b_max, b_min);

			double b_width_temp = b_max[0] - b_min[0];
			double b_height_temp = b_max[1] - b_min[1];
			double b_ratio_hw = b_height_temp / b_width_temp;
			double width_temp = std::sqrt(atlas_area / b_ratio_hw);
			double height_temp = width_temp * b_ratio_hw;

			int corner1;
			int corner2;
			int corner3;
			int corner4;
			OpenMesh::Vec3d start_point(b_min[0], b_max[1], 0);
			std::vector<double> start_length;
			double boundary_length = 0;
			for (int i = 0; i < cc_boundary_.size(); i++)
			{
				OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i]);
				OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				start_length.push_back((point_temp - start_point).norm());
				boundary_length += (point_temp1 - point_temp).norm();
			}
			auto smallest = std::min_element(start_length.begin(), start_length.end());
			corner1 = std::distance(start_length.begin(), smallest);
			/*int quad_vertex_num = n_boundary_edges / 4;
			corner2 = (corner1 + quad_vertex_num) % n_boundary_edges;
			corner3 = (corner2 + quad_vertex_num) % n_boundary_edges;
			corner4 = (corner3 + quad_vertex_num) % n_boundary_edges;*/

			double boundary_length1 = boundary_length * (height_temp / (2 * (height_temp + width_temp)));
			double boundary_length2 = boundary_length * (width_temp / (2 * (height_temp + width_temp)));
			double edge_length_temp = 0;
			for (int i = corner1; i < corner1 + cc_boundary_.size(); i++)
			{
				OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i%cc_boundary_.size()]);
				OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				edge_length_temp += (point_temp1 - point_temp).norm();
				if (edge_length_temp > boundary_length1)
				{
					corner2 = (i + 1) % cc_boundary_.size();
					break;
				}
			}
			edge_length_temp = 0;
			for (int i = corner2; i < corner2 + cc_boundary_.size(); i++)
			{
				OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i%cc_boundary_.size()]);
				OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				edge_length_temp += (point_temp1 - point_temp).norm();
				if (edge_length_temp > boundary_length2)
				{
					corner3 = (i + 1) % cc_boundary_.size();
					break;
				}
			}
			edge_length_temp = 0;
			for (int i = corner3; i < corner3 + cc_boundary_.size(); i++)
			{
				OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i%cc_boundary_.size()]);
				OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				edge_length_temp += (point_temp1 - point_temp).norm();
				if (edge_length_temp > boundary_length1)
				{
					corner4 = (i + 1) % cc_boundary_.size();
					break;
				}
			}
			quad_corner.push_back(corner1);
			quad_corner.push_back(corner2);
			quad_corner.push_back(corner3);
			quad_corner.push_back(corner4);
			quad_edge_length.push_back(height_temp);
			quad_edge_length.push_back(width_temp);
			quad_edge_length.push_back(height_temp);
			quad_edge_length.push_back(width_temp);

			OpenMesh::Vec3d center_point = (b_max + b_min) / 2;
			b_max = b_max - center_point;
			b_min = b_min - center_point;
			int width_int = std::lround(width_temp * (my_scale / 2));
			int height_int = std::lround(height_temp * (my_scale / 2));
			double ratio1 = 2 * width_int / (b_max[0] - b_min[0]);
			double ratio2 = 2 * height_int / (b_max[1] - b_min[1]);
			double ratio = ratio1 > ratio2 ? ratio2 : ratio1;
			for (auto v_h : atlas_temp1.vertices())
			{
				OpenMesh::Vec3d new_point = atlas_temp1.point(v_h) - center_point;
				new_point = new_point * ratio;
				atlas_temp1.set_point(v_h, new_point);
			}
			my_origin_atlas.push_back(atlas_temp1);
		}
	}
		
	if (corner_num > 4)
	{
		//std::cout << "case more corner" << std::endl;
		quad_corner.clear();
		bool is_first_valid = 1;
		if (my_para_atlas.size() <= 20)
		{
			//std::cout << "more case 1" << std::endl;
			std::vector<double> angle_plus;
			for (int i = 0; i < quad_corner_candidate.size(); i++)
			{
				int corner1 = quad_corner_candidate[i];
				int corner2 = quad_corner_candidate[(i + 1) % quad_corner_candidate.size()];
				angle_plus.push_back(angle_list[corner1] + angle_list[corner2]);
				if (corner2 > corner1)
				{
					double length0 = 0;
					for (int j = corner1; j < corner2; j++)
					{
						OpenMesh::Vec3d point_temp1 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[j%cc_boundary_.size()].idx()]));
						OpenMesh::Vec3d point_temp2 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[(j + 1) % cc_boundary_.size()].idx()]));
						length0 += (point_temp1 - point_temp2).norm();
					}
					quad_length_temp.push_back(length0);
				}
				else
				{
					double length0 = 0;
					for (int j = corner1; j < corner2 + cc_boundary_.size(); j++)
					{
						OpenMesh::Vec3d point_temp1 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[j%cc_boundary_.size()].idx()]));
						OpenMesh::Vec3d point_temp2 = new_mesh.point(new_mesh.vertex_handle(v_atlas2mesh[ii_][cc_boundary_[(j + 1) % cc_boundary_.size()].idx()]));
						length0 += (point_temp1 - point_temp2).norm();
					}
					quad_length_temp.push_back(length0);
				}
			}
			auto small_angle = std::min_element(angle_plus.begin(), angle_plus.end());
			int small_angle_p = std::distance(angle_plus.begin(), small_angle);
			double length1 = quad_length_temp[small_angle_p];
			std::vector<double> quad_energy;
			int oppo_edge_start = (small_angle_p + 2) % quad_length_temp.size();
			int oppo_edge_end = (small_angle_p + quad_length_temp.size() - 1) % quad_length_temp.size();
			if (oppo_edge_end > oppo_edge_start)
			{
				for (int i = oppo_edge_start; i < oppo_edge_end; i++)
				{
					double length3 = quad_length_temp[i];
					double length2 = quad_length_temp[(small_angle_p + 1) % quad_length_temp.size()];
					for (int j = oppo_edge_start; j < i; j++)
					{
						length2 += quad_length_temp[j];
					}
					double length4 = boundary_length - length1 - length2 - length3;
					double quad_energy_temp = (length1 - length3)*(length1 - length3) + (length2 - length4)*(length2 - length4);
					quad_energy.push_back(quad_energy_temp);
				}
			}
			else
			{
				for (int i = oppo_edge_start; i < oppo_edge_end + quad_length_temp.size(); i++)
				{
					double length3 = quad_length_temp[i%quad_length_temp.size()];
					double length2 = quad_length_temp[(small_angle_p + 1) % quad_length_temp.size()];
					for (int j = oppo_edge_start; j < i; j++)
					{
						length2 += quad_length_temp[j%quad_length_temp.size()];
					}
					double length4 = boundary_length - length1 - length2 - length3;
					double quad_energy_temp = (length1 - length3)*(length1 - length3) + (length2 - length4)*(length2 - length4);
					quad_energy.push_back(quad_energy_temp);
				}
			}

			auto small_energy = std::min_element(quad_energy.begin(), quad_energy.end());
			int small_energy_p = std::distance(quad_energy.begin(), small_energy);
			int oppo_edge_p = (small_angle_p + 2 + small_energy_p) % quad_length_temp.size();
			double fix_length = quad_length_temp[small_angle_p] > quad_length_temp[oppo_edge_p] ? quad_length_temp[small_angle_p] : quad_length_temp[oppo_edge_p];
			quad_corner.push_back(quad_corner_candidate[small_angle_p]);
			quad_corner.push_back(quad_corner_candidate[(small_angle_p + 1) % quad_corner_candidate.size()]);
			quad_corner.push_back(quad_corner_candidate[oppo_edge_p]);
			quad_corner.push_back(quad_corner_candidate[(oppo_edge_p + 1) % quad_corner_candidate.size()]);
			quad_edge_length.push_back(fix_length);
			quad_edge_length.push_back((boundary_length - 2 * fix_length) / 2);
			quad_edge_length.push_back(fix_length);
			quad_edge_length.push_back((boundary_length - 2 * fix_length) / 2);

			if ((boundary_length - 2 * fix_length) <= 0 || (fix_length / (boundary_length - 2 * fix_length)) >= 10 || (fix_length / (boundary_length - 2 * fix_length)) <= 0.1)
			{
				is_first_valid = 0;
			}
			else
			{
				Mesh atlas_temp1 = atlas_temp;
				build_rotation_message(atlas_temp);
				my_calc_global_rotation();
				for (auto v_h : atlas_temp1.vertices())
				{
					OpenMesh::Vec3d point_temp(my_uv_x[2 * v_h.idx()], my_uv_x[2 * v_h.idx() + 1], 0);
					atlas_temp1.set_point(v_h, point_temp);
				}
				OpenMesh::Vec3d b_max, b_min;
				find_BB(atlas_temp1, b_max, b_min);

				OpenMesh::Vec3d center_point = (b_max + b_min) / 2;
				b_max = b_max - center_point;
				b_min = b_min - center_point;
				int width_int = std::lround(fix_length * (my_scale / 2));
				int height_int = std::lround(((boundary_length - 2 * fix_length) / 2) * (my_scale / 2));
				double ratio1 = 2 * width_int / (b_max[0] - b_min[0]);
				double ratio2 = 2 * height_int / (b_max[1] - b_min[1]);
				double ratio = ratio1 > ratio2 ? ratio2 : ratio1;
				for (auto v_h : atlas_temp1.vertices())
				{
					OpenMesh::Vec3d new_point = atlas_temp1.point(v_h) - center_point;
					new_point = new_point * ratio;
					atlas_temp1.set_point(v_h, new_point);
				}
				my_origin_atlas.push_back(atlas_temp1);
			}
		}
		if (my_para_atlas.size() > 20 || !is_first_valid)
		{
			quad_corner.clear();
			quad_edge_length.clear();
			Mesh atlas_temp1 = atlas_temp;
			build_rotation_message(atlas_temp1);
			my_calc_global_rotation();
			for (auto v_h : atlas_temp1.vertices())
			{
				OpenMesh::Vec3d point_temp(my_uv_x[2 * v_h.idx()], my_uv_x[2 * v_h.idx() + 1], 0);
				atlas_temp1.set_point(v_h, point_temp);
			}
			OpenMesh::Vec3d b_max, b_min;
			find_BB(atlas_temp1, b_max, b_min);

			double b_width_temp = b_max[0] - b_min[0];
			double b_height_temp = b_max[1] - b_min[1];
			double width_temp = std::sqrt(atlas_area*b_width_temp / b_height_temp);
			double height_temp = width_temp * b_height_temp / b_width_temp;
			int corner1;
			int corner2;
			int corner3;
			int corner4;
			OpenMesh::Vec3d start_point(b_min[0], b_max[1], 0);
			std::vector<double> start_length;
			double boundary_length = 0;
			for (int i = 0; i < cc_boundary_.size(); i++)
			{
				OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i]);
				OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				start_length.push_back((point_temp - start_point).norm());
				boundary_length += (point_temp1 - point_temp).norm();
			}
			auto smallest = std::min_element(start_length.begin(), start_length.end());
			corner1 = std::distance(start_length.begin(), smallest);

			double boundary_length1 = boundary_length * (height_temp / (2 * (height_temp + width_temp)));
			double boundary_length2 = boundary_length * (width_temp / (2 * (height_temp + width_temp)));
			double edge_length_temp = 0;
			for (int i = corner1; i < corner1 + cc_boundary_.size(); i++)
			{
				OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i%cc_boundary_.size()]);
				OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				edge_length_temp += (point_temp1 - point_temp).norm();
				if (edge_length_temp > boundary_length1)
				{
					corner2 = (i + 1) % cc_boundary_.size();
					break;
				}
			}
			edge_length_temp = 0;
			for (int i = corner2; i < corner2 + cc_boundary_.size(); i++)
			{
				OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i%cc_boundary_.size()]);
				OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				edge_length_temp += (point_temp1 - point_temp).norm();
				if (edge_length_temp > boundary_length2)
				{
					corner3 = (i + 1) % cc_boundary_.size();
					break;
				}
			}
			edge_length_temp = 0;
			for (int i = corner3; i < corner3 + cc_boundary_.size(); i++)
			{
				OpenMesh::Vec3d point_temp = atlas_temp1.point(cc_boundary_[i%cc_boundary_.size()]);
				OpenMesh::Vec3d point_temp1 = atlas_temp1.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
				edge_length_temp += (point_temp1 - point_temp).norm();
				if (edge_length_temp > boundary_length1)
				{
					corner4 = (i + 1) % cc_boundary_.size();
					break;
				}
			}
			
			quad_corner.push_back(corner1);
			quad_corner.push_back(corner2);
			quad_corner.push_back(corner3);
			quad_corner.push_back(corner4);
			quad_edge_length.push_back(height_temp);
			quad_edge_length.push_back(width_temp);
			quad_edge_length.push_back(height_temp);
			quad_edge_length.push_back(width_temp);

			OpenMesh::Vec3d center_point = (b_max + b_min) / 2;
			b_max = b_max - center_point;
			b_min = b_min - center_point;
			int width_int = std::lround(width_temp * (my_scale / 2));
			int height_int = std::lround(height_temp * (my_scale / 2));
			double ratio1 = 2 * width_int / (b_max[0] - b_min[0]);
			double ratio2 = 2 * height_int / (b_max[1] - b_min[1]);
			double ratio = ratio1 > ratio2 ? ratio2 : ratio1;
			for (auto v_h : atlas_temp1.vertices())
			{
				OpenMesh::Vec3d new_point = atlas_temp1.point(v_h) - center_point;
				new_point = new_point * ratio;
				atlas_temp1.set_point(v_h, new_point);
			}
			my_origin_atlas.push_back(atlas_temp1);
		}
	}
		
	double ratio_hw = (quad_edge_length[0] + quad_edge_length[2]) / (quad_edge_length[1] + quad_edge_length[3]);
	
	double quad_width = std::sqrt(atlas_area / ratio_hw);
	double quad_height = quad_width * ratio_hw;
	
	std::vector<int> edge_vertex1, edge_vertex2, edge_vertex3, edge_vertex4;
	int edge_vertex1_begin = quad_corner[0];
	int edge_vertex1_end = quad_corner[1] > quad_corner[0] ? quad_corner[1] : quad_corner[1] + cc_boundary_.size();
	for (int i = edge_vertex1_begin; i < edge_vertex1_end; i++)
	{
		edge_vertex1.push_back(i%cc_boundary_.size());
	}
	int edge_vertex2_begin = quad_corner[1];
	int edge_vertex2_end = quad_corner[2] > quad_corner[1] ? quad_corner[2] : quad_corner[2] + cc_boundary_.size();
	for (int i = edge_vertex2_begin; i < edge_vertex2_end; i++)
	{
		edge_vertex2.push_back(i%cc_boundary_.size());
	}
	int edge_vertex3_begin = quad_corner[2];
	int edge_vertex3_end = quad_corner[3] > quad_corner[2] ? quad_corner[3] : quad_corner[3] + cc_boundary_.size();
	for (int i = edge_vertex3_begin; i < edge_vertex3_end; i++)
	{
		edge_vertex3.push_back(i%cc_boundary_.size());
	}
	int edge_vertex4_begin = quad_corner[3];
	int edge_vertex4_end = quad_corner[0] > quad_corner[3] ? quad_corner[0] : quad_corner[0] + cc_boundary_.size();
	for (int i = edge_vertex4_begin; i < edge_vertex4_end; i++)
	{
		edge_vertex4.push_back(i%cc_boundary_.size());
	}
	
	int quad_width_int = std::lround(quad_width * (my_scale/2));
	int quad_height_int = std::lround(quad_height * (my_scale / 2));
	if (quad_width_int == 0)
	{
		quad_width_int = 1;
	}
	if (quad_height_int == 0)
	{
		quad_height_int = 1;
	}
	double quad_width_half = quad_width_int;
	double quad_height_half = quad_height_int;
	OpenMesh::Vec3d quad_point1(-quad_width_half, quad_height_half, 0);
	OpenMesh::Vec3d quad_point2(-quad_width_half, -quad_height_half, 0);
	OpenMesh::Vec3d quad_point3(quad_width_half, -quad_height_half, 0);
	OpenMesh::Vec3d quad_point4(quad_width_half, quad_height_half, 0);
	
	new_bpoint_.resize(cc_boundary_.size());
	for (int i = 0; i < edge_vertex1.size(); i++)
	{
		OpenMesh::Vec3d new_bpoint_temp = quad_point1 + (quad_point2 - quad_point1) / edge_vertex1.size()*i;
		new_bpoint_[edge_vertex1[i]] = new_bpoint_temp;
	}
	for (int i = 0; i < edge_vertex2.size(); i++)
	{
		OpenMesh::Vec3d new_bpoint_temp = quad_point2 + (quad_point3 - quad_point2) / edge_vertex2.size()*i;
		new_bpoint_[edge_vertex2[i]] = new_bpoint_temp;
	}
	for (int i = 0; i < edge_vertex3.size(); i++)
	{
		OpenMesh::Vec3d new_bpoint_temp = quad_point3 + (quad_point4 - quad_point3) / edge_vertex3.size()*i;
		new_bpoint_[edge_vertex3[i]] = new_bpoint_temp;
	}
	for (int i = 0; i < edge_vertex4.size(); i++)
	{
		OpenMesh::Vec3d new_bpoint_temp = quad_point4 + (quad_point1 - quad_point4) / edge_vertex4.size()*i;
		new_bpoint_[edge_vertex4[i]] = new_bpoint_temp;
	}

	OpenMesh::Vec2i ret_min(-quad_width_int, -quad_height_int);
	OpenMesh::Vec2i ret_max(quad_width_int, quad_height_int);
	retangle_type ret_temp;
	ret_temp.min = ret_min;
	ret_temp.max = ret_max;
	my_ret_candidate.push_back(ret_temp);
}

void ParaQuadCutting::uniform_tutte()
{
	my_total_area = 0;
	for (int ii = 0; ii < my_para_atlas.size(); ii++)
	{
		Mesh& atlas_temp = my_para_atlas[ii];
		std::vector<OpenMesh::VertexHandle> cc_boundary;
		std::vector<OpenMesh::Vec3d> new_bpoint;
		set_boundary(ii, cc_boundary, new_bpoint);

		int num = atlas_temp.n_vertices();
		Eigen::SparseMatrix<double> A(num, num);
		std::vector<Eigen::Triplet<double>> triplets;
		Eigen::MatrixXd b;
		for (auto v_it = atlas_temp.vertices_begin(); v_it != atlas_temp.vertices_end(); v_it++)
		{
			auto vertex = v_it.handle();
			if (atlas_temp.is_boundary(vertex))
			{
				triplets.push_back(Eigen::Triplet<double>(vertex.idx(), vertex.idx(), 1));
			}
			else
			{
				triplets.push_back(Eigen::Triplet<double>(vertex.idx(), vertex.idx(), 1));
				int neighbornum = 0;
				for (auto v_it = atlas_temp.vv_begin(vertex); v_it != atlas_temp.vv_end(vertex); v_it++)
				{
					neighbornum++;
				}
				for (auto v_it = atlas_temp.vv_begin(vertex); v_it != atlas_temp.vv_end(vertex); v_it++)
				{
					auto nvertex = v_it.handle();
					triplets.push_back(Eigen::Triplet<double>(vertex.idx(), nvertex.idx(), -1.0 / neighbornum));
				}
			}
		}
		A.setFromTriplets(triplets.begin(), triplets.end());
		A.makeCompressed();
		Eigen::SparseLU<Eigen::SparseMatrix<double>> LU;
		LU.compute(A);

		b.setZero(num, 3);
		for (int i = 0; i < cc_boundary.size(); i++)
		{
			b(cc_boundary[i].idx(), 0) = new_bpoint[i][0];
			b(cc_boundary[i].idx(), 1) = new_bpoint[i][1];
			b(cc_boundary[i].idx(), 2) = new_bpoint[i][2];
		}

		Eigen::MatrixXd sol(num, 3);
		sol = LU.solve(b);
		
		for (auto v_h : atlas_temp.vertices())
		{
			OpenMesh::Vec3d newposition;
			newposition[0] = sol(v_h.idx(), 0);
			newposition[1] = sol(v_h.idx(), 1);
			newposition[2] = sol(v_h.idx(), 2);
			atlas_temp.set_point(v_h, newposition);
		}
	}
	para_atlas = my_para_atlas;
}

void ParaQuadCutting::set_boundary1(int ii_, std::vector<OpenMesh::VertexHandle> &cc_boundary_, std::vector<OpenMesh::Vec3d> &new_bpoint_)
{
	Mesh para_temp = para_atlas[ii_];
	Mesh mesh_temp = my_mesh_atlas[ii_];
	Eigen::MatrixXd v_pos;
	Eigen::MatrixXd uv_v_pos;
	Eigen::MatrixXi fv_id;
	Eigen::MatrixXi uv_fv_id;
	my_get_scaf_info(v_pos, uv_v_pos, fv_id, uv_fv_id, para_temp, mesh_temp);

	double chart_gap = 0;
	double pe_bound = 0.7;
	StateManager s_;
	s_.run_interface1(v_pos, uv_v_pos, fv_id, uv_fv_id, para_distortion * 4.0, chart_gap, pe_bound);

	Eigen::MatrixXd uv_v_pos_after = s_.get_uv();
	for (int j = 0; j < para_temp.n_vertices(); j++)
	{
		auto& pt = para_temp.point(para_temp.vertex_handle(j));
		pt[0] = uv_v_pos_after(j, 0);
		pt[1] = uv_v_pos_after(j, 1);
	}

	OpenMesh::Vec3d b_max;
	OpenMesh::Vec3d b_min;
	find_BB(para_temp, b_max, b_min);

	std::vector<int> f_atlas2mesh_temp = f_atlas2mesh[ii_];
	double para_area = 0;
	for (auto f_h : para_temp.faces())
	{
		OpenMesh::FaceHandle mesh_face_temp = new_mesh.face_handle(f_atlas2mesh_temp[f_h.idx()]);
		std::vector<OpenMesh::VertexHandle> face_vertex;
		for (auto cfv_it = new_mesh.cfv_begin(mesh_face_temp); cfv_it != new_mesh.cfv_end(mesh_face_temp); cfv_it++)
		{
			face_vertex.push_back(cfv_it.handle());
		}
		OpenMesh::Vec3d A = new_mesh.point(face_vertex[0]);
		OpenMesh::Vec3d B = new_mesh.point(face_vertex[1]);
		OpenMesh::Vec3d C = new_mesh.point(face_vertex[2]);
		//double area_temp = (A[0] * B[1] + B[0] * C[1] + C[0] * A[1] - A[0] * C[1] - B[0] * A[1] - C[0] * B[1]) / 2;
		double area_temp = std::fabs(((B - A) % (C - A)).norm() / 2);
		para_area += area_temp;
	}

	OpenMesh::HalfedgeHandle hedge_init;
	OpenMesh::HalfedgeHandle hedge_temp;
	OpenMesh::VertexHandle vertex_init;
	OpenMesh::VertexHandle vertex_temp;
	for (auto he_h : para_temp.halfedges())
	{
		if (para_temp.is_boundary(he_h))
		{
			hedge_init = he_h;
			break;
		}
	}
	vertex_init = para_temp.to_vertex_handle(hedge_init);
	hedge_temp = hedge_init;
	vertex_temp = vertex_init;
	do
	{
		cc_boundary_.push_back(vertex_temp);
		vertex_temp = para_temp.from_vertex_handle(hedge_temp);
		hedge_temp = para_temp.prev_halfedge_handle(hedge_temp);
	} while (vertex_temp != vertex_init);

	double ratio_hw = (b_max[1] - b_min[0]) / (b_max[0] - b_min[0]);
	double quad_width = std::sqrt(para_area / ratio_hw);
	double quad_height = quad_width * ratio_hw;

	int quad_width_int = std::lround(quad_width * (my_scale / 2));
	int quad_height_int = std::lround(quad_height * (my_scale / 2));
	double quad_width_half = quad_width_int;
	double quad_height_half = quad_height_int;
	OpenMesh::Vec3d quad_point1(-quad_width_half, quad_height_half, 0);
	OpenMesh::Vec3d quad_point2(-quad_width_half, -quad_height_half, 0);
	OpenMesh::Vec3d quad_point3(quad_width_half, -quad_height_half, 0);
	OpenMesh::Vec3d quad_point4(quad_width_half, quad_height_half, 0);

	int corner1;
	int corner2;
	int corner3;
	int corner4;
	OpenMesh::Vec3d start_point(b_min[0], b_max[1], 0);
	std::vector<double> start_length;
	double boundary_length = 0;
	for (int i = 0; i < cc_boundary_.size(); i++)
	{
		OpenMesh::Vec3d point_temp = para_temp.point(cc_boundary_[i]);
		OpenMesh::Vec3d point_temp1 = para_temp.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
		start_length.push_back((point_temp - start_point).norm());
		boundary_length += (point_temp1 - point_temp).norm();
	}
	auto smallest = std::min_element(start_length.begin(), start_length.end());
	corner1 = std::distance(start_length.begin(), smallest);

	double boundary_length1 = boundary_length * (quad_height_half / (2 * (quad_height_half + quad_width_half)));
	double boundary_length2 = boundary_length * (quad_width_half / (2 * (quad_height_half + quad_width_half)));
	std::vector<int> edge_vertex1, edge_vertex2, edge_vertex3, edge_vertex4;
	double edge_length_temp = 0;
	for (int i = corner1; i < corner1 + cc_boundary_.size(); i++)
	{
		OpenMesh::Vec3d point_temp = para_temp.point(cc_boundary_[i%cc_boundary_.size()]);
		OpenMesh::Vec3d point_temp1 = para_temp.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
		edge_length_temp += (point_temp1 - point_temp).norm();
		edge_vertex1.push_back(i%cc_boundary_.size());
		if (edge_length_temp > boundary_length1)
		{
			corner2 = (i + 1) % cc_boundary_.size();
			break;
		}
	}
	edge_length_temp = 0;
	for (int i = corner2; i < corner2 + cc_boundary_.size(); i++)
	{
		OpenMesh::Vec3d point_temp = para_temp.point(cc_boundary_[i%cc_boundary_.size()]);
		OpenMesh::Vec3d point_temp1 = para_temp.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
		edge_length_temp += (point_temp1 - point_temp).norm();
		edge_vertex2.push_back(i%cc_boundary_.size());
		if (edge_length_temp > boundary_length2)
		{
			corner3 = (i + 1) % cc_boundary_.size();
			break;
		}
	}
	edge_length_temp = 0;
	for (int i = corner3; i < corner3 + cc_boundary_.size(); i++)
	{
		OpenMesh::Vec3d point_temp = para_temp.point(cc_boundary_[i%cc_boundary_.size()]);
		OpenMesh::Vec3d point_temp1 = para_temp.point(cc_boundary_[(i + 1) % cc_boundary_.size()]);
		edge_length_temp += (point_temp1 - point_temp).norm();
		edge_vertex3.push_back(i%cc_boundary_.size());
		if (edge_length_temp > boundary_length1)
		{
			corner4 = (i + 1) % cc_boundary_.size();
			break;
		}
	}
	int edge_vertex4_start = corner4;
	int edge_vertex4_end = corner1 > corner4 ? corner1 : corner1 + cc_boundary_.size();
	for (int i = edge_vertex4_start; i < edge_vertex4_end; i++)
	{
		edge_vertex4.push_back(i%cc_boundary_.size());
	}
	
	new_bpoint_.resize(cc_boundary_.size());
	for (int i = 0; i < edge_vertex1.size(); i++)
	{
		OpenMesh::Vec3d new_bpoint_temp = quad_point1 + (quad_point2 - quad_point1) / edge_vertex1.size()*i;
		new_bpoint_[edge_vertex1[i]] = new_bpoint_temp;
	}
	for (int i = 0; i < edge_vertex2.size(); i++)
	{
		OpenMesh::Vec3d new_bpoint_temp = quad_point2 + (quad_point3 - quad_point2) / edge_vertex2.size()*i;
		new_bpoint_[edge_vertex2[i]] = new_bpoint_temp;
	}
	for (int i = 0; i < edge_vertex3.size(); i++)
	{
		OpenMesh::Vec3d new_bpoint_temp = quad_point3 + (quad_point4 - quad_point3) / edge_vertex3.size()*i;
		new_bpoint_[edge_vertex3[i]] = new_bpoint_temp;
	}
	for (int i = 0; i < edge_vertex4.size(); i++)
	{
		OpenMesh::Vec3d new_bpoint_temp = quad_point4 + (quad_point1 - quad_point4) / edge_vertex4.size()*i;
		new_bpoint_[edge_vertex4[i]] = new_bpoint_temp;
	}

	OpenMesh::Vec2i ret_min(-quad_width_int, -quad_height_int);
	OpenMesh::Vec2i ret_max(quad_width_int, quad_height_int);
	retangle_type ret_temp;
	ret_temp.min = ret_min;
	ret_temp.max = ret_max;
	my_ret_candidate.push_back(ret_temp);
}

void ParaQuadCutting::uniform_tutte1(bool is_reconstruct_)
{
	my_origin_atlas.clear();
	for (int ii = 0; ii < para_atlas.size(); ii++)
	{
		Mesh& atlas_temp = para_atlas[ii];
		std::vector<OpenMesh::VertexHandle> cc_boundary;
		std::vector<OpenMesh::Vec3d> new_bpoint;
		if (!is_reconstruct_)
		{
			set_boundary(ii, cc_boundary, new_bpoint);
		}
		else
		{
			set_boundary1(ii, cc_boundary, new_bpoint);
		}

		int num = atlas_temp.n_vertices();
		Eigen::SparseMatrix<double> A(num, num);
		std::vector<Eigen::Triplet<double>> triplets;
		Eigen::MatrixXd b;
		for (auto v_it = atlas_temp.vertices_begin(); v_it != atlas_temp.vertices_end(); v_it++)
		{
			auto vertex = v_it.handle();
			if (atlas_temp.is_boundary(vertex))
			{
				triplets.push_back(Eigen::Triplet<double>(vertex.idx(), vertex.idx(), 1));
			}
			else
			{
				triplets.push_back(Eigen::Triplet<double>(vertex.idx(), vertex.idx(), 1));
				int neighbornum = 0;
				for (auto v_it = atlas_temp.vv_begin(vertex); v_it != atlas_temp.vv_end(vertex); v_it++)
				{
					neighbornum++;
				}
				for (auto v_it = atlas_temp.vv_begin(vertex); v_it != atlas_temp.vv_end(vertex); v_it++)
				{
					auto nvertex = v_it.handle();
					triplets.push_back(Eigen::Triplet<double>(vertex.idx(), nvertex.idx(), -1.0 / neighbornum));
				}
			}
		}
		A.setFromTriplets(triplets.begin(), triplets.end());
		A.makeCompressed();
		Eigen::SparseLU<Eigen::SparseMatrix<double>> LU;
		LU.compute(A);

		b.setZero(num, 3);
		for (int i = 0; i < cc_boundary.size(); i++)
		{
			b(cc_boundary[i].idx(), 0) = new_bpoint[i][0];
			b(cc_boundary[i].idx(), 1) = new_bpoint[i][1];
			b(cc_boundary[i].idx(), 2) = new_bpoint[i][2];
		}

		Eigen::MatrixXd sol(num, 3);
		sol = LU.solve(b);

		for (auto v_h : atlas_temp.vertices())
		{
			OpenMesh::Vec3d newposition;
			newposition[0] = sol(v_h.idx(), 0);
			newposition[1] = sol(v_h.idx(), 1);
			newposition[2] = sol(v_h.idx(), 2);
			atlas_temp.set_point(v_h, newposition);
		}


	}
}

double ParaQuadCutting::my_retangle_packing(std::vector<retangle_type> my_ret_candidate_)
{
	my_cf_candidate.clear();
	my_ct_candidate.clear();
	int n_charts = my_ret_candidate_.size();
	my_cf_candidate.assign(n_charts, 0);
	my_ct_candidate.assign(n_charts, OpenMesh::Vec2d(0.0, 0.0));

	constexpr bool allow_flip = true;
	const auto runtime_flipping_mode = rectpack2D::flipping_option::ENABLED;
	using spaces_type = rectpack2D::empty_spaces<allow_flip, rectpack2D::default_empty_spaces>;
	using rect_type = rectpack2D::output_rect_t<spaces_type>;

	auto report_successful = [](rect_type&)
	{
		return rectpack2D::callback_result::CONTINUE_PACKING;
	};

	auto report_unsuccessful = [](rect_type&)
	{
		return rectpack2D::callback_result::ABORT_PACKING;
	};

	const auto max_side = 8192;
	const auto discard_step = 1;

	std::vector<rect_type> rectangles;
	std::vector<OpenMesh::Vec2i> chart_offset(n_charts);
	std::vector<OpenMesh::Vec2i> chart_p_min(n_charts);

	double area_charts = 0;
	for (int i = 0; i < n_charts; i++)
	{
		area_charts += my_ret_candidate_[i].area();
		auto p_min = my_ret_candidate_[i].min;
		auto range = my_ret_candidate_[i].max - my_ret_candidate_[i].min;
	
		OpenMesh::Vec2i offset(5, 5);

		range += 2 * offset;

		rectangles.emplace_back(p_min[0], p_min[1], range[0], range[1], false);
		chart_offset[i] = offset;

		chart_p_min[i] = p_min;
	}

	auto report_result = [&]()
	{
		std::cout << "-------------------------------------" << std::endl;
		for (int i = 0; i < rectangles.size(); i++)
		{
			const auto& r = rectangles[i];
			std::cout << i << " " << (r.flipped ? "T" : "F") << " " << r.x + chart_offset[i][0] << " " << r.y + chart_offset[i][1] << " " << r.w - 2 * chart_offset[i][0] << " " << r.h - 2 * chart_offset[i][1] << std::endl;
		}
	};

	//report_result();
	const auto result_size = rectpack2D::find_best_packing<spaces_type>(rectangles, rectpack2D::make_finder_input(max_side, discard_step, report_successful, report_unsuccessful, runtime_flipping_mode));

	//std::cout << "-------------------------------------\nResultant bin: " << result_size.w << " " << result_size.h << std::endl;
	//report_result();

	retangle_type bb;
	for (int i = 0; i < n_charts; i++)
	{
		OpenMesh::Vec2i trans_i;
		OpenMesh::Vec2i packing_p_min, packing_p_max;
		const auto& r = rectangles[i];
		if (r.flipped)
		{
			OpenMesh::Vec2i flipped_p_min(-chart_p_min[i][1] - r.w + 2 * chart_offset[i][1], chart_p_min[i][0]);
			packing_p_min = OpenMesh::Vec2i(r.x + chart_offset[i][1], r.y + chart_offset[i][0]);
			packing_p_max = OpenMesh::Vec2i(r.x + r.w - chart_offset[i][1], r.y + r.h - chart_offset[i][0]);
			trans_i = packing_p_min - flipped_p_min;
		}
		else
		{
			packing_p_min = OpenMesh::Vec2i(r.x + chart_offset[i][0], r.y + chart_offset[i][1]);
			packing_p_max = OpenMesh::Vec2i(r.x + r.w - chart_offset[i][0], r.y + r.h - chart_offset[i][1]);
			trans_i = packing_p_min - chart_p_min[i];
		}
		if (r.flipped)
		{
			my_cf_candidate[i] = 1;
		}
		else
		{
			my_cf_candidate[i] = 0;
		}
		//my_cf_candidate[i] = r.flipped;

		my_ct_candidate[i][0] = double(trans_i[0]);
		my_ct_candidate[i][1] = double(trans_i[1]);

		bb.min = OpenMesh::minimize(bb.min, packing_p_min);

		bb.max = OpenMesh::maximize(bb.max, packing_p_max);
	}
	//std::cout << "bb max===============" << bb.max[0] << " " << bb.min[0] << std::endl;
	//	if (bb.width() > 3 * bb.height() || bb.height() > 3 * bb.width()) isgood = false;

	return double(area_charts) / double(bb.area());
}

void ParaQuadCutting::find_new_cut(bool my_switch)
{
	std::vector<Mesh> my_para_atlas0 = my_para_atlas;
	int count_num = 0;
	std::vector<double> ratio_list;
	std::vector<int> ratio_flag;
	std::vector<OpenMesh::Vec3d> B_max_list;
	std::vector<OpenMesh::Vec3d> B_min_list;
	for (int i = 0; i < my_para_atlas0.size(); i++)
	{
		OpenMesh::Vec3d b_max, b_min;
		find_BB(my_para_atlas0[i], b_max, b_min);
		double edge_length1 = b_max[1] - b_min[1];
		double edge_length2 = b_max[0] - b_min[0];
		double ratio_temp = edge_length1 > edge_length2 ? (edge_length1 / edge_length2) : (edge_length2 / edge_length1);
		//double ratio_temp = edge_length1 > edge_length2 ? edge_length1 : edge_length2;
		double flag_temp;
		if (!my_switch)
		{
			flag_temp = edge_length1 > edge_length2 ? 1 : 0;
		}
		else
		{
			flag_temp = edge_length1 > edge_length2 ? 0 : 1;
		}

		ratio_list.push_back(ratio_temp);
		ratio_flag.push_back(flag_temp);
		B_max_list.push_back(b_max);
		B_min_list.push_back(b_min);
	}
	std::vector<double>::iterator biggest;
	int cut_id;
	if (!my_switch)
	{
		biggest = std::max_element(ratio_list.begin(), ratio_list.end());
		cut_id = std::distance(ratio_list.begin(), biggest);
	}
	else
	{
		biggest = std::min_element(ratio_list.begin(), ratio_list.end());
		cut_id = std::distance(ratio_list.begin(), biggest);
	}

	OpenMesh::Vec3d middle_point;
	OpenMesh::Vec3d middle_point1;
	if (ratio_flag[cut_id] == 0)
	{
		middle_point[0] = (B_max_list[cut_id][0] + B_min_list[cut_id][0]) / 2;
		middle_point[1] = B_max_list[cut_id][1];
		middle_point[2] = 0;
		middle_point1[0] = (B_max_list[cut_id][0] + B_min_list[cut_id][0]) / 2;
		middle_point1[1] = B_min_list[cut_id][1];
		middle_point1[2] = 0;
	}
	if (ratio_flag[cut_id] == 1)
	{
		middle_point[0] = B_max_list[cut_id][0];
		middle_point[1] = (B_max_list[cut_id][1] + B_min_list[cut_id][1]) / 2;
		middle_point[2] = 0;
		middle_point1[0] = B_min_list[cut_id][0];
		middle_point1[1] = (B_max_list[cut_id][1] + B_min_list[cut_id][1]) / 2;
		middle_point1[2] = 0;
	}

	std::vector<OpenMesh::Vec3d> intersect_point;
	std::vector<double> intersect_distance;
	OpenMesh::Vec3d center_point = (middle_point + middle_point1) / 2;
	for (auto he_h : my_para_atlas0[cut_id].halfedges())
	{
		if (my_para_atlas0[cut_id].is_boundary(he_h))
		{
			OpenMesh::VertexHandle fvertex = my_para_atlas0[cut_id].from_vertex_handle(he_h);
			OpenMesh::VertexHandle tvertex = my_para_atlas0[cut_id].to_vertex_handle(he_h);
			OpenMesh::Vec3d fpoint = my_para_atlas0[cut_id].point(fvertex);
			OpenMesh::Vec3d tpoint = my_para_atlas0[cut_id].point(tvertex);
			OpenMesh::Vec3d vec0 = middle_point1 - middle_point;
			OpenMesh::Vec3d vec1 = fpoint - middle_point;
			OpenMesh::Vec3d vec2 = tpoint - middle_point;
			bool is_intersect_temp = 0;
			if ((vec0%vec1)[2] * (vec0%vec2)[2] <= 0)
			{
				is_intersect_temp = 1;
			}
			if (is_intersect_temp)
			{
				intersect_point.push_back(fpoint);
				intersect_distance.push_back((fpoint - center_point).norm());
			}
		}
	}
	if (intersect_distance.size() >= 2)
	{
		auto smallest_d = std::min_element(intersect_distance.begin(), intersect_distance.end());
		int smallest_p = std::distance(intersect_distance.begin(), smallest_d);
		OpenMesh::Vec3d new_segment1_p = intersect_point[smallest_p];
		intersect_point.erase(intersect_point.begin() + smallest_p);
		intersect_distance.erase(intersect_distance.begin() + smallest_p);
		OpenMesh::Vec3d direction1 = (new_segment1_p - center_point).normalize();
		OpenMesh::Vec3d new_segment2_p;
		while (intersect_point.size() > 0)
		{
			auto smallest_d1 = std::min_element(intersect_distance.begin(), intersect_distance.end());
			int smallest_p1 = std::distance(intersect_distance.begin(), smallest_d1);
			OpenMesh::Vec3d new_segment2_p_temp = intersect_point[smallest_p1];
			OpenMesh::Vec3d direction2 = (new_segment2_p_temp - center_point).normalize();
			if ((direction1 | direction2) < 0)
			{
				new_segment2_p = new_segment2_p_temp;
				break;
			}
			intersect_point.erase(intersect_point.begin() + smallest_p1);
			intersect_distance.erase(intersect_distance.begin() + smallest_p1);
		}

		std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d> new_segment_pair;
		new_segment_pair.first = new_segment1_p;
		new_segment_pair.second = new_segment2_p;
		my_cut_segments_new.push_back(new_segment_pair);
	}
	else
	{
		return;
	}
}

void ParaQuadCutting::my_packing()
{
	std::vector<int> nv_offset(para_atlas.size());

	nv_offset[0] = 0;
	for (int i = 0; i < para_atlas.size() - 1; i++)
	{
		nv_offset[i + 1] = nv_offset[i] + para_atlas[i].n_vertices();
	}

	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
	OpenMesh::HPropHandleT<int> hvt_index;

	new_mesh.add_property(mvt_list, "mvt_list");
	new_mesh.add_property(hvt_index, "hvt_index");

	new_mesh.property(mvt_list).resize(nv_offset.back() + para_atlas.back().n_vertices());

	for (int i = 0; i < para_atlas.size(); i++)
	{
		for (int j = 0; j < para_atlas[i].n_vertices(); j++)
		{
			const auto& point_v = para_atlas[i].point(para_atlas[i].vertex_handle(j));		
			new_mesh.property(mvt_list)[nv_offset[i] + j][0] = point_v[0];
			new_mesh.property(mvt_list)[nv_offset[i] + j][1] = point_v[1];
		}
	}

	for (int i = 0; i < new_mesh.n_faces(); i++)
	{
		auto f_h = new_mesh.face_handle(i);
		int layer_id = face_layer[i];

		for (auto fh_h = new_mesh.cfh_begin(f_h); fh_h != new_mesh.cfh_end(f_h); fh_h++)
		{
			auto h_para = get_mesh2para(*fh_h);
			if (!h_para.is_valid())
			{
				std::cout << new_mesh.from_vertex_handle(*fh_h).idx() << " " << new_mesh.to_vertex_handle(*fh_h).idx() << " " << i << std::endl;
				continue;
			}

			new_mesh.property(hvt_index, *fh_h) = v_para2layer[layer_id][new_para.to_vertex_handle(h_para).idx()] + nv_offset[layer_id];
		}
	}

	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list_para;
	OpenMesh::HPropHandleT<int> hvt_index_para;

	new_para.add_property(mvt_list_para, "mvt_list");
	new_para.add_property(hvt_index_para, "hvt_index");

	new_para.property(mvt_list_para) = new_mesh.property(mvt_list);
	for (int i = 0; i < new_para.n_faces(); i++)
	{
		auto f_h = new_para.face_handle(i);
		int layer_id = face_layer[i];

		for (auto fh_h : new_para.fh_range(f_h))
		{
			new_para.property(hvt_index_para, fh_h) = v_para2layer[layer_id][new_para.to_vertex_handle(fh_h).idx()] + nv_offset[layer_id];
		}
	}
}

void ParaQuadCutting::find_BB(Mesh mesh_, OpenMesh::Vec3d &b_max_, OpenMesh::Vec3d &b_min_)
{
	OpenMesh::Vec3d b_max;
	OpenMesh::Vec3d b_min;
	b_max = b_min = mesh_.point(mesh_.vertex_handle(0));
	for (auto v_h : mesh_.vertices())
	{
		const auto& p = mesh_.point(v_h);
		b_max = b_max.maximize(p);
		b_min = b_min.minimize(p);
	}
	b_max_ = b_max;
	b_min_ = b_min;
}

void ParaQuadCutting::my_get_scaf_info(Eigen::MatrixXd& v_pos, Eigen::MatrixXd& uv_v_pos, Eigen::MatrixXi& fv_id, Eigen::MatrixXi& uv_fv_id, Mesh para_, Mesh mesh_)
{
	v_pos.resize(mesh_.n_vertices(), 3);
	uv_v_pos.resize(para_.n_vertices(), 2);
	fv_id.resize(mesh_.n_faces(), 3);
	uv_fv_id.resize(para_.n_faces(), 3);

	for (int i = 0; i < mesh_.n_vertices(); i++)
	{
		const auto& pt = mesh_.point(mesh_.vertex_handle(i));
		v_pos(i, 0) = pt[0];
		v_pos(i, 1) = pt[1];
		v_pos(i, 2) = pt[2];
	}

	for (int i = 0; i < para_.n_vertices(); i++)
	{
		const auto& pt = para_.point(para_.vertex_handle(i));
		uv_v_pos(i, 0) = pt[0];
		uv_v_pos(i, 1) = pt[1];
	}

	for (int i = 0; i < mesh_.n_faces(); i++)
	{
		int j = 0;
		for (auto vf : mesh_.fv_range(mesh_.face_handle(i))) fv_id(i, j++) = vf.idx();
	}

	for (int i = 0; i < para_.n_faces(); i++)
	{
		int j = 0;
		for (auto vf : para_.fv_range(para_.face_handle(i))) uv_fv_id(i, j++) = vf.idx();
	}
}

void my_angle_func(const size_t N, const std::vector<double>& x, double& f, std::vector<double>& g, void* user_supply)
{
	auto ptr_this = static_cast<ParaQuadCutting*>(user_supply);
	ptr_this->my_rotating_angle_evalfunc(x, f, g);
};

void ParaQuadCutting::build_rotation_message(Mesh mesh_)
{
	my_uv_x.clear();
	my_boundary_h_svec.clear();
	my_boundary_h_len0.clear();
	n_vertices = mesh_.n_vertices();
	my_uv_x.resize(2 * n_vertices);
	for (auto v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); v_it++)
	{
		my_uv_x[2 * v_it.handle().idx()] = mesh_.point(v_it.handle())[0];
		my_uv_x[2 * v_it.handle().idx() + 1] = mesh_.point(v_it.handle())[1];
	}

	std::vector<std::pair<int, int>> boundary_h_vert;
	std::vector<OpenMesh::HalfedgeHandle> he_boundary;
	OpenMesh::HalfedgeHandle hedge_init;
	OpenMesh::HalfedgeHandle hedge_temp;
	for (auto he_h : mesh_.halfedges())
	{
		if (mesh_.is_boundary(he_h))
		{
			hedge_init = he_h;
			break;
		}
	}
	hedge_temp = hedge_init;
	do
	{
		he_boundary.push_back(hedge_temp);
		hedge_temp = mesh_.next_halfedge_handle(hedge_temp);
	} while (hedge_temp != hedge_init);

	n_boundary_edges = he_boundary.size();
	//my_boundary_h_svec.resize(n_boundary_edges);
	//my_boundary_h_len0.resize(n_boundary_edges);
	double total_boundary_length = 0;
	for (int i = 0; i < n_boundary_edges; i++)
	{
		OpenMesh::HalfedgeHandle he_temp = he_boundary[i];
		OpenMesh::Vec3d fpoint = mesh_.point(mesh_.from_vertex_handle(he_temp));
		OpenMesh::Vec3d tpoint = mesh_.point(mesh_.to_vertex_handle(he_temp));
		OpenMesh::Vec3d vec_temp3 = tpoint - fpoint;
		OpenMesh::Vec2d vec_temp2(vec_temp3[0], vec_temp3[1]);
		my_boundary_h_svec.push_back(vec_temp2);
		my_boundary_h_len0.push_back(vec_temp2.norm());
		total_boundary_length += vec_temp2.norm();
		std::pair<int, int> pair_temp;
		pair_temp.first = mesh_.from_vertex_handle(he_temp).idx();
		pair_temp.second = mesh_.to_vertex_handle(he_temp).idx();
		boundary_h_vert.push_back(pair_temp);
	}
	
	double kernel_width = 4 * (double)n_boundary_edges / 1000.0;
	kernel_width = std::min(std::max(0.0, kernel_width), 20.0);
	double avg_len = total_boundary_length / n_boundary_edges;
	double gaussian_delta = kernel_width * avg_len;
	double gaussian_thres = 3.0 * gaussian_delta;
	double delta_sqr = 2.0 * gaussian_delta * gaussian_delta;
	for (int i = 0; i < 3; i++)
	{
		std::vector<OpenMesh::Vec2d> boundary_svec(n_boundary_edges);
		for (int j = 0; j < n_boundary_edges; j++)
		{
			OpenMesh::Vec2d smooth_vec = my_boundary_h_svec[j];
			for (int step : {1, -1})
			{
				int cur_id = j;
				double dist = my_boundary_h_len0[cur_id] / 2.0;
				while (dist < gaussian_thres)
				{
					int x = cur_id + step;
					cur_id = (x < 0) ? (x + n_boundary_edges) : (x >= n_boundary_edges ? x - n_boundary_edges : x);
					//cur_id = boundary(cur_id + step);
					int vert_h = (step == 1) ? boundary_h_vert[cur_id].first : boundary_h_vert[cur_id].second;

					double alpha = std::min((gaussian_thres - dist) / my_boundary_h_len0[cur_id], 1.0);
					double dist_b = dist + my_boundary_h_len0[cur_id] * alpha / 2.0;

					smooth_vec += ig::FastNegExp1(dist_b * dist_b / delta_sqr) * my_boundary_h_svec[cur_id] * alpha;

					dist += my_boundary_h_len0[cur_id];
					
				}
			}
	
			if (OpenMesh::dot(smooth_vec, my_boundary_h_svec[j]) < 0.0)
			{
				boundary_svec[j] = my_boundary_h_svec[j];
			}
			else
			{
				int x = j + 1;
				int k = (x < 0) ? (x + n_boundary_edges) : (x >= n_boundary_edges ? x - n_boundary_edges : x);
				//int k = boundary(j + 1);
				double d_angle = CommonFunctions::vec_angle_atan2(my_boundary_h_svec[j], smooth_vec);

				//boundary_v_sangle[j] += d_angle;
				//boundary_v_sangle[k] -= d_angle;
				boundary_svec[j] = smooth_vec.normalized() * my_boundary_h_svec[j].norm();
			}
		}
		my_boundary_h_svec = std::move(boundary_svec);
	}

	for (int j = 0; j < n_boundary_edges; j++) my_boundary_h_svec[j].normalize();
}

void ParaQuadCutting::my_calc_global_rotation()
{
	HLBFGS angle_solver;
	angle_solver.set_number_of_variables(1);
	angle_solver.set_verbose(false);
	angle_solver.set_func_callback(my_angle_func, 0, 0, 0, 0);

	double global_theta = 0.0;
	double theta[1] = { global_theta };
	angle_solver.optimize_without_constraints(theta, 100, this);
	global_theta = theta[0];

	OpenMesh::Vec2d chart_center(0.0, 0.0);
	for (int i = 0; i < n_vertices; i++)
	{
		chart_center[0] += my_uv_x[2 * i + 0];
		chart_center[1] += my_uv_x[2 * i + 1];
	}
	chart_center /= n_vertices;

	double cos_theta = std::cos(global_theta);
	double sin_theta = std::sin(global_theta);

	for (int i = 0; i < n_vertices; i++)
	{
		OpenMesh::Vec2d point(my_uv_x[2 * i + 0], my_uv_x[2 * i + 1]);
		point -= chart_center;
		my_rotate_uv(point, cos_theta, sin_theta);
		point += chart_center;

		my_uv_x[2 * i + 0] = point[0];
		my_uv_x[2 * i + 1] = point[1];
	}

	for (int i = 0; i < n_boundary_edges; i++)
	{
		my_rotate_uv(my_boundary_h_svec[i], cos_theta, sin_theta);
	}
}

void ParaQuadCutting::my_rotating_angle_evalfunc(const std::vector<double>& x, double& f, std::vector<double>& g)
{
	double cos_theta_t = cos(x[0]);
	double sin_theta_t = sin(x[0]);

	f = 0.0;
	g[0] = 0.0;

	for (int i = 0; i < n_boundary_edges; i++)
	{
		OpenMesh::Vec2d h_vec = my_boundary_h_svec[i];
		double length = my_boundary_h_len0[i];

		my_rotate_uv(h_vec, cos_theta_t, sin_theta_t);

		f += length * h_vec[0] * h_vec[0] * h_vec[1] * h_vec[1];
		g[0] += 2.0 * length * h_vec[0] * h_vec[1] * (h_vec[1] * h_vec[1] - h_vec[0] * h_vec[0]);
	}
}
