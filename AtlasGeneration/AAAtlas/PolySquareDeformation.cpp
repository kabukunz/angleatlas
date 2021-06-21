#include <iostream>
#include <queue>
#include "PolySquareDeformation.h"
#include "ChartDeformation.h"
#include "MeshViewer/Mesh_doubleIO.h"

#include <fstream>
#include <QTime>
#include <QDir>

OpenMesh::Vec3uc tag_colors[4] =
{
	{ 255,  41,  30 },
	{ 229, 198,  56 },
	{ 6, 216,  20 },
	{ 0, 237, 255 }
};

PolySquareDeformation::PolySquareDeformation(const Mesh& _mesh, Mesh& _para, std::map<int, int>& _vk)
	: mesh(_mesh), para(_para), viewer_vk(_vk)
{
	vertex_uv_vec.clear();

	std::vector<OpenMesh::Vec2d> vertex_uv(_para.n_vertices());

	for (auto v_h : _para.vertices())
	{
		vertex_uv[v_h.idx()][0] = _para.point(v_h)[0];
		vertex_uv[v_h.idx()][1] = _para.point(v_h)[1];
	}

	vertex_uv_vec.emplace_back(std::move(vertex_uv));
}

PolySquareDeformation::~PolySquareDeformation()
{
}

bool PolySquareDeformation::calc(double kernel, double goal, double exp, std::string str_)
{
	kernel_width = kernel;
//	goal_factor = goal;
	goal_length = goal;
	energy_exp_factor = exp;
	amips_exp = (exp > 0.0);

	QTime timer;
	timer.start();

	if (!load_face_info()) return false;

	build_boundary_cycles();
	if (atlas.size() != para.n_vertices() + para.n_faces() - para.n_edges())
	{
		std::cout << "Charts have holes!" << std::endl;
		return false;
	}

	build_charts_info();

	set_uv();
	for (auto& chart : atlas) chart.deformation(str_);
	get_uv();

	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "Total Time " << timer.elapsed() / 1000.0 << "s" << std::endl;

 	check_flip(vertex_uv_vec.back());
// 	color_edges();
 
 	gather_vk();
	return true;
}

void PolySquareDeformation::gather_vk()
{
	viewer_vk.clear();
	for (int i = 0; i < atlas.size(); i++)
	{
		int euler = 0;
		for (int j = 0; j < atlas[i].n_boundary_edges; j++)
		{
			int cur = atlas[i].boundary_v_k[j];
			int v0 = para.from_vertex_handle(para.halfedge_handle(atlas[i].boundary_h_meshid[j])).idx();
			if (cur == 0) continue;
			viewer_vk.emplace(v0, cur);
			euler += cur;
		}
		if (euler != 4) std::cout << "Euler Chart " << euler << std::endl;
	}
}

bool PolySquareDeformation::load_face_info()
{
	mesh_face_info.resize(para.n_faces());

	for (auto f_h : mesh.faces())
	{
		auto fn = para.normal(f_h);
		bool z_pos = para.normal(f_h)[2] >= 0;
		auto& f_info = mesh_face_info[f_h.idx()];

		f_info.normal_towards = z_pos ? 1 : 0;

		OpenMesh::Vec2d uv_p1, uv_p2;
		calc_face_vectors(f_h.idx(), uv_p1, uv_p2, vertex_uv_vec[0]);

		if (uv_p1[0] * uv_p2[1] - uv_p1[1] * uv_p2[0] < 1.0e-30)
		{
			std::cout << "Input is incorrect! " << f_info.det_p << std::endl;
			return false;
		}

		OpenMesh::Vec3d mesh_p[3];
		auto fh_iter = mesh.cfh_begin(f_h);
		for (int i = 0; i < 3; i++, fh_iter++)
		{
			mesh_p[i] = mesh.point(mesh.to_vertex_handle(*fh_iter));
		}

		mesh_p[1] -= mesh_p[0];
		mesh_p[2] -= mesh_p[0];

		f_info.l2_p1 = mesh_p[1].sqrnorm();
		f_info.l2_p2 = mesh_p[2].sqrnorm();
		f_info.dot_p = OpenMesh::dot(mesh_p[1], mesh_p[2]);
		f_info.det_p = OpenMesh::cross(mesh_p[1], mesh_p[2]).norm();
	}

	return true;
}

void PolySquareDeformation::calc_face_vectors(int f, OpenMesh::Vec2d& p1, OpenMesh::Vec2d& p2, const std::vector<OpenMesh::Vec2d>& points)
{
	auto& f_info = mesh_face_info[f];

	auto f_h = para.face_handle(f);

	int vid[3];
	OpenMesh::Vec2d vec[3];

	auto fh_iter = para.cfh_begin(f_h);
	for (int i = 0; i < 3; i++, fh_iter++)
	{
		vid[i] = para.to_vertex_handle(*fh_iter).idx();
		vec[i] = points[vid[i]];
	}

	f_info.uv0 = vid[0];
	f_info.uv1 = vid[1];
	f_info.uv2 = vid[2];

	p1 = vec[1] - vec[0];
	p2 = vec[2] - vec[0];
}

void PolySquareDeformation::build_boundary_cycles()
{
	std::vector<bool> e_visited(para.n_edges(), false);

	for (auto e_h : para.edges())
	{
		if (!para.is_boundary(e_h) || e_visited[e_h.idx()]) continue;

		auto h0 = para.is_boundary(para.halfedge_handle(e_h, 0)) ? para.halfedge_handle(e_h, 0) : para.halfedge_handle(e_h, 1);
		int cycle_id = atlas.size();

		atlas.emplace_back(*this);
		auto& chart = atlas.back();

		auto h_iter = h0;
		do
		{
			h_boundary[h_iter.idx()] = std::make_pair(cycle_id, int(chart.boundary_h_meshid.size()));

			chart.boundary_h_meshid.emplace_back(h_iter.idx());
			chart.boundary_h_len0.emplace_back(para.calc_edge_length(h_iter));

			e_visited[h_iter.idx() / 2] = true;

			h_iter = para.next_halfedge_handle(h_iter);
		} while (h_iter != h0);
	}
}

void PolySquareDeformation::build_charts_info()
{
	v_chart.assign(para.n_vertices(), { -1, -1 });
	f_chart.assign(para.n_faces(), -1);

	for (int i = 0; i < atlas.size(); i++)
	{
		atlas[i].build_chart();
		auto f0 = para.opposite_face_handle(para.halfedge_handle(atlas[i].boundary_h_meshid[0]));

		std::queue<int> f_bfs;
		f_bfs.push(f0.idx());
		while (!f_bfs.empty())
		{
			auto f_cur = para.face_handle(f_bfs.front());
			f_bfs.pop();

			if (f_chart[f_cur.idx()] != -1) continue;

			for (auto fv_h : para.fv_range(f_cur))
			{
				v_chart[fv_h.idx()].first = i;
			}
			f_chart[f_cur.idx()] = i;

			for (auto fh_h : para.fh_range(f_cur))
			{
				auto f_next = para.opposite_face_handle(fh_h);
				if (f_next.is_valid() && f_chart[f_next.idx()] == -1)
				{
					f_bfs.push(f_next.idx());
				}
			}
		}
	}

	for (int i = 0; i < para.n_vertices(); i++)
	{
		v_chart[i].second = atlas[v_chart[i].first].n_vertices;
		atlas[v_chart[i].first].n_vertices++;
	}
	for (int i = 0; i < atlas.size(); i++)
	{
		atlas[i].build_vertices();
	}

	for (auto f_h : para.faces())
	{
		atlas[f_chart[f_h.idx()]].n_faces++;
		atlas[f_chart[f_h.idx()]].total_area += face_area(f_h.idx());
	}
	for (int i = 0; i < atlas.size(); i++)
	{
		atlas[i].mesh_faces.reserve(atlas[i].n_faces);
	}
	for (int i = 0; i < para.n_faces(); i++)
	{
		atlas[f_chart[i]].mesh_faces.push_back(i);
	}
	
	for (int i = 0; i < atlas.size(); i++)
	{
		atlas[i].goal_length = goal_length;
	}
}

void PolySquareDeformation::set_uv()
{
	const auto& vertex_uv = vertex_uv_vec[0];
	for (int i = 0; i < para.n_vertices(); i++)
	{
		int cid = v_chart[i].first;
		int vid = v_chart[i].second;

		atlas[cid].uv_x[2 * vid + 0] = vertex_uv[i][0];
		atlas[cid].uv_x[2 * vid + 1] = vertex_uv[i][1];
	}
}

void PolySquareDeformation::get_uv()
{
	size_t n_phase = 1;
	for (const auto& chart : atlas)
	{
		n_phase = std::max(n_phase, chart.uv_x_vec.size() + 1);
	}
	
	vertex_uv_vec.resize(n_phase);
	for (int i = 1; i < vertex_uv_vec.size(); i++)
	{
		vertex_uv_vec[i].resize(para.n_vertices());
		for (int j = 0; j < para.n_vertices(); j++)
		{
			int cid = v_chart[j].first;
			int vid = v_chart[j].second;
			
			int vec_id;
			int vec_size = atlas[cid].uv_x_vec.size();

			if (i < vec_size - 2) 
				vec_id = i - 1;
			else if (i >= n_phase - 3) 
				vec_id = vec_size + i - n_phase;
			else 
				vec_id = vec_size - 4;

			vertex_uv_vec[i][j][0] = atlas[cid].uv_x_vec[vec_id][2 * vid + 0];
			vertex_uv_vec[i][j][1] = atlas[cid].uv_x_vec[vec_id][2 * vid + 1];
		}
	}
}

void PolySquareDeformation::check_flip(const std::vector<OpenMesh::Vec2d>& vertex_uv)
{
	int n_flip = 0;
	for (auto f_h : mesh.faces())
	{
		auto& f_info = mesh_face_info[f_h.idx()];
		bool z_pos = (f_info.normal_towards == 1);

		f_info.normal_towards = z_pos ? 1 : 0;

		OpenMesh::Vec2d p1, p2;
		calc_face_vectors(f_h.idx(), p1, p2, vertex_uv);

		double det_p = p1[0] * p2[1] - p1[1] * p2[0];

		if (z_pos != (det_p > 0))
		{
			if (n_flip < 20) std::cout << "Flip at " << f_h.idx() << std::endl;
			n_flip++;
		}
	}

	if (n_flip > 20) std::cout << "Flips " << n_flip << std::endl;
}

void PolySquareDeformation::set_vertices(const std::vector<OpenMesh::Vec2d>& vertex_uv)
{
	for (auto v_h : para.vertices())
	{
		para.point(v_h)[0] = vertex_uv[v_h.idx()][0];
		para.point(v_h)[1] = vertex_uv[v_h.idx()][1];
	}
}

void PolySquareDeformation::color_edges()
{
	for (int i = 0; i < atlas.size(); i++)
	{
		auto& chart = atlas[i];

		for (int j = 0; j < chart.boundary_h_meshid.size(); j++)
		{
			para.set_color(para.edge_handle(chart.boundary_h_meshid[j] / 2), tag_colors[chart.boundary_h_tag[j]]);
		}
	}
}