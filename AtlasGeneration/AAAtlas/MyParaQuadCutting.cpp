#include "MyParaQuadCutting.h"
#include "MeshViewer/Mesh_doubleIO.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <queue>

#include <Eigen/Eigen>
#include <Eigen/CholmodSupport>

#include <QTime>

MyParaQuadCutting::MyParaQuadCutting(Mesh& _mesh, const char* path, const char* file)
	:origin_mesh(_mesh)
{
	str_path = std::string(path);
	str_file = std::string(file);
	init();
}

void MyParaQuadCutting::init()
{

}

bool MyParaQuadCutting::get_para_mesh()
{
	//origin_para.clear();
	//new_para.clear();

	//OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
	//OpenMesh::HPropHandleT<int> hvt_index;

	//if (!origin_mesh.get_property_handle(mvt_list, "mvt_list") || !origin_mesh.get_property_handle(hvt_index, "hvt_index"))
	//{
	//	std::cout << "Texture data is invalid." << std::endl;
	//	return false;
	//}

	//for (int i = 0; i < origin_mesh.property(mvt_list).size(); i++)
	//{
	//	origin_para.add_vertex(Mesh::Point(origin_mesh.property(mvt_list)[i][0], origin_mesh.property(mvt_list)[i][1], 0.0));
	//	new_para.add_vertex(Mesh::Point(origin_mesh.property(mvt_list)[i][0], origin_mesh.property(mvt_list)[i][1], 0.0));
	//}

	//for (int i = 0; i < origin_mesh.n_faces(); i++)
	//{
	//	auto f_h = origin_mesh.face_handle(i);
	//	std::vector<OpenMesh::VertexHandle> para_f;
	//	para_f.reserve(3);

	//	for (auto fh_h = origin_mesh.cfh_begin(f_h); fh_h != origin_mesh.cfh_end(f_h); fh_h++)
	//	{
	//		para_f.push_back(origin_para.vertex_handle(origin_mesh.property(hvt_index, *fh_h)));
	//	}

	//	origin_para.add_face(para_f);
	//	new_para.add_face(para_f);
	//}

	//origin_h_mesh2para.resize(origin_mesh.n_halfedges(), -1);
	//new_mesh.add_property(h_mesh2para);
	//new_para.add_property(h_para2mesh);
	//for (int i = 0; i < origin_mesh.n_faces(); i++)
	//{
	//	auto f0 = origin_mesh.face_handle(i);
	//	auto f1 = new_mesh.face_handle(i);
	//	auto f2 = new_para.face_handle(i);

	//	std::map<int, int> v_para2mesh;
	//	for (auto f0_h : origin_mesh.fh_range(f0))
	//	{
	//		v_para2mesh[origin_mesh.property(hvt_index, f0_h)] = origin_mesh.to_vertex_handle(f0_h).idx();
	//	}

	//	std::map<int, int> mesh_v2h;
	//	for (auto f1_h : new_mesh.fh_range(f1))
	//	{
	//		mesh_v2h[new_mesh.to_vertex_handle(f1_h).idx()] = f1_h.idx();
	//	}

	//	for (auto f2_h : new_para.fh_range(f2))
	//	{
	//		auto f1_h = new_mesh.halfedge_handle(mesh_v2h[v_para2mesh[new_para.to_vertex_handle(f2_h).idx()]]);

	//		new_mesh.property(h_mesh2para, f1_h) = std::make_pair(new_para.from_vertex_handle(f2_h).idx(), new_para.to_vertex_handle(f2_h).idx());
	//		new_para.property(h_para2mesh, f2_h) = std::make_pair(new_mesh.from_vertex_handle(f1_h).idx(), new_mesh.to_vertex_handle(f1_h).idx());

	//		origin_h_mesh2para[f1_h.idx()] = f2_h.idx();
	//	}
	//}

	//cut_length = 0.0;
	//OpenMesh::EPropHandleT<bool> e_oncut;
	//origin_mesh.get_property_handle(e_oncut, "e_oncut");
	//for (auto e_h : origin_mesh.edges())
	//{
	//	if (origin_mesh.is_boundary(e_h))
	//	{
	//		cut_length += origin_mesh.calc_edge_length(e_h);
	//		continue;
	//	}

	//	int to0 = origin_mesh.property(hvt_index, origin_mesh.halfedge_handle(e_h, 0));
	//	int to1 = origin_mesh.property(hvt_index, origin_mesh.halfedge_handle(e_h, 1));
	//	int from0 = origin_mesh.property(hvt_index, origin_mesh.prev_halfedge_handle(origin_mesh.halfedge_handle(e_h, 0)));
	//	int from1 = origin_mesh.property(hvt_index, origin_mesh.prev_halfedge_handle(origin_mesh.halfedge_handle(e_h, 1)));

	//	bool is_cut = ((to0 != from1) || (to1 != from0));
	//	origin_mesh.property(e_oncut, e_h) = is_cut;
	//	if (is_cut) cut_length += origin_mesh.calc_edge_length(e_h) * 2.0;
	//}

	//origin_mesh.request_vertex_normals();

	//new_para.request_edge_status();
	//new_para.request_face_status();

	//new_mesh.request_edge_status();
	//new_mesh.request_face_status();

	//origin_para.request_face_normals();
	//origin_para.update_face_normals();

	//origin_para.request_edge_colors();
	//for (auto e_h : origin_para.edges())
	//{
	//	origin_para.set_color(e_h, origin_para.is_boundary(e_h) ? boundary_edge_color : interior_edge_color);
	//}
	////	std::cout << "Euler of Para : " << origin_para.n_vertices() + origin_para.n_faces() - origin_para.n_edges() << std::endl;

	//BB_Max = BB_Min = origin_para.point(origin_para.vertex_handle(0));
	//for (auto v_h : origin_para.vertices())
	//{
	//	const auto& p = origin_para.point(v_h);
	//	BB_Max = BB_Max.maximize(p);
	//	BB_Min = BB_Min.minimize(p);
	//}
	////	std::cout << "BB of Para : X [" << BB_Min[0] << ", " << BB_Max[0] << "] Y [" << BB_Min[1] << ", " << BB_Max[1] << "]" << std::endl;
	return true;
}

void MyParaQuadCutting::read_cut()
{
	std::string instr = str_file.substr(0, str_file.find_last_of('.')) + ".txt";
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
		cut_segments.push_back(cut_segments_temp);
	}
}

bool MyParaQuadCutting::is_point_on_edge(OpenMesh::Vec3d p_, OpenMesh::EdgeHandle edge_, Mesh &mesh_)
{
	OpenMesh::HalfedgeHandle he = mesh_.halfedge_handle(edge_, 0);
	OpenMesh::VertexHandle fromv = mesh_.from_vertex_handle(he);
	OpenMesh::VertexHandle tov = mesh_.to_vertex_handle(he);
	OpenMesh::Vec3d q1 = mesh_.point(fromv);
	OpenMesh::Vec3d q2 = mesh_.point(tov);
	if ((p_[0] - q1[0])*(q2[1] - q1[1]) == (q2[0] - q1[0])*(p_[1] - q1[1])
		&& std::min(q1[0], q2[0]) <= p_[0] && p_[0] <= std::max(q1[0], q2[0])
		&& std::min(q1[1], q2[1]) <= p_[1] && p_[1] <= std::max(q1[1], q2[1]))
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool MyParaQuadCutting::is_point_in_face(OpenMesh::Vec3d p_, OpenMesh::FaceHandle face_, Mesh &mesh_)
{
	std::vector<OpenMesh::Vec3d> face_vertex;
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

	return u + v < 1;
}

OpenMesh::HalfedgeHandle MyParaQuadCutting::find_split_uv(OpenMesh::Vec3d p_, OpenMesh::FaceHandle face_, Mesh &mesh_,double &u, double &v)
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

double MyParaQuadCutting::find_split_t(OpenMesh::Vec3d p_, OpenMesh::HalfedgeHandle he_)
{
	const auto& uv0 = new_para.point(new_para.from_vertex_handle(he_));
	const auto& uv1 = new_para.point(new_para.to_vertex_handle(he_));

	auto v_e = uv1 - uv0;
	auto v_p = p_ - uv0;

	return OpenMesh::dot(v_e, v_p) / v_e.sqrnorm();
}

void MyParaQuadCutting::find_replace(OpenMesh::HalfedgeHandle& h0_p, OpenMesh::VertexHandle v_hf_p, OpenMesh::VertexHandle v_ht_p, OpenMesh::HalfedgeHandle& h0_m, OpenMesh::VertexHandle v_hf_m, OpenMesh::VertexHandle v_ht_m)
{
	if (!new_para.status(h0_p).deleted()) return;

	auto new_h0_p = new_para.find_halfedge(v_hf_p, v_ht_p);
	auto new_h0_m = new_mesh.find_halfedge(v_hf_m, v_ht_m);

	//new_para.property(e_segment, new_para.edge_handle(new_h0_p)) = new_para.property(e_segment, new_para.edge_handle(h0_p));

	new_para.property(h_para2mesh, new_h0_p) = new_para.property(h_para2mesh, h0_p);
	new_mesh.property(h_mesh2para, new_h0_m) = new_mesh.property(h_mesh2para, h0_m);

	h0_p = new_h0_p;
	h0_m = new_h0_m;
}

OpenMesh::VertexHandle MyParaQuadCutting::split_edge(OpenMesh::HalfedgeHandle h_para, double t)
{
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

	for (int i = 0; i < 2; i++)
	{
		if (is_bh_mesh[i]) continue;

		h_split[i] = get_mesh2para(h_split_mesh[i]);

		//int es_seg = new_para.property(e_segment, new_para.edge_handle(h_split[i]));

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

		//new_para.property(e_segment, new_para.edge_handle(new_para.next_halfedge_handle(h_face[0]))) = -1;

		new_para.property(h_para2mesh, new_para.next_halfedge_handle(h_face[0])) = std::make_pair(v_ht_mesh[0].idx(), v_new_mesh.idx());
		new_para.property(h_para2mesh, new_para.prev_halfedge_handle(h_face[1])) = std::make_pair(v_new_mesh.idx(), v_hf_mesh[1].idx());

		new_mesh.property(h_mesh2para, new_mesh.next_halfedge_handle(h_face_mesh[0])) = std::make_pair(v_ht[0].idx(), v_new[i].idx());
		new_mesh.property(h_mesh2para, new_mesh.prev_halfedge_handle(h_face_mesh[1])) = std::make_pair(v_new[i].idx(), v_hf[1].idx());

		//new_para.property(e_segment, new_para.edge_handle(new_para.prev_halfedge_handle(h_face[0]))) = es_seg;
		//new_para.property(e_segment, new_para.edge_handle(new_para.next_halfedge_handle(h_face[1]))) = es_seg;

		new_para.property(h_para2mesh, new_para.prev_halfedge_handle(h_face[0])) = std::make_pair(v_new_mesh.idx(), v_hf_mesh[0].idx());
		new_para.property(h_para2mesh, new_para.next_halfedge_handle(h_face[1])) = std::make_pair(v_ht_mesh[1].idx(), v_new_mesh.idx());

		new_mesh.property(h_mesh2para, new_mesh.prev_halfedge_handle(h_face_mesh[0])) = std::make_pair(v_new[i].idx(), v_hf[0].idx());
		new_mesh.property(h_mesh2para, new_mesh.next_halfedge_handle(h_face_mesh[1])) = std::make_pair(v_ht[1].idx(), v_new[i].idx());
	}

	return v_new[0];
}

OpenMesh::VertexHandle MyParaQuadCutting::split_face(OpenMesh::FaceHandle f_para, double u, double v)
{
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

	/*for (auto ve : new_para.ve_range(v_new))
	{
		new_para.property(e_segment, ve) = -1;
	}*/

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

void MyParaQuadCutting::trans_textured(Mesh& tar)
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

void MyParaQuadCutting::cut_segment()
{
	new_para.add_property(e_cut);
	for (auto e_h : new_mesh.edges())
	{
		new_para.property(e_cut, e_h) = 0;
	}
	for (int i = 0; i < cut_segments.size(); i++)
	{
		OpenMesh::Vec3d segment1 = cut_segments[i].first;
		OpenMesh::Vec3d segment2 = cut_segments[i].second;
		process_point_in_face(segment1);
		process_point_in_face(segment2);
		OpenMesh::EdgeHandle edge_init;
		for (auto e_h : new_para.edges())
		{
			bool edge_init_flag = is_point_on_edge(segment1, e_h, new_para);
			if (edge_init_flag)
			{
				edge_init = e_h;
				break;
			}
		}

		OpenMesh::HalfedgeHandle he_temp = new_para.halfedge_handle(edge_init, 0);
		OpenMesh::Vec3d fpoint = new_para.point(new_para.from_vertex_handle(he_temp));
		OpenMesh::Vec3d tpoint = new_para.point(new_para.to_vertex_handle(he_temp));
		double t_temp = (segment1 - fpoint).norm() / (tpoint - fpoint).norm();
		OpenMesh::VertexHandle new_vertex_pre_temp = split_edge(he_temp, t_temp);
		for (auto vv_h : new_para.vv_range(new_vertex_pre_temp))
		{
			OpenMesh::HalfedgeHandle new_he_temp = new_para.find_halfedge(new_vertex_pre_temp, vv_h);
			OpenMesh::EdgeHandle  new_edge_temp = new_para.edge_handle(new_he_temp.idx() / 2);
			new_para.property(e_cut, new_edge_temp) = 0;
		}
		OpenMesh::VertexHandle new_vertex_temp = new_vertex_pre_temp;
		he_temp = find_next_point(new_vertex_temp, segment1, segment2, t_temp);

		OpenMesh::Vec3d point_temp = segment1;
		bool cut_over = 0;
		while (!cut_over)
		{
			new_vertex_temp = split_edge(he_temp, t_temp);
			mark_e_cut(new_vertex_pre_temp, new_vertex_temp);
			new_vertex_pre_temp = new_vertex_temp;
			he_temp = find_next_point(new_vertex_temp, segment1, segment2, t_temp);
			fpoint = new_para.point(new_para.from_vertex_handle(he_temp));
			tpoint = new_para.point(new_para.to_vertex_handle(he_temp));
			point_temp = (1 - t_temp)*fpoint + t_temp * tpoint;
			OpenMesh::Vec3d new_point_temp = new_para.point(new_vertex_temp);
			if ((segment1 - point_temp).norm() < (segment1 - new_point_temp).norm())
			{
				cut_over = 1;
			}
		}
	}
}

bool MyParaQuadCutting::is_intersect(OpenMesh::Vec3d p1_, OpenMesh::Vec3d p2_, OpenMesh::Vec3d p3_, OpenMesh::Vec3d p4_)
{
	bool flag1, flag2;
	OpenMesh::Vec3d vec12 = p2_ - p1_;
	OpenMesh::Vec3d vec13 = p3_ - p1_;
	OpenMesh::Vec3d vec14 = p4_ - p1_;
	if (((vec13%vec12) | (vec14%vec12)) <= 0)
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
	if (((vec31%vec34) | (vec32%vec34)) <= 0)
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

OpenMesh::Vec3d MyParaQuadCutting::compute_intersect_point(OpenMesh::Vec3d p1_, OpenMesh::Vec3d p2_, OpenMesh::Vec3d p3_, OpenMesh::Vec3d p4_)
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

OpenMesh::HalfedgeHandle MyParaQuadCutting::find_next_point(OpenMesh::VertexHandle vertex_, OpenMesh::Vec3d segment1_, OpenMesh::Vec3d segment2_, double &t_)
{
	std::vector<OpenMesh::HalfedgeHandle> he_list;
	std::vector<double> t_list;
	std::vector<double> length_list;
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
		t_list.push_back(t_temp);
		double length_temp = (segment1_ - new_point_temp).norm();
		length_list.push_back(length_temp);
	}
	std::vector<double>::iterator biggest = std::max_element(length_list.begin(), length_list.end());
	int biggest_position = std::distance(length_list.begin(), biggest);
	t_ = t_list[biggest_position];
	return he_list[biggest_position];
}

void MyParaQuadCutting::process_point_in_face(OpenMesh::Vec3d p_)
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
			OpenMesh::Vec3d new_point_temp = (1 - t_t)*new_para.point(new_para.from_vertex_handle(face_he_temp)) + t_t * new_para.point(new_para.to_vertex_handle(face_he_temp));
			update_cut_segments(p_, new_point_temp);
		}
	}
}

void MyParaQuadCutting::update_cut_segments(OpenMesh::Vec3d p_old_, OpenMesh::Vec3d p_new_)
{
	for (int i = 0; i < cut_segments.size(); i++)
	{
		OpenMesh::Vec3d segment1 = cut_segments[i].first;
		OpenMesh::Vec3d segment2 = cut_segments[i].second;
		if (segment1 == p_old_)
		{
			cut_segments[i].first = p_new_;
		}
		if (segment2 == p_old_)
		{
			cut_segments[i].second = p_new_;
		}
	}
}

void MyParaQuadCutting::mark_e_cut(OpenMesh::VertexHandle v_pre_, OpenMesh::VertexHandle v_now_)
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
			new_para.property(e_cut, new_edge_temp) = 0;
		}
	}
}

void MyParaQuadCutting::para_decomposition()
{
	for (auto e_h : new_para.edges())
	{
		if (new_para.is_boundary(e_h))
		{
			new_para.property(e_cut, e_h) = 1;
		}
	}

	std::vector<int> is_check_face;
	is_check_face.resize(new_para.n_faces());
	while (std::accumulate(is_check_face.begin(), is_check_face.end(), 0) != new_para.n_faces())
	{
		Mesh atlas_temp;
		std::vector<OpenMesh::FaceHandle> atlas_face;
		std::vector<int> f_atlas2mesh_temp;
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
		while (!myqueue.empty())
		{
			OpenMesh::FaceHandle face_temp = myqueue.front();
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
			myqueue.pop();
		}
		std::vector<int> is_check_vertex;
		is_check_vertex.resize(new_para.n_vertices());
		for (int i = 0; i < atlas_face.size(); i++)
		{
			for (auto fv_h : new_para.fv_range(atlas_face[i]))
			{
				is_check_vertex[fv_h.idx()] = 1;
			}
		}
		std::vector<OpenMesh::VertexHandle> v2newv;
		v2newv.resize(new_para.n_vertices());
		for (int i = 0; i < is_check_vertex.size(); i++)
		{
			if (is_check_vertex[i] == 1)
			{
				OpenMesh::VertexHandle newv = atlas_temp.add_vertex(new_para.point(new_para.vertex_handle(i)));
				v2newv[i] = newv;
			}
		}
		f_atlas2mesh_temp.resize(atlas_face.size());
		for (int i = 0; i < atlas_face.size(); i++)
		{
			std::vector<int> face_vertex;
			for (auto cfv = new_para.cfv_begin(atlas_face[i]); cfv != new_para.cfv_end(atlas_face[i]); cfv++)
			{
				face_vertex.push_back(cfv.handle().idx());
			}
			OpenMesh::FaceHandle newf = atlas_temp.add_face(v2newv[face_vertex[0]], v2newv[face_vertex[1]], v2newv[face_vertex[2]]);
			f_atlas2mesh_temp[newf.idx()] = atlas_face[i].idx();
		}
		my_para_atlas.push_back(atlas_temp);
		f_atlas2mesh.push_back(f_atlas2mesh_temp);
	}
}