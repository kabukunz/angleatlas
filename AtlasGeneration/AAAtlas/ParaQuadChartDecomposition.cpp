#include "ParaQuadChartDecomposition.h"

#include <queue>
#include <forward_list>
#include <iostream>
#include <ctime>

#include "finders_interface.h"

#include "MeshViewer/Mesh_doubleIO.h"
#include "Common/CommonFunctions.h"

#ifdef MY_DEBUG
#define MY_DOUBT(cond, msg) if (cond) std::cout << msg << std::endl;
#else
#define MY_DOUBT(cond, msg)
#endif

OpenMesh::Vec3d tag_direction[4] =
{
	{ 1.0, 0.0, 0.0 },
	{ 0.0, 1.0, 0.0 },
	{-1.0, 0.0, 0.0 },
	{ 0.0,-1.0, 0.0 }
};

ParaQuadChartDecomposition::ParaQuadChartDecomposition(const Mesh& IGPP_mesh)
	: mesh(IGPP_mesh)
{
	viewer = std::make_unique<Viewer2D_Segments>();
	viewer->hide();
}

bool ParaQuadChartDecomposition::range_test(int seg_id, int val)
{
	auto seg_range = get_seg_range(seg_id);

	return (val > std::min(seg_range.first, seg_range.second)) && (val < std::max(seg_range.first, seg_range.second));
}

std::pair<int, int> ParaQuadChartDecomposition::get_seg_range(int seg_id)
{
	int t = segments[seg_id].tag & 1;
	int v0 = mesh.point(mesh.vertex_handle(segments[seg_id].vert0))[t];
	int v1 = (segments[seg_id].tag & 4) ? segments[segments[seg_id].vert1].coord : mesh.point(mesh.vertex_handle(segments[seg_id].vert1))[t];

	return std::make_pair(v0, v1);
}

void ParaQuadChartDecomposition::decomposition()
{
	std::srand(std::time(0));

	segment_detective();
	quad_composition();

	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "Inner intersections : " << seg_intersections.size() << ", Quad Vertices : " << quad_vertices.size() << std::endl;
	std::cout << "Boundary Segments : " << n_bsegs << ", Inner Vertices : " << segments.size() - n_bsegs << std::endl;

	show_quad_decomposition();

	motorcycle();

	chart_decomposition();

	cut_preprocess();
	extract_cut();
}

inline int ParaQuadChartDecomposition::chart_boundary_cycle(int cur, int delta)
{
	int cur_chart = seg_chart[bh_segment[cur]];
	int n_bh_chart = accu_bh_chart[cur_chart + 1] - accu_bh_chart[cur_chart];
	int bh0_chart = accu_bh_chart[cur_chart];
	return CommonFunctions::period_id(cur - bh0_chart + delta, n_bh_chart) + bh0_chart;
}

inline int ParaQuadChartDecomposition::chart_segment_cycle(int cur, int delta)
{
	int cur_chart = seg_chart[cur];
	int n_bseg_chart = accu_bseg_chart[cur_chart + 1] - accu_bseg_chart[cur_chart];
	int bseg0_chart = accu_bseg_chart[cur_chart];
	return CommonFunctions::period_id(cur - bseg0_chart + delta, n_bseg_chart) + bseg0_chart;
}

void ParaQuadChartDecomposition::show_quad_decomposition()
{
	for (int i = 0; i < segments.size(); i++)
	{
		auto p0 = quad_vertices[seg_nodes[i][0]];
		auto p1 = quad_vertices[seg_nodes[i].back()];

		viewer->add_segment({ p0[0], p0[1] }, { p1[0], p1[1] });
		
		for (int j = 0; j < seg_nodes[i].size(); j++)
		{
			auto pp = quad_vertices[seg_nodes[i][j]];
			viewer->add_point(pp[0], pp[1]);
		}
	}

	viewer->add_mesh(mesh);

	viewer->show();
	viewer->adjust_range();
}

void ParaQuadChartDecomposition::segment_detective()
{
	std::set<int> add_c;
//	std::set<int> add_c = sample_boundary_vertices();

	int n_inner_segs = 0;
	vert_corner.clear();
	quad_vertices.clear();
	accu_bh_chart.assign(1, 0);
	accu_bseg_chart.assign(1, 0);
	std::vector<bool> h_visited(mesh.n_halfedges(), false);
	for (auto h : mesh.halfedges())
	{
		if (!mesh.is_boundary(h) || h_visited[h.idx()]) continue;
		
		int euler = 0;
		double z_normal = 0.0;
		auto bh_iter = h;
		do
		{
			h_visited[bh_iter.idx()] = true;

			double angle_v = 0;
			int to_v = mesh.to_vertex_handle(bh_iter).idx();
			bh_iter = mesh.opposite_halfedge_handle(bh_iter);

			while (!mesh.is_boundary(bh_iter))
			{
				auto vec0 = mesh.calc_edge_vector(bh_iter);
				bh_iter = mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(bh_iter));
				auto vec1 = mesh.calc_edge_vector(bh_iter);

				z_normal += vec0[0] * vec1[1] - vec0[1] * vec1[0];
				angle_v += CommonFunctions::vec_angle_acos(vec0, vec1);              //ÇóÄÚ½Ç
			}

			int vk = std::lround(angle_v / M_PI_2);
			if (vk != 2 || add_c.count(to_v) == 1)
			{
				vert_corner.push_back(to_v);
				n_inner_segs += vk - 1;

				v_in_quad[to_v] = quad_vertices.size();
				const auto& uv_bv = mesh.point(mesh.vertex_handle(to_v));
				quad_vertices.emplace_back(std::lround(uv_bv[0]), std::lround(uv_bv[1]));

//				std::cout << to_v << " " << angle_v / M_PI_2 << " " << vk << std::endl;
			}
			euler += 2 - vk;
			vert_k[to_v] = vk;
		} while (bh_iter != h);

		accu_bh_chart.push_back(vert_k.size());
		accu_bseg_chart.push_back(vert_corner.size());
		MY_DOUBT(euler != 4, "Euler Error");
		MY_DOUBT(z_normal <= 0, "Normal Error");
	}

	int n_charts = accu_bseg_chart.size() - 1;
	seg_chart.resize(vert_corner.size());
	boundary_v.reserve(vert_k.size());
	boundary_h.reserve(vert_k.size());
	for (int i = 0; i < n_charts; i++)
	{
		auto bv0 = mesh.vertex_handle(vert_corner[accu_bseg_chart[i]]);
		OpenMesh::HalfedgeHandle bh_iter;
		for (auto voh : mesh.voh_range(bv0))
		{
			bh_iter = voh;
			if (mesh.is_boundary(voh)) break;
		}

		for (int j = accu_bh_chart[i]; j < accu_bh_chart[i + 1]; j++)
		{
			auto bv_iter = mesh.from_vertex_handle(bh_iter);
			boundary_v_index[bv_iter.idx()] = boundary_v.size();
			boundary_h_index[bh_iter.idx()] = boundary_h.size();

			boundary_v.push_back(bv_iter.idx());
			boundary_h.push_back(bh_iter.idx());

			bh_iter = mesh.next_halfedge_handle(bh_iter);
		}

		for (int j = accu_bseg_chart[i]; j < accu_bseg_chart[i + 1]; j++)
		{
			seg_chart[j] = i;
		}
	}
	n_bsegs = vert_corner.size();

	MY_DOUBT(n_bsegs - n_inner_segs != 4 * n_charts, "Segs Num Error");

	segments.clear();
	segments.reserve(n_bsegs + n_inner_segs);
	segments.resize(n_bsegs);

	seg_nodes.clear();
	seg_nodes.reserve(n_bsegs + n_inner_segs);
	seg_nodes.resize(n_bsegs);

	triangle_segs.clear();
	triangle_segs.resize(mesh.n_faces());

	bh_segment.resize(boundary_h.size());

	for (int i = 0; i < vert_corner.size(); i++)
	{
		segments[i].vert0 = vert_corner[i];
		segments[i].vert1 = vert_corner[chart_segment_cycle(i, 1)];

		auto vec_seg = mesh.point(mesh.vertex_handle(segments[i].vert1)) - mesh.point(mesh.vertex_handle(segments[i].vert0));
		segments[i].tag = CommonFunctions::get_tag(vec_seg);
		segments[i].coord = std::lround(mesh.point(mesh.vertex_handle(segments[i].vert0))[(segments[i].tag & 1) ^ 1]);

		seg_nodes[i] = { v_in_quad[segments[i].vert0], v_in_quad[segments[i].vert1] };

		for (auto vf_h : mesh.vf_range(mesh.vertex_handle(vert_corner[i])))
		{
			triangle_corner[vf_h.idx()].push_back(v_in_quad[vert_corner[i]]);
		}
	}

	int cur_seg = 0;
	for (int i = 0; i < boundary_h.size(); i++)
	{
		if (boundary_v[i] == segments[cur_seg].vert1 || i == accu_bh_chart[seg_chart[cur_seg] + 1]) cur_seg++;
	
		triangle_segs[mesh.opposite_face_handle(mesh.halfedge_handle(boundary_h[i])).idx()].insert(cur_seg);
		bh_segment[i] = cur_seg;
	}
}

void ParaQuadChartDecomposition::quad_composition()
{
	cornor_k_seg.resize(vert_corner.size());
	for (int i = 0; i < vert_corner.size(); i++)
	{
		//if (vert_k[vert_corner[i]] <= 2) continue;
		
		cornor_k_seg[i].assign(vert_k[vert_corner[i]] + 1, -1);
		cornor_k_seg[i][0] = bh_segment[chart_boundary_cycle(boundary_v_index[vert_corner[i]], -1)];
		cornor_k_seg[i][vert_k[vert_corner[i]]] = bh_segment[boundary_v_index[vert_corner[i]]];
	}
	
	for (int vert_id : vert_corner)
	{
		//if (vert_k[vert_id] <= 2) continue;
		
		tracing_from_vert(vert_id);
	}

	for (int i = n_bsegs; i < segments.size(); i++)
	{
		if (segments[i].tag & 4)
		{
			auto uv1 = mesh.point(mesh.vertex_handle(segments[i].vert0));
			uv1[segments[i].tag & 1] = segments[segments[i].vert1].coord;
			
			seg_nodes[i].push_back(quad_vertices.size());
			seg_nodes[segments[i].vert1].push_back(quad_vertices.size());
			quad_vertices.emplace_back(std::lround(uv1[0]), std::lround(uv1[1]));
		}
	}
	
	for (const auto& inter : seg_intersections)
	{
		int new_quad_v = quad_vertices.size();
		quad_vertices.emplace_back(segments[inter.first].coord, segments[inter.second].coord);
		
		seg_nodes[inter.first].push_back(new_quad_v);
		seg_nodes[inter.second].push_back(new_quad_v);
	}

	for (int i = 0; i < segments.size(); i++)
	{
		switch (segments[i].tag & 3)
		{
		case 0:
			std::sort(seg_nodes[i].begin(), seg_nodes[i].end(), [&](int a, int b){return quad_vertices[a][0] < quad_vertices[b][0]; });
			break;
		case 1:
			std::sort(seg_nodes[i].begin(), seg_nodes[i].end(), [&](int a, int b){return quad_vertices[a][1] < quad_vertices[b][1]; });
			break;
		case 2:
			std::sort(seg_nodes[i].begin(), seg_nodes[i].end(), [&](int a, int b){return quad_vertices[a][0] > quad_vertices[b][0]; });
			break;
		case 3:
			std::sort(seg_nodes[i].begin(), seg_nodes[i].end(), [&](int a, int b){return quad_vertices[a][1] > quad_vertices[b][1]; });
			break;
		default:
			MY_DOUBT(true, "Tag Error");
		}
	}
	
	std::vector<bool> multi_next(quad_vertices.size(), false);
	quad_v_neighbor.resize(quad_vertices.size());
	for (int vert_id : vert_corner)
	{
		if (vert_k[vert_id] >= 4)
		{
			int cornor_id = v_in_quad[vert_id];
			int new_size = (vert_k[vert_id] / 4 + 2) * 4;
			quad_v_neighbor[cornor_id].next_v.assign(new_size, -1);
			quad_v_neighbor[cornor_id].next_seg.assign(new_size, -1);
			multi_next[cornor_id] = true;

			auto k_segs = cornor_k_seg[cornor_id];

			int k0_h = chart_boundary_cycle(boundary_v_index[vert_id], -1);
			int tag0 = CommonFunctions::get_tag(-mesh.calc_edge_vector(mesh.halfedge_handle(boundary_h[k0_h])));

			for (int i = 0; i < k_segs.size(); i++)
			{
				quad_v_neighbor[cornor_id].next_seg[tag0 + i] = k_segs[i];
				quad_v_neighbor[cornor_id].next_v[tag0 + i] = (cornor_id == seg_nodes[k_segs[i]][0]) ? seg_nodes[k_segs[i]][1] : *(seg_nodes[k_segs[i]].end() - 2);
			} 
		}
	}

	for (int i = 0; i < segments.size(); i++)
	{
		uint t0 = segments[i].tag & 3;
		uint t1 = t0 ^ 2;

		for (int j = 0; j < seg_nodes[i].size() - 1; j++)
		{
			if (!multi_next[seg_nodes[i][j]])
			{
				quad_v_neighbor[seg_nodes[i][j]].next_v[t0] = seg_nodes[i][j + 1];
				quad_v_neighbor[seg_nodes[i][j]].next_seg[t0] = i;
			}

			if (!multi_next[seg_nodes[i][j + 1]])
			{
				quad_v_neighbor[seg_nodes[i][j + 1]].next_v[t1] = seg_nodes[i][j];
				quad_v_neighbor[seg_nodes[i][j + 1]].next_seg[t1] = i;
			}
		}
	}

	for (int i = 0; i < quad_vertices.size(); i++)
	{
		quad_mesh.add_vertex(Mesh::Point(quad_vertices[i][0], quad_vertices[i][1], 0.0));
	}

	for (int i = 0; i < segments.size(); i++)
	{
		const int tag_step = 1;
		if ((segments[i].tag & 1) || (i < n_bsegs && ((segments[i].tag + tag_step) & 3) == 3)) continue;

		std::vector<int> quad_tl = seg_nodes[i];
		if ((segments[i].tag & 3) == 2) std::reverse(quad_tl.begin(), quad_tl.end());
		
		for (int j = 0; j < quad_tl.size() - 1; j++)
		{
			int cur_seg = i;
			OpenMesh::VertexHandle quad_v[4];

			quad_v[0] = quad_mesh.vertex_handle(quad_tl[j]);
 			int cur_tag = 0;

			for (int k = 1; k < 4; k++)
			{
				if (multi_next[quad_v[k - 1].idx()])
				{
					for (int l = 0; l < quad_v_neighbor[quad_v[k - 1].idx()].next_v.size() / 4; l++)
					{
						if (quad_v_neighbor[quad_v[k - 1].idx()].next_seg[l * 4 + cur_tag] == cur_seg)
						{
							quad_v[k] = quad_mesh.vertex_handle(quad_v_neighbor[quad_v[k - 1].idx()].next_v[l * 4 + cur_tag]);
							break;
						}
					}
				}
				else
				{
					quad_v[k] = quad_mesh.vertex_handle(quad_v_neighbor[quad_v[k - 1].idx()].next_v[cur_tag]);
				}

				if (multi_next[quad_v[k].idx()])
				{
					for (int l = 0; l < quad_v_neighbor[quad_v[k].idx()].next_v.size() / 4; l++)
					{
						if (quad_v_neighbor[quad_v[k].idx()].next_seg[l * 4 + (cur_tag ^ 2)] == cur_seg)
						{
							cur_seg = quad_v_neighbor[quad_v[k].idx()].next_seg[l * 4 + (cur_tag ^ 2) + 1];
							break;
						}
					}
					cur_tag = (cur_tag - 1) & 3;
				}
				else
				{
					cur_tag = (cur_tag - 1) & 3;
					cur_seg = quad_v_neighbor[quad_v[k].idx()].next_seg[cur_tag];
				}
			}

			quad_mesh.add_face(quad_v[0], quad_v[3], quad_v[2], quad_v[1]);
			quads.emplace_back(quad_vertices[quad_v[3].idx()], quad_vertices[quad_v[1].idx()]);
		}
	}
}

void ParaQuadChartDecomposition::tracing_from_vert(int vert_id)
{
	int bv_id = boundary_v_index[vert_id]; 
	int chart_id = seg_chart[bh_segment[bv_id]];
	auto h_iter = mesh.opposite_halfedge_handle(mesh.halfedge_handle(boundary_h[chart_boundary_cycle(bv_id, -1)]));
	auto uv0 = mesh.point(mesh.vertex_handle(vert_id));

	int k_cur = 1;
	const int tag_step = 1;
	int tag_iter = (CommonFunctions::get_tag(mesh.calc_edge_vector(h_iter)) + tag_step) & 3;
	
	while (h_iter.idx() != boundary_h[bv_id])
	{
		int gtag = (tag_iter + tag_step) & 3;

		auto v1 = mesh.point(mesh.to_vertex_handle(h_iter)) - uv0;
		auto v2 = mesh.point(mesh.from_vertex_handle(mesh.prev_halfedge_handle(h_iter))) - uv0;
		
		if ((!(gtag & 2) && (v1[gtag & 1] <= 0) && (v2[gtag & 1] >= 0)) || ((gtag & 2) && (v1[gtag & 1] >= 0) && (v2[gtag & 1] <= 0)))
		{
			if (cornor_k_seg[v_in_quad[vert_id]][k_cur] < 0)
			{
				//tracing
				int tracing_seg_id = segments.size();
				seg_chart.push_back(chart_id);
				seg_nodes.emplace_back();
				cornor_k_seg[v_in_quad[vert_id]][k_cur] = tracing_seg_id;
				std::set<int> intersections;

				double tracing_val = uv0[gtag & 1];
				auto tracing_h = mesh.next_halfedge_handle(h_iter);

				double p1_t = mesh.point(mesh.from_vertex_handle(tracing_h))[gtag & 1];
				double p2_t = mesh.point(mesh.to_vertex_handle(tracing_h))[gtag & 1];

				MY_DOUBT((p1_t - tracing_val) * (p2_t - tracing_val) > 0, "Tracing Error 1, v " << vert_id << ", k " << k_cur);

				tracing_h = mesh.opposite_halfedge_handle(tracing_h);
				while (true)
				{
					tracing_h = mesh.opposite_halfedge_handle(tracing_h);
					auto tracing_f = mesh.face_handle(tracing_h);

					int intersect_bseg = -1;
					for (int seg_id : triangle_segs[tracing_f.idx()])
					{
						if ((segments[seg_id].tag & 1) != (gtag & 1)) continue;
						if (!range_test(seg_id, tracing_val)) continue;
						if ((!(tag_iter & 2) && (uv0[tag_iter & 1] >= segments[seg_id].coord)) || ((tag_iter & 2) && (uv0[tag_iter & 1] <= segments[seg_id].coord))) continue;
						
						if (seg_id < n_bsegs && (segments[seg_id].vert0 != vert_id) && (segments[seg_id].vert1 != vert_id))
						{
							intersect_bseg = seg_id;
						}
						else if (seg_id >= n_bsegs)
						{
							intersections.insert(seg_id);
						}
					}

					bool arrival = false;
					int arrival_quad_v = -1;
					auto tc = triangle_corner.find(tracing_f.idx());
					if (tc != triangle_corner.end())
					{
						for (int qv : tc->second)
						{
							if (qv != v_in_quad[vert_id] && std::lround(tracing_val) == quad_vertices[qv][gtag & 1])
							{
								arrival_quad_v = qv;
								break;
							}
						}

						if (arrival_quad_v != -1)
						{
							//std::cout << uv0 << " -> " << quad_vertices[arrival_quad_v] << std::endl;
							arrival = true;
							seg_nodes[tracing_seg_id].push_back(arrival_quad_v);
							segments.emplace_back(segment_type{ vert_id, vert_corner[arrival_quad_v], std::lround(tracing_val), tag_iter });

							int end_bv = boundary_v_index[vert_corner[arrival_quad_v]]; 
							auto end_h_iter = mesh.opposite_halfedge_handle(mesh.halfedge_handle(boundary_h[chart_boundary_cycle(end_bv, -1)]));

							double end_angle = 0.0;
							while (mesh.face_handle(end_h_iter) != tracing_f)
							{
								auto vec0 = mesh.calc_edge_vector(end_h_iter);
								end_h_iter = mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(end_h_iter));
								auto vec1 = mesh.calc_edge_vector(end_h_iter);

								end_angle += CommonFunctions::vec_angle_atan2(vec0, vec1);
							}

							end_angle += CommonFunctions::vec_angle_atan2(mesh.calc_edge_vector(end_h_iter), tag_direction[tag_iter ^ 2]);
							end_angle = (tag_step == 1) ? end_angle : -end_angle;

							cornor_k_seg[arrival_quad_v][std::lround(end_angle / M_PI_2)] = tracing_seg_id;
						}
					}
					if (arrival_quad_v == -1 && intersect_bseg >= 0)
					{
						arrival = true;
						segments.emplace_back(segment_type{ vert_id, intersect_bseg, std::lround(tracing_val), tag_iter + 4 });
					}
					triangle_segs[tracing_f.idx()].insert(tracing_seg_id);
					if (arrival) break;

					if (tracing_h == mesh.next_halfedge_handle(h_iter)) continue;

					double c0 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(tracing_h)))[gtag & 1];
					double c1 = mesh.point(mesh.from_vertex_handle(tracing_h))[gtag & 1];
					double c2 = mesh.point(mesh.to_vertex_handle(tracing_h))[gtag & 1];

					MY_DOUBT((c1 - tracing_val) * (c2 - tracing_val) > 0, "Tracing Error 2, v " << vert_id << ", k " << k_cur);

					auto t1t0 = mesh.prev_halfedge_handle(tracing_h);
					auto t2t0 = mesh.next_halfedge_handle(tracing_h);
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
						double d_t1 = std::abs(mesh.point(mesh.from_vertex_handle(tracing_h))[(gtag & 1) ^ 1] - uv0[(gtag & 1) ^ 1]);
						double d_t2 = std::abs(mesh.point(mesh.to_vertex_handle(tracing_h))[(gtag & 1) ^ 1] - uv0[(gtag & 1) ^ 1]);

						tracing_h = (d_t1 <= d_t2) ? t2t0 : t1t0;
					}
				}

				auto range_pair = get_seg_range(tracing_seg_id);
				auto range_min = std::min(range_pair.first, range_pair.second);
				auto range_max = std::max(range_pair.first, range_pair.second);

				for (int seg_id : intersections)
				{
					if (segments[seg_id].coord >= range_min || segments[seg_id].coord <= range_max)
					{
						seg_intersections.emplace((gtag & 1) ? seg_id : tracing_seg_id, (gtag & 1) ? tracing_seg_id : seg_id);
					}
				}

				seg_nodes[tracing_seg_id].push_back(v_in_quad[vert_id]);
			}

			k_cur++;
			tag_iter = (tag_iter + tag_step) & 3;

			if (k_cur == vert_k[vert_id]) break;
		}
		else
		{
			h_iter = mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(h_iter));
		}
	}
}

void ParaQuadChartDecomposition::motorcycle()
{
	double weight_len = 0.1;
	quad_mesh_perimeter = 0.0;
	for (auto e_h : quad_mesh.edges())
	{
		if (quad_mesh.is_boundary(e_h)) quad_mesh_perimeter += quad_mesh.calc_edge_length(e_h);
	}
	
	double cur_eff = 0.0, cur_eff_bad = 0.0;
	double max_scr = -std::numeric_limits<double>::max();
	double max_scr_bad = -std::numeric_limits<double>::max();
	std::vector<bool> quad_edge_inner_label_bad;
	std::vector<std::vector<int>> chart_faces_bad;

	std::vector<bool> chart_flipped_bad;
	std::vector<OpenMesh::Vec2d> chart_translation_bad;

	std::vector<double> effs;
	for (int i = 0; i < 500; i++)
	{
		double new_cut_len = 0.0;
		motorcycle_canididate();
		bool isgood = retangle_decomposition(new_cut_len);
		double eff = retangle_packing(isgood);

		effs.push_back(eff);
		double score = eff - weight_len * new_cut_len / quad_mesh_perimeter;

		if (isgood)
		{
			if (score <= max_scr) continue;

			max_scr = score;
			cur_eff = eff;
			rects_candidate.clear();
			quad_edge_inner_label = std::move(qeil_candidate);
			chart_faces = std::move(chart_faces_candidate);
			chart_flipped = std::move(cf_candidate);
			chart_translation = std::move(ct_candidate);
		}
		else
		{
			if (score <= max_scr_bad) continue;

			max_scr_bad = score;
			cur_eff_bad = eff;
			rects_candidate.clear();
			quad_edge_inner_label_bad = std::move(qeil_candidate);
			chart_faces_bad = std::move(chart_faces_candidate);
			chart_flipped_bad = std::move(cf_candidate);
			chart_translation_bad = std::move(ct_candidate);
		}
	}

	std::cout << "PE Good " << cur_eff << ", Bad " << cur_eff_bad << std::endl;
	if (cur_eff == 0.0)
	{
		max_scr = max_scr_bad;
		cur_eff = cur_eff_bad;
		rects_candidate.clear();
		quad_edge_inner_label = std::move(quad_edge_inner_label_bad);
		chart_faces = std::move(chart_faces_bad);
		chart_flipped = std::move(chart_flipped_bad);
		chart_translation = std::move(chart_translation_bad);
	}

	std::cout << "Packing Efficiency" << std::endl;
	std::cout << "Max " << *std::max_element(effs.begin(), effs.end()) << std::endl;
	std::cout << "Min " << *std::min_element(effs.begin(), effs.end()) << std::endl;
	double eff_avg = std::accumulate(effs.begin(), effs.end(), 0.0) / effs.size();
	double eff_dev = 0.0;
	std::for_each(effs.begin(), effs.end(), [&](double x) {eff_dev += (x - eff_avg) * (x - eff_avg); });
	eff_dev = std::sqrt(eff_dev / effs.size());
	std::cout << "Avg " << eff_avg << std::endl;
	std::cout << "Dev " << eff_dev << std::endl;

	std::cout << "Res " << cur_eff << ", " << (cur_eff - max_scr) / weight_len << std::endl;
}

void ParaQuadChartDecomposition::motorcycle_canididate()
{
	qeil_candidate.resize(quad_mesh.n_edges());
	for (auto e_h : quad_mesh.edges())
	{
		qeil_candidate[e_h.idx()] = !quad_mesh.is_boundary(e_h);
	}
	
	struct motor_status
	{
		int seg_id;
		int tag;
		int cur_vert;
		int next_vert;
	};

	std::vector<motor_status> motor_list;
	std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> motorcycle_queue;

	auto motor_next_step = [&](int motor_id, int cur_dist, int next_v = -1)
	{
		auto& motor_path = motor_list[motor_id];
		int next_vert = (next_v == -1) ? quad_v_neighbor[motor_path.cur_vert].next_v[motor_path.tag] : next_v;
		int dist_coord = motor_path.tag & 1;
		int step_dist = std::abs(quad_vertices[next_vert][dist_coord] - quad_vertices[motor_path.cur_vert][dist_coord]);
		motor_path.next_vert = next_vert;
		motorcycle_queue.emplace(step_dist + cur_dist, motor_id);
	};

	for (int vert_id : vert_corner)
	{
		int quad_v = v_in_quad[vert_id];
		int k0 = (cornor_k_seg[quad_v].size() & 1) ? 2 : ((std::rand() & 1) + 1);

		for (int i = k0; i < cornor_k_seg[quad_v].size(); i += 2)
		{
			int seg_id = cornor_k_seg[quad_v][i];
			int tag = segments[seg_id].tag & 3;
			tag = (vert_id == segments[seg_id].vert0) ? tag : (tag ^ 2);
			MY_DOUBT(vert_id != segments[seg_id].vert0 && vert_id != segments[seg_id].vert1, "Seg End Point Error " << vert_id << " " << segments[seg_id].vert0 << " " << segments[seg_id].vert1);

			int next_v = -1;
			for (int l = 0; l < quad_v_neighbor[quad_v].next_v.size() / 4; l++)
			{
				if (quad_v_neighbor[quad_v].next_seg[l * 4 + tag] == seg_id)
				{
					next_v = quad_v_neighbor[quad_v].next_v[l * 4 + tag];
					break;
				}
			}
			
			MY_DOUBT(next_v == -1, "CANNOT FIND NEXT V");
			motor_list.emplace_back(motor_status{ seg_id, tag, quad_v, -1 });
			motor_next_step(motor_list.size() - 1, 0, next_v);
		}
	}

	std::vector<bool> motor_continue(motor_list.size(), true);
	std::vector<std::vector<int>> quad_v_adj(quad_vertices.size(), std::vector<int>(4, -1));

	while (!motorcycle_queue.empty())
	{
		auto next_step = motorcycle_queue.top();
		motorcycle_queue.pop();

		int motor_id = next_step.second;
		if (!motor_continue[motor_id]) continue;
		
		auto& motor_path = motor_list[motor_id];

		qeil_candidate[quad_mesh.find_halfedge(quad_mesh.vertex_handle(motor_path.cur_vert), quad_mesh.vertex_handle(motor_path.next_vert)).idx() / 2] = false;

		motor_path.cur_vert = motor_path.next_vert;
		int cur_dist = next_step.first;
		quad_v_adj[motor_path.next_vert][motor_path.tag] = motor_id;

		if (quad_mesh.is_boundary(quad_mesh.vertex_handle(motor_path.next_vert)))
		{
			motor_continue[motor_id] = false;
		}
		else if (quad_v_adj[motor_path.next_vert][(motor_path.tag + 2) & 3] >= 0)
		{
			motor_continue[motor_id] = false;
			motor_continue[quad_v_adj[motor_path.next_vert][(motor_path.tag + 2) & 3]] = false;
		}
		else if (quad_v_adj[motor_path.next_vert][(motor_path.tag + 1) & 3] >= 0 || quad_v_adj[motor_path.next_vert][(motor_path.tag - 1) & 3] >= 0)
		{
			motor_continue[motor_id] = false;
		}

		if (!motor_continue[motor_id]) continue;
		motor_next_step(motor_id, cur_dist);
	}
}

bool ParaQuadChartDecomposition::retangle_decomposition(double& cut_len)
{
	int n_charts = 0;
	std::vector<int> face_chart(quad_mesh.n_faces(), -1);

	for (auto f_h : quad_mesh.faces())
	{
		if (face_chart[f_h.idx()] >= 0) continue;

		int cur_chart = (n_charts++);
		std::queue<int> face_queue;
		face_queue.push(f_h.idx());

		while (!face_queue.empty())
		{
			int fid_cur = face_queue.front();
			face_queue.pop();

			if (face_chart[fid_cur] >= 0) continue;

			face_chart[fid_cur] = cur_chart;

			for (auto fh : quad_mesh.fh_range(quad_mesh.face_handle(fid_cur)))
			{
				if (qeil_candidate[quad_mesh.edge_handle(fh).idx()] && face_chart[quad_mesh.opposite_face_handle(fh).idx()] < 0)
				{
					face_queue.push(quad_mesh.opposite_face_handle(fh).idx());
				}
			}
		}
	}

	chart_faces_candidate.clear();
	chart_faces_candidate.resize(n_charts);
	rects_candidate.clear();
	rects_candidate.resize(n_charts);
	
	for (auto f_h : quad_mesh.faces())
	{
		int chart_id = face_chart[f_h.idx()];
		chart_faces_candidate[chart_id].push_back(f_h.idx());

		for (auto fv_h : quad_mesh.fv_range(f_h))
		{
			const auto& vec3d_p = quad_mesh.point(fv_h);
			OpenMesh::Vec2i vec2i_p(std::lround(vec3d_p[0]), std::lround(vec3d_p[1]));

			rects_candidate[chart_id].min = OpenMesh::minimize(rects_candidate[chart_id].min, vec2i_p);
			rects_candidate[chart_id].max = OpenMesh::maximize(rects_candidate[chart_id].max, vec2i_p);
		}
	}

	bool isgood = true;
	int new_cut_length = 0;
	for (int i = 0; i < n_charts; i++)
	{
		OpenMesh::Vec2i rect_diff = rects_candidate[i].max - rects_candidate[i].min;
		int rect_min = rect_diff.min();
		int rect_max = rect_diff.max();

		if (rect_min < 50 && (double)rect_max / (double)rect_min >= 15.0) isgood = false;

		new_cut_length += 2 * (rect_min + rect_max);
	}
	cut_len = (double)new_cut_length;
	return isgood;
}

double ParaQuadChartDecomposition::retangle_packing(bool& isgood)
{
	int n_charts = chart_faces_candidate.size();
	cf_candidate.assign(n_charts, false);
	ct_candidate.assign(n_charts, OpenMesh::Vec2d(0.0, 0.0));

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

	int area_charts = 0;
	for (int i = 0; i < n_charts; i++)
	{
		area_charts += rects_candidate[i].area();
		auto p_min = rects_candidate[i].min;
		auto range = rects_candidate[i].max - rects_candidate[i].min;

		OpenMesh::Vec2i offset(8, 8);
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
		cf_candidate[i] = r.flipped;

		ct_candidate[i][0] = double(trans_i[0]);
		ct_candidate[i][1] = double(trans_i[1]);

		bb.min = OpenMesh::minimize(bb.min, packing_p_min);

		bb.max = OpenMesh::maximize(bb.max, packing_p_max);
	}

//	if (bb.width() > 3 * bb.height() || bb.height() > 3 * bb.width()) isgood = false;

	return double(area_charts) / double(bb.area());
}

void ParaQuadChartDecomposition::cut_preprocess()
{
	cutting_vertex_valence.assign(quad_mesh.n_vertices(), 0);
	for (auto e_h : quad_mesh.edges())
	{
		if (!quad_edge_inner_label[e_h.idx()])
		{
			auto h_h = quad_mesh.halfedge_handle(e_h, 0);
			int tag = CommonFunctions::get_tag(quad_mesh.calc_edge_vector(h_h));

			int v0 = quad_mesh.from_vertex_handle(h_h).idx();
			int v1 = quad_mesh.to_vertex_handle(h_h).idx();

			for (int l = 0; l < quad_v_neighbor[v0].next_v.size() / 4; l++)
			{
				if (quad_v_neighbor[v0].next_v[l * 4 + tag] == v1)
				{
					cutting_edge_seg[e_h.idx()] = quad_v_neighbor[v0].next_seg[l * 4 + tag];
					break;
				}
			}

			cutting_vertex_valence[v0]++;
			cutting_vertex_valence[v1]++;
		}
	}

	for (int i = 0; i < cutting_vertex_valence.size(); i++)
	{
		MY_DOUBT(cutting_vertex_valence[i] < 2 && cutting_vertex_valence[i] != 0, "CVV Error " << i << " "
			<< cutting_vertex_valence[i] << " " << (quad_mesh.is_boundary(quad_mesh.vertex_handle(i)) ? "T" : "F"));

		if (!(quad_mesh.is_boundary(quad_mesh.vertex_handle(i)) && cutting_vertex_valence[i] > 2 && i >= vert_corner.size())) continue;
		
		for (int j = 0; j < quad_v_neighbor[i].next_seg.size(); j++)
		{
			int j_seg = quad_v_neighbor[i].next_seg[j];
			if (j_seg < 0 || j_seg >= n_bsegs) continue;
			
			boundary_cutting_nodes.emplace_back(i, j_seg);
			break;
		}
	}
}

void ParaQuadChartDecomposition::chart_decomposition()
{
	int n_charts = chart_faces.size();
	atlas.resize(n_charts);
	v_chart2quad.resize(n_charts);

	std::vector<bool> face_visited(quad_mesh.n_faces(), false);

	for (int i = 0; i < n_charts; i++)
	{
		auto& comp = atlas[i];
		auto& v_comp2quad = v_chart2quad[i];
		std::vector<int> vertex_index(quad_mesh.n_vertices(), -1);

		for (int fid : chart_faces[i])
		{
			auto f_h = quad_mesh.face_handle(fid);

			std::vector<OpenMesh::VertexHandle> v_f;
			for (auto fv : quad_mesh.fv_range(f_h))
			{
				if (vertex_index[fv.idx()] == -1)
				{
					vertex_index[fv.idx()] = comp.n_vertices();
					comp.add_vertex(quad_mesh.point(fv));
					v_comp2quad.push_back(fv.idx());
				}
				v_f.emplace_back(vertex_index[fv.idx()]);
			}
			comp.add_face(v_f);
		}
	}

	std::cout << "Components " << n_charts << std::endl;

	int n_verts = 0;
	for (int i = 0; i < n_charts; i++)
	{
		n_verts += atlas[i].n_vertices();
	}

	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
	OpenMesh::HPropHandleT<int> hvt_index;

	quad_mesh.add_property(mvt_list, "mvt_list");
	quad_mesh.add_property(hvt_index, "hvt_index");

	quad_mesh.property(mvt_list).clear();
	quad_mesh.property(mvt_list).reserve(n_verts);

	int v_offset = 0;
	for (int i = 0; i < n_charts; i++)
	{
		OpenMesh::Vec3d trans(chart_translation[i][0], chart_translation[i][1], 0.0);
		for (const auto v_h : atlas[i].vertices())
		{
			OpenMesh::Vec3d point_v = atlas[i].point(v_h);
			if (chart_flipped[i])
			{
				point_v[0] = -atlas[i].point(v_h)[1];
				point_v[1] = atlas[i].point(v_h)[0];
			}
			point_v += trans;
			packing_mesh.add_vertex(point_v);
			quad_mesh.property(mvt_list).emplace_back(point_v[0], point_v[1]);
		}

		for (const auto f_h : atlas[i].faces())
		{
			std::vector<OpenMesh::VertexHandle> v_fcur;
			for (auto fv : atlas[i].fv_range(f_h))
			{
				v_fcur.emplace_back(packing_mesh.vertex_handle(fv.idx() + v_offset));
			}
			packing_mesh.add_face(v_fcur);

			for (auto fh : atlas[i].fh_range(f_h))
			{
				int v0 = atlas[i].from_vertex_handle(fh).idx();
				int v1 = atlas[i].to_vertex_handle(fh).idx();

				auto quad_h = quad_mesh.find_halfedge(quad_mesh.vertex_handle(v_chart2quad[i][v0]), quad_mesh.vertex_handle(v_chart2quad[i][v1]));
				quad_mesh.property(hvt_index, quad_h) = v1 + v_offset;
			}
		}

		v_offset += atlas[i].n_vertices();
	}
}

void ParaQuadChartDecomposition::extract_cut()
{
	std::vector<bool> already_cut(quad_mesh.n_edges());

	for (auto e_h : quad_mesh.edges())
	{
		already_cut[e_h.idx()] = quad_mesh.is_boundary(e_h);
	}

	for (auto e_h : quad_mesh.edges())
	{
		if (!quad_edge_inner_label[e_h.idx()] && !already_cut[e_h.idx()])
		{
			auto h0 = quad_mesh.halfedge_handle(e_h, 0);

			OpenMesh::HalfedgeHandle h_iter;
			std::vector<int> h_before, h_after;
			std::vector<std::pair<int, int>> v_before, v_after;

			h_iter = h0;
			while (!already_cut[h_iter.idx() / 2])
			{
				h_before.push_back(h_iter.idx());
				h_iter = quad_mesh.prev_halfedge_handle(h_iter);

				int v_k = 1;
				while (quad_edge_inner_label[h_iter.idx() / 2])
				{
					v_k++;
					h_iter = quad_mesh.prev_halfedge_handle(quad_mesh.opposite_halfedge_handle(h_iter));
				}

				int to_v = quad_mesh.to_vertex_handle(h_iter).idx();
				if (v_k != 2 || cutting_vertex_valence[to_v] != 2)
				{
					v_before.emplace_back(to_v, v_k);
				}
			}
			h_before.push_back(h_iter.idx());

			h_iter = h0;
			while (!already_cut[h_iter.idx() / 2])
			{
				h_after.push_back(h_iter.idx());
				h_iter = quad_mesh.next_halfedge_handle(h_iter);

				int v_k = 1;
				while (quad_edge_inner_label[h_iter.idx() / 2])
				{
					v_k++;
					h_iter = quad_mesh.next_halfedge_handle(quad_mesh.opposite_halfedge_handle(h_iter));
				}

				int from_v = quad_mesh.from_vertex_handle(h_iter).idx();
				if (v_k != 2 || cutting_vertex_valence[from_v] != 2)
				{
					v_after.emplace_back(from_v, v_k);
				}
			}

			std::vector<int> h_path = h_before;
			std::reverse(h_path.begin(), h_path.end());

			h_path.insert(h_path.end(), h_after.begin() + 1, h_after.end());

			std::vector<std::pair<int, int>> v_path = v_before;
			std::reverse(v_path.begin(), v_path.end());
			v_path.insert(v_path.end(), v_after.begin(), v_after.end());

			for (int h : h_path)
			{
				already_cut[h / 2] = true;
			}
			
			cut_path_seg.emplace_back();
			int v_idx = 1;
			for (int i = 1; i < h_path.size(); i++)
			{
				auto h = quad_mesh.halfedge_handle(h_path[i]);
				if (quad_mesh.to_vertex_handle(h).idx() == v_path[v_idx].first)
				{
					cut_path_seg.back().push_back(cutting_edge_seg[h.idx() / 2]);
					v_idx++;
				}
			}

			MY_DOUBT(v_idx != v_path.size(), "CPC Size Error.");

			cut_path_h.emplace_back(std::move(h_path));
			cut_path_v.emplace_back(std::move(v_path));
		}
	}

	MY_DOUBT(cut_path_h.size() != cut_path_v.size() || cut_path_h.size() != cut_path_seg.size(), "Path Size Error.");
}

std::set<int> ParaQuadChartDecomposition::sample_boundary_vertices()
{
	std::vector<int> boundary_vertices;
	for (auto h : mesh.halfedges())
	{
		if (!mesh.is_boundary(h)) continue;

		boundary_vertices.push_back(mesh.to_vertex_handle(h).idx());
	}

	std::set<int> res;

	while (res.size() < 15)
	{
		int new_res = std::rand() % boundary_vertices.size();
		res.insert(boundary_vertices[new_res]);
	}

	return res;
}
