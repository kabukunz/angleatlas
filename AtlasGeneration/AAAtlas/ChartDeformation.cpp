#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
#include "ChartDeformation.h"

#include "Common\CommonFunctions.h"
#include "Optimization\Numeric\FastMath.h"
#include "Mosek\ConvexQuadOptimization.h"

#include "ChartUntangle.h"

#include <QTime>

ChartDeformation::ChartDeformation(PolySquareDeformation& _parent)
	: parent(_parent)
{
}

ChartDeformation::~ChartDeformation()
{
}

void ChartDeformation::build_chart()
{
	const auto& para = parent.para;
	n_boundary_edges = boundary_h_meshid.size();

	boundary_h_vert.resize(n_boundary_edges);
	boundary_h_tag.resize(n_boundary_edges);
	boundary_h_svec.resize(n_boundary_edges);

	double total_blen = 0.0;
	for (double len : boundary_h_len0)
	{
		total_blen += len;
	}
	total_area = 0.0;
	total_boundary_length = total_blen;
	n_vertices = 0;
	n_edges = 0;
	n_faces = 0;
}

void ChartDeformation::build_vertices()
{
	const auto& para = parent.para;
	const auto& v_chart = parent.v_chart;

	uv_x.resize(n_vertices * 2);
	for (int j = 0; j < n_boundary_edges; j++)
	{
		auto h_h = para.halfedge_handle(boundary_h_meshid[j]);

		boundary_h_vert[j].first = v_chart[para.from_vertex_handle(h_h).idx()].second;
		boundary_h_vert[j].second = v_chart[para.to_vertex_handle(h_h).idx()].second;
	}

	boundary_v_k.assign(n_boundary_edges, 0);
}

void ChartDeformation::build_faces()
{
	const auto& v_chart = parent.v_chart;

	double total_uv_area = 0.0;
	chart_face_info.reserve(mesh_faces.size());
	for (int fid : mesh_faces)
	{
		chart_face_info.push_back(parent.mesh_face_info[fid]);

		auto& finfo = chart_face_info.back();
		finfo.uv0 = v_chart[finfo.uv0].second;
		finfo.uv1 = v_chart[finfo.uv1].second;
		finfo.uv2 = v_chart[finfo.uv2].second;

		auto p1 = get_vec(finfo.uv0, finfo.uv1, uv_x);
		auto p2 = get_vec(finfo.uv0, finfo.uv2, uv_x);

		total_uv_area += 0.5 * (p1[0] * p2[1] - p1[1] * p2[0]);
	}

	z_pos = (chart_face_info[0].normal_towards == 1);

	std::cout << "==============scale:" << std::sqrt(total_uv_area / total_area) << std::endl;
	scale_finfo(std::sqrt(total_uv_area / total_area));
}

double ChartDeformation::interior_angle(int next_boundary_h, const std::vector<double>& uv)
{
	const auto& para = parent.para;
	const auto& v_chart = parent.v_chart;

	auto h_h = para.halfedge_handle(boundary_h_meshid[next_boundary_h]);

	double angle = 0.0;
	int v0 = boundary_h_vert[next_boundary_h].first;
	auto h_iter = para.opposite_halfedge_handle(h_h);
	while (!para.is_boundary(h_iter))
	{
		int v1 = v_chart[para.from_vertex_handle(h_iter).idx()].second;
		h_iter = para.next_halfedge_handle(h_iter);
		int v2 = v_chart[para.to_vertex_handle(h_iter).idx()].second;

		OpenMesh::Vec2d p1 = get_vec(v0, v1, uv);
		OpenMesh::Vec2d p2 = get_vec(v0, v2, uv);
		angle += CommonFunctions::vec_angle_atan2(p1, p2);
		
		h_iter = para.opposite_halfedge_handle(h_iter);
	}

	return std::abs(angle);
}

ChartDeformation::corner ChartDeformation::next_corner(corner it)
{
	if (it->status & 4) return it;
	auto ith = boundary_corners.begin();
	auto itt = boundary_corners.end();

	do 
	{
		++it;
		if (it == itt) it = ith;
	} while (it->status & 4);

	return it;
}

ChartDeformation::corner ChartDeformation::prev_corner(corner it)
{
	if (it->status & 4) return it;
	auto ith = boundary_corners.begin();
	auto itt = boundary_corners.end();

	do
	{
		if (it == ith) it = itt;
		--it;
	} while (it->status & 4);

	return it;
}

void ChartDeformation::calc_seglength(corner it, const std::vector<double>& bh_len)
{
	it->seg_length = 0.0;
	auto next_it = next_corner(it);
	for (int bh = it->bvid; bh != next_it->bvid; bh = boundary(bh + 1))
	{
		it->seg_length += bh_len[bh];
	}
}

void ChartDeformation::update_segs_len(int uv_id)
{
	std::vector<double> bh_len;
	if (uv_id >= 0)
	{
		const auto& uv = (uv_id < uv_x_vec.size()) ? uv_x_vec[uv_id] : uv_x;

		bh_len.resize(n_boundary_edges);
		for (int i = 0; i < n_boundary_edges; i++)
		{
			int v0 = boundary_h_vert[i].first;
			int v1 = boundary_h_vert[i].second;

			bh_len[i] = get_vec(v0, v1, uv).norm();
		}
	}

	const std::vector<double>& bh_len_ref = (uv_id == -1) ? boundary_h_len0 : bh_len;
	for (auto it = boundary_corners.begin(); it != boundary_corners.end(); it++)
	{
		calc_seglength(it, bh_len_ref);
	}
}

void ChartDeformation::deformation(std::string str1_)
{
	QTime timer;
	kernel_width = parent.kernel_width * (double)n_boundary_edges / 1000.0;
	kernel_width = std::min(std::max(min_sigma, kernel_width), max_sigma);

	align_alpha = 0.3;
	amips_alpha = 0.5;
	energy_amips_threshold = 0;

	auto align_energy = parent.amips_exp ? &ChartDeformation::calc_align_energy<true> : &ChartDeformation::calc_align_energy<false>;
	auto align_deformation = parent.amips_exp ? &ChartDeformation::calc_align_deformation<true> : &ChartDeformation::calc_align_deformation<false>;
	auto inner_deformation = parent.amips_exp ? &ChartDeformation::calc_inner_deformation<true> : &ChartDeformation::calc_inner_deformation<false>;
	auto final_deformation = parent.amips_exp ? &ChartDeformation::calc_final_deformation<true> : &ChartDeformation::calc_final_deformation<false>;
	auto align_energy1 = parent.amips_exp ? &ChartDeformation::calc_align_energy1<true> : &ChartDeformation::calc_align_energy1<false>;
	auto align_deformation1 = parent.amips_exp ? &ChartDeformation::calc_align_deformation1<true> : &ChartDeformation::calc_align_deformation1<false>;
	auto line_energy = parent.amips_exp ? &ChartDeformation::calc_line_energy<true> : &ChartDeformation::calc_line_energy<false>;
	auto line_deformation = parent.amips_exp ? &ChartDeformation::calc_line_deformation<true> : &ChartDeformation::calc_line_deformation<false>;
	auto output_energy = [&](bool changed = false)
	{
		auto default_precision = std::cout.precision();
		std::cout << std::setprecision(4) << std::scientific << std::setw(14) << energy_lambda << std::setw(15) << energy_align;
		if (energy_amips >= 1.0e4) std::cout << std::setw(15) << energy_amips;
		else std::cout << std::setprecision(8) << std::fixed << std::setw(15) << energy_amips;
		std::cout << std::setprecision(8) << std::fixed << std::setw(15) << max_angle_align;
		std::cout << std::setprecision(3) << std::fixed << std::setw(10) << timer.elapsed() / 1000.0;
		std::cout << std::setw(10) << (changed ? "True " : "False");
		std::cout << std::resetiosflags(std::ios::floatfield) << std::setprecision(default_precision) << std::endl;
	};

	build_faces();
	build_scaffold(parent.para);
	//scaffold

	timer.start();
	calc_boundary_directions();     
	calc_global_rotation();         
	uvplus2uv();
	uv_x_vec.push_back(uv_x);       

	find_corners();
//	tag is invalid from here
	modify_short_segments();
	conflicting_segments();
	move_corners();

	tag_from_segments();
//	tag is valid again

	hessian_preparation();
	energy_lambda = 1.0;
	(this->*align_energy)();
	std::cout << std::setw(4) << "Iter" << std::setw(18) << "Lambda"
		<< std::setw(15) << "E_align" << std::setw(15) << "E_iso" 
		<< std::setw(15) << "Max Diff" << std::setw(10) << "Time" << std::setw(10) << "Changed" << std::endl;
	std::cout << std::setw(4) << 0 << "    ";
	output_energy();
	energy_lambda = 1.0 / (energy_align + 1e-8);

	energy_amips_threshold = energy_amips + 0.1;
	std::vector<double> uv_x_plus_pre = uv_x_plus;
	double max_angle_align_pre = 0;
	int axis_iter_num = 5;
	for (int i = 0; i < axis_iter_num; i++)
	{
		std::cout << std::setw(4) << i + 1;

		(this->*align_deformation)(500);
		std::cout << "-A";
		
		uvplus2uv();
		uv_x_vec.push_back(uv_x);
		CM_inner_deformation(5);

		std::cout << "-I";
		
		uvplus2uv();
		uv_x_vec.push_back(uv_x);
		(this->*align_energy)();
		
		bool seg_changed = modify_short_segments();
		output_energy(seg_changed);

		if (energy_amips > energy_amips_threshold)
		{
			if (i > 0)
			{
				uv_x_plus = uv_x_plus_pre;
				uvplus2uv();
				int vec_num = uv_x_vec.size();
				uv_x_vec.erase(uv_x_vec.begin() + vec_num - 1);
				uv_x_vec.erase(uv_x_vec.begin() + vec_num - 2);
				update_segs_len(uv_x_vec.size());
				break;
			}
			else
			{
				break;
			}
		}
		if ((!seg_changed && max_angle_align < 0.3) || std::fabs(max_angle_align_pre - max_angle_align) < 1e-6) break;

		update_segs_len(uv_x_vec.size());            
		tag_from_segments();

		energy_lambda *= seg_changed ? 3.0 : 6.0;
		uv_x_plus_pre = uv_x_plus;
		max_angle_align_pre = max_angle_align;
		if (i == 4 && max_angle_align > 1)
		{
			axis_iter_num = 10;
		}
	}
	if (energy_amips > 10)
	{
		uv_x_plus = uv_x_plus_pre;
		uvplus2uv();
	}

	/////////////////////////////////////////////////////////////////角点对齐
	std::vector<OpenMesh::VertexHandle> v2newv;
	v2newv.resize(n_vertices);
	for (int i = 0; i < n_vertices; i++)
	{
		OpenMesh::Vec3d point_temp(uv_x[2 * i], uv_x[2 * i + 1], 0);
		auto newv = mesh_temp.add_vertex(point_temp);
		v2newv[i] = newv;
	}
	for (auto finfo : chart_face_info)
	{
		mesh_temp.add_face(v2newv[finfo.uv0], v2newv[finfo.uv1], v2newv[finfo.uv2]);
	}
	
	//记录Corner点
	for (auto it = boundary_corners.begin(); it != boundary_corners.end(); it++)
	{
		int v0 = boundary_h_vert[it->bvid].first;
		origin_corner_list.push_back(v0);
	}

	build_scaffold(mesh_temp);
	angle_lambda = 0.6;
	for (int iter = 0; iter < 8; iter++)
	{
		std::cout << "=================================================iter " << iter << std::endl;
		get_align_direction_corner();
		uv_x_vec.clear();
		uv_x_vec.push_back(uv_x);
		timer.start();
		calc_boundary_directions();

		hessian_preparation();
		energy_lambda = 1.0;
		(this->*align_energy1)();
		std::cout << std::setw(4) << "Iter" << std::setw(18) << "Lambda"
			<< std::setw(15) << "E_align" << std::setw(15) << "E_iso"
			<< std::setw(15) << "Max Diff" << std::setw(10) << "Time" << std::setw(10) << "Changed" << std::endl;
		std::cout << std::setw(4) << 0 << "    ";
		output_energy();
		energy_lambda = 1.0 / (energy_align + 1e-8);
		angle_lambda = 0.6 + iter*0.1;
		if (angle_lambda > 1)
		{
			angle_lambda = 1;
		}

		uv_x_plus_pre = uv_x_plus;
		max_angle_align_pre = 0;
		for (int i = 0; i < 15; i++)
		{
			std::cout << std::setw(4) << i + 1;

			(this->*align_deformation1)(500);
			std::cout << "-A";
			uvplus2uv();
			uv_x_vec.push_back(uv_x);
			CM_inner_deformation(5);
			std::cout << "-I";
			uvplus2uv();
			uv_x_vec.push_back(uv_x);
			(this->*align_energy1)();
			bool seg_changed = false;
			output_energy(seg_changed);
			if (energy_amips > energy_amips_threshold)
			{
				if (i > 0)
				{
					uv_x_plus = uv_x_plus_pre;
					uvplus2uv();
					int vec_num = uv_x_vec.size();
					uv_x_vec.erase(uv_x_vec.begin() + vec_num - 1);
					uv_x_vec.erase(uv_x_vec.begin() + vec_num - 2);
					update_segs_len(uv_x_vec.size());
					break;
				}
				else
				{
					break;
				}
			}
			if (max_angle_align < 0.03 || std::fabs(max_angle_align_pre - max_angle_align) < 1e-6) break;
			energy_lambda *= 6.0;
			update_segs_len(uv_x_vec.size());            //这里用到uv_x了，注意      不需要修改
			uv_x_plus_pre = uv_x_plus;
			max_angle_align_pre = max_angle_align;
		}
		if (energy_amips > 10)
		{
			uv_x_plus = uv_x_plus_pre;
			uvplus2uv();
		}
	}

	double max_angle_angle_last = max_angle_align;

	bool is_bad_result = 0;
	is_bad_result = 1;
	angle_lambda = 0;               
///////////////////////////////////////////////////////////////////Lloyd迭代
	get_align_direction();
	uv_x_vec.clear();
	uv_x_vec.push_back(uv_x);
	timer.start();
	calc_boundary_directions();

	hessian_preparation();
	energy_lambda = 1.0;
	(this->*align_energy1)();
	is_bad_result = 0;
	if (max_angle_angle_last >= 0.1 && max_angle_align <3)
	{
		is_bad_result = 1;
		std::cout << std::setw(4) << "Iter" << std::setw(18) << "Lambda"
			<< std::setw(15) << "E_align" << std::setw(15) << "E_iso"
			<< std::setw(15) << "Max Diff" << std::setw(10) << "Time" << std::setw(10) << "Changed" << std::endl;
		std::cout << std::setw(4) << 0 << "    ";
		output_energy();
		energy_lambda = 1.0 / (energy_align + 1e-8);

		uv_x_plus_pre = uv_x_plus;
		max_angle_align_pre = 0;
		for (int i = 0; i < 15; i++)
		{
			std::cout << std::setw(4) << i + 1;

			(this->*align_deformation1)(500);
			std::cout << "-A";
			//修改
			uvplus2uv();
			uv_x_vec.push_back(uv_x);
			CM_inner_deformation(5);

			std::cout << "-I";
			//修改
			uvplus2uv();
			uv_x_vec.push_back(uv_x);
			(this->*align_energy1)();


			bool seg_changed = false;
			output_energy(seg_changed);
			if (energy_amips > energy_amips_threshold)
			{
				if (i > 0)
				{
					uv_x_plus = uv_x_plus_pre;
					uvplus2uv();
					int vec_num = uv_x_vec.size();
					uv_x_vec.erase(uv_x_vec.begin() + vec_num - 1);
					uv_x_vec.erase(uv_x_vec.begin() + vec_num - 2);
					update_segs_len(uv_x_vec.size());
					break;
				}
				else
				{
					break;
				}
			}
			if (max_angle_align < /*0.3*/0.03 || std::fabs(max_angle_align_pre - max_angle_align) < 1e-6) break;

			energy_lambda *= 6.0;

			update_segs_len(uv_x_vec.size());            //这里用到uv_x了，注意      不需要修改
			uv_x_plus_pre = uv_x_plus;
			max_angle_align_pre = max_angle_align;
		}
		if (energy_amips > 10)
		{
			uv_x_plus = uv_x_plus_pre;
			uvplus2uv();
		}
		

		Mesh mesh_temp1;
		std::vector<OpenMesh::VertexHandle> v2newv1;
		v2newv1.resize(n_vertices);
		for (int i = 0; i < n_vertices; i++)
		{
			OpenMesh::Vec3d point_temp(uv_x[2 * i], uv_x[2 * i + 1], 0);
			auto newv = mesh_temp1.add_vertex(point_temp);
			v2newv1[i] = newv;
		}
		for (auto finfo : chart_face_info)
		{
			mesh_temp1.add_face(v2newv1[finfo.uv0], v2newv1[finfo.uv1], v2newv1[finfo.uv2]);
		}
	}

	////////////////////////////////////////////////////////////////////降低扭曲
	std::cout << "===========================================================" << std::endl;
	construct_line();
	timer.start();
	calc_boundary_directions();

	hessian_preparation();
	energy_lambda = 1.0;
	(this->*line_energy)();
	std::cout << std::setw(4) << "Iter" << std::setw(18) << "Lambda"
		<< std::setw(15) << "E_align" << std::setw(15) << "E_iso"
		<< std::setw(15) << "Max Diff" << std::setw(10) << "Time" << std::setw(10) << "Changed" << std::endl;
	std::cout << std::setw(4) << 0 << "    ";
	output_energy();

	double pre_energy_amips = 0;
	uv_x_plus_pre = uv_x_plus;
	line_lambda = 10;
	for (int i = 0; i < 20; i++)
	{
		std::cout << std::setw(4) << i + 1;
		(this->*line_deformation)(500);
		std::cout << "-A";
		uvplus2uv();
		uv_x_vec.push_back(uv_x);
		CM_inner_deformation(5);
		std::cout << "-I";
		uvplus2uv();
		uv_x_vec.push_back(uv_x);
		(this->*line_energy)();
		bool seg_changed = false;
		output_energy(seg_changed);
		if (std::fabs(energy_amips - pre_energy_amips)/ energy_amips < 1e-6)
		{
			break;
		}
		pre_energy_amips = energy_amips;
		uv_x_plus_pre = uv_x_plus;
	}

	//////////////////////////////////////////////////////////////////再次角点对齐

	std::vector<std::vector<double>> uv_x_vec_lloyd = uv_x_vec;
	bool first_wrong = 0;
	uv_x_vec.clear();
	uv_x_vec.push_back(uv_x);
	timer.start();
	calc_boundary_directions();

	hessian_preparation();
	energy_lambda = 1.0;
	(this->*align_energy1)();
	std::cout << std::setw(4) << "Iter" << std::setw(18) << "Lambda"
		<< std::setw(15) << "E_align" << std::setw(15) << "E_iso"
		<< std::setw(15) << "Max Diff" << std::setw(10) << "Time" << std::setw(10) << "Changed" << std::endl;
	std::cout << std::setw(4) << 0 << "    ";
	output_energy();
	energy_lambda = 1.0 / (energy_align + 1e-8);

	uv_x_plus_pre = uv_x_plus;
	max_angle_align_pre = 0;
	for (int i = 0; i < 15; i++)
	{
		std::cout << std::setw(4) << i + 1;

		(this->*align_deformation1)(500);
		std::cout << "-A";
		uvplus2uv();
		uv_x_vec.push_back(uv_x);
		CM_inner_deformation(5);
		std::cout << "-I";
		uvplus2uv();
		uv_x_vec.push_back(uv_x);
		(this->*align_energy1)();
		bool seg_changed = false;
		output_energy(seg_changed);
		if (energy_amips > energy_amips_threshold)
		{
			if (i > 0)
			{
				uv_x_plus = uv_x_plus_pre;
				uvplus2uv();
				int vec_num = uv_x_vec.size();
				uv_x_vec.erase(uv_x_vec.begin() + vec_num - 1);
				uv_x_vec.erase(uv_x_vec.begin() + vec_num - 2);
				update_segs_len(uv_x_vec.size());
				break;
			}
			else
			{
				first_wrong = 1;
				break;
			}
		}
		if (max_angle_align < 0.03 || std::fabs(max_angle_align_pre - max_angle_align) < 1e-6) break;
		energy_lambda *= 6.0;
		update_segs_len(uv_x_vec.size());            
		uv_x_plus_pre = uv_x_plus;
		max_angle_align_pre = max_angle_align;
	}
	if (energy_amips > 10 && first_wrong)
	{
		uv_x_plus = uv_x_plus_pre;
		uvplus2uv();
		uv_x_vec = uv_x_vec_lloyd;
	}
	uv_x_vec.push_back(uv_x);
	uv_x_vec.push_back(uv_x);
	uv_x_vec.push_back(uv_x);

}

bool ChartDeformation::modify_short_segments(double thres_factor /*= 1.0*/)
{
	double avg_len0 = total_boundary_length / n_boundary_edges;
	double length_thres = thres_factor * kernel_width * avg_len0;
	struct queue_entry
	{
		corner it;
		double len;
		double score;
	};
	struct queue_comparer
	{
		constexpr bool operator ()(const queue_entry& a, const queue_entry& b) const
		{
			return (a.score == b.score) ? (a.len > b.len) : (a.score < b.score);
		}
	};
	std::priority_queue<queue_entry, std::vector<queue_entry>, queue_comparer> short_segments;
	auto angle_score = [&](int bv, int vk)
	{
		double t = (double)(2 - vk) - boundary_v_sangle[bv];
		return std::abs(t);
	};
	auto add_segment = [&](const corner& it)
	{
		if (it->seg_length >= length_thres) return;
		int vp = it->bvid;
		int vn = next_corner(it)->bvid;
		int vk0_p = boundary_v_k[vp];
		int vk0_n = boundary_v_k[vn];
		if ((vk0_p > 0) == (vk0_n > 0)) return;

		bool new_at_p = (std::abs(vk0_p) > std::abs(vk0_n));
		int vk1_p = new_at_p ? (vk0_p + vk0_n) : 0;
		int vk1_n = new_at_p ? 0 : (vk0_p + vk0_n);

		double score0 = angle_score(vp, vk0_p) + angle_score(vn, vk0_n);
		double score1 = angle_score(vp, vk1_p) + angle_score(vn, vk1_n);

		short_segments.emplace(queue_entry{ it, it->seg_length, score0 - score1 - 2.0 * it->seg_length / avg_len0 / kernel_width });
	};

	for (auto it = boundary_corners.begin(); it != boundary_corners.end(); it++) add_segment(it);

	while (!short_segments.empty())
	{
		auto it = short_segments.top().it;
		double len0 = short_segments.top().len;
		short_segments.pop();

		if ((it->status & 4) || it->seg_length != len0) continue;

		auto cp = prev_corner(it);
		auto cn = next_corner(it);
		int vk_p = boundary_v_k[it->bvid];
		int vk_n = boundary_v_k[cn->bvid];

		if ((vk_p > 0) == (vk_n > 0)) continue;

		if (vk_p + vk_n == 0)
		{
			boundary_v_k[it->bvid] = 0;
			boundary_v_k[cn->bvid] = 0;
			cp->seg_length += it->seg_length + cn->seg_length;

			it->status |= 4;
			cn->status |= 4;

			add_segment(cp);
		}
		else if (std::abs(vk_p) > std::abs(vk_n))
		{
			boundary_v_k[it->bvid] = vk_p + vk_n;
			boundary_v_k[cn->bvid] = 0;

			it->seg_length += cn->seg_length;

			it->status = ((it->status & 4) | (cn->status & 3));
			cn->status |= 4;
			add_segment(it);
		}
		else if (std::abs(vk_p) < std::abs(vk_n))
		{
			boundary_v_k[it->bvid] = 0;
			boundary_v_k[cn->bvid] = vk_p + vk_n;

			cp->seg_length += it->seg_length;

			it->status |= 4;
			add_segment(cp);
		}
	}

	for (auto it = boundary_corners.begin(); it != boundary_corners.end(); it++)
	{
		if ((it->status & 4) || it->seg_length >= length_thres) continue;

		auto cp = prev_corner(it);
		auto cn = next_corner(it);
		int vk_p = boundary_v_k[it->bvid];
		int vk_n = boundary_v_k[cn->bvid];

		if ((vk_p >= 0) || (vk_n >= 0)) continue;

		int vk_new = vk_p + vk_n;
		double score_p = std::abs(boundary_v_angle0[it->bvid] - (double)(2 - vk_new) * M_PI_2);
		double score_n = std::abs(boundary_v_angle0[cn->bvid] - (double)(2 - vk_new) * M_PI_2);
		if (score_p < score_n)
		{
			boundary_v_k[it->bvid] = vk_new;
			boundary_v_k[cn->bvid] = 0;

			it->seg_length += cn->seg_length;

			it->status = ((it->status & 4) | (cn->status & 3));
			cn->status |= 4;
		}
		else
		{
			boundary_v_k[it->bvid] = 0;
			boundary_v_k[cn->bvid] = vk_new;

			cp->seg_length += it->seg_length;

			it->status |= 4;
		}
	}

	bool has_seg_deleted = false;
	for (auto it = boundary_corners.begin(); it != boundary_corners.end();)
	{
		if (it->status & 4)
		{
			it = boundary_corners.erase(it);
			has_seg_deleted = true;
		}
		else ++it;
	}

	return has_seg_deleted;
}

bool ChartDeformation::conflicting_segments()
{
	bool no_neg_quad = true;

	double length_thres = kernel_width * total_boundary_length / n_boundary_edges;

	for (auto it = boundary_corners.begin(); it != boundary_corners.end(); it++)
	{
		if ((it->status & 4) || it->seg_length >= length_thres) continue;

		auto cp = prev_corner(it);
		auto cn = next_corner(it);
		int vk_p = boundary_v_k[it->bvid];
		int vk_n = boundary_v_k[cn->bvid];

		if ((vk_p != 1) || (vk_n != 1)) continue;

		if (cp->seg_length > cn->seg_length)
		{
			double thres = std::min(length_thres, cp->seg_length / 2.0);
			double dist = it->seg_length;

			int cur_id = it->bvid;
			while (dist < thres)
			{
				cur_id = boundary(cur_id - 1);
				dist += boundary_h_len0[cur_id];
			}
			boundary_v_k[it->bvid] = 0;
			it->bvid = cur_id;
			boundary_v_k[it->bvid] = 1;

			calc_seglength(cp, boundary_h_len0);
			calc_seglength(it, boundary_h_len0);
		}
		else
		{
			double thres = std::min(length_thres, cn->seg_length / 2.0);
			double dist = it->seg_length;

			int cur_id = cn->bvid;
			while (dist < thres)
			{
				dist += boundary_h_len0[cur_id];
				cur_id = boundary(cur_id + 1);
			}
			boundary_v_k[cn->bvid] = 0;
			cn->bvid = cur_id;
			boundary_v_k[cn->bvid] = 1;

			calc_seglength(cn, boundary_h_len0);
			calc_seglength(it, boundary_h_len0);
		}
	}

	for (auto it = boundary_corners.begin(); it != boundary_corners.end(); it++)
	{
		int k0 = boundary_v_k[it->bvid];
		no_neg_quad = no_neg_quad && (k0 <= 2);
		if (k0 <= 1 || k0 > 2) continue;
		
		auto cp = prev_corner(it);
		auto cn = next_corner(it);
		if (cp->seg_length > it->seg_length)
		{
			double thres = std::min(length_thres, cp->seg_length / 2.0);
			double dist = 0.0;
			auto new_corner = boundary_corners.emplace(it);
			new_corner->status = (cp->status - 1) & 3;

			int cur_id = it->bvid;
			while (dist < thres)
			{
				cur_id = boundary(cur_id - 1);
				dist += boundary_h_len0[cur_id];
			}
			new_corner->bvid = cur_id;

			calc_seglength(cp, boundary_h_len0);
			calc_seglength(new_corner, boundary_h_len0);

			boundary_v_k[it->bvid] = 1;
			boundary_v_k[new_corner->bvid] = 1;
		}
		else
		{
			double thres = std::min(length_thres, it->seg_length / 2.0);
			double dist = 0.0;
			auto new_corner = boundary_corners.emplace(it);
			new_corner->bvid = it->bvid;
			new_corner->status = (cp->status - 1) & 3;

			int cur_id = it->bvid;
			while (dist < thres)
			{
				dist += boundary_h_len0[cur_id];
				cur_id = boundary(cur_id + 1);
			}
			it->bvid = cur_id;

			calc_seglength(it, boundary_h_len0);
			calc_seglength(new_corner, boundary_h_len0);

			boundary_v_k[it->bvid] = 1;
			boundary_v_k[new_corner->bvid] = 1;
		}
	}

	if (!no_neg_quad) std::cout << "Unhandled Conflicts" << std::endl;
	return no_neg_quad;
}

void ChartDeformation::tag_from_segments()
{
	for (auto it = boundary_corners.begin(); it != boundary_corners.end(); it++)
	{
		auto next_it = next_corner(it);

		for (int bh = it->bvid; bh != next_it->bvid; bh = boundary(bh + 1))
		{
			boundary_h_tag[bh] = (it->status & 3);
		}
	}
}

void ChartDeformation::built_segments()
{
	segments.reserve(boundary_corners.size());
	for (auto cur_it = boundary_corners.begin(); cur_it != boundary_corners.end(); cur_it++)
	{
		segments.emplace_back();

		segments.back().begin = cur_it->bvid;
		segments.back().end = next_corner(cur_it)->bvid;
		segments.back().tag = boundary_h_tag[cur_it->bvid];
		segments.back().length0 = cur_it->seg_length;
		segments.back().size = boundary(segments.back().end - segments.back().begin);
	}
}

int ChartDeformation::boundary(int bid)
{
	return CommonFunctions::period_id(bid, n_boundary_edges);
}

void ChartDeformation::calc_boundary_directions()
{
	fixed_uv_x.clear();

	boundary_h_svec.resize(n_boundary_edges);
	boundary_v_angle0.resize(n_boundary_edges);
	for (int j = 0; j < n_boundary_edges; j++)
	{
		int v0 = boundary_h_vert[j].first;
		int v1 = boundary_h_vert[j].second;

		boundary_h_svec[j] = get_vec(v0, v1, uv_x);
		boundary_v_angle0[j] = interior_angle(j, uv_x);

		fixed_uv_x.insert(2 * v0 + 0);
		fixed_uv_x.insert(2 * v0 + 1);
	}
	update_uv2kkt();

	boundary_v_sangle = boundary_v_angle0;
	double avg_len = total_boundary_length / n_boundary_edges;
	double gaussian_delta = kernel_width * avg_len;
	double gaussian_thres = 3.0 * gaussian_delta;
	double delta_sqr = 2.0 * gaussian_delta * gaussian_delta;

	for (int i = 0; i < 3; i++)
	{
		std::vector<OpenMesh::Vec2d> boundary_svec(n_boundary_edges);
		for (int j = 0; j < n_boundary_edges; j++)
		{
			OpenMesh::Vec2d smooth_vec = boundary_h_svec[j];
			for (int step : {1, -1})
			{
				int cur_id = j;
				double dist = boundary_h_len0[cur_id] / 2.0;
				while (dist < gaussian_thres)
				{
					cur_id = boundary(cur_id + step);
					int vert_h = (step == 1) ? boundary_h_vert[cur_id].first : boundary_h_vert[cur_id].second;

					double alpha = std::min((gaussian_thres - dist) / boundary_h_len0[cur_id], 1.0);
					double dist_b = dist + boundary_h_len0[cur_id] * alpha / 2.0;

					smooth_vec += ig::FastNegExp1(dist_b * dist_b / delta_sqr) * boundary_h_svec[cur_id] * alpha;

					dist += boundary_h_len0[cur_id];
				}
			}

			if (OpenMesh::dot(smooth_vec, boundary_h_svec[j]) < 0.0)
			{
				boundary_svec[j] = boundary_h_svec[j];
			}
			else
			{
				int k = boundary(j + 1);
				double d_angle = CommonFunctions::vec_angle_atan2(boundary_h_svec[j], smooth_vec);

				boundary_v_sangle[j] += d_angle;
				boundary_v_sangle[k] -= d_angle;
				boundary_svec[j] = smooth_vec.normalized() * boundary_h_svec[j].norm();
			}
		}
		boundary_h_svec = std::move(boundary_svec);
	}

	for (int j = 0; j < n_boundary_edges; j++) boundary_h_svec[j].normalize();
// 
// 	std::cout << "SUM ANGLE " << 2 * n_boundary_edges - std::accumulate(boundary_v_angle0.begin(), boundary_v_angle0.end(), 0.0) / M_PI_2 << std::endl;
// 	std::cout << "SUM SMOOTHED ANGLE " << 2 * n_boundary_edges - std::accumulate(boundary_v_sangle.begin(), boundary_v_sangle.end(), 0.0) / M_PI_2 << std::endl;
}

void ChartDeformation::find_corners()
{
	std::vector<double> v_angle2 = boundary_v_sangle;
	for (int j = 0; j < n_boundary_edges; j++)
	{
		auto angle = std::atan2(boundary_h_svec[j][1], boundary_h_svec[j][0]);
		int quad_angle = std::lround(angle / M_PI_2);
		angle -= (double)quad_angle * M_PI_2;

		int k = boundary(j + 1);

		v_angle2[j] -= angle;
		v_angle2[k] += angle;

		boundary_h_tag[j] = (quad_angle & 3);
	}

	for (int j = 0; j < n_boundary_edges; j++)
	{
		boundary_v_k[j] = 2 - std::lround(v_angle2[j] / M_PI_2);
	}

	boundary_corners.clear();
	for (int j = 0; j < n_boundary_edges; j++)
	{
		if (boundary_v_k[j] != 0)
		{
			boundary_corners.emplace_back();
			boundary_corners.back().bvid = j;
			boundary_corners.back().status = boundary_h_tag[j];

			//std::cout << boundary_h_vert[j].first << ", " << boundary_h_vert[j].second << ", " << boundary_v_k[j] << ", " << boundary_h_tag[j] << std::endl;
		}
	}

	update_segs_len();
}

void ChartDeformation::move_corners()
{
	double dist_thres = 2.0 * kernel_width * total_boundary_length / n_boundary_edges;

	auto k_diff = [&](int bvid, int k0) {return boundary_v_angle0[bvid] / M_PI_2 - (double)(2 - k0); };

	auto moving_condition_pos = [&](int bvid, int k0) { return boundary_v_angle0[bvid] / M_PI_2 > 1.5; };
	auto moving_condition_neg = [&](int bvid, int k0) { return std::abs(k_diff(bvid, k0)) > 0.5; };
	auto candidate_condition_pos = [&](int bvid, int k0) { return boundary_v_angle0[bvid] / M_PI_2 < 1.4; };
	auto candidate_condition_neg = [&](int bvid, int k0) { return std::abs(k_diff(bvid, k0)) < 0.4; };

	for (auto it = boundary_corners.begin(); it != boundary_corners.end(); it++)
	{
		int new_id = -1;
		auto cp = prev_corner(it);
		double min_dist = std::numeric_limits<double>::max();

		int k0 = boundary_v_k[it->bvid];
		if (boundary_v_angle0[it->bvid] > 1.75 * M_PI) continue;
		if ((k0 > 0 && !moving_condition_pos(it->bvid, k0)) || (k0 < 0 && !moving_condition_neg(it->bvid, k0))) continue;

		for (int step : {1, -1})
		{
			int cur_id = it->bvid;
			double dist = 0.0;
			double thres0 = (step == 1) ? (it->seg_length - dist_thres * 0.5) : (cp->seg_length - dist_thres * 0.5);
			thres0 = std::min(thres0, dist_thres);
			while (dist < thres0)
			{
				cur_id = boundary(cur_id + step);
				if (boundary_v_k[cur_id] != 0) break;
				if ((k0 > 0 && candidate_condition_pos(cur_id, k0)) || (k0 < 0 && candidate_condition_neg(cur_id, k0))) break;
				dist += boundary_h_len0[(step == 1) ? boundary(cur_id - 1) : cur_id];
			}

			if (boundary_v_k[cur_id] == 0 && dist < thres0 && dist < min_dist)
			{
				min_dist = dist;
				new_id = cur_id;
			}
		}

		if (new_id != -1)
		{
			boundary_v_k[new_id] = boundary_v_k[it->bvid];
			boundary_v_k[it->bvid] = 0;
			it->bvid = new_id;

			calc_seglength(cp, boundary_h_len0);
			calc_seglength(it, boundary_h_len0);
		}
	}
}

void ChartDeformation::segment_flattening()
{
	auto vf = [&](int bh) {return boundary_h_vert[bh].first; };
	auto vt = [&](int bh) {return boundary_h_vert[bh].second; };
	for (int i = 0; i < 2 * n_vertices; i++)
	{
		uv_x[i] /= goal_length;
	}

	std::vector<std::vector<double>> seg_pos(segments.size());
	for (int i = 0; i < segments.size(); i++)
	{
		auto& seg = segments[i];
		seg_pos[i].reserve(seg.size);
		double len_accu = 0.0;
		double w_coord = 0.0;
		int coord_tag = ((seg.tag & 1) ^ 1);
		for (int bh = seg.begin; bh != seg.end; bh = boundary(bh + 1))
		{
			int v0 = vf(bh);
			int v1 = vt(bh);

			w_coord += (uv_x[2 * v0 + coord_tag] + uv_x[2 * v1 + coord_tag]) / 2.0 * boundary_h_len0[bh];

			len_accu += get_vec(v0, v1, uv_x).norm();
			seg_pos[i].emplace_back(len_accu);
		}
		seg.coordinate = w_coord / seg.length0;
		for (double& pos : seg_pos[i]) pos /= len_accu;
	}

	std::vector<double> sol_x;
	mosek_flattening(sol_x);

	auto set_rounding = [&](int id, double val)
	{
		uv_x[id] = val;
		fixed_uv_x.insert(id);
	};

	for (double& val : sol_x) val = std::round(val);
	
	for (int i = 0; i < segments.size(); i++)
	{
		auto& seg = segments[i];

		if (boundary_v_k[seg.begin] != 0)
		{
			int v0 = vf(seg.begin);
			set_rounding(2 * v0 + 0, sol_x[2 * i + 0]);
			set_rounding(2 * v0 + 1, sol_x[2 * i + 1]);
		}

		int coord_tag = ((seg.tag & 1) ^ 1);
		seg.coordinate = sol_x[2 * i + coord_tag];

		for (int bh = seg.begin; bh != seg.end; bh = boundary(bh + 1))
		{
			set_rounding(2 * vf(bh) + coord_tag, seg.coordinate);
		}
	}

	for (int i = 0; i < segments.size(); i++)
	{
		auto& seg = segments[i];
		int tag = seg.tag & 1;
		double seg_length = uv_x[2 * vf(seg.end) + tag] - uv_x[2 * vf(seg.begin) + tag];
		double uv_begin = uv_x[2 * vf(seg.begin) + tag];

		int j = 0;
		int bh_end = boundary(seg.end - 1);
		for (int bh = seg.begin; bh != bh_end; bh = boundary(bh + 1))
		{
			uv_x[2 * vt(bh) + tag] = uv_begin + seg_pos[i][j++] * seg_length;
		}
	}

	update_uv2kkt();
}

bool ChartDeformation::mosek_flattening(std::vector<double>& sol_x)
{
	/*---------------------------------------------------------------

		min 1/2*x^T*Q*x + c^T*x
		s.t. blc <= A*x <= buc

		x: vector x, s.t. blx <= x <= bux
			bkx: bound type of each xi
				MSK_BK_LO: blxi <= xi <= +Inf (lower bound)
				MSK_BK_UP: -Inf <= xi <= buxi (upper bound)
				MSK_BK_FX: blxi == xi == buxi (fixed)
				MSK_BK_FR: -Inf <= xi <= +Inf (free)
				MSK_BK_RA: blxi <= xi <= buxi (ranged)
		c: vector c
		Q: sparse representation
			Q(qsubi, qsubj) = qval

		constraints:
		bkc, blc, buc: similar to bkx, blx, bux
		A: sparse representation, different from Q

	---------------------------------------------------------------*/

	std::map<int, int> bv2seg;
	for (int i = 0; i < segments.size(); i++) bv2seg[boundary_h_vert[segments[i].begin].first] = i;

	int n_segs = segments.size();

	std::vector<double> c(n_segs * 2, 0.0);
	
	std::vector<MSKidxt> qsubi(n_segs * 2);
	std::vector<MSKidxt> qsubj(n_segs * 2);
	std::vector<double> qval(n_segs * 2);
	
	std::vector<MSKboundkeye> bkx(n_segs * 2, MSK_BK_FR);
	std::vector<double> blx(n_segs * 2, -MSK_INFINITY);
	std::vector<double> bux(n_segs * 2, +MSK_INFINITY);
	
	std::vector<MSKboundkeye> bkc;
	std::vector<double> blc;
	std::vector<double> buc;
	
	std::vector<MSKlidxt> aptrb(n_segs * 2 + 1);
	std::vector<MSKidxt> asub;
	std::vector<double> aval;

	double lambda = 0.001 * total_boundary_length / n_boundary_edges;
	for (int i = 0; i < n_segs; i++)
	{
		qsubi[2 * i + 0] = 2 * i + 0;
		qsubj[2 * i + 0] = 2 * i + 0;

		qsubi[2 * i + 1] = 2 * i + 1;
		qsubj[2 * i + 1] = 2 * i + 1;

		int tag0 = segments[i].tag & 1;
		int tag1 = tag0 ^ 1;

		qval[2 * i + tag0] = lambda;
		qval[2 * i + tag1] = segments[i].length0;

		c[2 * i + tag0] = -lambda * uv_x[2 * boundary_h_vert[segments[i].begin].first + tag0];
		c[2 * i + tag1] = -segments[i].length0 * segments[i].coordinate;
	}

	int n_constraints = n_segs * 2;

	bkc.resize(n_constraints);
	blc.resize(n_constraints);
	buc.resize(n_constraints);

	std::vector<Eigen::Triplet<int>> constraints_sparse;
	constraints_sparse.reserve(n_constraints * 2);
	for (int i = 0; i < n_segs; i++)
	{
		int tag = segments[i].tag;
		int coord_tag = (tag & 1) ^ 1;
		int j = CommonFunctions::period_id(i + 1, segments.size());

		constraints_sparse.emplace_back(i, 2 * i + coord_tag, +1);
		constraints_sparse.emplace_back(i, 2 * j + coord_tag, -1);

		bkc[i] = MSK_BK_FX;
		blc[i] = 0.0;
		buc[i] = 0.0;

		int i_sign = (tag & 2) ? +1 : -1;
		constraints_sparse.emplace_back(i + n_segs, 2 * i + (tag & 1), +i_sign);
		constraints_sparse.emplace_back(i + n_segs, 2 * j + (tag & 1), -i_sign);

		double len_seg = get_vec(boundary_h_vert[segments[i].begin].first, boundary_h_vert[segments[j].begin].first, uv_x).norm();
		bkc[i + n_segs] = MSK_BK_LO;
		blc[i + n_segs] = std::max(2.0, 0.9 * len_seg);
		buc[i + n_segs] = +MSK_INFINITY;
	}

	Eigen::SparseMatrix<int> constraints_mat(n_constraints, n_segs * 2);
	constraints_mat.setFromTriplets(constraints_sparse.begin(), constraints_sparse.end());
	constraints_mat.makeCompressed();

	asub.reserve(n_constraints * 2);
	aval.reserve(n_constraints * 2);

	aptrb[0] = 0;
	for (int i = 0; i < n_segs * 2; i++)
	{
		aptrb[i + 1] = aptrb[i];
		for (Eigen::SparseMatrix<int>::InnerIterator it(constraints_mat, i); it; ++it)
		{
			asub.push_back(it.row());
			aval.push_back(it.value());

			aptrb[i + 1]++;
		}
	}

	sol_x.resize(n_segs * 2);
	return solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, sol_x);
}

void ChartDeformation::polysquare_post_deformation()
{
	uv_x_vec.push_back(uv_x);
	uv_x_vec.push_back(uv_x);
	uv_x_vec.push_back(uv_x);

	auto get_tag = [](const OpenMesh::Vec2d& vec)
	{
		uint b = (abs(vec[0]) <= abs(vec[1]));
		uint a = (vec[b] <= 0.0);

		return (a << 1) + b;
	};

	fixed_uv_x.clear();
	for (int j = 0; j < n_boundary_edges; j++)
	{
		int v_p = boundary_h_vert[boundary(j - 1)].first;
		int v_c = boundary_h_vert[j].first;
		int v_n = boundary_h_vert[j].second;

		int tag_p = get_tag(get_vec(v_p, v_c, uv_x));
		int tag_n = get_tag(get_vec(v_c, v_n, uv_x));
		
		if (tag_p != tag_n)
		{
			fixed_uv_x.insert(2 * v_c + 0);
			fixed_uv_x.insert(2 * v_c + 1);
		}
		else
		{
			int tag_uv = ((tag_p & 1) ^ 1);
			fixed_uv_x.insert(2 * v_c + tag_uv);
		}
	}
	update_uv2kkt();

	ChartDeformation::calc_final_deformation<true>(1000);
	for (int i : fixed_uv_x) uv_x[i] = std::round(uv_x[i]);
	std::cout << std::endl;

	uv_x_vec.push_back(uv_x);
	uv_x_vec.push_back(uv_x);
	uv_x_vec.push_back(uv_x);
}

void ChartDeformation::hessian_preparation()
{
	H_det.setZero();
	H_det(0, 4) = 1.0;
	H_det(0, 5) = -1.0;
	H_det(1, 3) = -1.0;
	H_det(1, 5) = 1.0;
	H_det(2, 3) = 1.0;
	H_det(2, 4) = -1.0;

	H_det(3, 1) = -1.0;
	H_det(3, 2) = 1.0;
	H_det(4, 0) = 1.0;
	H_det(4, 2) = -1.0;
	H_det(5, 0) = -1.0;
	H_det(5, 1) = 1.0;
}

void ChartDeformation::update_uv2kkt()
{
	uv2kkt.assign(2 * n_vertices_plus, -2);
	for (int i : fixed_uv_x) uv2kkt[i] = -1;

	int n_kkt = 0;
	for (int i = 0; i < 2 * n_vertices_plus; i++)
	{
		if (uv2kkt[i] == -2) uv2kkt[i] = (n_kkt++);
	}
}

void ChartDeformation::update_face_info()
{
// 	for (auto& finfo : chart_face_info)
// 	{
// 		int v0 = finfo.uv0;
// 		int v1 = finfo.uv1;
// 		int v2 = finfo.uv2;
// 
// 		auto p1 = get_vec(v0, v1, uv_x);
// 		auto p2 = get_vec(v0, v2, uv_x);
// 
// 		finfo.l2_p1 = p1.sqrnorm();
// 		finfo.l2_p2 = p2.sqrnorm();
// 		finfo.dot_p = OpenMesh::dot(p1, p2);
// 		finfo.det_p = p1[0] * p2[1] - p1[1] * p2[0];
// 	}
}

void ChartDeformation::scale_finfo(double scale)
{
	double scale2 = scale * scale;
	for (auto& finfo : chart_face_info)
	{
		finfo.l2_p1 *= scale2;
		finfo.l2_p2 *= scale2;
		finfo.dot_p *= scale2;
		finfo.det_p *= scale2;
	}

// 	total_boundary_length *= scale;
 	total_area *= scale2;
}

void ChartDeformation::check_vk()
{
	for (int j = 0; j < n_boundary_edges; j++)
	{
		auto v_h = parent.para.from_vertex_handle(parent.para.halfedge_handle(boundary_h_meshid[j]));
		double angle = interior_angle(j, uv_x);
		int quad0 = 2 - boundary_v_k[j];
		int quad1 = std::lround(angle / M_PI_2);

		if (quad0 != quad1) std::cout << "K error at " << v_h.idx() << ", " << quad1 << "/" << quad0 << std::endl;
	}
}

std::vector<OpenMesh::Vec3d> ChartDeformation::build_frame(Mesh &mesh_)
{
	std::vector<OpenMesh::Vec3d> frame_points;

	auto para_temp = mesh_;
	OpenMesh::Vec3d bmax_temp, bmin_temp;
	bmax_temp = bmin_temp = para_temp.point(para_temp.vertex_handle(0));
	for (auto v_h : para_temp.vertices())
	{
		const auto& p = para_temp.point(v_h);
		bmax_temp = bmax_temp.maximize(p);
		bmin_temp = bmin_temp.minimize(p);
	}

	OpenMesh::Vec3d half_size(bmax_temp[0] - bmin_temp[0], bmax_temp[1] - bmin_temp[1], 0);
	bmax_temp = bmax_temp + half_size/2;
	bmin_temp = bmin_temp - half_size/2;

	OpenMesh::Vec3d original_point(bmin_temp[0], bmax_temp[1], 0);
	int frame_size = 20;
	double rect_lenx = (bmax_temp[0] - bmin_temp[0])/frame_size;
	double rect_leny = (bmax_temp[1] - bmin_temp[1])/frame_size;

	frame_points.push_back(original_point);
	for (int i = 1; i < frame_size + 1; i++)
	{
		OpenMesh::Vec3d point_temp(original_point[0] + rect_lenx * i, original_point[1], 0);
		frame_points.push_back(point_temp);
	}
	original_point = frame_points[frame_points.size() - 1];
	for (int i = 1; i < frame_size + 1; i++)
	{
		OpenMesh::Vec3d point_temp(original_point[0], original_point[1] - rect_leny * i, 0);
		frame_points.push_back(point_temp);
	}
	original_point = frame_points[frame_points.size() - 1];
	for (int i = 1; i < frame_size + 1; i++)
	{
		OpenMesh::Vec3d point_temp(original_point[0] - rect_lenx * i, original_point[1], 0);
		frame_points.push_back(point_temp);
	}
	original_point = frame_points[frame_points.size() - 1];
	for (int i = 1; i < frame_size; i++)
	{
		OpenMesh::Vec3d point_temp(original_point[0], original_point[1] + rect_leny * i, 0);
		frame_points.push_back(point_temp);
	}

	return frame_points;
}

void ChartDeformation::build_boundary(Mesh &mesh_, std::vector<OpenMesh::Vec3d> &boundary_, std::vector<int> &boundary_size_, std::vector<int> &boundary_id_)
{
	auto para_temp = mesh_;
	std::vector<int> all_boundary;
	int boundary_num = 0;
	all_boundary.resize(para_temp.n_vertices());
	for (auto v_h : para_temp.vertices())
	{
		if (para_temp.is_boundary(v_h))
		{
			all_boundary[v_h.idx()] = 1;
			boundary_num++;
		}
	}

	OpenMesh::HalfedgeHandle hedge_init;
	OpenMesh::HalfedgeHandle hedge_temp;
	OpenMesh::VertexHandle vertex_init;
	OpenMesh::VertexHandle vertex_temp;
	int allboundary_size = 0;
	while (std::accumulate(all_boundary.begin(), all_boundary.end(), 0) != 0)
	{
		for (auto he_h : para_temp.halfedges())
		{
			OpenMesh::VertexHandle from_v = para_temp.from_vertex_handle(he_h);
			if (para_temp.is_boundary(he_h) && all_boundary[from_v.idx()] == 1)
			{
				hedge_init = he_h;
				vertex_init = from_v;
				break;
			}
		}

		int size_count = 0;
		Mesh::Point point_temp;
		hedge_temp = hedge_init;
		vertex_temp = vertex_init;
		boundary_.resize(boundary_num);
		do
		{
			point_temp = para_temp.point(vertex_temp);
			boundary_[allboundary_size] = point_temp;
			all_boundary[vertex_temp.idx()] = 0;
			boundary_id_.push_back(vertex_temp.idx());
			vertex_temp = para_temp.to_vertex_handle(hedge_temp);
			hedge_temp = para_temp.next_halfedge_handle(hedge_temp);
			size_count++;
			allboundary_size++;
		} while (vertex_temp != vertex_init);
		boundary_size_.push_back(size_count);
	}
}

void ChartDeformation::build_scaffold(Mesh &mesh_)
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi E;
	Eigen::MatrixXd H;
	Eigen::MatrixXd NewV;
	Eigen::MatrixXi NewF;

	auto para_temp = mesh_;
	std::vector<OpenMesh::Vec3d> frame_points;
	std::vector<OpenMesh::Vec3d> boundary_points;
	std::vector<int> boundary_size;
	std::vector<int> boundary_id;
	frame_points = build_frame(mesh_);
	build_boundary(mesh_, boundary_points, boundary_size, boundary_id);

	int count_num = 0;
	for (auto v_h : para_temp.vertices())
	{
		if (para_temp.is_boundary(v_h))
		{
			count_num++;
		}
	}
	/*std::cout << "boundary number" << boundary_id.size() << " // " << count_num << std::endl;
	std::cout << "frame number" << frame_points.size() << std::endl;*/

	std::vector<OpenMesh::Vec3d> V_points;
	std::vector<int> V_sizes;
	for (int i = 0; i < boundary_points.size(); i++)
	{
		V_points.push_back(boundary_points[i]);
	}
	for (int i = 0; i < frame_points.size(); i++)
	{
		V_points.push_back(frame_points[i]);
	}
	for (int i = 0; i < boundary_size.size(); i++)
	{
		V_sizes.push_back(boundary_size[i]);
	}
	V_sizes.push_back(frame_points.size());

	//构建V
	V.resize(V_points.size(), 2);
	for (int i = 0; i < V_points.size(); i++)
	{
		V(i, 0) = V_points[i][0];
		V(i, 1) = V_points[i][1];
	}

	//构建E
	E.resize(V.rows(), 2);
	for (int i = 0; i < E.rows(); i++)
	{
		E.row(i) << i, i + 1;
	}
	int acc_bs = 0;
	for (auto bs : V_sizes)
	{
		E(acc_bs + bs - 1, 1) = acc_bs;
		acc_bs += bs;
	}

	//构建H
	OpenMesh::Vec3d center_point(0, 0, 0);
	OpenMesh::FaceHandle face0 = para_temp.face_handle(0);
	int neighbor_num = 0;
	for (auto fv_h : para_temp.fv_range(face0))
	{
		center_point += para_temp.point(fv_h);
		neighbor_num++;
	}
	center_point = center_point / neighbor_num;
	H.resize(1, 2);
	H(0, 0) = center_point[0];
	H(0, 1) = center_point[1];

	triangulate(V, E, H, NewV, NewF);

	uv_x_plus.clear();
	uv_x_plus.resize(2*(para_temp.n_vertices() + NewV.rows() - boundary_id.size()));
	for (int i = 0; i < uv_x.size(); i++)
	{
		uv_x_plus[i] = uv_x[i];
	}

	std::vector<int> newv2v;
	newv2v.resize(NewV.rows());
	for (int i = 0; i < boundary_id.size(); i++)
	{
		newv2v[i] = boundary_id[i];
	}
	for (int i = boundary_id.size(); i < NewV.rows(); i++)
	{
		int plus_id = para_temp.n_vertices() - boundary_id.size() + i;
		newv2v[i] = plus_id;
		uv_x_plus[2 * plus_id] = NewV(i, 0);
		uv_x_plus[2 * plus_id + 1] = NewV(i, 1);
	}

	std::vector<OpenMesh::Vec3d> point_plus;
	for (int i = 0; i < uv_x_plus.size() / 2; i++)
	{
		OpenMesh::Vec3d point_temp(uv_x_plus[2 * i], uv_x_plus[2 * i + 1], 0);
		point_plus.push_back(point_temp);
	}

	Eigen::MatrixXi NewF1;
	NewF1.resize(NewF.rows(), 3);
	for (int i = 0; i < NewF.rows(); i++)
	{
		NewF1(i, 0) = newv2v[NewF(i, 0)];
		NewF1(i, 1) = newv2v[NewF(i, 1)];
		NewF1(i, 2) = newv2v[NewF(i, 2)];
	}

	std::vector<PolySquareDeformation::face_info> chart_face_info_temp;
	chart_face_info_temp.resize(NewF.rows());
	for (int i = 0; i < NewF.rows(); i++)
	{
		auto &f_info = chart_face_info_temp[i];
		f_info.uv0 = NewF1(i, 0);
		f_info.uv1 = NewF1(i, 1);
		f_info.uv2 = NewF1(i, 2);

		OpenMesh::Vec3d vec1, vec2;
		vec1 = point_plus[f_info.uv1] - point_plus[f_info.uv0];
		vec2 = point_plus[f_info.uv2] - point_plus[f_info.uv0];

		OpenMesh::Vec3d normal_temp = OpenMesh::cross(vec1, vec2);
		bool z_pos = normal_temp[2] >= 0;
		f_info.normal_towards = z_pos ? 1 : 0;

		f_info.l2_p1 = vec1.sqrnorm();
		f_info.l2_p2 = vec2.sqrnorm();
		f_info.dot_p = OpenMesh::dot(vec1, vec2);
		f_info.det_p = OpenMesh::cross(vec1, vec2).norm();
	}

	chart_face_info_plus.clear();
	for (auto f_info : chart_face_info)
	{
		chart_face_info_plus.push_back(f_info);
	}
	for (auto f_info : chart_face_info_temp)
	{
		chart_face_info_plus.push_back(f_info);
	}

	//构建total_area_plus
	total_area_plus = 0;
	for (auto f_info : chart_face_info_plus)
	{
		total_area_plus += std::abs(f_info.det_p) / 2.0;
	}

	//构建n_vertices_plus
	n_vertices_plus = uv_x_plus.size() / 2;
}

void ChartDeformation::uvplus2uv()
{
	for (int i = 0; i < n_vertices; i++)
	{
		uv_x[2 * i] = uv_x_plus[2 * i];
		uv_x[2 * i + 1] = uv_x_plus[2 * i + 1];
	}
}

void ChartDeformation::uv2uvplus()
{
	for (int i = 0; i < n_vertices; i++)
	{
		uv_x_plus[2 * i] = uv_x[2 * i];
		uv_x_plus[2 * i + 1] = uv_x[2 * i + 1];
	}
}

void ChartDeformation::get_align_direction_corner()
{
	boundary_cce_list.clear();
	boundary_cc_list.clear();
	corner_vid_list.clear();
	segment_list.clear();
	align_direction.clear();
	Mesh mesh_temp1;
	std::vector<OpenMesh::VertexHandle> v2newv;
	v2newv.resize(n_vertices);
	for (int i = 0; i < n_vertices; i++)
	{
		OpenMesh::Vec3d point_temp(uv_x[2 * i], uv_x[2 * i + 1], 0);
		auto newv = mesh_temp1.add_vertex(point_temp);
		v2newv[i] = newv;
	}
	for (auto finfo : chart_face_info)
	{
		mesh_temp1.add_face(v2newv[finfo.uv0], v2newv[finfo.uv1], v2newv[finfo.uv2]);
	}
	mesh_temp = mesh_temp1;
	//OpenMesh::IO::write_mesh(mesh_temp, ".\\simple_result\\middle_result_p.obj");

	std::vector<OpenMesh::VertexHandle> boundary_list;
	OpenMesh::HalfedgeHandle hedge_init;
	OpenMesh::HalfedgeHandle hedge_temp;
	OpenMesh::VertexHandle vertex_init;
	OpenMesh::VertexHandle vertex_temp;

	for (auto he_h : mesh_temp.halfedges())
	{
		if (mesh_temp.is_boundary(he_h))
		{
			hedge_init = he_h;
			break;
		}
	}
	vertex_init = mesh_temp.to_vertex_handle(hedge_init);

	hedge_temp = hedge_init;
	vertex_temp = vertex_init;
	do
	{
		boundary_list.push_back(vertex_temp);
		boundary_cc_list.push_back(vertex_temp.idx());
		boundary_cce_list.push_back(mesh_temp.edge_handle(hedge_temp).idx());
		vertex_temp = mesh_temp.from_vertex_handle(hedge_temp);
		hedge_temp = mesh_temp.prev_halfedge_handle(hedge_temp);
	} while (vertex_temp != vertex_init);

	std::vector<OpenMesh::VertexHandle> corner_list;
	for (int i = origin_corner_list.size() - 1; i >= 0; i--)
	{
		corner_list.push_back(v2newv[origin_corner_list[i]]);
		corner_vid_list.push_back(v2newv[origin_corner_list[i]].idx());
	}
	//std::cout << "corner corner_list.............." << corner_list.size() << std::endl;
	construct_segment();
	std::vector<double> align_direction1;
	align_direction1.resize(mesh_temp.n_vertices());
	for (int i = 0; i < segment_list.size(); i++)
	{
		OpenMesh::Vec2d proxy_temp = get_proxy(segment_list[i]);
		OpenMesh::Vec2d proxy_direction(-proxy_temp[1], proxy_temp[0]);
		double direction_temp = polar_angle(proxy_direction);
		for (int k = 0; k < segment_list[i].size(); k++)
		{
			OpenMesh::HalfedgeHandle he_temp = mesh_temp.halfedge_handle(mesh_temp.edge_handle(segment_list[i][k]), 0);
			he_temp = mesh_temp.is_boundary(he_temp) ? he_temp : mesh_temp.opposite_halfedge_handle(he_temp);
			OpenMesh::VertexHandle fromv = mesh_temp.from_vertex_handle(he_temp);
			align_direction1[fromv.idx()] = direction_temp;
		}
	}

	align_direction.resize(n_vertices);
	for (int i = 0; i < n_boundary_edges; i++)
	{
		int v0 = boundary_h_vert[i].first;
		align_direction[v0] = align_direction1[v2newv[v0].idx()];
	}
}

void ChartDeformation::get_align_direction()
{
	boundary_cce_list.clear();
	boundary_cc_list.clear();
	corner_vid_list.clear();
	segment_list.clear();
	Mesh mesh_temp1;
	double Fmax = M_PI / 9.0;
	std::vector<OpenMesh::VertexHandle> v2newv;
	v2newv.resize(n_vertices);
	for (int i = 0; i < n_vertices; i++)
	{
		OpenMesh::Vec3d point_temp(uv_x[2 * i], uv_x[2 * i + 1], 0);
		auto newv = mesh_temp1.add_vertex(point_temp);
		v2newv[i] = newv;
	}
	for (auto finfo : chart_face_info)
	{
		mesh_temp1.add_face(v2newv[finfo.uv0], v2newv[finfo.uv1], v2newv[finfo.uv2]);
	}
	mesh_temp = mesh_temp1;
	//OpenMesh::IO::write_mesh(mesh_temp1, ".\\simple_result\\middle_result_p2.obj");

	std::vector<OpenMesh::VertexHandle> boundary_list;
	OpenMesh::HalfedgeHandle hedge_init;
	OpenMesh::HalfedgeHandle hedge_temp;
	OpenMesh::VertexHandle vertex_init;
	OpenMesh::VertexHandle vertex_temp;

	for (auto he_h : mesh_temp.halfedges())
	{
		if (mesh_temp.is_boundary(he_h))
		{
			hedge_init = he_h;
			break;
		}
	}
	vertex_init = mesh_temp.to_vertex_handle(hedge_init);

	hedge_temp = hedge_init;
	vertex_temp = vertex_init;
	do
	{
		boundary_list.push_back(vertex_temp);
		boundary_cc_list.push_back(vertex_temp.idx());
		boundary_cce_list.push_back(mesh_temp.edge_handle(hedge_temp).idx());
		vertex_temp = mesh_temp.from_vertex_handle(hedge_temp);
		hedge_temp = mesh_temp.prev_halfedge_handle(hedge_temp);
	} while (vertex_temp != vertex_init);

	//std::vector<OpenMesh::VertexHandle> corner_list;
	//for (int i = 0; i < boundary_list.size(); i++)
	//{
	//	boundary_cc_list.push_back(boundary_list[i].idx());
	//	OpenMesh::Vec3d vec1 = mesh_temp.point(boundary_list[(i + 1) % boundary_list.size()]) - mesh_temp.point(boundary_list[i]);
	//	OpenMesh::Vec3d vec2 = mesh_temp.point(boundary_list[(i + boundary_list.size() - 1) % boundary_list.size()]) - mesh_temp.point(boundary_list[i]);
	//	double theta = std::acos((vec1 | vec2) / (vec1.norm()*vec2.norm()));
	//	if (theta < 7.0 / 9.0 * M_PI)
	//	{
	//		corner_list.push_back(boundary_list[i]);
	//		//corner_vid_list.push_back(boundary_list[i].idx());
	//	}
	//}

	std::vector<OpenMesh::VertexHandle> corner_list;
	for (int i = origin_corner_list.size() - 1; i >= 0; i--)
	{
		corner_list.push_back(v2newv[origin_corner_list[i]]);
		corner_vid_list.push_back(v2newv[origin_corner_list[i]].idx());
	}

	//std::cout << "step1==============" << std::endl;
	std::vector<double> align_direction1 = lloyd_iteration();
	std::vector<OpenMesh::VertexHandle> corner_list1;
	std::vector<int> corner_position1;
	double boundary_length = 0;
	for (int i = 0; i < boundary_list.size(); i++)
	{
		OpenMesh::Vec3d point1 = mesh_temp.point(boundary_list[i]);
		OpenMesh::Vec3d point2 = mesh_temp.point(boundary_list[(i+1)%boundary_list.size()]);
		boundary_length += (point1 - point2).norm();
		if (align_direction1[boundary_list[i].idx()] != align_direction1[boundary_list[(i + 1) % boundary_list.size()].idx()])
		{
			corner_list1.push_back(boundary_list[i]);
			corner_position1.push_back(i);
		}
	}
	
	align_direction.clear();
	align_direction.resize(n_vertices);
	for (int i = 0; i < n_boundary_edges; i++)
	{
		int v0 = boundary_h_vert[i].first;
		align_direction[v0] = align_direction1[v2newv[v0].idx()];
	}
	int count_temp = 0;
	for (int i = 0; i < n_boundary_edges; i++)
	{
		int v0 = boundary_h_vert[i].first;
		int v1 = boundary_h_vert[i].second;
		if (align_direction[v0] != align_direction[v1])
		{
			count_temp++;
		}
	}
	//std::cout << "align_direction size: " << count_temp <<std::endl;
}

void ChartDeformation::write_txt(const std::string &outstr)
{
	Mesh mesh_temp1;
	std::vector<OpenMesh::VertexHandle> v2newv;
	v2newv.resize(n_vertices);
	for (int i = 0; i < n_vertices; i++)
	{
		OpenMesh::Vec3d point_temp(uv_x[2 * i], uv_x[2 * i + 1], 0);
		auto newv = mesh_temp1.add_vertex(point_temp);
		v2newv[i] = newv;
	}
	for (auto finfo : chart_face_info)
	{
		mesh_temp1.add_face(v2newv[finfo.uv0], v2newv[finfo.uv1], v2newv[finfo.uv2]);
	}
	//OpenMesh::IO::write_mesh(mesh_temp1, ".\\simple_result\\middle_result_nn.obj");

	std::vector<OpenMesh::VertexHandle> boundary_list;
	OpenMesh::HalfedgeHandle hedge_init;
	OpenMesh::HalfedgeHandle hedge_temp;
	OpenMesh::VertexHandle vertex_init;
	OpenMesh::VertexHandle vertex_temp;

	for (auto he_h : mesh_temp1.halfedges())
	{
		if (mesh_temp1.is_boundary(he_h))
		{
			hedge_init = he_h;
			break;
		}
	}
	vertex_init = mesh_temp1.to_vertex_handle(hedge_init);

	hedge_temp = hedge_init;
	vertex_temp = vertex_init;
	do
	{
		boundary_list.push_back(vertex_temp);
		vertex_temp = mesh_temp1.from_vertex_handle(hedge_temp);
		hedge_temp = mesh_temp1.prev_halfedge_handle(hedge_temp);
	} while (vertex_temp != vertex_init);

	std::vector<OpenMesh::VertexHandle> corner_list;
	for (int i = 0; i < boundary_list.size(); i++)
	{
		OpenMesh::Vec3d vec1 = mesh_temp1.point(boundary_list[(i + 1) % boundary_list.size()]) - mesh_temp1.point(boundary_list[i]);
		OpenMesh::Vec3d vec2 = mesh_temp1.point(boundary_list[(i + boundary_list.size() - 1) % boundary_list.size()]) - mesh_temp1.point(boundary_list[i]);
		double x = (vec1 | vec2) / (vec1.norm()*vec2.norm());
		if (x >= 1)
		{
			x = 1;
		}
		if (x <= -1)
		{
			x = -1;
		}
		double theta = std::acos(x);
		if (theta <= 35.0 / 36.0 * M_PI)
		{
			corner_list.push_back(boundary_list[i]);
		}
	}

	std::vector<OpenMesh::VertexHandle> corner_list1;
	for (int i = 0; i < boundary_list.size(); i++)
	{
		if (align_direction[boundary_list[i].idx()] != align_direction[boundary_list[(i + 1) % boundary_list.size()].idx()])
		{
			corner_list1.push_back(boundary_list[i]);
		}
	}
	
	bool is_flip = 0;
	for (int i = 0; i < corner_list1.size(); i++)
	{
		OpenMesh::Vec3d point1 = mesh_temp1.point(corner_list1[i]);
		OpenMesh::Vec3d point2 = mesh_temp1.point(corner_list1[(i + 1) % corner_list1.size()]);
		for (int j = i + 2; j < i + corner_list1.size() - 2; j++)
		{
			OpenMesh::Vec3d point_temp1 = mesh_temp1.point(corner_list1[j%corner_list1.size()]);
			OpenMesh::Vec3d point_temp2 = mesh_temp1.point(corner_list1[(j + 1) % corner_list1.size()]);
			is_flip = detect_intersect(point1, point2, point_temp1, point_temp2);
			if (is_flip)
			{
				break;
			}
		}
		if (is_flip)
		{
			break;
		}
	}
	//std::cout << "is_flip: " << is_flip << std::endl;

	//std::cout << "corner_list size:" << corner_list.size() << "  corner_list1 size:" << corner_list1.size() << std::endl;

	std::vector<OpenMesh::VertexHandle> corner_list_temp;
	if (corner_list.size() <= corner_list1.size())
	{
		corner_list_temp = corner_list1;
	}
	else
	{
		corner_list_temp = corner_list;
	}
	std::ofstream of_txt(outstr, std::ios::trunc);
	of_txt << corner_list_temp.size() << std::endl;
	for (int i = 0; i < corner_list_temp.size(); i++)
	{
		of_txt << mesh_temp1.point(corner_list_temp[i])[0] << " " << mesh_temp1.point(corner_list_temp[i])[1] << std::endl;
	}
	of_txt.close();
}

void ChartDeformation::write_obj(const std::string & outstr)
{
	Mesh mesh_temp1;
	std::vector<OpenMesh::VertexHandle> v2newv;
	v2newv.resize(n_vertices);
	for (int i = 0; i < n_vertices; i++)
	{
		OpenMesh::Vec3d point_temp(uv_x[2 * i], uv_x[2 * i + 1], 0);
		auto newv = mesh_temp1.add_vertex(point_temp);
		v2newv[i] = newv;
	}
	for (auto finfo : chart_face_info)
	{
		mesh_temp1.add_face(v2newv[finfo.uv0], v2newv[finfo.uv1], v2newv[finfo.uv2]);
	}
	//OpenMesh::IO::write_mesh(mesh_temp1, ".\\simple_result\\middle_result_nn.obj");

	std::vector<std::vector<int>> face_vertex_list;
	std::vector<OpenMesh::Vec3d> vertex_point_list;
	face_vertex_list.resize(mesh_temp1.n_faces());
	vertex_point_list.resize(mesh_temp1.n_vertices());
	for (auto f_h : mesh_temp1.faces())
	{
		for (auto fv_h : mesh_temp1.fv_range(f_h))
		{
			face_vertex_list[f_h.idx()].push_back(fv_h.idx());
		}
	}
	for (auto v_h : mesh_temp1.vertices())
	{
		OpenMesh::Vec3d point_temp = mesh_temp1.point(v_h);
		vertex_point_list[v_h.idx()] = point_temp;
	}

	std::ofstream of_obj(outstr, std::ios::trunc);

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

bool ChartDeformation::detect_intersect(OpenMesh::Vec3d start_point1, OpenMesh::Vec3d end_point1, OpenMesh::Vec3d start_point2, OpenMesh::Vec3d end_point2)
{
	bool intersect_flag = 0;
	OpenMesh::Vec3d vec10 = end_point1 - start_point1;
	OpenMesh::Vec3d vec11 = start_point2 - start_point1;
	OpenMesh::Vec3d vec12 = end_point2 - start_point1;
	double n1 = (vec10%vec11)[2];
	double n2 = (vec10%vec12)[2];
	bool intersect_flag1 = n1 * n2 >= 0 ? 1 : 0;
	OpenMesh::Vec3d vec20 = end_point2 - start_point2;
	OpenMesh::Vec3d vec21 = start_point1 - start_point2;
	OpenMesh::Vec3d vec22 = end_point1 - start_point2;
	double n3 = (vec20%vec21)[2];
	double n4 = (vec20%vec22)[2];
	bool intersect_flag2 = n3 * n4 >= 0 ? 1 : 0;
	if (!intersect_flag1 && !intersect_flag2)
	{
		intersect_flag = 1;
	}
	return intersect_flag;
}

void ChartDeformation::construct_segment()
{
	for (int i = 0; i < corner_vid_list.size(); i++)
	{
		std::vector<int> segment_temp;
		int segment_start = corner_vid_list[i];
		int segment_end = corner_vid_list[(i + 1) % corner_vid_list.size()];
		std::vector<int>::iterator it_start = std::find(boundary_cc_list.begin(), boundary_cc_list.end(), segment_start);
		std::vector<int>::iterator it_end = std::find(boundary_cc_list.begin(), boundary_cc_list.end(), segment_end);
		int it_start_position = std::distance(boundary_cc_list.begin(), it_start);
		int it_end_position = std::distance(boundary_cc_list.begin(), it_end);
		if (it_end_position > it_start_position)
		{
			for (int j = it_start_position; j < it_end_position; j++)
			{
				segment_temp.push_back(boundary_cce_list[j]);
			}
			segment_list.push_back(segment_temp);
		}
		else
		{
			for (int j = it_start_position; j < it_end_position + boundary_cce_list.size(); j++)
			{
				segment_temp.push_back(boundary_cce_list[j%boundary_cce_list.size()]);
			}
			segment_list.push_back(segment_temp);
		}
	}
}

double ChartDeformation::polar_angle(OpenMesh::Vec2d vec2)
{
	OpenMesh::Vec2d x_axis(1, 0);
	vec2 = vec2 / vec2.norm();
	double theta = std::atan2(vec2[1], vec2[0]);
	return theta;
}

double ChartDeformation::get_error(OpenMesh::Vec2d e1, OpenMesh::Vec2d e2)
{
	double x = (e1 | e2) / (e1.norm()*e2.norm());
	if (x < -1)
	{
		x = -1;
	}
	if (x > 1)
	{
		x = 1;
	}
	return std::acos(x);
}

OpenMesh::Vec2d ChartDeformation::get_normal(int edge_)
{
	OpenMesh::EdgeHandle edge_temp = mesh_temp.edge_handle(edge_);
	OpenMesh::HalfedgeHandle he_temp = mesh_temp.halfedge_handle(edge_temp,0);
	if (mesh_temp.is_boundary(he_temp))
	{
		he_temp = mesh_temp.opposite_halfedge_handle(he_temp);
	}
	OpenMesh::VertexHandle from_vertex = mesh_temp.from_vertex_handle(he_temp);
	OpenMesh::VertexHandle to_vertex = mesh_temp.to_vertex_handle(he_temp);
	OpenMesh::Vec3d vec3 = mesh_temp.point(to_vertex) - mesh_temp.point(from_vertex);
	OpenMesh::Vec2d normal_temp;
	normal_temp[0] = -vec3[1];
	normal_temp[1] = vec3[0];
	normal_temp = normal_temp / normal_temp.norm();
	return normal_temp;
}

OpenMesh::Vec2d ChartDeformation::get_proxy(std::vector<int> edge_part)
{
	OpenMesh::Vec2d normal_temp(0,0);
	for (int i = 0; i < edge_part.size(); i++)
	{
		OpenMesh::Vec2d edge_normal = get_normal(edge_part[i]);
		double edge_length = mesh_temp.calc_edge_length(mesh_temp.edge_handle(edge_part[i]));
		normal_temp += edge_length * edge_normal;
	}
	normal_temp = normal_temp / normal_temp.norm();
	return normal_temp;
}

int ChartDeformation::get_seed(std::vector<int> edge_part, OpenMesh::Vec2d proxy_)
{
	std::vector <double> error_list;
	for (int i = 0; i < edge_part.size(); i++)
	{
		OpenMesh::Vec2d edge_normal = get_normal(edge_part[i]);
		double edge_error = get_error(edge_normal , proxy_);
		error_list.push_back(edge_error);
	}
	int min_position = std::min_element(error_list.begin(), error_list.end()) - error_list.begin();
	return edge_part[min_position];
}

void ChartDeformation::merge_hole(std::vector<std::vector<int>> &segment_part_, std::vector<int> hole_)
{
	OpenMesh::EdgeHandle hole_edge_start = mesh_temp.edge_handle(hole_[0]);
	OpenMesh::EdgeHandle hole_edge_end = mesh_temp.edge_handle(hole_[hole_.size() - 1]);
	std::vector<int> segment_id;
	std::vector<int> is_insert;
	for (int i = 0; i < segment_part_.size(); i++)
	{
		OpenMesh::EdgeHandle segment_part_start = mesh_temp.edge_handle(segment_part_[i][0]);
		OpenMesh::EdgeHandle segment_part_end = mesh_temp.edge_handle(segment_part_[i][segment_part_[i].size() - 1]);
		OpenMesh::HalfedgeHandle start_he = mesh_temp.halfedge_handle(hole_edge_start, 0);
		OpenMesh::HalfedgeHandle end_he = mesh_temp.halfedge_handle(hole_edge_end, 0);
		start_he = mesh_temp.is_boundary(start_he)? start_he : mesh_temp.opposite_halfedge_handle(start_he);
		end_he = mesh_temp.is_boundary(end_he) ? end_he : mesh_temp.opposite_halfedge_handle(end_he);
		if (mesh_temp.edge_handle(mesh_temp.next_halfedge_handle(start_he)) == segment_part_end)
		{
			segment_id.push_back(i);
			is_insert.push_back(0);
		}
		if (mesh_temp.edge_handle(mesh_temp.prev_halfedge_handle(end_he)) == segment_part_start)
		{
			segment_id.push_back(i);
			is_insert.push_back(1);
		}
	}
	if (segment_id.size() == 1)
	{
		if (is_insert[0] == 1)
		{
			segment_part_[segment_id[0]].insert(segment_part_[segment_id[0]].begin(), hole_.begin(), hole_.end());
		}
		else
		{
			for (int i = 0; i < hole_.size(); i++)
			{
				segment_part_[segment_id[0]].push_back(hole_[i]);
			}
		}
	}
	if (segment_id.size() == 2)
	{
		OpenMesh::Vec2d hole_proxy = get_proxy(hole_);
		OpenMesh::Vec2d segment_part_proxy1 = get_proxy(segment_part_[segment_id[0]]);
		OpenMesh::Vec2d segment_part_proxy2 = get_proxy(segment_part_[segment_id[1]]);
		int id_temp = get_error(segment_part_proxy1 , hole_proxy) > get_error(segment_part_proxy2 , hole_proxy) ? 1 : 0;
		if (is_insert[id_temp] == 1)
		{
			segment_part_[segment_id[id_temp]].insert(segment_part_[segment_id[id_temp]].begin(), hole_.begin(), hole_.end());
		}
		else
		{
			for (int i = 0; i < hole_.size(); i++)
			{
				segment_part_[segment_id[id_temp]].push_back(hole_[i]);
			}
		}
	}
}

void ChartDeformation::merge_segment_part(std::vector<std::vector<int>> &segment_part_, int merge_id_)
{
	std::vector<int> segment_part_temp = segment_part_[merge_id_ + 1];
	for (int i = 0; i < segment_part_temp.size(); i++)
	{
		segment_part_[merge_id_].push_back(segment_part_temp[i]);
	}
	segment_part_.erase(segment_part_.begin() + merge_id_ + 1);
}

std::vector<int> ChartDeformation::get_new_segmentpart(std::vector<int> hole_)
{
	double Fmax = M_PI / 9.0;
	int seed = hole_[hole_.size() / 2];
	//double proxy = polar_angle(get_normal(seed));
	OpenMesh::Vec2d proxy = get_normal(seed);
	std::vector<int> segment_part;
	segment_part.push_back(seed);
	std::vector<int> is_edge_checked;
	is_edge_checked.resize(mesh_temp.n_edges());
	for (int i = 0; i < hole_.size(); i++)
	{
		is_edge_checked[hole_[i]] = 1;
	}
	is_edge_checked[seed] = 0;
	int candidate_first, candidate_last;
	candidate_first = candidate_last = hole_.size() / 2;
	double segment_part_length = 0;
	double segment_length = 0;
	for (auto i : hole_)
	{
		segment_length += mesh_temp.calc_edge_length(mesh_temp.edge_handle(i));
	}

	bool is_still_growth = 1;
	while (is_still_growth)
	{
		is_still_growth = 0;
		int candidate_first0 = candidate_first - 1 >= 0 ? candidate_first - 1 : candidate_first;
		int candidate_last0 = candidate_last + 1 <= hole_.size() - 1 ? candidate_last + 1 : candidate_last;
		OpenMesh::Vec2d edge_proxy1 = get_normal(candidate_first0);
		OpenMesh::Vec2d edge_proxy2 = get_normal(candidate_last0);
		double edge_error1 = get_error(edge_proxy1 , proxy);
		double edge_error2 = get_error(edge_proxy2 , proxy);
		if (!is_edge_checked[hole_[candidate_first0]])
		{
			edge_error1 = 1e3;
		}
		if (!is_edge_checked[hole_[candidate_last0]])
		{
			edge_error2 = 1e3;
		}
		
		if (edge_error1 < edge_error2 && edge_error1 < Fmax)
		{
			segment_part.insert(segment_part.begin(), hole_[candidate_first0]);
			candidate_first = candidate_first0;
			is_edge_checked[hole_[candidate_first0]] = 0;
			is_still_growth = 1;
		}
		if (edge_error1 >= edge_error2 && edge_error2 < Fmax)
		{
			segment_part.push_back(hole_[candidate_last0]);
			candidate_last = candidate_last0;
			is_edge_checked[hole_[candidate_last0]] = 0;
			is_still_growth = 1;
		}
	}
	for (auto i : segment_part)
	{
		segment_part_length += mesh_temp.calc_edge_length(mesh_temp.edge_handle(i));
	}

	bool is_convergence = 0;
	int iter_num = 0;
	while (!is_convergence)
	{
		seed = get_seed(segment_part, get_proxy(segment_part));
		proxy = get_normal(seed);
		segment_part.clear();
		segment_part.push_back(seed);
		for (int i = 0; i < hole_.size(); i++)
		{
			is_edge_checked[hole_[i]] = 1;
		}
		is_edge_checked[seed] = 0;
		int seed_position = std::distance(hole_.begin(),std::find(hole_.begin(), hole_.end(), seed));
		candidate_first = candidate_last = seed_position;

		bool is_still_growth1 = 1;
		while (is_still_growth1)
		{
			is_still_growth1 = 0;
			int candidate_first0 = candidate_first - 1 >= 0 ? candidate_first - 1 : candidate_first;
			int candidate_last0 = candidate_last + 1 <= hole_.size() - 1 ? candidate_last + 1 : candidate_last;
			OpenMesh::Vec2d edge_proxy1 = get_normal(candidate_first0);
			OpenMesh::Vec2d edge_proxy2 = get_normal(candidate_last0);
			double edge_error1 = get_error(edge_proxy1 , proxy);
			double edge_error2 = get_error(edge_proxy2 , proxy);
			if (!is_edge_checked[hole_[candidate_first0]])
			{
				edge_error1 = 1e3;
			}
			if (!is_edge_checked[hole_[candidate_last0]])
			{
				edge_error2 = 1e3;
			}
			if (edge_error1 < edge_error2 && edge_error1 < Fmax)
			{
				segment_part.insert(segment_part.begin(), hole_[candidate_first0]);
				candidate_first = candidate_first0;
				is_edge_checked[hole_[candidate_first]] = 0;
				is_still_growth1 = 1;
			}
			if (edge_error1 >= edge_error2 && edge_error2 < Fmax)
			{
				segment_part.push_back(hole_[candidate_last0]);
				candidate_last = candidate_last0;
				is_edge_checked[hole_[candidate_last]] = 0;
				is_still_growth1 = 1;
			}
		}
		double segment_part_length1 = 0;
		for (auto i : segment_part)
		{
			segment_part_length1 += mesh_temp.calc_edge_length(mesh_temp.edge_handle(i));
		}
		if (std::fabs(segment_part_length1 - segment_part_length) / segment_part_length1 < 0.05)
		{
			is_convergence = 1;
		}
		iter_num++;
		if (iter_num > 9)
		{
			break;
		}
	}

	return segment_part;
}

std::vector<double> ChartDeformation::lloyd_iteration()
{
	std::vector<double> align_direction_;
	align_direction_.resize(mesh_temp.n_vertices());
	double Fmax = M_PI / 9.0;
	construct_segment();
	//std::cout << "step2==============" << std::endl;
	//std::cout << "segment_list size"<<segment_list.size() << std::endl;
	for (int i = 0; i < segment_list.size(); i++)
	{
		//初始所需变量
		int seed1 = segment_list[i][0];
		int seed2 = segment_list[i][segment_list[i].size() - 1];
		OpenMesh::Vec2d proxy1 = get_normal(seed1);
		OpenMesh::Vec2d proxy2 = get_normal(seed2);
		std::vector<std::vector<int>> segment_part;
		segment_part.resize(2);
		segment_part[0].push_back(seed1);
		segment_part[1].push_back(seed2);
		//标记遍历边
		std::vector<int> is_edge_checked;
		is_edge_checked.resize(mesh_temp.n_edges());
		for (int j = 0; j < segment_list[i].size(); j++)
		{
			is_edge_checked[segment_list[i][j]] = 1;
		}
		is_edge_checked[seed1] = 0;
		is_edge_checked[seed2] = 0;
		int candidate1_first, candidate1_last, candidate2_first, candidate2_last;
		candidate1_first = candidate1_last = 0;
		candidate2_first = candidate2_last = segment_list[i].size() - 1;
		std::vector<double> segment_part_length;
		segment_part_length.resize(segment_part.size());
		double segment_length = 0;
		for (auto j : segment_list[i])
		{
			segment_length += mesh_temp.calc_edge_length(mesh_temp.edge_handle(j));
		}

		bool is_still_growth = 1;
		while (is_still_growth)
		{
			is_still_growth = 0;
			int candidate1_last0 = candidate1_last + 1 <= segment_list[i].size() - 1 ? candidate1_last + 1 : candidate1_last ;
			int candidate2_first0 = candidate2_first - 1 >= 0 ? candidate2_first - 1 : candidate2_first;
			OpenMesh::Vec2d edge_proxy1 = get_normal(candidate1_last0);
			OpenMesh::Vec2d edge_proxy2 = get_normal(candidate2_first0);
			double edge_error1 = get_error(edge_proxy1 , proxy1);
			double edge_error2 = get_error(edge_proxy2 , proxy2);
			if (!is_edge_checked[segment_list[i][candidate1_last0]])
			{
				edge_error1 = 1e3;
			}
			if (!is_edge_checked[segment_list[i][candidate2_first0]])
			{
				edge_error2 = 1e3;
			}
			if (edge_error1 < edge_error2 && edge_error1 < Fmax)
			{
				segment_part[0].push_back(segment_list[i][candidate1_last0]);
				candidate1_last = candidate1_last0;
				is_edge_checked[segment_list[i][candidate1_last]] = 0;
				is_still_growth = 1;
			}
			if (edge_error1 >= edge_error2 && edge_error2 < Fmax)
			{
				segment_part[1].insert(segment_part[1].begin(),segment_list[i][candidate2_first0]);
				candidate2_first = candidate2_first0;
				is_edge_checked[segment_list[i][candidate2_first]] = 0;
				is_still_growth = 1;
			}
		}
		for (auto j : segment_part[0])
		{
			segment_part_length[0] += mesh_temp.calc_edge_length(mesh_temp.edge_handle(j));
		}
		for (auto j : segment_part[1])
		{
			segment_part_length[1] += mesh_temp.calc_edge_length(mesh_temp.edge_handle(j));
		}

		//std::cout << "step2.5===============" << std::endl;

		bool is_convergence = 0;
		int iter_num = 0;
		while (!is_convergence)
		{
			seed1 = get_seed(segment_part[0], get_proxy(segment_part[0]));
			seed2 = get_seed(segment_part[1], get_proxy(segment_part[1]));
			proxy1 = get_normal(seed1);
			proxy2 = get_normal(seed2);
			segment_part.clear();
			segment_part.resize(2);
			segment_part[0].push_back(seed1);
			segment_part[1].push_back(seed2);
			//标记遍历边
			//is_edge_checked.resize(boundary_cce_list.size());
			for (int j = 0; j < segment_list[i].size(); j++)
			{
				is_edge_checked[segment_list[i][j]] = 1;
			}
			is_edge_checked[seed1] = 0;
			is_edge_checked[seed2] = 0;
			int seed1_position = std::distance(segment_list[i].begin(),std::find(segment_list[i].begin(), segment_list[i].end(), seed1));
			int seed2_position = std::distance(segment_list[i].begin(),std::find(segment_list[i].begin(), segment_list[i].end(), seed2));
			candidate1_first = candidate1_last = seed1_position;
			candidate2_first = candidate2_last = seed2_position;

			bool is_still_growth1 = 1;
			while (is_still_growth1)
			{
				is_still_growth1 = 0;
				int candidate1_first0 = candidate1_first - 1 >= 0 ? candidate1_first - 1 : candidate1_first;
				int candidate1_last0 = candidate1_last + 1 <= segment_list[i].size() - 1 ? candidate1_last + 1 : candidate1_last;
				int candidate2_first0 = candidate2_first - 1 >= 0 ? candidate2_first - 1 : candidate2_first;
				int candidate2_last0 = candidate2_last + 1 <= segment_list[i].size() - 1 ? candidate2_last + 1 : candidate2_last;
				OpenMesh::Vec2d candidate1_first_proxy = get_normal(candidate1_first0);
				OpenMesh::Vec2d candidate1_last_proxy = get_normal(candidate1_last0);
				OpenMesh::Vec2d candidate2_first_proxy = get_normal(candidate2_first0);
				OpenMesh::Vec2d candidate2_last_proxy = get_normal(candidate2_last0);
				double edge_error1 = get_error(candidate1_first_proxy , proxy1);
				double edge_error2 = get_error(candidate1_last_proxy , proxy1);
				double edge_error3 = get_error(candidate2_first_proxy , proxy2);
				double edge_error4 = get_error(candidate2_last_proxy , proxy2);
				if (!is_edge_checked[segment_list[i][candidate1_first0]])
				{
					edge_error1 = 1e3;
				}
				if (!is_edge_checked[segment_list[i][candidate1_last0]])
				{
					edge_error2 = 1e3;
				}
				if (!is_edge_checked[segment_list[i][candidate2_first0]])
				{
					edge_error3 = 1e3;
				}
				if (!is_edge_checked[segment_list[i][candidate2_last0]])
				{
					edge_error4 = 1e3;
				}
				std::vector<double> error_list;
				error_list.push_back(edge_error1);
				error_list.push_back(edge_error2);
				error_list.push_back(edge_error3);
				error_list.push_back(edge_error4);
				int min_proxy_position = std::min_element(error_list.begin(), error_list.end()) - error_list.begin();
				if (error_list[min_proxy_position] < Fmax) /////////////////////here is a threshold
				{
					switch (min_proxy_position)
					{
					case 0:
						segment_part[0].insert(segment_part[0].begin(), segment_list[i][candidate1_first0]);
						candidate1_first = candidate1_first0;
						is_edge_checked[segment_list[i][candidate1_first]] = 0;
					case 1:
						segment_part[0].push_back(segment_list[i][candidate1_last0]);
						candidate1_last = candidate1_last0;
						is_edge_checked[segment_list[i][candidate1_last]] = 0;
					case 2:
						segment_part[1].insert(segment_part[1].begin(), segment_list[i][candidate2_first0]);
						candidate2_first = candidate2_first0;
						is_edge_checked[segment_list[i][candidate2_first]] = 0;
					case 3:
						segment_part[1].push_back(segment_list[i][candidate2_last0]);
						candidate2_last = candidate2_last0;
						is_edge_checked[segment_list[i][candidate2_last]] = 0;
					}
					is_still_growth1 = 1;
				}
			}
			double segment_part_length1 = 0;
			double segment_part_length2 = 0;
			for (auto j : segment_part[0])
			{
				segment_part_length1 += mesh_temp.calc_edge_length(mesh_temp.edge_handle(j));
			}
			for (auto j : segment_part[1])
			{
				segment_part_length2 += mesh_temp.calc_edge_length(mesh_temp.edge_handle(j));
			}
			if (std::fabs(segment_part_length1 - segment_part_length[0]) / segment_part_length1 < 0.05 && std::fabs(segment_part_length2 - segment_part_length[1]) / segment_part_length2 < 0.05)
			{
				is_convergence = 1;
			}
			iter_num++;
			if (iter_num > 9)
			{
				break;
			}
		}

		//std::cout << "step3==============" << std::endl;

		//filing holes
		while (std::accumulate(is_edge_checked.begin(), is_edge_checked.end(), 0) > 0)
		{
			std::vector<std::vector<int>> segment_holes;
			for (int j = 0; j < segment_list[i].size(); j++)
			{
				if (is_edge_checked[segment_list[i][j]] == 0)
				{
					continue;
				}
				else
				{
					std::vector<int> hole_temp;
					hole_temp.push_back(segment_list[i][j]);
					if (j + 1 == segment_list[i].size())
					{
						segment_holes.push_back(hole_temp);
						break;
					}
					while (is_edge_checked[segment_list[i][j + 1]] == 1)
					{
						j++;
						hole_temp.push_back(segment_list[i][j]);
						if (j + 1 == segment_list[i].size())
						{
							break;
						}
					}
					segment_holes.push_back(hole_temp);
				}
			}

			for (int j = 0; j < segment_holes.size(); j++)
			{
				double hole_length = 0;
				for (auto k : segment_holes[j])
				{
					hole_length += mesh_temp.calc_edge_length(mesh_temp.edge_handle(k));
				}
				if (hole_length/segment_length < 0.1)
				{
					merge_hole(segment_part, segment_holes[j]);
					for (int k = 0; k < segment_holes[j].size(); k++)
					{
						is_edge_checked[segment_holes[j][k]] = 0;
					}
				}
				else
				{
					std::vector<int> segment_part_new = get_new_segmentpart(segment_holes[j]);
					int position0 = std::distance(segment_list[i].begin(),std::find(segment_list[i].begin(), segment_list[i].end(), segment_part_new[0]));
					int insert_position = segment_part.size();
					for (int k = 0; k < segment_part.size(); k++)
					{
						int position1 = std::distance(segment_list[i].begin(),std::find(segment_list[i].begin(), segment_list[i].end(), segment_part[k][0]));
						if (position0 < position1)
						{
							insert_position = k;
							break;
						}
					}
					if (insert_position == segment_part.size())
					{
						segment_part.push_back(segment_part_new);
					}
					else
					{
						segment_part.insert(segment_part.begin() + insert_position, segment_part_new);
					}
					for (int k = 0; k < segment_part_new.size(); k++)
					{
						is_edge_checked[segment_part_new[k]] = 0;
					}
				}
			}
		}

		//std::cout << "step4==============" << std::endl;
		//modify segment_part
		bool is_exit_small = 1;
		while (is_exit_small)
		{
			is_exit_small = 0;
			for (int j = 0; j < segment_part.size(); j++)
			{
				double part_length = 0;
				for (auto k : segment_part[j])
				{
					part_length += mesh_temp.calc_edge_length(mesh_temp.edge_handle(k));
				}

				if (part_length / segment_length < 0.1)
				{
					merge_hole(segment_part, segment_part[j]);
					segment_part.erase(segment_part.begin() + j);
					is_exit_small = 1;
					break;
				}
			}
		}

		////////////////////gai
		bool is_can_merge = 1;
		while (is_can_merge)
		{
			if (segment_part.size() == 1)
			{
				break;
			}
			is_can_merge = 0;
			std::vector<double> error_list;
			for (int j = 0; j < segment_part.size() - 1; j++)
			{
				OpenMesh::Vec2d proxy_temp1 = get_proxy(segment_part[j]);
				OpenMesh::Vec2d proxy_temp2 = get_proxy(segment_part[j + 1]);
				error_list.push_back(get_error(proxy_temp1, proxy_temp2));
			}

			int merge_id;
			bool is_id_intersect = 1;
			while (is_id_intersect)
			{
				merge_id = std::min_element(error_list.begin(), error_list.end()) - error_list.begin();
				if (error_list[merge_id] == 1e3)
				{
					is_id_intersect = 0;
					break;
				}
				OpenMesh::EdgeHandle start_edge = mesh_temp.edge_handle(segment_part[merge_id][0]);
				OpenMesh::EdgeHandle end_edge = mesh_temp.edge_handle(segment_part[merge_id][segment_part[merge_id].size() - 1]);
				OpenMesh::HalfedgeHandle start_he = mesh_temp.halfedge_handle(start_edge, 0);
				OpenMesh::HalfedgeHandle end_he = mesh_temp.halfedge_handle(end_edge, 0);
				start_he = mesh_temp.is_boundary(start_he) ? start_he : mesh_temp.opposite_halfedge_handle(start_he);
				end_he = mesh_temp.is_boundary(end_he) ? end_he : mesh_temp.opposite_halfedge_handle(end_he);
				OpenMesh::Vec3d start_point = mesh_temp.point(mesh_temp.to_vertex_handle(start_he));
				OpenMesh::Vec3d end_point = mesh_temp.point(mesh_temp.from_vertex_handle(end_he));
				for (int k = i + 2; k < segment_list.size() + i-1; k++)
				{
					OpenMesh::EdgeHandle segmentedge0 = mesh_temp.edge_handle(segment_list[k%segment_list.size()][0]);
					OpenMesh::EdgeHandle segmentedge1 = mesh_temp.edge_handle(segment_list[k%segment_list.size()][segment_list[k%segment_list.size()].size() - 1]);
					OpenMesh::HalfedgeHandle segmenthe0 = mesh_temp.halfedge_handle(segmentedge0, 0);
					OpenMesh::HalfedgeHandle segmenthe1 = mesh_temp.halfedge_handle(segmentedge1, 0);
					segmenthe0 = mesh_temp.is_boundary(segmenthe0) ? segmenthe0 : mesh_temp.opposite_halfedge_handle(segmenthe0);
					segmenthe1 = mesh_temp.is_boundary(segmenthe1) ? segmenthe1 : mesh_temp.opposite_halfedge_handle(segmenthe1);
					OpenMesh::Vec3d point_start_temp = mesh_temp.point(mesh_temp.to_vertex_handle(segmenthe0));
					OpenMesh::Vec3d point_end_temp = mesh_temp.point(mesh_temp.from_vertex_handle(segmenthe1));
					is_id_intersect = detect_intersect(start_point, end_point, point_start_temp, point_end_temp);
					if (is_id_intersect)
					{
						break;
					}
				}
				if (is_id_intersect)
				{
					error_list[merge_id] = 1e3;
				}
			}
			if (error_list[merge_id] < 0.5*Fmax)
			{
				merge_segment_part(segment_part, merge_id);
				is_can_merge = 1;
			}
		}

		//std::cout << "step4.5==============" << std::endl;
		/*bool is_flip = 1;
		double merge_angle = Fmax;
		int iter_num1 = 0;
		while (is_flip)
		{
			iter_num1++;
			merge_angle = merge_angle / 2.0;
			is_flip = 0;
			std::vector<std::vector<int>> segment_part_temp = segment_part;
			bool is_can_merge = 1;
			while (is_can_merge)
			{
				if (segment_part_temp.size() == 1)
				{
					break;
				}
				is_can_merge = 0;
				std::vector<double> error_list;
				for (int j = 0; j < segment_part_temp.size() - 1; j++)
				{
					OpenMesh::Vec2d proxy_temp1 = get_proxy(segment_part_temp[j]);
					OpenMesh::Vec2d proxy_temp2 = get_proxy(segment_part_temp[j + 1]);
					error_list.push_back(get_error(proxy_temp1, proxy_temp2));
				}
				int merge_id = std::min_element(error_list.begin(), error_list.end()) - error_list.begin();
				if (error_list[merge_id] < merge_angle)
				{
					merge_segment_part(segment_part_temp, merge_id);
					is_can_merge = 1;
				}
			}

			for (int j = 0; j < segment_part_temp.size(); j++)
			{
				OpenMesh::EdgeHandle partedge0 = mesh_temp.edge_handle(segment_part_temp[j][0]);
				OpenMesh::EdgeHandle partedge1 = mesh_temp.edge_handle(segment_part_temp[j][segment_part_temp[j].size() - 1]);
				OpenMesh::HalfedgeHandle parthe0 = mesh_temp.halfedge_handle(partedge0, 0);
				OpenMesh::HalfedgeHandle parthe1 = mesh_temp.halfedge_handle(partedge1, 0);
				parthe0 = mesh_temp.is_boundary(parthe0) ? parthe0 : mesh_temp.opposite_halfedge_handle(parthe0);
				parthe1 = mesh_temp.is_boundary(parthe1) ? parthe1 : mesh_temp.opposite_halfedge_handle(parthe1);
				OpenMesh::Vec3d point_start = mesh_temp.point(mesh_temp.to_vertex_handle(parthe0));
				OpenMesh::Vec3d point_end = mesh_temp.point(mesh_temp.from_vertex_handle(parthe1));
				for (int k = i + 1; k < segment_list.size() + i; k++)
				{
					OpenMesh::EdgeHandle segmentedge0 = mesh_temp.edge_handle(segment_list[k%segment_list.size()][0]);
					OpenMesh::EdgeHandle segmentedge1 = mesh_temp.edge_handle(segment_list[k%segment_list.size()][segment_list[k%segment_list.size()].size() - 1]);
					OpenMesh::HalfedgeHandle segmenthe0 = mesh_temp.halfedge_handle(segmentedge0, 0);
					OpenMesh::HalfedgeHandle segmenthe1 = mesh_temp.halfedge_handle(segmentedge1, 0);
					segmenthe0 = mesh_temp.is_boundary(segmenthe0) ? segmenthe0 : mesh_temp.opposite_halfedge_handle(segmenthe0);
					segmenthe1 = mesh_temp.is_boundary(segmenthe1) ? segmenthe1 : mesh_temp.opposite_halfedge_handle(segmenthe1);
					OpenMesh::Vec3d point_start_temp = mesh_temp.point(mesh_temp.to_vertex_handle(segmenthe0));
					OpenMesh::Vec3d point_end_temp = mesh_temp.point(mesh_temp.from_vertex_handle(segmenthe1));
					is_flip = detect_intersect(point_start, point_end, point_start_temp, point_end_temp);
					if (is_flip)
					{
						break;
					}
				}
				if (is_flip)
				{
					std::cout << "???????????????" << std::endl;
					break;
				}
			}

			if (!is_flip)
			{
				segment_part = segment_part_temp;
			}
			if(iter_num1 > 5)
			{
				std::cout << "segment can not merge" << std::endl;
				break;
			}
		}*/

		//std::cout << "step5==============" << std::endl;
		for (int j = 0; j < segment_part.size(); j++)
		{
			OpenMesh::Vec2d proxy_temp = get_proxy(segment_part[j]);
			OpenMesh::Vec2d proxy_direction(-proxy_temp[1], proxy_temp[0]);
			double direction_temp = polar_angle(proxy_direction);
			for (int k = 0; k < segment_part[j].size(); k++)
			{
				OpenMesh::HalfedgeHandle he_temp = mesh_temp.halfedge_handle(mesh_temp.edge_handle(segment_part[j][k]), 0);
				he_temp = mesh_temp.is_boundary(he_temp) ? he_temp : mesh_temp.opposite_halfedge_handle(he_temp);
				OpenMesh::VertexHandle fromv = mesh_temp.from_vertex_handle(he_temp);
				align_direction_[fromv.idx()] = direction_temp;
			}
		}
		//std::cout << "step6==============" << std::endl;
	}
	return align_direction_;
}

void ChartDeformation::construct_line()
{
	std::vector<int> corner_list;
	for (int i = 0; i < n_boundary_edges; i++)
	{
		int v0 = boundary_h_vert[i].first;
		int v1 = boundary_h_vert[i].second;
		if (align_direction[v0] != align_direction[v1])
		{
			corner_list.push_back(v1);
		}
	}
	//std::cout << "??last corner_list size: " << corner_list.size() << std::endl;
	line_list.resize(corner_list.size());
	std::vector<int> boundary_c_vid;
	for (int i = 0; i < n_boundary_edges; i++)
	{
		int v0 = boundary_h_vert[i].first;
		boundary_c_vid.push_back(v0);
	}
	for (int i = 0; i < corner_list.size(); i++)
	{
		int v0 = corner_list[i];
		int v1 = corner_list[(i + 1) % corner_list.size()];
		int position0 = std::distance(boundary_c_vid.begin(), std::find(boundary_c_vid.begin(), boundary_c_vid.end(), v0));
		int position1 = std::distance(boundary_c_vid.begin(), std::find(boundary_c_vid.begin(), boundary_c_vid.end(), v1));
		if (position1 > position0)
		{
			for (int j = position0; j <= position1; j++)
			{
				line_list[i].push_back(boundary_c_vid[j]);
			}
		}
		else
		{
			for (int j = position0; j <= position1 + boundary_c_vid.size(); j++)
			{
				line_list[i].push_back(boundary_c_vid[j%boundary_c_vid.size()]);
			}
		}
	}


	/*line_list.resize(origin_corner_list.size());
	std::vector<int> boundary_c_vid;
	for (int i = 0; i < n_boundary_edges; i++)
	{
		int v0 = boundary_h_vert[i].first;
		boundary_c_vid.push_back(v0);
	}
	for (int i = 0; i < origin_corner_list.size(); i++)
	{
		int v0 = origin_corner_list[i];
		int v1 = origin_corner_list[(i + 1) % origin_corner_list.size()];
		int position0 = std::distance(boundary_c_vid.begin(), std::find(boundary_c_vid.begin(), boundary_c_vid.end(), v0));
		int position1 = std::distance(boundary_c_vid.begin(), std::find(boundary_c_vid.begin(), boundary_c_vid.end(), v1));
		if (position1 > position0)
		{
			for (int j = position0; j <= position1; j++)
			{
				line_list[i].push_back(boundary_c_vid[j]);
			}
		}
		else
		{
			for (int j = position0; j <= position1 + boundary_c_vid.size(); j++)
			{
				line_list[i].push_back(boundary_c_vid[j%boundary_c_vid.size()]);
			}
		}
	}*/
}