#pragma once

#include <set>

#include "MeshViewer/MeshDefinition.h"
#include "viewer2d_segments.h"

#include <memory>

#define MY_DEBUG

class ParaQuadChartDecomposition
{
public:
	ParaQuadChartDecomposition(const Mesh& IGPP_mesh);
	~ParaQuadChartDecomposition() = default;

	struct segment_type
	{
		int vert0, vert1, coord, tag;
	};

	void decomposition();

	const Mesh& get_quad_mesh() const { return quad_mesh; };
	const std::vector<Mesh>& get_atlas() const { return atlas; }
	const Mesh& get_packing_mesh() const { return packing_mesh; };

	const std::vector<int>& get_boundary_h() const { return boundary_h; }
	const std::vector<int>& get_boundary_v() const { return boundary_v; }

	const std::map<int, int>& get_boundary_h_index() const { return boundary_h_index; }
	const std::map<int, int>& get_boundary_v_index() const { return boundary_v_index; }

	const std::vector<int>& get_bh_segment() const { return bh_segment; }
	const std::vector<int>& get_vert_corner() const { return vert_corner; }

	const std::vector<std::pair<int, int>>& get_boundary_cutting_nodes() const { return boundary_cutting_nodes; }
	const std::vector<segment_type>& get_segments() const { return segments; }

	const std::vector<std::vector<int>>& get_cut_path_h() const { return cut_path_h; }
	const std::vector<std::vector<int>>& get_cut_path_seg() const { return cut_path_seg; }
	const std::vector<std::vector<std::pair<int, int>>>& get_cut_path_v() const { return cut_path_v; }

	const std::vector<std::vector<int>>& get_seg_nodes() const { return seg_nodes; }
	const std::vector<OpenMesh::Vec2i>& get_quad_vertices() const { return quad_vertices; }

	const std::vector<bool>& get_chart_flipped() const { return chart_flipped; }
	const std::vector<OpenMesh::Vec2d>& get_chart_translation() const { return chart_translation; }

	const std::vector<std::vector<int>>& get_v_chart2quad() const { return v_chart2quad; }

	const std::map<int, int>& get_cutting_edge_seg() const { return cutting_edge_seg; }

	bool range_test(int seg_id, int val);

private:
	const Mesh& mesh;
	Mesh quad_mesh, packing_mesh;

	std::vector<Mesh> atlas;
	std::vector<std::vector<int>> v_chart2quad;

	std::map<int, int> vert_k;
	
	std::vector<int> boundary_h;
	std::vector<int> boundary_v;

	std::map<int, int> boundary_h_index;
	std::map<int, int> boundary_v_index;

	std::vector<int> bh_segment, seg_chart;
	std::vector<int> accu_bh_chart, accu_bseg_chart;

	inline int chart_boundary_cycle(int cur, int delta);
	inline int chart_segment_cycle(int cur, int delta);
	
	std::unique_ptr<Viewer2D_Segments> viewer;
	void show_quad_decomposition();

	std::vector<int> vert_corner;
	int n_bsegs;
	std::vector<segment_type> segments;
	std::vector<std::set<int>> triangle_segs;
	std::map<int, std::vector<int>> triangle_corner;
	void segment_detective();

	std::set<std::pair<int, int>> seg_intersections;
	void quad_composition();

	void tracing_from_vert(int vert_id);

	std::map<int, int> v_in_quad;
	std::vector<OpenMesh::Vec2i> quad_vertices;
	std::vector<std::vector<int>> seg_nodes;

	std::vector<std::vector<int>> cornor_k_seg;

	struct quad_neighbor_type 
	{
		quad_neighbor_type() : next_v(4, -1), next_seg(4, -1){};

		std::vector<int> next_v, next_seg;
	};
	std::vector<quad_neighbor_type> quad_v_neighbor;

	struct retangle_type
	{
		retangle_type() :min(std::numeric_limits<int>::max(), std::numeric_limits<int>::max()), max(std::numeric_limits<int>::min(), std::numeric_limits<int>::min()) {};
		retangle_type(const OpenMesh::Vec2i& _min, const OpenMesh::Vec2i& _max) :min(_min), max(_max) {};

		OpenMesh::Vec2i min, max;

		inline int width() const { return max[0] - min[0]; };
		inline int height() const { return max[1] - min[1]; };
		inline int area() const { return width() * height(); };
	};

	std::vector<retangle_type> quads;
	std::vector<std::pair<int, int>> overlap_quads;

	std::vector<retangle_type> rects_candidate;
	std::vector<bool> quad_edge_inner_label, qeil_candidate;
	std::vector<std::vector<int>> chart_faces, chart_faces_candidate;

	std::vector<bool> chart_flipped, cf_candidate;
	std::vector<OpenMesh::Vec2d> chart_translation, ct_candidate;

	std::map<int, int> cutting_edge_seg;
	std::vector<int> cutting_vertex_valence;
	std::vector<std::pair<int, int>> boundary_cutting_nodes;

	double quad_mesh_perimeter;

	void motorcycle();
	void motorcycle_canididate();
	
	bool retangle_decomposition(double& new_cut_length);
	double retangle_packing(bool& isgood);
	void chart_decomposition();

	void cut_preprocess();
	std::vector<std::vector<int>> cut_path_h, cut_path_seg;
	std::vector<std::vector<std::pair<int, int>>> cut_path_v;
	void extract_cut();

	std::set<int> sample_boundary_vertices();

	inline std::pair<int, int> get_seg_range(int seg_id);
};

