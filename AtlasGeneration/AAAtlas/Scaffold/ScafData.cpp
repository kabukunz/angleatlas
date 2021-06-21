
#include "ScafData.h"

using namespace std;
using namespace Eigen;

void ScafData::update_scaffold()
{
  mv_num = m_V.rows();
  mf_num = m_T.rows();

  v_num = w_uv.rows();
  sf_num = s_T.rows();

  sv_num = v_num - mv_num;
  f_num = sf_num + mf_num;

  s_M = Eigen::VectorXd::Constant(sf_num, scaffold_factor);
}


void ScafData::mesh_improve(bool in_packing = false, double bigger_factor) {
  if (dim == 2) {

    MatrixXd m_uv = w_uv.topRows(mv_num);
    MatrixXd V_bnd;
    V_bnd.resize(internal_bnd.size(), 2);
    for (int i = 0; i < internal_bnd.size(); i++) // redoing step 1.
    {
      V_bnd.row(i) = m_uv.row(internal_bnd(i));
    }

    if(rect_frame_V.size() == 0) {
      Matrix2d ob;// = rect_corners;
      {
//      double scaf_range = 3;
        Eigen::Array2d scaf_range(bigger_factor, bigger_factor);
        ob.row(0) = initbox.uv_mid.array() + scaf_range * ((initbox.uv_min - initbox.uv_mid).array());
        ob.row(1) = initbox.uv_mid.array() + scaf_range * ((initbox.uv_max - initbox.uv_mid).array());
      }
      Vector2d rect_len;
      rect_len << ob(1, 0) - ob(0, 0), ob(1, 1) - ob(0, 1);

	  area_bbox = rect_len(0)*rect_len(1);
      int frame_points = 5;

      if(in_packing) {
        frame_points = 20;
      }
      
      rect_frame_V.resize(4 * frame_points, 2);
      for (int i = 0; i < frame_points; i++) {
        // 0,0;0,1
        rect_frame_V.row(i) << ob(0, 0), ob(0, 1) + i * rect_len(1) / frame_points;
        // 0,0;1,1
        rect_frame_V.row(i + frame_points)
            << ob(0, 0) + i * rect_len(0) / frame_points, ob(1, 1);
        // 1,0;1,1
        rect_frame_V.row(i + 2 * frame_points) << ob(1, 0), ob(1, 1)
            - i * rect_len(1) /
                frame_points;
        // 1,0;0,1
        rect_frame_V.row(i + 3 * frame_points)
            << ob(1, 0) - i * rect_len(0) / frame_points, ob(0, 1);
        // 0,0;0,1
      }
      frame_ids = Eigen::VectorXi::LinSpaced(rect_frame_V.rows(), mv_num,
      mv_num + rect_frame_V.rows());
    }

    // Concatenate Vert and Edge
    MatrixXd V;
    MatrixXi E;
	{
		V.resize(V_bnd.rows() + rect_frame_V.rows(), V_bnd.cols());
		V << V_bnd, rect_frame_V;
	}
    E.resize(V.rows(), 2);
    for (int i = 0; i < E.rows(); i++)
      E.row(i) << i, i + 1;
    int acc_bs = 0;
    for(auto bs:bnd_sizes) {
      E(acc_bs + bs - 1,1 ) = acc_bs;
      acc_bs += bs;
    }
    E(V.rows() - 1, 1) = acc_bs;
    assert(acc_bs== internal_bnd.size());

    MatrixXd H = MatrixXd::Zero(component_sizes.size(), 2);
    {
      int hole_f = 0;
      int hole_i = 0;
      for (auto cs:component_sizes) {
        for (int i = 0; i < 3; i++)
          H.row(hole_i) += m_uv.row(m_T(hole_f, i)); // redoing step 2
        hole_f += cs;
        hole_i++;
      }
    }
    H /= 3.;

    MatrixXd uv2;
	triangulate(V, E, H, uv2, s_T);

    auto bnd_n = internal_bnd.size();

    for (auto i = 0; i < s_T.rows(); i++)
      for (auto j = 0; j < s_T.cols(); j++) {
        auto &x = s_T(i, j);
        if (x < bnd_n) x = internal_bnd(x);
        else x += m_uv.rows() - bnd_n;
      }

	{
		surface_F.resize(m_T.rows() + s_T.rows(), 3);
		surface_F << m_T, s_T;
	}
	w_uv.conservativeResize(m_uv.rows() - bnd_n + uv2.rows(), 2);
    w_uv.bottomRows(uv2.rows() - bnd_n) = uv2.bottomRows(-bnd_n + uv2.rows());
  }
  neibor_1.clear();
  update_scaffold();
}


void ScafData::add_new_chart(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& UV_V, const Eigen::MatrixXi& UV_F)
{
	using namespace std;
	using namespace Eigen;

	VectorXd M;
	doublearea(V, F, M);

	Eigen::MatrixXd uv_init;

	std::vector<std::vector<int>> all_bnds;
	boundary_loop(F, all_bnds);
	int num_holes = all_bnds.size() - 1;

	uv_init = UV_V;

	mesh_measure += M.sum() / 2;
	//std::cout << "Mesh Measure" << M.sum() / 2 << std::endl;

	component_sizes.push_back(F.rows());

	if (mv_num == 0)
	{
		w_uv = uv_init;
	}
	else
	{
		MatrixXd m_uv = w_uv.topRows(mv_num);
		{
			w_uv.resize(m_uv.rows() + uv_init.rows(), 2);
			w_uv << m_uv, uv_init;
		}
	}


	m_M.conservativeResize(mf_num + M.size());
	m_M.bottomRows(M.size()) = M / 2;

	for (auto cur_bnd : all_bnds)
	{
		internal_bnd.conservativeResize(internal_bnd.size() + cur_bnd.size());
		internal_bnd.bottomRows(cur_bnd.size()) = Map<ArrayXi>(cur_bnd.data(), cur_bnd.size()) + mv_num;
		bnd_sizes.push_back(cur_bnd.size());
	}

	m_T.conservativeResize(mf_num + F.rows(), 3);
	m_T.bottomRows(F.rows()) = F.array() + mv_num;
	mf_num += F.rows();

	m_V.conservativeResize(mv_num + V.rows(), 3);
	m_V.bottomRows(V.rows()) = V;
	mv_num += V.rows();

	//rect_frame_V = MatrixXd();
	//mesh_improve(true);
}

ScafData::ScafData() {
  dim = 2;
  mv_num = 0;
  mf_num = 0;
  sf_num= 0;
  sv_num = 0;
  mesh_measure = 0;
}

void ScafData::gen_bbox(const Eigen::MatrixXd & uv, double dx, int frame_points, Eigen::MatrixXd & boundv_cc, Eigen::MatrixXd & boundv_c)
{
	Eigen::VectorXd uv_max = uv.colwise().maxCoeff();
	Eigen::VectorXd uv_min = uv.colwise().minCoeff();
	//Eigen::VectorXd uv_mid = (uv_max + uv_min) / 2.;
	Eigen::Vector2d delta;
	delta << dx, dx;
	Matrix2d ob;

	ob.row(0) = uv_min - delta;
	ob.row(1) = uv_max + delta;
	Vector2d rect_len;
	rect_len << ob(1, 0) - ob(0, 0), ob(1, 1) - ob(0, 1);

	boundv_c.resize(4 * frame_points, 2);
	for (int i = 0; i < frame_points; i++) {
		// 0,0;0,1
		boundv_c.row(i) << ob(0, 0), ob(0, 1) + i * rect_len(1) / frame_points;
		// 0,0;1,1
		boundv_c.row(i + frame_points)
			<< ob(0, 0) + i * rect_len(0) / frame_points, ob(1, 1);
		// 1,0;1,1
		boundv_c.row(i + 2 * frame_points) << ob(1, 0), ob(1, 1)
			- i * rect_len(1) /
			frame_points;
		// 1,0;0,1
		boundv_c.row(i + 3 * frame_points)
			<< ob(1, 0) - i * rect_len(0) / frame_points, ob(0, 1);
		// 0,0;0,1
	}
	boundv_cc.resize(4 * frame_points, 2);
	for (size_t i = 0; i < boundv_cc.rows(); i++)
	{
		boundv_cc.row(i) = boundv_c.row(4 * frame_points - i - 1);
	}


}


void ScafData::mesh_improve_inner(bool in_packing, double bigger_factor)
{
	MatrixXd m_uv = w_uv.topRows(mv_num);
	MatrixXd V_bnd;
	V_bnd.resize(internal_bnd.size(), 2);
	for (int i = 0; i < internal_bnd.size(); i++) // redoing step 1.
	{
		V_bnd.row(i) = m_uv.row(internal_bnd(i));
	}
	if (rect_frame_V.size() == 0) {
		Matrix2d ob;// = rect_corners;
		{
			Eigen::Array2d scaf_range(bigger_factor, bigger_factor);
			ob.row(0) = initbox.uv_mid.array() + scaf_range * ((initbox.uv_min - initbox.uv_mid).array());
			ob.row(1) = initbox.uv_mid.array() + scaf_range * ((initbox.uv_max - initbox.uv_mid).array());
		}
		Vector2d rect_len;
		rect_len << ob(1, 0) - ob(0, 0), ob(1, 1) - ob(0, 1);

		area_bbox = rect_len(0)*rect_len(1);
		int frame_points = 5;

		if (in_packing) {
			frame_points = 20;
		}

		rect_frame_V.resize(4 * frame_points, 2);
		for (int i = 0; i < frame_points; i++) {
			// 0,0;0,1
			rect_frame_V.row(i) << ob(0, 0), ob(0, 1) + i * rect_len(1) / frame_points;
			// 0,0;1,1
			rect_frame_V.row(i + frame_points)
				<< ob(0, 0) + i * rect_len(0) / frame_points, ob(1, 1);
			// 1,0;1,1
			rect_frame_V.row(i + 2 * frame_points) << ob(1, 0), ob(1, 1)
				- i * rect_len(1) /
				frame_points;
			// 1,0;0,1
			rect_frame_V.row(i + 3 * frame_points)
				<< ob(1, 0) - i * rect_len(0) / frame_points, ob(0, 1);
			// 0,0;0,1
		}
		frame_ids = Eigen::VectorXi::LinSpaced(rect_frame_V.rows(), mv_num,
			mv_num + rect_frame_V.rows());
	}

	vector< Eigen::MatrixXd> boundv_cc;
	int components_size = separated_V_UV.size();
	boundv_cc.resize(components_size);

	int inner_bound_pnum = 20;
	for (size_t i = 0; i < components_size; i++)
	{
		boundv_cc[i].resize(4 * inner_bound_pnum, 2);
		for (size_t j = 0; j < 4 * inner_bound_pnum; j++)
		{
			boundv_cc[i].row(j) = w_uv.row(inner_bound_ids(i * 4 * inner_bound_pnum + j));
		}
	}


	Eigen::MatrixXd inner_frame_V;
	inner_frame_V.resize(components_size * 4 * inner_bound_pnum, 2);
	for (size_t i = 0; i < components_size; i++)
	{
		inner_frame_V.block(i * 4 * inner_bound_pnum, 0, 4 * inner_bound_pnum, 2) = boundv_cc[i];
	}

	MatrixXd H = MatrixXd::Zero(component_sizes.size(), 2);
	{
		int hole_f = 0;
		int hole_i = 0;
		for (auto cs : component_sizes) {
			for (int i = 0; i < 3; i++)
				H.row(hole_i) += m_uv.row(m_T(hole_f, i)); // redoing step 2
			hole_f += cs;
			hole_i++;
		}
	}
	H /= 3.;


	Eigen::MatrixXi s_T1;
	MatrixXd uv21;

	{
		Eigen::MatrixXd V = rect_frame_V;
		V.conservativeResize(rect_frame_V.rows() + components_size * 4 * inner_bound_pnum, 2);
		V.bottomRows(components_size * 4 * inner_bound_pnum) = inner_frame_V;

		//std::vector<int> bnd_sizes_copy = bnd_sizes;	
		std::vector<int> bnd_sizes_copy;
		bnd_sizes_copy.push_back(rect_frame_V.rows());
		for (size_t i = 0; i < components_size; i++)
		{
			bnd_sizes_copy.push_back(4 * inner_bound_pnum);
		}

		MatrixXi E;
		E.resize(V.rows(), 2);
		for (int i = 0; i < E.rows(); i++)
			E.row(i) << i, i + 1;
		int acc_bs = 0;
		for (auto bs : bnd_sizes_copy) {
			E(acc_bs + bs - 1, 1) = acc_bs;
			acc_bs += bs;
		}
		triangulate(V, E, H, uv21, s_T1);
	}

	for (auto i = 0; i < s_T1.rows(); i++)
		for (auto j = 0; j < s_T1.cols(); j++) {
			auto &x = s_T1(i, j);
			x += m_uv.rows();
		}


	int vn_1 = V_bnd.rows() + inner_frame_V.rows();
	int fn_1 = m_T.rows() + s_T1.rows();
	neibor_1.clear();

	Eigen::MatrixXi s_T2 = gz.g_T;
	MatrixXd uv22= w_uv.bottomRows(gz.gv_inter_num);

	gz.gf_num_start = m_T.rows() + s_T1.rows();

	int gv_inter_start_new = m_uv.rows() + uv21.rows();

	for (auto i = 0; i < s_T2.rows(); i++)
	{
		neibor_1.insert(fn_1 + i);
		for (auto j = 0; j < s_T2.cols(); j++) {
			auto &x = s_T2(i, j);
			if (x >=gz.gv_inter_start)
			{
				x += gv_inter_start_new-gz.gv_inter_start;
			}
		}
	}

	{
		s_T.resize(s_T1.rows() + s_T2.rows(), 3);
		s_T << s_T1, s_T2;
		surface_F.resize(m_T.rows() + s_T.rows(), 3);
		surface_F << m_T, s_T;
	}

	w_uv.conservativeResize(m_uv.rows() + uv21.rows() + gz.gv_inter_num, 2);
	w_uv.block(m_uv.rows(), 0, uv21.rows(), 2) = uv21;
	w_uv.bottomRows(gz.gv_inter_num) = uv22;

	update_scaffold();
}

void ScafData::mesh_improve_inner_init(bool in_packing, double bigger_factor)
{
	MatrixXd m_uv = w_uv.topRows(mv_num);
	MatrixXd V_bnd;
	V_bnd.resize(internal_bnd.size(), 2);
	for (int i = 0; i < internal_bnd.size(); i++) // redoing step 1.
	{
		V_bnd.row(i) = m_uv.row(internal_bnd(i));
	}

	if (rect_frame_V.size() == 0) {
		Matrix2d ob;// = rect_corners;
		{
			Eigen::Array2d scaf_range(bigger_factor, bigger_factor);
			ob.row(0) = initbox.uv_mid.array() + scaf_range * ((initbox.uv_min - initbox.uv_mid).array());
			ob.row(1) = initbox.uv_mid.array() + scaf_range * ((initbox.uv_max - initbox.uv_mid).array());
		}
		Vector2d rect_len;
		rect_len << ob(1, 0) - ob(0, 0), ob(1, 1) - ob(0, 1);

		area_bbox = rect_len(0)*rect_len(1);
		int frame_points = 5;

		if (in_packing) {
			frame_points = 20;
		}

		rect_frame_V.resize(4 * frame_points, 2);
		for (int i = 0; i < frame_points; i++) {
			// 0,0;0,1
			rect_frame_V.row(i) << ob(0, 0), ob(0, 1) + i * rect_len(1) / frame_points;
			// 0,0;1,1
			rect_frame_V.row(i + frame_points)
				<< ob(0, 0) + i * rect_len(0) / frame_points, ob(1, 1);
			// 1,0;1,1
			rect_frame_V.row(i + 2 * frame_points) << ob(1, 0), ob(1, 1)
				- i * rect_len(1) /
				frame_points;
			// 1,0;0,1
			rect_frame_V.row(i + 3 * frame_points)
				<< ob(1, 0) - i * rect_len(0) / frame_points, ob(0, 1);
			// 0,0;0,1
		}
		frame_ids = Eigen::VectorXi::LinSpaced(rect_frame_V.rows(), mv_num,
			mv_num + rect_frame_V.rows());
	}


	vector< Eigen::MatrixXd> boundv_cc, boundv_c;
	int components_size = separated_V_UV.size();
	boundv_c.resize(components_size);
	boundv_cc.resize(components_size);

	int inner_bound_pnum = 20;
	for (size_t i = 0; i < components_size; i++)
	{
		gen_bbox(separated_V_UV[i], initbox.dx, inner_bound_pnum, boundv_cc[i], boundv_c[i]);
	}

	Eigen::MatrixXd inner_frame_V;
	inner_frame_V.resize(components_size * 4 * inner_bound_pnum, 2);
	for (size_t i = 0; i < components_size; i++)
	{
		inner_frame_V.block(i * 4 * inner_bound_pnum, 0, 4 * inner_bound_pnum, 2) = boundv_cc[i];
	}


	{
		int inner_b_id_start = m_uv.rows() + rect_frame_V.rows();
		inner_bound_ids.resize(components_size * 4 * inner_bound_pnum);
		for (int i = 0; i < components_size * 4 * inner_bound_pnum; i++)
		{
			inner_bound_ids(i) = inner_b_id_start + i;
		}
	}

	MatrixXd H = MatrixXd::Zero(component_sizes.size(), 2);
	{
		int hole_f = 0;
		int hole_i = 0;
		for (auto cs : component_sizes) {
			for (int i = 0; i < 3; i++)
				H.row(hole_i) += m_uv.row(m_T(hole_f, i)); // redoing step 2
			hole_f += cs;
			hole_i++;
		}
	}
	H /= 3.;


	Eigen::MatrixXi s_T1, s_T2;
	MatrixXd uv21, uv22;

	{

		Eigen::MatrixXd V = rect_frame_V;
		V.conservativeResize(rect_frame_V.rows() + components_size * 4 * inner_bound_pnum, 2);
		V.bottomRows(components_size * 4 * inner_bound_pnum) = inner_frame_V;

		//std::vector<int> bnd_sizes_copy = bnd_sizes;	
		std::vector<int> bnd_sizes_copy;
		bnd_sizes_copy.push_back(rect_frame_V.rows());
		for (size_t i = 0; i < components_size; i++)
		{
			bnd_sizes_copy.push_back(4 * inner_bound_pnum);
		}

		MatrixXi E;
		E.resize(V.rows(), 2);
		for (int i = 0; i < E.rows(); i++)
			E.row(i) << i, i + 1;
		int acc_bs = 0;
		for (auto bs : bnd_sizes_copy) {
			E(acc_bs + bs - 1, 1) = acc_bs;
			acc_bs += bs;
		}
		triangulate(V, E, H, uv21, s_T1);
	}


	{
		Eigen::MatrixXd V;
		MatrixXi E;

		{
			V.resize(V_bnd.rows() + inner_frame_V.rows(), 2);
			V << V_bnd, inner_frame_V;
		}

		std::vector<int> bnd_sizes_copy = bnd_sizes;
		for (size_t i = 0; i < components_size; i++)
		{
			bnd_sizes_copy.push_back(4 * inner_bound_pnum);
		}
		E.resize(V.rows(), 2);
		for (int i = 0; i < E.rows(); i++)
			E.row(i) << i, i + 1;
		int acc_bs = 0;
		for (auto bs : bnd_sizes_copy) {
			E(acc_bs + bs - 1, 1) = acc_bs;
			acc_bs += bs;
		}

		triangulate(V, E, H, uv22, s_T2);
	}


	for (auto i = 0; i < s_T1.rows(); i++)
		for (auto j = 0; j < s_T1.cols(); j++) {
			auto &x = s_T1(i, j);
			x += m_uv.rows();
		}


	int vn_1 = V_bnd.rows() + inner_frame_V.rows();
	int fn_1 = m_T.rows() + s_T1.rows();
	neibor_1.clear();


	gz.g_F = s_T2;
	gz.g_V.resize(uv22.rows(), 3);
	gz.g_V.setZero();
	gz.g_V.block(0, 0, uv22.rows(), 2) = uv22;

	for (auto i = 0; i < s_T2.rows(); i++)
	{
		neibor_1.insert(fn_1 + i);
		for (auto j = 0; j < s_T2.cols(); j++) {
			auto &x = s_T2(i, j);
			if (x < V_bnd.rows())
			{
				x = internal_bnd(x);
			}
			else if (x < vn_1)
			{
				x += m_uv.rows() + rect_frame_V.rows() - V_bnd.rows();
			}
			else
			{
				x += m_uv.rows() + uv21.rows() - vn_1;
			}

		}
	}

	{
		s_T.resize(s_T1.rows() + s_T2.rows(), 3);
		s_T << s_T1, s_T2;
		surface_F.resize(m_T.rows() + s_T.rows(), 3);
		surface_F << m_T, s_T;
	}
	gz.gv_inter_start = m_uv.rows() + uv21.rows();
	gz.g_T = s_T2;
	gz.gf_num_start = m_T.rows()+s_T1.rows();
	gz.gv_inter_num = uv22.rows() - vn_1;

	gz.g_M.resize(gz.g_F.rows());
	double gM_sum = 0.0;
	for (size_t i = 0; i < gz.g_F.rows(); i++)
	{
		int f0 = gz.g_F(i, 0);
		int f1 = gz.g_F(i, 1);
		int f2 = gz.g_F(i, 2);

		Vector2d x_(gz.g_V(f1, 0) - gz.g_V(f0, 0), gz.g_V(f1, 1) - gz.g_V(f0, 1));
		Vector2d l_(gz.g_V(f2, 0) - gz.g_V(f0, 0), gz.g_V(f2, 1) - gz.g_V(f0, 1));

		double area_tri = 0.5*abs(x_(0)*l_(1) - x_(1)*l_(0));
		gz.g_M(i) = area_tri;
		gM_sum += area_tri;
	}
	gz.g_M /= gM_sum;

	w_uv.conservativeResize(m_uv.rows() + uv21.rows() + uv22.rows() - vn_1, 2);
	w_uv.block(m_uv.rows(), 0, uv21.rows(), 2) = uv21;
	w_uv.bottomRows(uv22.rows() - vn_1) = uv22.bottomRows(uv22.rows() - vn_1);

	update_scaffold();
}
