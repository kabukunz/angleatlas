
#ifndef SCAFFOLD_TEST_SCAFDATA_H
#define SCAFFOLD_TEST_SCAFDATA_H

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include "ScafCommon.h"
#include "TriangleInterface.h"

struct ScafData {
  // dimension for domain of parameterization/deformation
 public:
  ScafData();


  void add_new_chart(const Eigen::MatrixXd&, const Eigen::MatrixXi&,
	  const Eigen::MatrixXd&, const Eigen::MatrixXi&);

  void mesh_improve(bool,double bigger_factor=3.0);

  void update_scaffold();
  void mesh_improve_inner(bool in_packing = false, double bigger_factor=3.0);
  void mesh_improve_inner_init(bool in_packing = false, double bigger_factor = 3.0);

  void gen_bbox(const Eigen::MatrixXd& uv, double dx, int frame_points, Eigen::MatrixXd & boundv_cc, Eigen::MatrixXd & boundv_c);

  double scaffold_factor = 10;


// Output
  double energy; // objective value

// INTERNAL
  long mv_num, mf_num;
  long sv_num, sf_num;
  Eigen::MatrixXd m_V; // input initial mesh V
  Eigen::MatrixXi m_T; // input initial mesh F/T

  Eigen::MatrixXd w_uv; // whole domain uv: mesh + free vertices
  Eigen::MatrixXi s_T; // scaffold domain tets: scaffold tets
  Eigen::MatrixXi w_T;

  Eigen::VectorXd m_M; // mesh area or volume
  Eigen::VectorXd s_M; // scaffold area or volume
  Eigen::VectorXd w_M; // area/volume weights for whole
  double mesh_measure; // area or volume
  double avg_edge_length;
  long v_num;
  long f_num;
  double proximal_p = 1e-8; //unused


  Eigen::MatrixXd V, UV_V;
  Eigen::MatrixXi F, UV_F;
  //my
  Eigen::MatrixXd F_N;

  std::vector<Eigen::MatrixXd> separated_V_UV;
  std::vector<Eigen::VectorXi> UV_VI;


  double atlas_factor;
  double bigger_factor_init = 1.005;

  Eigen::VectorXi frame_ids;

 public: // public for ser
  // caching
  Eigen::VectorXi internal_bnd;
  Eigen::MatrixXd rect_frame_V;

  // multi-chart support
  std::vector<int> component_sizes;
  std::vector<int> bnd_sizes;

  //3D
 public:
  Eigen::MatrixXi surface_F;

  int dim; // dimension for ambient space. Same for mesh/scaf

  // flow arap
  int inner_scaf_tets = 0; 
          
//  std::string weight_file;
 // std::vector<std::vector<double>> weight_l_c;
//  std::vector<double> local_control;

  double area_bbox;
  struct InitBox
  {
	  Eigen::VectorXd uv_max;
	  Eigen::VectorXd uv_min;
	  Eigen::VectorXd uv_mid;
	  double dx;
  };
  InitBox initbox;

  Eigen::VectorXi inner_bound_ids;
  std::set<int> neibor_1;
  double gap_distance = 0.0;
  bool is_gap = false;
  bool is_toBound = false;
  double PE_BOUND=0.8;
  double interZone_weight = 0.1;

  struct GapZone
  {
	  Eigen::MatrixXi g_T;
	  int gf_num_start;
	  int gv_inter_num;
	  int gv_inter_start;
	  Eigen::MatrixXd g_V;
	  Eigen::MatrixXi g_F;
	  Eigen::VectorXd g_M;
  };
  GapZone gz;
};


  #endif //SCAFFOLD_TEST_SCAFDATA_H
