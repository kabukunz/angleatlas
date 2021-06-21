#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math.h>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "ScafData.h"
#include "PardisoSolver.h"
#include <time.h>

using namespace Eigen;
using namespace std;



double get_smallest_pos_quad_zero(double a, double b, double c);
class Parafun
{

public:
	Parafun(ScafData & data) :d_(data) {
		pardiso = NULL;
		convgence_con_rate = 1e-6;
		MAX_ITER_NUM = 100;
		Plus_area_term = false;
		Plus_area_term_search = false;

		Plus_shrink_term = false;
		shrink_weight=1.0;
		free_scaffold_weight=10;
		shrink_constant=1.21;
		area_term_weight = 1e-10;
		packing_efficiency = d_.PE_BOUND;
		inter_weight = d_.interZone_weight;

		if (!d_.is_gap)
			d_.gz.gf_num_start = d_.f_num;

		init_once();

	};
	~Parafun();

	void after_mesh_improve();

	void init_once();

	//void BPE();
	//void calc_gradient_norm(const VectorXd &x);

	void init();
	void init_area();

	void run_bpe();
	double perform_iteration_cm(bool use_CM);
	//double perform_iteration_cm(bool is_ip_convrate, bool is_slim_convrate);
	void id_vs_index();

	void Update_source_same_t();
	//void Update_source_same_t();

	void Pre_calculate();
	void Pre_calculate_once();

	void CM();
	void SLIM();
	//pp
	//void CM(bool is_interp = false);
	//void SLIM(bool is_interp = false);

	void Energy(const VectorXd &x, double &energy,double step=0.0,const Vector4d& quad=Vector4d(0,0,1,0));
	void Energysource(double step = 0.0);
	double compute_energy(const Eigen::MatrixXd & x, bool whole=false);
	void adjust_scaf_weight(double new_weight);
	void handle_mintri();

	double newton_equation(const double & a, const double & b, const double & K);

	void backtracking_line_search(const VectorXd &x, const VectorXd &d, const VectorXd &negetive_grad, double &alpha);

	void local_coordinate_inverse(int i, double &p00, double &p01, double &p10, double &p11);
	void local_coordinate_inverse_scaf(int i, double &p00, double &p01, double &p10, double &p11);

	void local_coordinate_inverse_gap(int i, double &p00, double &p01, double &p10, double &p11);

	void max_step(const VectorXd &xx, const VectorXd &dd, double &step);

	void descent_direction_modified(const VectorXd &pos, VectorXd &d);
	void gradient_modified(const VectorXd &pos, vector<double> &b);

	void out_dis(const string& file_str);

	ScafData &d_;

	double area_threshold;
	vector<int> var_ids;
	vector<int> id2index;

	double convgence_con_rate;
	double time_consumption;
	int MAX_ITER_NUM;

	struct Packinfo
	{
		double distortion;
		double pe;
		Eigen::MatrixXd w_uv;
		int v_num;
		Eigen::MatrixXi surface_F;
	};

	double inter_weight;
	double inter_weight_1;

	Packinfo best_cache;
	//total points number(mesh+scaffold)
	int total_num;

	//mesh faces number
	int mf_num;
	//mesh points number
	int mv_num;
	//total faces number
	int F_N;
	//variable points number(mesh+free scaffold point)
	int V_N;
	int dim = 2;
	std::vector<std::set<int>> VV_cache;

	vector<double> area;
	vector<double> area_scaf;
	vector<double> area_src;


	vector<double> source_p00;
	vector<double> source_p01;
	vector<double> source_p10;
	vector<double> source_p11;

	//pp
	/*double Intp_T_Min;
	double changetocm_flag;
	vector<double> update_p00;
	vector<double> update_p01;
	vector<double> update_p10;
	vector<double> update_p11;*/

	vector<int> F0;
	vector<int> F1;
	vector<int> F2;

	PardisoSolver* pardiso;
	vector<int> pardiso_ia;
	vector<int> pardiso_ja;
	vector<double> pardiso_a;
	vector<double> pardiso_b;

	double energy_uniform;
	double energy_area;

	double bound_distortion_K;
	VectorXd position_of_mesh;

	vector<int> id_h00; vector<int> id_h01; vector<int> id_h02; vector<int> id_h03; vector<int> id_h04; vector<int> id_h05;
	vector<int> id_h11; vector<int> id_h12; vector<int> id_h13; vector<int> id_h14; vector<int> id_h15;
	vector<int> id_h22; vector<int> id_h23; vector<int> id_h24; vector<int> id_h25;
	vector<int> id_h33; vector<int> id_h34; vector<int> id_h35;
	vector<int> id_h44; vector<int> id_h45;
	vector<int> id_h55;

	// variable and function about packing area constraint;
	bool Plus_area_term;
	double area_term_weight;
	double packing_efficiency;
	double area_bbox;
	double Cur_gap;
	double Cur_pe;
	Vector3d quad_equation;
	bool Plus_area_term_search;
	//my
	double my_packing_efficiency;
	bool is_to_bound;

	Eigen::Matrix2d bbox;
	vector<vector<int>> bnds;
	vector<vector<int>> id_bnd_h;

	void init_boundary();
	void calc_area_H();
	double calc_efficiency_area(const VectorXd &x);
	double calc_Curgap(const VectorXd &x);
	void calc_areaTerm_hessian_gradient();
	void areaCondition(const VectorXd &xx, const VectorXd &dd, Vector3d& quad_equ);
	double calc_cur_bbox(const VectorXd &x);

	//variable and function about shrink the scaffold area

	bool Plus_shrink_term;
	double shrink_weight;
	double free_scaffold_weight;
	double shrink_constant;

	set<int> inner_scaffold;

	void separate_scaffold();

};

