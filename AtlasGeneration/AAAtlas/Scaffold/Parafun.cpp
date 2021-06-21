#include "Parafun.h"

void Parafun::after_mesh_improve()
{
	id_vs_index();
	total_num = d_.v_num;
	F_N = d_.f_num;
	V_N = d_.v_num-d_.frame_ids.size();

	position_of_mesh.resize(2 * total_num);
	for (size_t i = 0; i < dim; i++)
	{
		position_of_mesh.block(i*total_num, 0, total_num, 1) = d_.w_uv.col(i);
	}
	if (pardiso)
	{
		delete pardiso;
		pardiso = NULL;
	}
	if (!d_.is_gap)
		d_.gz.gf_num_start = d_.f_num;
	init();
}

void Parafun::init_once()
{
	mf_num=d_.mf_num;
	mv_num=d_.mv_num;
	total_num = d_.v_num;
	position_of_mesh.resize(2 * total_num);
	for (size_t i = 0; i < dim; i++)
	{
		position_of_mesh.block(i*total_num, 0, total_num, 1) = d_.w_uv.col(i);
	}
	F_N = d_.f_num;
	V_N = d_.v_num - d_.frame_ids.size();
	F0.resize(F_N);
	F1.resize(F_N);
	F2.resize(F_N);
	for (int i = 0; i < F_N; ++i)
	{
		F0[i] = d_.surface_F(i, 0);
		F1[i] = d_.surface_F(i, 1);
		F2[i] = d_.surface_F(i, 2);
	}
	init_boundary();
	//area 
	{
		area.resize(F_N);
		area_scaf.resize(F_N - mf_num);
		area_src.resize(mf_num);
		for (int i = 0; i < mf_num; ++i)
		{
			area_src[i] = d_.m_M(i);// *d_.local_control[i];
			area[i] = d_.m_M(i);//*d_.local_control[i];
		}
		for (int i = mf_num; i < F_N; ++i)
		{
			if (d_.neibor_1.count(i) > 0)
			{
				area[i] = d_.gz.g_M(i - d_.gz.gf_num_start);
				area_scaf[i - mf_num] = area[i];
			}
			else
			{

				area_scaf[i - mf_num] = d_.s_M(i - mf_num);
				area[i] = d_.s_M(i - mf_num);
			}
		}

	}

	//id2index,var_ids
	{
		auto& fr_ids = d_.frame_ids;
		var_ids.clear();
		{
			var_ids.reserve(total_num);
			int assign = 0, i = 0;
			for (int get = 0; i < total_num && get < fr_ids.size(); i++)
			{
				if (fr_ids(get) == i)
					get++;
				else
					var_ids.push_back(i);
			}
			while (i < total_num)
			{
				var_ids.push_back(i);
				i++;
			}
		}

		id2index.resize(total_num, -1);
		for (size_t i = 0; i < var_ids.size(); i++)
		{
			id2index[var_ids[i]] = i;
		}
	}
	handle_mintri();
	Pre_calculate_once();

	area_bbox = calc_cur_bbox(position_of_mesh);
	Cur_pe = calc_efficiency_area(position_of_mesh) / area_bbox;

	{
		energy_area = compute_energy(d_.w_uv, false);
		best_cache.distortion = energy_area / d_.mesh_measure - 4.0;
		best_cache.pe = Cur_pe;
		best_cache.w_uv = d_.w_uv;
		best_cache.v_num = total_num;
		best_cache.surface_F = d_.surface_F;
	}

// 	printf("Setting info==================================================\n");
// 	printf("area_bbox: %f; init_distortion: %f; init_packing_efficiency: %f; PE_BOUND=%.4f\n", area_bbox, best_cache.distortion, Cur_pe, packing_efficiency);
// 	printf("Plus_area_term=%d; area_term_weight=%.20f", Plus_area_term, area_term_weight);
// 	printf("Plus_area_term_search=%d\n", Plus_area_term_search);
// 	printf("Plus_shrink_term=%d; shrink_weight=%f; free_scaffold_weight=%f\n", Plus_shrink_term, shrink_weight, free_scaffold_weight);
// 	printf("is_interZone=%d; interZone_weight=%.2f\n", d_.is_gap, inter_weight);
// 	printf("Setting info==================================================\n");

	std::cout << "Scaf Distortion " << best_cache.distortion * 0.25 + 1.0 << std::endl;
	std::cout << "PE Input " << Cur_pe << std::endl;
	std::cout << "PE Bound " << packing_efficiency << std::endl;
}

Parafun::~Parafun()
{
}


void Parafun::init()
{

	F0.resize(F_N);
	F1.resize(F_N);
	F2.resize(F_N);
	for (int i = mf_num; i < F_N; ++i)
	{
		F0[i] = d_.surface_F(i, 0);
		F1[i] = d_.surface_F(i, 1);
		F2[i] = d_.surface_F(i, 2);
	}

	handle_mintri();
	init_area();
	Pre_calculate();
}

void Parafun::init_boundary()
{
	bnds.resize(d_.bnd_sizes.size());
	int start = 0;
	int cur_b_size;
	for (size_t i = 0; i < d_.bnd_sizes.size(); i++)
	{
		cur_b_size = d_.bnd_sizes[i];
		bnds[i].resize(cur_b_size);
		for (size_t j = 0; j < cur_b_size; j++)
		{
			bnds[i][j] = d_.internal_bnd(start + j);
		}
		start += cur_b_size;
	}
}

void Parafun::init_area()
{
	area.resize(F_N);
	area_scaf.resize(F_N - mf_num);


	for (int i = mf_num; i < F_N; ++i)
	{
		if (d_.neibor_1.count(i) > 0)
		{
			area[i] = d_.gz.g_M(i - d_.gz.gf_num_start);
			area_scaf[i - mf_num] = area[i];
		}
		else
		{
			area_scaf[i - mf_num] = d_.s_M(i - mf_num);
			area[i] = d_.s_M(i - mf_num);
		}
	}

	double area_inter_factor = inter_weight_1 * inter_weight;
	for (auto&var: d_.neibor_1)
	{
		area[var] *= area_inter_factor;
		area_scaf[var - mf_num] = area[var];
	}

}

void Parafun::run_bpe()
{
	init();
	//BPE();
}
 
double Parafun::perform_iteration_cm(bool use_CM)
{
	if (pardiso != NULL)
	{
		delete pardiso;
		pardiso = NULL;
	}
	pardiso = new PardisoSolver();
	pardiso->ia = pardiso_ia;
	pardiso->ja = pardiso_ja;
	pardiso->a.resize(pardiso_ja.size());
	pardiso->nnz = pardiso_ja.size();
	pardiso->num = 2 * V_N;
	pardiso->pardiso_init();
	//pp
	//Update_source_same_t();

	long time_beg, time_end;
	time_beg = clock();

	//pp
	/*std::cout << "Intp_T_Min: " << Intp_T_Min << std::endl;
	
	bool is_interp = is_ip_convrate && (Intp_T_Min < 0.999);
	bool is_slim = is_interp&& is_slim_convrate && (changetocm_flag < 0.99);
	std::cout << "is_slim: " << is_slim << std::endl;*/
	if(use_CM)
		CM();
	else
		SLIM();

	//pp
	/*if (!is_interp)
	{
		for (int i = 0; i < d_.m_T.rows(); ++i)
		{
			update_p00[i] = source_p00[i];
			update_p01[i] = source_p01[i];
			update_p10[i] = source_p10[i];
			update_p11[i] = source_p11[i];
		}
	}
	if (is_slim)
	{
		SLIM(is_interp);
	}
	else
	{
		CM(is_interp);
	}*/

	d_.w_uv = Map<Matrix<double, -1, -1, Eigen::ColMajor>>(position_of_mesh.data(), total_num, dim);

	time_end = clock();
	time_consumption = (time_end - time_beg) / 1000.0;

	area_bbox = calc_cur_bbox(position_of_mesh);
	Cur_pe = calc_efficiency_area(position_of_mesh) / area_bbox;

	if (Cur_pe >= packing_efficiency)
	{
		if (best_cache.distortion > energy_area / d_.mesh_measure-4.0)
		{
			best_cache.distortion = energy_area / d_.mesh_measure-4.0;
			best_cache.pe = Cur_pe;
			best_cache.w_uv = d_.w_uv;
			best_cache.v_num = total_num;
			best_cache.surface_F = d_.surface_F;
			//printf("Record Cur_pe: %f; Cur_dis: %f\n", Cur_pe, best_cache.distortion);
		}
		is_to_bound = 1;
	}
	if (!is_to_bound && (Cur_pe >= my_packing_efficiency))
	{
		if (best_cache.distortion > energy_area / d_.mesh_measure - 4.0)
		{
			best_cache.distortion = energy_area / d_.mesh_measure - 4.0;
			best_cache.pe = Cur_pe;
			best_cache.w_uv = d_.w_uv;
			best_cache.v_num = total_num;
			best_cache.surface_F = d_.surface_F;
			//printf("Record Cur_pe: %f; Cur_dis: %f\n", Cur_pe, best_cache.distortion);
		}
	}
	printf("PE: %f, ED: %f\n", Cur_pe, energy_area / d_.mesh_measure * 0.25);

	delete pardiso;
	pardiso = NULL;
	return energy_area/d_.mesh_measure;
}

void Parafun::id_vs_index()
{
	var_ids.resize(d_.v_num-d_.frame_ids.size());
	id2index.resize(d_.v_num);
	if (d_.v_num > total_num)
	{
		for (size_t i = 0; i < d_.v_num- total_num; i++)
		{
			id2index[i+total_num] = V_N + i;
			var_ids[i + V_N] = total_num + i;
		}
	}

}


void Parafun::Pre_calculate()
{
	source_p00.resize(F_N);
	source_p01.resize(F_N);
	source_p10.resize(F_N);
	source_p11.resize(F_N);

	for (int i = mf_num; i < d_.gz.gf_num_start; ++i)
	{
		double p00, p01, p10, p11;
		local_coordinate_inverse_scaf(i, p00, p01, p10, p11);

		source_p00[i] = p00;
		source_p01[i] = p01;
		source_p10[i] = p10;
		source_p11[i] = p11;
	}
	for (int i = d_.gz.gf_num_start; i < F_N; ++i)
	{
		double p00, p01, p10, p11;
		local_coordinate_inverse_gap(i - d_.gz.gf_num_start, p00, p01, p10, p11);
		source_p00[i] = p00;
		source_p01[i] = p01;
		source_p10[i] = p10;
		source_p11[i] = p11;
	}

	pardiso_ia.clear(); pardiso_ia.reserve(2 * V_N + 1);
	pardiso_ja.clear(); pardiso_ja.reserve(8 * V_N);

	typedef Triplet<int> T;
	std::vector<T> tripletlist;

	std::vector<std::set<int>> VV_tmp;
	VV_tmp.resize(V_N);
	for (size_t i = 0; i < mv_num; i++)
	{
		VV_tmp[i] = VV_cache[i];
	}

	if (Plus_area_term)
	{
		//for (size_t i = 0; i < bnds.size(); i++)
		//{
		//	for (size_t j = 0; j < bnds[i].size(); j++)
		//	{
		//		bnds[i][j] = id2index[bnds[i][j]];
		//	}
		//}

		for (size_t i = 0; i < bnds.size(); i++)
		{
			for (size_t j = 0; j < bnds[i].size(); j++)
			{
				for (auto&var : bnds[i])
				{
					//VV_tmp[id2index[bnds[i][j]]].insert(id2index[var]);
					VV_tmp[bnds[i][j]].insert(var);
				}
			}
		}

	}

	vector<int> s_vid;
	for (size_t i = 0; i < d_.s_T.rows(); i++)
	{
		s_vid.clear();
		for (size_t j = 0; j < d_.s_T.cols(); j++)
		{
			int s_id = id2index[d_.s_T(i, j)];
			if (s_id !=-1)
				s_vid.push_back(s_id);
		}
		if (s_vid.size() <= 1)
			continue;
		if (s_vid.size() == 2)
		{
			VV_tmp[s_vid[0]].insert(s_vid[1]);
			VV_tmp[s_vid[1]].insert(s_vid[0]);
		}
		else
		{
			VV_tmp[s_vid[0]].insert(s_vid[1]);
			VV_tmp[s_vid[0]].insert(s_vid[2]);

			VV_tmp[s_vid[1]].insert(s_vid[0]);
			VV_tmp[s_vid[1]].insert(s_vid[2]);

			VV_tmp[s_vid[2]].insert(s_vid[0]);
			VV_tmp[s_vid[2]].insert(s_vid[1]);
		}
	}

	for (int i = 0; i < V_N; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		VV_tmp[i].insert(i);
		vector<int> row_id;
		for (auto&var : VV_tmp[i])
		{
			row_id.push_back(var);
		}

		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k]);
			tripletlist.push_back(T(i, row_id[k], dd));
			++dd;
		}
		for (int k = 0; k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k]+V_N);
			tripletlist.push_back(T(i, row_id[k]+V_N, dd));
			++dd;
		}
	}
	for (int i = V_N; i < 2*V_N; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		vector<int> row_id;
		for (auto&var : VV_tmp[i-V_N])
		{
			row_id.push_back(var);
		}
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i-V_N);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k] + V_N);
			tripletlist.push_back(T(i, row_id[k] + V_N, dd));
			++dd;
		}
	}

	SparseMatrix<int> find_id_in_rows;
	find_id_in_rows.resize(2 * V_N, 2 * V_N);
	find_id_in_rows.setFromTriplets(tripletlist.begin(), tripletlist.end());

	pardiso_ia.push_back(pardiso_ja.size());

	id_h00.resize(F_N, -1); id_h01.resize(F_N, -1); id_h02.resize(F_N, -1); id_h03.resize(F_N, -1); id_h04.resize(F_N, -1); id_h05.resize(F_N, -1);
	id_h11.resize(F_N, -1); id_h12.resize(F_N, -1); id_h13.resize(F_N, -1); id_h14.resize(F_N, -1); id_h15.resize(F_N, -1);
	id_h22.resize(F_N, -1); id_h23.resize(F_N, -1); id_h24.resize(F_N, -1); id_h25.resize(F_N, -1);
	id_h33.resize(F_N, -1); id_h34.resize(F_N, -1); id_h35.resize(F_N, -1);
	id_h44.resize(F_N, -1); id_h45.resize(F_N, -1);
	id_h55.resize(F_N, -1);
	
	for (int i = 0; i < d_.m_T.rows(); i++)
	{
		int f0 = id2index[F0[i]]; int f1 = id2index[F1[i]]; int f2 = id2index[F2[i]]; int f3 = f0 + V_N; int f4 = f1 + V_N; int f5 = f2 + V_N;

		int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
		int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
		int min12 = min(f1, f2); int max12 = f1 + f2 - min12;

		id_h00[i] = pardiso_ia[f0]; id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01); id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
		id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3); id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4); id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);

		id_h11[i] = pardiso_ia[f1]; id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
		id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3); id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4); id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);

		id_h22[i] = pardiso_ia[f2];
		id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3); id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4); id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);

		id_h33[i] = pardiso_ia[f3]; id_h34[i] = pardiso_ia[min01 + V_N] + find_id_in_rows.coeff(min01 + V_N, max01 + V_N); id_h35[i] = pardiso_ia[min02 + V_N] + find_id_in_rows.coeff(min02 + V_N, max02 + V_N);

		id_h44[i] = pardiso_ia[f4]; id_h45[i] = pardiso_ia[min12 + V_N] + find_id_in_rows.coeff(min12 + V_N, max12 + V_N);

		id_h55[i] = pardiso_ia[f5];

	}

	for (int i = d_.m_T.rows(); i < F_N; i++)
	{
		int f0 = id2index[F0[i]]; int f1 = id2index[F1[i]]; int f2 = id2index[F2[i]];
		int f3 = f0 + V_N; int f4 = f1 + V_N; int f5 = f2 + V_N;
		if (f0 != -1)
		{
			id_h00[i] = pardiso_ia[f0];
			id_h33[i] = pardiso_ia[f3];
			id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3);
		}
		if (f1 != -1)
		{
			id_h11[i] = pardiso_ia[f1];
			id_h44[i] = pardiso_ia[f4];
			id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4);
		}
		if (f2 != -1)
		{
			id_h22[i] = pardiso_ia[f2];
			id_h55[i] = pardiso_ia[f5];
			id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);
		}

		if (f1 != -1 && f2 != -1)
		{
			int min12 = min(f1, f2); int max12 = f1 + f2 - min12;
			id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
			id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);
			id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4);
			id_h45[i] = pardiso_ia[min12 + V_N] + find_id_in_rows.coeff(min12 + V_N, max12 + V_N);
		}
		if (f0 != -1 && f2 != -1)
		{
			int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
			id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
			id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);
			id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3);
			id_h35[i] = pardiso_ia[min02 + V_N] + find_id_in_rows.coeff(min02 + V_N, max02 + V_N);
		}
		if (f1 != -1 && f0 != -1)
		{
			int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
			id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01);
			id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4);
			id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3);
			id_h34[i] = pardiso_ia[min01 + V_N] + find_id_in_rows.coeff(min01 + V_N, max01 + V_N);
		}

	}


	if (Plus_area_term)
	{
		id_bnd_h.resize(bnds.size());
		for (size_t i = 0; i < bnds.size(); i++)
		{
			int bnds_i_size = bnds[i].size();
			id_bnd_h[i].resize(bnds_i_size*(2 * bnds_i_size + 1));
			int index = 0;
			for (size_t j = 0; j < bnds_i_size; j++)
			{
				id_bnd_h[i][index] = pardiso_ia[bnds[i][j]];
				index++;
				for (size_t k = j + 1; k < bnds_i_size; k++)
				{
					int v1 = bnds[i][j];
					int v2 = bnds[i][k];
					if (v1 > v2) swap(v1, v2);
					id_bnd_h[i][index] = pardiso_ia[v1] + find_id_in_rows.coeff(v1, v2);
					index++;
				}
				for (size_t k = 0; k < bnds_i_size; k++)
				{
					int v1 = bnds[i][j];
					int v2 = bnds[i][k] + V_N;
					id_bnd_h[i][index] = pardiso_ia[v1] + find_id_in_rows.coeff(v1, v2);
					index++;
				}
			}

			for (size_t j = 0; j < bnds_i_size; j++)
			{
				id_bnd_h[i][index] = pardiso_ia[bnds[i][j]];
				index++;
				for (size_t k = j + 1; k < bnds_i_size; k++)
				{
					int v1 = bnds[i][j];
					int v2 = bnds[i][k];
					if (v1 > v2) swap(v1, v2);
					id_bnd_h[i][index] = pardiso_ia[v1 + V_N] + find_id_in_rows.coeff(v1 + V_N, v2 + V_N);
					index++;
				}
			}
		}

	}

}

void Parafun::Pre_calculate_once()
{
	source_p00.resize(F_N);
	source_p01.resize(F_N);
	source_p10.resize(F_N);
	source_p11.resize(F_N);

	for (int i = 0; i < d_.m_T.rows(); ++i)
	{
		double p00, p01, p10, p11;
		local_coordinate_inverse(i, p00, p01, p10, p11);

		source_p00[i] = p00;
		source_p01[i] = p01;
		source_p10[i] = p10;
		source_p11[i] = p11;
	}
	for (int i = mf_num; i < d_.gz.gf_num_start; ++i)
	{
		double p00, p01, p10, p11;
		local_coordinate_inverse_scaf(i, p00, p01, p10, p11);

		source_p00[i] = p00;
		source_p01[i] = p01;
		source_p10[i] = p10;
		source_p11[i] = p11;
	}
	for (int i = d_.gz.gf_num_start; i < F_N; ++i)
	{
		double p00, p01, p10, p11;
		local_coordinate_inverse_gap(i- d_.gz.gf_num_start, p00, p01, p10, p11);
		source_p00[i] = p00;
		source_p01[i] = p01;
		source_p10[i] = p10;
		source_p11[i] = p11;
	}

	pardiso_ia.clear(); pardiso_ia.reserve(2 * V_N + 1);
	pardiso_ja.clear(); pardiso_ja.reserve(8 * V_N);

	typedef Triplet<int> T;
	std::vector<T> tripletlist;

	std::vector<std::set<int>> VV_tmp;
	VV_tmp.resize(V_N);
	VV_cache.resize(mv_num);
	for (size_t i = 0; i < d_.m_T.rows(); i++)
	{
		int vid[3];

		for (size_t j = 0; j < d_.m_T.cols(); j++)
		{
			vid[j] = id2index[d_.m_T(i, j)];
		}
		VV_tmp[vid[0]].insert(vid[1]);
		VV_tmp[vid[0]].insert(vid[2]);

		VV_tmp[vid[1]].insert(vid[0]);
		VV_tmp[vid[1]].insert(vid[2]);

		VV_tmp[vid[2]].insert(vid[0]);
		VV_tmp[vid[2]].insert(vid[1]);

		VV_cache[vid[0]].insert(vid[1]);
		VV_cache[vid[0]].insert(vid[2]);

		VV_cache[vid[1]].insert(vid[0]);
		VV_cache[vid[1]].insert(vid[2]);

		VV_cache[vid[2]].insert(vid[0]);
		VV_cache[vid[2]].insert(vid[1]);
	}

	if (Plus_area_term)
	{
		for (size_t i = 0; i < bnds.size(); i++)
		{
			for (size_t j = 0; j < bnds[i].size(); j++)
			{
				bnds[i][j] = id2index[bnds[i][j]];
			}
		}

		for (size_t i = 0; i < bnds.size(); i++)
		{
			for (size_t j = 0; j < bnds[i].size(); j++)
			{
				for (auto&var : bnds[i])
				{
					//VV_tmp[id2index[bnds[i][j]]].insert(id2index[var]);
					VV_tmp[bnds[i][j]].insert(var);
				}
			}
		}

	}

	vector<int> s_vid;
	for (size_t i = 0; i < d_.s_T.rows(); i++)
	{
		s_vid.clear();
		for (size_t j = 0; j < d_.s_T.cols(); j++)
		{
			int s_id = id2index[d_.s_T(i, j)];
			if (s_id != -1)
				s_vid.push_back(s_id);
		}
		if (s_vid.size() <= 1)
			continue;
		if (s_vid.size() == 2)
		{
			VV_tmp[s_vid[0]].insert(s_vid[1]);
			VV_tmp[s_vid[1]].insert(s_vid[0]);
		}
		else
		{
			VV_tmp[s_vid[0]].insert(s_vid[1]);
			VV_tmp[s_vid[0]].insert(s_vid[2]);

			VV_tmp[s_vid[1]].insert(s_vid[0]);
			VV_tmp[s_vid[1]].insert(s_vid[2]);

			VV_tmp[s_vid[2]].insert(s_vid[0]);
			VV_tmp[s_vid[2]].insert(s_vid[1]);
		}
	}

	for (int i = 0; i < V_N; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		VV_tmp[i].insert(i);
		vector<int> row_id;
		for (auto&var : VV_tmp[i])
		{
			row_id.push_back(var);
		}

		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k]);
			tripletlist.push_back(T(i, row_id[k], dd));
			++dd;
		}
		for (int k = 0; k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k] + V_N);
			tripletlist.push_back(T(i, row_id[k] + V_N, dd));
			++dd;
		}
	}
	for (int i = V_N; i < 2 * V_N; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		vector<int> row_id;
		for (auto&var : VV_tmp[i - V_N])
		{
			row_id.push_back(var);
		}
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i - V_N);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k] + V_N);
			tripletlist.push_back(T(i, row_id[k] + V_N, dd));
			++dd;
		}
	}

	SparseMatrix<int> find_id_in_rows;
	find_id_in_rows.resize(2 * V_N, 2 * V_N);
	find_id_in_rows.setFromTriplets(tripletlist.begin(), tripletlist.end());

	pardiso_ia.push_back(pardiso_ja.size());

	id_h00.resize(F_N, -1); id_h01.resize(F_N, -1); id_h02.resize(F_N, -1); id_h03.resize(F_N, -1); id_h04.resize(F_N, -1); id_h05.resize(F_N, -1);
	id_h11.resize(F_N, -1); id_h12.resize(F_N, -1); id_h13.resize(F_N, -1); id_h14.resize(F_N, -1); id_h15.resize(F_N, -1);
	id_h22.resize(F_N, -1); id_h23.resize(F_N, -1); id_h24.resize(F_N, -1); id_h25.resize(F_N, -1);
	id_h33.resize(F_N, -1); id_h34.resize(F_N, -1); id_h35.resize(F_N, -1);
	id_h44.resize(F_N, -1); id_h45.resize(F_N, -1);
	id_h55.resize(F_N, -1);

	for (int i = 0; i < d_.m_T.rows(); i++)
	{
		int f0 = id2index[F0[i]]; int f1 = id2index[F1[i]]; int f2 = id2index[F2[i]]; int f3 = f0 + V_N; int f4 = f1 + V_N; int f5 = f2 + V_N;

		int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
		int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
		int min12 = min(f1, f2); int max12 = f1 + f2 - min12;

		id_h00[i] = pardiso_ia[f0]; id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01); id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
		id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3); id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4); id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);

		id_h11[i] = pardiso_ia[f1]; id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
		id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3); id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4); id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);

		id_h22[i] = pardiso_ia[f2];
		id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3); id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4); id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);

		id_h33[i] = pardiso_ia[f3]; id_h34[i] = pardiso_ia[min01 + V_N] + find_id_in_rows.coeff(min01 + V_N, max01 + V_N); id_h35[i] = pardiso_ia[min02 + V_N] + find_id_in_rows.coeff(min02 + V_N, max02 + V_N);

		id_h44[i] = pardiso_ia[f4]; id_h45[i] = pardiso_ia[min12 + V_N] + find_id_in_rows.coeff(min12 + V_N, max12 + V_N);

		id_h55[i] = pardiso_ia[f5];

	}

	for (int i = d_.m_T.rows(); i < F_N; i++)
	{
		int f0 = id2index[F0[i]]; int f1 = id2index[F1[i]]; int f2 = id2index[F2[i]];
		int f3 = f0 + V_N; int f4 = f1 + V_N; int f5 = f2 + V_N;
		if (f0 != -1)
		{
			id_h00[i] = pardiso_ia[f0];
			id_h33[i] = pardiso_ia[f3];
			id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3);
		}
		if (f1 != -1)
		{
			id_h11[i] = pardiso_ia[f1];
			id_h44[i] = pardiso_ia[f4];
			id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4);
		}
		if (f2 != -1)
		{
			id_h22[i] = pardiso_ia[f2];
			id_h55[i] = pardiso_ia[f5];
			id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);
		}

		if (f1 != -1 && f2 != -1)
		{
			int min12 = min(f1, f2); int max12 = f1 + f2 - min12;
			id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
			id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);
			id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4);
			id_h45[i] = pardiso_ia[min12 + V_N] + find_id_in_rows.coeff(min12 + V_N, max12 + V_N);
		}
		if (f0 != -1 && f2 != -1)
		{
			int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
			id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
			id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);
			id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3);
			id_h35[i] = pardiso_ia[min02 + V_N] + find_id_in_rows.coeff(min02 + V_N, max02 + V_N);
		}
		if (f1 != -1 && f0 != -1)
		{
			int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
			id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01);
			id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4);
			id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3);
			id_h34[i] = pardiso_ia[min01 + V_N] + find_id_in_rows.coeff(min01 + V_N, max01 + V_N);
		}

	}


	if (Plus_area_term)
	{
		id_bnd_h.resize(bnds.size());
		for (size_t i = 0; i < bnds.size(); i++)
		{
			int bnds_i_size = bnds[i].size();
			id_bnd_h[i].resize(bnds_i_size*(2 * bnds_i_size + 1));
			int index = 0;
			for (size_t j = 0; j < bnds_i_size; j++)
			{
				id_bnd_h[i][index] = pardiso_ia[bnds[i][j]];
				index++;
				for (size_t k = j + 1; k < bnds_i_size; k++)
				{
					int v1 = bnds[i][j];
					int v2 = bnds[i][k];
					if (v1 > v2) swap(v1, v2);
					id_bnd_h[i][index] = pardiso_ia[v1] + find_id_in_rows.coeff(v1, v2);
					index++;
				}
				for (size_t k = 0; k < bnds_i_size; k++)
				{
					int v1 = bnds[i][j];
					int v2 = bnds[i][k] + V_N;
					id_bnd_h[i][index] = pardiso_ia[v1] + find_id_in_rows.coeff(v1, v2);
					index++;
				}
			}

			for (size_t j = 0; j < bnds_i_size; j++)
			{
				id_bnd_h[i][index] = pardiso_ia[bnds[i][j]];
				index++;
				for (size_t k = j + 1; k < bnds_i_size; k++)
				{
					int v1 = bnds[i][j];
					int v2 = bnds[i][k];
					if (v1 > v2) swap(v1, v2);
					id_bnd_h[i][index] = pardiso_ia[v1 + V_N] + find_id_in_rows.coeff(v1 + V_N, v2 + V_N);
					index++;
				}
			}
		}

	}
}

//void Parafun::Update_source_same_t()
//{
//	double t_min = 1;
//	int geqK = 0;
//
//	int update_fn = d_.m_T.rows();
//	vector<double> all_s0; all_s0.resize(update_fn);
//	vector<double> all_s1; all_s1.resize(update_fn);
//
//	vector<double> all_w00; all_w00.resize(update_fn);
//	vector<double> all_w01; all_w01.resize(update_fn);
//	vector<double> all_w10; all_w10.resize(update_fn);
//	vector<double> all_w11; all_w11.resize(update_fn);
//
//
//	int f0, f1, f2;
//	double x0, y0, x1, y1, x2, y2;
//	double det;
//	double E_d;
//	double tt;
//	double new_sig0, new_sig1;
//	double j00, j01, j10, j11;
//	double p00, p01, p10, p11;
//	double q00, q01, q10, q11;
//
//	double *position = position_of_mesh.data();
//
//	for (int i = 0; i < update_fn; ++i)
//	{
//		f0 = F0[i];
//		f1 = F1[i];
//		f2 = F2[i];
//
//		x0 = position[f0];
//		y0 = position[f0 + total_num];
//
//		x1 = position[f1];
//		y1 = position[f1 + total_num];
//
//		x2 = position[f2];
//		y2 = position[f2 + total_num];
//
//		q00 = x1 - x0; q01 = x2 - x0;
//		q10 = y1 - y0; q11 = y2 - y0;
//
//		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];
//
//		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;
//
//
//		det = j00 * j11 - j01 * j10;
//		E_d = (1 + 1 / (det*det)) * (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
//
//		double alpha_0 = j00 + j11; double alpha_1 = j10 - j01;
//		double beta_0 = j00 - j11; double beta_1 = j10 + j01;
//
//		double alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
//		double beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1 * beta_1);
//
//		double sig0 = alpha_norm + beta_norm;
//		double sig1 = alpha_norm - beta_norm;
//		all_s0[i] = sig0;
//		all_s1[i] = sig1;
//
//		if (beta_norm < 1e-15)
//		{
//			all_w00[i] = 0.0;
//			all_w01[i] = 0.0;
//			all_w10[i] = 0.0;
//			all_w11[i] = 0.0;
//		}
//		else
//		{
//			double temp = 1 / (sig1*sig1 - sig0 * sig0);
//			all_w00[i] = temp * (j00*j00 + j10 * j10 - 0.5*(sig0*sig0 + sig1 * sig1));
//			all_w01[i] = temp * (j00*j01 + j10 * j11);
//			all_w10[i] = temp * (j01*j00 + j11 * j10);
//			all_w11[i] = temp * (j01*j01 + j11 * j11 - 0.5*(sig0*sig0 + sig1 * sig1));
//		}
//
//
//
//		if (E_d <= bound_distortion_K)
//		{
//			geqK++;
//		}
//		else
//		{
//			tt = newton_equation(sig0, sig1, bound_distortion_K);
//			if (tt < t_min)
//			{
//				t_min = tt;
//			}
//		}
//	}
//
//	changetocm_flag = (double)geqK / update_fn;
//
//	update_p00 = source_p00;
//	update_p01 = source_p01;
//	update_p10 = source_p10;
//	update_p11 = source_p11;
//
//	for (int i = 0; i < update_fn; ++i)
//	{
//		double sig0 = all_s0[i];
//		double sig1 = all_s1[i];
//
//		new_sig0 = pow(sig0, t_min - 1);
//		new_sig1 = pow(sig1, t_min - 1);
//
//		double delta_new = new_sig1 - new_sig0;
//		double plus_new = 0.5*(new_sig1 + new_sig0);
//
//		double w00 = delta_new * all_w00[i] + plus_new;
//		double w01 = delta_new * all_w01[i];
//		double w10 = delta_new * all_w10[i];
//		double w11 = delta_new * all_w11[i] + plus_new;
//
//		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];
//
//		update_p00[i] = p00 * w00 + p01 * w10;
//		update_p01[i] = p00 * w01 + p01 * w11;
//		update_p10[i] = p10 * w00 + p11 * w10;
//		update_p11[i] = p10 * w01 + p11 * w11;
//	}
//
//	Intp_T_Min = t_min;
//}

void Parafun::SLIM()
{
	cout << "perform_iteration_slim begin-----" << endl;

	double area_now;
	int f0, f1, f2;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double x0, y0, x1, y1, x2, y2;

	double alpha_norm, beta_norm;

	double alpha_0, alpha_1, beta_0, beta_1;

	double sig0, sig1;

	double det, tr;
	double r0, r1, r2, r3;
	double d00, d01, d02,
		d10, d11, d12;

	double new_sig0, new_sig1;
	double temp;
	double w00, w01, w10, w11;
	double p1, p2, p3, w1, w2, w3;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double *position = position_of_mesh.data();

	int nnz = pardiso_ja.size();
	pardiso_a.clear(); pardiso_b.clear();
	pardiso_a.resize(nnz, 0.0);
	pardiso_b.resize(2 * V_N, 0.0);
	int src_t_num = d_.m_T.rows();
	for (int i = 0; i < src_t_num; ++i)
	{
		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];


		x0 = position[f0];
		y0 = position[f0 + total_num];

		x1 = position[f1];
		y1 = position[f1 + total_num];

		x2 = position[f2];
		y2 = position[f2 + total_num];

		f0 = id2index[f0];
		f1 = id2index[f1];
		f2 = id2index[f2];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11; beta_1 = j10 + j01;

		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1*alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1*beta_1);

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0*sig0) + 1 / (sig0*sig0*sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1*sig1) + 1 / (sig1*sig1*sig1));

		if (beta_norm < 1e-10)
		{
			temp = 0;
		}
		else
		{
			temp = (new_sig1 - new_sig0) / (sig1*sig1 - sig0 * sig0);
		}

		w00 = temp*(j00*j00 + j01*j01 - 0.5*(sig0*sig0 + sig1*sig1)) + 0.5*(new_sig0 + new_sig1);
		w01 = temp*(j00*j10 + j01*j11);
		w10 = temp*(j10*j00 + j11*j01);
		w11 = temp*(j10*j10 + j11*j11 - 0.5*(sig0*sig0 + sig1*sig1)) + 0.5*(new_sig0 + new_sig1);

		p1 = p00*p00 + p01*p01; p2 = p00*p10 + p01*p11; p3 = p10*p10 + p11*p11;
		w1 = w00*w00 + w10*w10; w2 = w00*w01 + w10*w11; w3 = w01*w01 + w11*w11;

		//area_now *= 2;

		h00 = area_now *(p1 + p2 + p2 + p3)*w1; h01 = -area_now *(p1 + p2)*w1; h02 = -area_now *(p2 + p3)*w1; h03 = area_now *(p1 + p2 + p2 + p3)*w2; h04 = -area_now *(p1 + p2)*w2; h05 = -area_now *(p2 + p3)*w2;
		h11 = area_now *p1*w1;                  h12 = area_now *p2*w1;    	 h13 = -area_now *(p1 + p2)*w2; h14 = area_now *p1*w2;                  h15 = area_now *p2*w2;
		h22 = area_now *p3*w1;                  h23 = -area_now *(p2 + p3)*w2; h24 = area_now *p2*w2;         h25 = area_now *p3*w2;
		h33 = area_now *(p1 + p2 + p2 + p3)*w3; h34 = -area_now *(p1 + p2)*w3; h35 = -area_now *(p2 + p3)*w3;
		h44 = area_now *p1*w3;                  h45 = area_now *p2*w3;
		h55 = area_now *p3*w3;


		det = j00*j11 - j01*j10;
		tr = (j00*j00 + j01*j01 + j10*j10 + j11*j11);

		d00 = -p00 - p10; d01 = p00; d02 = p10;
		d10 = -p01 - p11; d11 = p01; d12 = p11;

		r0 = area_now * ((1 + 1 / (det*det))*j00 - tr*j11 / (det*det*det));
		r1 = area_now * ((1 + 1 / (det*det))*j01 + tr*j10 / (det*det*det));
		r2 = area_now * ((1 + 1 / (det*det))*j10 + tr*j01 / (det*det*det));
		r3 = area_now * ((1 + 1 / (det*det))*j11 - tr*j00 / (det*det*det));


		pardiso_b[f0] -= r0*d00 + r1*d10;
		pardiso_b[f1] -= r0*d01 + r1*d11;
		pardiso_b[f2] -= r0*d02 + r1*d12;
		pardiso_b[f0 + V_N] -= r2*d00 + r3*d10;
		pardiso_b[f1 + V_N] -= r2*d01 + r3*d11;
		pardiso_b[f2 + V_N] -= r2*d02 + r3*d12;

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;

	}

	for (int i = src_t_num; i < F_N; ++i)
	{
		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];


		x0 = position[f0];
		y0 = position[f0 + total_num];

		x1 = position[f1];
		y1 = position[f1 + total_num];

		x2 = position[f2];
		y2 = position[f2 + total_num];

		f0 = id2index[f0];
		f1 = id2index[f1];
		f2 = id2index[f2];
		if (f0 == -1&& f1 == -1&& f2 == -1)
			continue;

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11; beta_1 = j10 + j01;

		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1 * beta_1);

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0*sig0) + 1 / (sig0*sig0*sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1*sig1) + 1 / (sig1*sig1*sig1));

		if (abs(sig1 - sig0) < 1e-10)
		{
			temp = 0;
		}
		else
		{
			temp = (new_sig1 - new_sig0) / (sig1*sig1 - sig0 * sig0);
		}

		w00 = temp * (j00*j00 + j01 * j01 - 0.5*(sig0*sig0 + sig1 * sig1)) + 0.5*(new_sig0 + new_sig1);
		w01 = temp * (j00*j10 + j01 * j11);
		w10 = temp * (j10*j00 + j11 * j01);
		w11 = temp * (j10*j10 + j11 * j11 - 0.5*(sig0*sig0 + sig1 * sig1)) + 0.5*(new_sig0 + new_sig1);

		p1 = p00 * p00 + p01 * p01; p2 = p00 * p10 + p01 * p11; p3 = p10 * p10 + p11 * p11;
		w1 = w00 * w00 + w10 * w10; w2 = w00 * w01 + w10 * w11; w3 = w01 * w01 + w11 * w11;

		//area_now *= 2;

		h00 = area_now * (p1 + p2 + p2 + p3)*w1; h01 = -area_now * (p1 + p2)*w1; h02 = -area_now * (p2 + p3)*w1; h03 = area_now * (p1 + p2 + p2 + p3)*w2; h04 = -area_now * (p1 + p2)*w2; h05 = -area_now * (p2 + p3)*w2;
		h11 = area_now * p1*w1;                  h12 = area_now * p2*w1;    	 h13 = -area_now * (p1 + p2)*w2; h14 = area_now * p1*w2;                  h15 = area_now * p2*w2;
		h22 = area_now * p3*w1;                  h23 = -area_now * (p2 + p3)*w2; h24 = area_now * p2*w2;         h25 = area_now * p3*w2;
		h33 = area_now * (p1 + p2 + p2 + p3)*w3; h34 = -area_now * (p1 + p2)*w3; h35 = -area_now * (p2 + p3)*w3;
		h44 = area_now * p1*w3;                  h45 = area_now * p2*w3;
		h55 = area_now * p3*w3;


		det = j00 * j11 - j01 * j10;
		tr = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);

		d00 = -p00 - p10; d01 = p00; d02 = p10;
		d10 = -p01 - p11; d11 = p01; d12 = p11;

		r0 = area_now * ((1 + 1 / (det*det))*j00 - tr * j11 / (det*det*det));
		r1 = area_now * ((1 + 1 / (det*det))*j01 + tr * j10 / (det*det*det));
		r2 = area_now * ((1 + 1 / (det*det))*j10 + tr * j01 / (det*det*det));
		r3 = area_now * ((1 + 1 / (det*det))*j11 - tr * j00 / (det*det*det));

		if (f0 != -1)
		{
			pardiso_b[f0] -= r0 * d00 + r1 * d10;
			pardiso_b[f0 + V_N] -= r2 * d00 + r3 * d10;
			pardiso_a[id_h00[i]] += h00;
			pardiso_a[id_h33[i]] += h33;
			pardiso_a[id_h03[i]] += h03;
		}
		if (f1 != -1)
		{
			pardiso_b[f1] -= r0 * d01 + r1 * d11;
			pardiso_b[f1 + V_N] -= r2 * d01 + r3 * d11;
			pardiso_a[id_h11[i]] += h11;
			pardiso_a[id_h44[i]] += h44;
			pardiso_a[id_h14[i]] += h14;
		}
		if (f2 != -1)
		{
			pardiso_b[f2] -= r0 * d02 + r1 * d12;
			pardiso_b[f2 + V_N] -= r2 * d02 + r3 * d12;
			pardiso_a[id_h22[i]] += h22;
			pardiso_a[id_h55[i]] += h55;
			pardiso_a[id_h25[i]] += h25;
		}

		if (f0 != -1 && f1 != -1)
		{
			pardiso_a[id_h01[i]] += h01;
			pardiso_a[id_h04[i]] += h04;
			pardiso_a[id_h13[i]] += h13;
			pardiso_a[id_h34[i]] += h34;
		}
		if (f0 != -1 && f2 != -1)
		{
			pardiso_a[id_h02[i]] += h02;
			pardiso_a[id_h05[i]] += h05;
			pardiso_a[id_h23[i]] += h23;
			pardiso_a[id_h35[i]] += h35;
		}
		if (f2 != -1 && f1 != -1)
		{
			pardiso_a[id_h12[i]] += h12;
			pardiso_a[id_h15[i]] += h15;
			pardiso_a[id_h24[i]] += h24;
			pardiso_a[id_h45[i]] += h45;
		}

	}
	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;

	pardiso->factorize();
	pardiso->pardiso_solver();

	vector<double> result_d = pardiso->result;

	VectorXd negative_grad(2 * total_num), d(2 * total_num);
	negative_grad.setZero();
	d.setZero();

	for (int i = 0; i < V_N; i++)
	{
		negative_grad(var_ids[i]) = pardiso_b[i];
		negative_grad(var_ids[i] + total_num) = pardiso_b[i + V_N];
		d(var_ids[i]) = result_d[i];
		d(var_ids[i] + total_num) = result_d[i + V_N];
	}


	double temp_t;
 	max_step(position_of_mesh, d, temp_t);

	double alpha = min(1.0, 0.8 * temp_t);
	backtracking_line_search(position_of_mesh, d, negative_grad, alpha);
	position_of_mesh += alpha * d;
	cout << "slim step length : " << alpha << endl;

	Energysource();
}
//void Parafun::SLIM(bool is_interp)
//{
//	cout << "perform_iteration_slim begin-----" << is_interp /*<< "---" << is_remesh*/ << endl;
//
//	double area_now;
//	int f0, f1, f2;
//	double j00, j01, j10, j11;
//	double p00, p01, p10, p11;
//	double q00, q01, q10, q11;
//
//	double x0, y0, x1, y1, x2, y2;
//
//	double alpha_norm, beta_norm;
//
//	double alpha_0, alpha_1, beta_0, beta_1;
//
//	double sig0, sig1;
//
//	double det, tr;
//	double r0, r1, r2, r3;
//	double d00, d01, d02,
//		d10, d11, d12;
//
//	double new_sig0, new_sig1;
//	double temp;
//	double w00, w01, w10, w11;
//	double p1, p2, p3, w1, w2, w3;
//
//	double h00, h01, h02, h03, h04, h05,
//		h11, h12, h13, h14, h15,
//		h22, h23, h24, h25,
//		h33, h34, h35,
//		h44, h45,
//		h55;
//	double *position = position_of_mesh.data();
//
//	int nnz = pardiso_ja.size();
//	pardiso_a.clear(); pardiso_b.clear();
//	pardiso_a.resize(nnz, 0.0);
//	pardiso_b.resize(2 * V_N, 0.0);
//
//	double* tmp_p00;
//	double* tmp_p01;
//	double* tmp_p10;
//	double* tmp_p11;
//
//	if (is_interp)
//	{
//		tmp_p00 = update_p00.data();
//		tmp_p01 = update_p01.data();
//		tmp_p10 = update_p10.data();
//		tmp_p11 = update_p11.data();
//	}
//	else
//	{
//		tmp_p00 = source_p00.data();
//		tmp_p01 = source_p01.data();
//		tmp_p10 = source_p10.data();
//		tmp_p11 = source_p11.data();
//	}
//
//	int src_t_num = d_.m_T.rows();
//	for (int i = 0; i < src_t_num; ++i)
//	{
//		area_now = area[i];
//		f0 = F0[i];
//		f1 = F1[i];
//		f2 = F2[i];
//
//
//		x0 = position[f0];
//		y0 = position[f0 + total_num];
//
//		x1 = position[f1];
//		y1 = position[f1 + total_num];
//
//		x2 = position[f2];
//		y2 = position[f2 + total_num];
//
//		f0 = id2index[f0];
//		f1 = id2index[f1];
//		f2 = id2index[f2];
//
//		q00 = x1 - x0; q01 = x2 - x0;
//		q10 = y1 - y0; q11 = y2 - y0;
//
//		p00 = tmp_p00[i]; p01 = tmp_p01[i]; p10 = tmp_p10[i]; p11 = tmp_p11[i];
//
//		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;
//
//		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
//		beta_0 = j00 - j11; beta_1 = j10 + j01;
//
//		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
//		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1 * beta_1);
//
//		sig0 = alpha_norm + beta_norm;
//		sig1 = alpha_norm - beta_norm;
//
//		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0*sig0) + 1 / (sig0*sig0*sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1*sig1) + 1 / (sig1*sig1*sig1));
//
//		if (beta_norm < 1e-6)
//		{
//			temp = 0;
//		}
//		else
//		{
//			temp = (new_sig1 - new_sig0) / (sig1*sig1 - sig0 * sig0);
//		}
//
//		w00 = temp * (j00*j00 + j01 * j01 - 0.5*(sig0*sig0 + sig1 * sig1)) + 0.5*(new_sig0 + new_sig1);
//		w01 = temp * (j00*j10 + j01 * j11);
//		w10 = temp * (j10*j00 + j11 * j01);
//		w11 = temp * (j10*j10 + j11 * j11 - 0.5*(sig0*sig0 + sig1 * sig1)) + 0.5*(new_sig0 + new_sig1);
//
//		p1 = p00 * p00 + p01 * p01; p2 = p00 * p10 + p01 * p11; p3 = p10 * p10 + p11 * p11;
//		w1 = w00 * w00 + w10 * w10; w2 = w00 * w01 + w10 * w11; w3 = w01 * w01 + w11 * w11;
//
//		//area_now *= 2;
//
//		h00 = area_now * (p1 + p2 + p2 + p3)*w1; h01 = -area_now * (p1 + p2)*w1; h02 = -area_now * (p2 + p3)*w1; h03 = area_now * (p1 + p2 + p2 + p3)*w2; h04 = -area_now * (p1 + p2)*w2; h05 = -area_now * (p2 + p3)*w2;
//		h11 = area_now * p1*w1;                  h12 = area_now * p2*w1;    	 h13 = -area_now * (p1 + p2)*w2; h14 = area_now * p1*w2;                  h15 = area_now * p2*w2;
//		h22 = area_now * p3*w1;                  h23 = -area_now * (p2 + p3)*w2; h24 = area_now * p2*w2;         h25 = area_now * p3*w2;
//		h33 = area_now * (p1 + p2 + p2 + p3)*w3; h34 = -area_now * (p1 + p2)*w3; h35 = -area_now * (p2 + p3)*w3;
//		h44 = area_now * p1*w3;                  h45 = area_now * p2*w3;
//		h55 = area_now * p3*w3;
//
//
//		det = j00 * j11 - j01 * j10;
//		tr = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
//
//		d00 = -p00 - p10; d01 = p00; d02 = p10;
//		d10 = -p01 - p11; d11 = p01; d12 = p11;
//
//		r0 = area_now * ((1 + 1 / (det*det))*j00 - tr * j11 / (det*det*det));
//		r1 = area_now * ((1 + 1 / (det*det))*j01 + tr * j10 / (det*det*det));
//		r2 = area_now * ((1 + 1 / (det*det))*j10 + tr * j01 / (det*det*det));
//		r3 = area_now * ((1 + 1 / (det*det))*j11 - tr * j00 / (det*det*det));
//
//
//		pardiso_b[f0] -= r0 * d00 + r1 * d10;
//		pardiso_b[f1] -= r0 * d01 + r1 * d11;
//		pardiso_b[f2] -= r0 * d02 + r1 * d12;
//		pardiso_b[f0 + V_N] -= r2 * d00 + r3 * d10;
//		pardiso_b[f1 + V_N] -= r2 * d01 + r3 * d11;
//		pardiso_b[f2 + V_N] -= r2 * d02 + r3 * d12;
//
//		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
//		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
//		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
//		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
//		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
//		pardiso_a[id_h55[i]] += h55;
//
//	}
//
//	for (int i = src_t_num; i < F_N; ++i)
//	{
//		area_now = area[i];
//		f0 = F0[i];
//		f1 = F1[i];
//		f2 = F2[i];
//
//
//		x0 = position[f0];
//		y0 = position[f0 + total_num];
//
//		x1 = position[f1];
//		y1 = position[f1 + total_num];
//
//		x2 = position[f2];
//		y2 = position[f2 + total_num];
//
//		f0 = id2index[f0];
//		f1 = id2index[f1];
//		f2 = id2index[f2];
//		if (f0 == -1 && f1 == -1 && f2 == -1)
//			continue;
//
//		q00 = x1 - x0; q01 = x2 - x0;
//		q10 = y1 - y0; q11 = y2 - y0;
//
//		p00 = tmp_p00[i]; p01 = tmp_p01[i]; p10 = tmp_p10[i]; p11 = tmp_p11[i];
//
//		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;
//
//		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
//		beta_0 = j00 - j11; beta_1 = j10 + j01;
//
//		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
//		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1 * beta_1);
//
//		sig0 = alpha_norm + beta_norm;
//		sig1 = alpha_norm - beta_norm;
//
//		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0*sig0) + 1 / (sig0*sig0*sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1*sig1) + 1 / (sig1*sig1*sig1));
//
//		if (abs(sig1 - sig0) < 1e-10)
//		{
//			temp = 0;
//		}
//		else
//		{
//			temp = (new_sig1 - new_sig0) / (sig1*sig1 - sig0 * sig0);
//		}
//
//		w00 = temp * (j00*j00 + j01 * j01 - 0.5*(sig0*sig0 + sig1 * sig1)) + 0.5*(new_sig0 + new_sig1);
//		w01 = temp * (j00*j10 + j01 * j11);
//		w10 = temp * (j10*j00 + j11 * j01);
//		w11 = temp * (j10*j10 + j11 * j11 - 0.5*(sig0*sig0 + sig1 * sig1)) + 0.5*(new_sig0 + new_sig1);
//
//		p1 = p00 * p00 + p01 * p01; p2 = p00 * p10 + p01 * p11; p3 = p10 * p10 + p11 * p11;
//		w1 = w00 * w00 + w10 * w10; w2 = w00 * w01 + w10 * w11; w3 = w01 * w01 + w11 * w11;
//
//		//area_now *= 2;
//
//		h00 = area_now * (p1 + p2 + p2 + p3)*w1; h01 = -area_now * (p1 + p2)*w1; h02 = -area_now * (p2 + p3)*w1; h03 = area_now * (p1 + p2 + p2 + p3)*w2; h04 = -area_now * (p1 + p2)*w2; h05 = -area_now * (p2 + p3)*w2;
//		h11 = area_now * p1*w1;                  h12 = area_now * p2*w1;    	 h13 = -area_now * (p1 + p2)*w2; h14 = area_now * p1*w2;                  h15 = area_now * p2*w2;
//		h22 = area_now * p3*w1;                  h23 = -area_now * (p2 + p3)*w2; h24 = area_now * p2*w2;         h25 = area_now * p3*w2;
//		h33 = area_now * (p1 + p2 + p2 + p3)*w3; h34 = -area_now * (p1 + p2)*w3; h35 = -area_now * (p2 + p3)*w3;
//		h44 = area_now * p1*w3;                  h45 = area_now * p2*w3;
//		h55 = area_now * p3*w3;
//
//
//		det = j00 * j11 - j01 * j10;
//		tr = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
//
//		d00 = -p00 - p10; d01 = p00; d02 = p10;
//		d10 = -p01 - p11; d11 = p01; d12 = p11;
//
//		r0 = area_now * ((1 + 1 / (det*det))*j00 - tr * j11 / (det*det*det));
//		r1 = area_now * ((1 + 1 / (det*det))*j01 + tr * j10 / (det*det*det));
//		r2 = area_now * ((1 + 1 / (det*det))*j10 + tr * j01 / (det*det*det));
//		r3 = area_now * ((1 + 1 / (det*det))*j11 - tr * j00 / (det*det*det));
//
//		if (f0 != -1)
//		{
//			pardiso_b[f0] -= r0 * d00 + r1 * d10;
//			pardiso_b[f0 + V_N] -= r2 * d00 + r3 * d10;
//			pardiso_a[id_h00[i]] += h00;
//			pardiso_a[id_h33[i]] += h33;
//			pardiso_a[id_h03[i]] += h03;
//		}
//		if (f1 != -1)
//		{
//			pardiso_b[f1] -= r0 * d01 + r1 * d11;
//			pardiso_b[f1 + V_N] -= r2 * d01 + r3 * d11;
//			pardiso_a[id_h11[i]] += h11;
//			pardiso_a[id_h44[i]] += h44;
//			pardiso_a[id_h14[i]] += h14;
//		}
//		if (f2 != -1)
//		{
//			pardiso_b[f2] -= r0 * d02 + r1 * d12;
//			pardiso_b[f2 + V_N] -= r2 * d02 + r3 * d12;
//			pardiso_a[id_h22[i]] += h22;
//			pardiso_a[id_h55[i]] += h55;
//			pardiso_a[id_h25[i]] += h25;
//		}
//
//		if (f0 != -1 && f1 != -1)
//		{
//			pardiso_a[id_h01[i]] += h01;
//			pardiso_a[id_h04[i]] += h04;
//			pardiso_a[id_h13[i]] += h13;
//			pardiso_a[id_h34[i]] += h34;
//		}
//		if (f0 != -1 && f2 != -1)
//		{
//			pardiso_a[id_h02[i]] += h02;
//			pardiso_a[id_h05[i]] += h05;
//			pardiso_a[id_h23[i]] += h23;
//			pardiso_a[id_h35[i]] += h35;
//		}
//		if (f2 != -1 && f1 != -1)
//		{
//			pardiso_a[id_h12[i]] += h12;
//			pardiso_a[id_h15[i]] += h15;
//			pardiso_a[id_h24[i]] += h24;
//			pardiso_a[id_h45[i]] += h45;
//		}
//
//	}
//	pardiso->a = pardiso_a;
//	pardiso->rhs = pardiso_b;
//
//	pardiso->factorize();
//	pardiso->pardiso_solver();
//
//	vector<double> result_d = pardiso->result;
//
//	VectorXd negative_grad(2 * total_num), d(2 * total_num);
//	negative_grad.setZero();
//	d.setZero();
//
//	for (int i = 0; i < V_N; i++)
//	{
//		negative_grad(var_ids[i]) = pardiso_b[i];
//		negative_grad(var_ids[i] + total_num) = pardiso_b[i + V_N];
//		d(var_ids[i]) = result_d[i];
//		d(var_ids[i] + total_num) = result_d[i + V_N];
//	}
//
//
//	double temp_t;
//	max_step(position_of_mesh, d, temp_t);
//
//	double alpha = min(1.0, 0.8 * temp_t);
//	backtracking_line_search(position_of_mesh, d, negative_grad, alpha);
//	position_of_mesh += alpha * d;
//	cout << "slim step length : " << alpha << endl;
//
//	Energysource();
//}


void Parafun::CM()
{
	//std::cout << "perform_iteration_cm begin-----" << endl;
	double area_now;
	int f0, f1, f2;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double x0, y0, x1, y1, x2, y2;

	double hi_0, hi_1;

	double alpha_0, alpha_1, beta_0, beta_1;

	double s1, s2, sig0, sig1;

	double alpha_norm, beta_norm;
	double h_u, h_v, walpha, wbeta;

	double a1x0, a1x1, a1x2, a1x3, a1x4, a1x5,
		a2x0, a2x1, a2x2, a2x3, a2x4, a2x5;

	double aa, bb;
	double uu, vv, uv;
	double u, v;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;

	double *position = position_of_mesh.data();
	int nnz = pardiso_ja.size();
	pardiso_a.clear(); pardiso_b.clear();
	pardiso_a.resize(nnz, 0.0);
	pardiso_b.resize(2 * V_N, 0.0);

	if (Plus_area_term)
	{
		calc_areaTerm_hessian_gradient();
	}
	if (Plus_shrink_term)
	{
		calc_cur_bbox(position_of_mesh);
		separate_scaffold();
	}

	for (int i = 0; i < d_.m_T.rows(); i++)
	{
		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = position[f0];
		y0 = position[f0 + total_num];

		x1 = position[f1];
		y1 = position[f1 + total_num];

		x2 = position[f2];
		y2 = position[f2 + total_num];

		f0 = id2index[f0];
		f1 = id2index[f1];
		f2 = id2index[f2];


		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;
		//pp
		//p00 = update_p00[i]; p01 = update_p01[i]; p10 = update_p10[i]; p11 = update_p11[i];
		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11;  beta_1 = j10 + j01;

		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1*alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1*beta_1);

		s1 = (p00)*(p00 + p10) + (p01)*(p01 + p11);
		s2 = (p10)*(p00 + p10) + (p11)*(p01 + p11);

		double h1 = p00*p00 + p01*p01;
		double h2 = p00*p10 + p01*p11;
		double h3 = p10*p10 + p11*p11;
		double h4 = p00*p11 - p01*p10;

		a1x0 = alpha_0*(-p00 - p10) + alpha_1*(p01 + p11);  a1x1 = alpha_0*p00 - alpha_1*p01; a1x2 = alpha_0*p10 - alpha_1*p11;
		a1x3 = alpha_0*(-p01 - p11) + alpha_1*(-p00 - p10); a1x4 = alpha_0*p01 + alpha_1*p00; a1x5 = alpha_0*p11 + alpha_1*p10;

		a2x0 = beta_0*(-p00 - p10) + beta_1*(-p01 - p11);   a2x1 = beta_0*p00 + beta_1*p01;   a2x2 = beta_0*p10 + beta_1*p11;
		a2x3 = beta_0*(p01 + p11) + beta_1*(-p00 - p10);    a2x4 = -beta_0*p01 + beta_1*p00;  a2x5 = -beta_0*p11 + beta_1*p10;

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		hi_0 = 2 + 6 * 1 / (sig0*sig0*sig0*sig0); hi_1 = 2 + 6 * 1 / (sig1*sig1*sig1*sig1);

		aa = 0.25 / alpha_norm; bb = 0.25 / beta_norm;

		uu = aa * aa*(area_now*hi_0 + area_now * hi_1);
		vv = bb * bb*(area_now*hi_0 + area_now * hi_1);
		uv = aa * bb*(area_now*hi_0 - area_now * hi_1);
		h_u = area_now * (2 * sig0 - 2 * 1 / (sig0*sig0*sig0));
		h_v = area_now * (2 * sig1 - 2 * 1 / (sig1*sig1*sig1));
		walpha = h_u + h_v;
		wbeta = h_u - h_v;
		double hwa1 = (walpha * 0.25 / alpha_norm); double hwa2 = -(walpha * 0.25*0.25 / (alpha_norm*alpha_norm*alpha_norm));
		double hwb1 = (wbeta * 0.25 / beta_norm); double hwb2 = -(wbeta *0.25*0.25 / (beta_norm*beta_norm*beta_norm));

		h00 = uu * a1x0*a1x0 + vv * a2x0*a2x0 + uv * a1x0*a2x0 + uv * a2x0*a1x0; h01 = uu * a1x0*a1x1 + vv * a2x0*a2x1 + uv * a1x0*a2x1 + uv * a2x0*a1x1; h02 = uu * a1x0*a1x2 + vv * a2x0*a2x2 + uv * a1x0*a2x2 + uv * a2x0*a1x2; h03 = uu * a1x0*a1x3 + vv * a2x0*a2x3 + uv * a1x0*a2x3 + uv * a2x0*a1x3; h04 = uu * a1x0*a1x4 + vv * a2x0*a2x4 + uv * a1x0*a2x4 + uv * a2x0*a1x4; h05 = uu * a1x0*a1x5 + vv * a2x0*a2x5 + uv * a1x0*a2x5 + uv * a2x0*a1x5;

		h11 = uu * a1x1*a1x1 + vv * a2x1*a2x1 + uv * a1x1*a2x1 + uv * a2x1*a1x1; h12 = uu * a1x1*a1x2 + vv * a2x1*a2x2 + uv * a1x1*a2x2 + uv * a2x1*a1x2; h13 = uu * a1x1*a1x3 + vv * a2x1*a2x3 + uv * a1x1*a2x3 + uv * a2x1*a1x3; h14 = uu * a1x1*a1x4 + vv * a2x1*a2x4 + uv * a1x1*a2x4 + uv * a2x1*a1x4; h15 = uu * a1x1*a1x5 + vv * a2x1*a2x5 + uv * a1x1*a2x5 + uv * a2x1*a1x5;

		h22 = uu * a1x2*a1x2 + vv * a2x2*a2x2 + uv * a1x2*a2x2 + uv * a2x2*a1x2; h23 = uu * a1x2*a1x3 + vv * a2x2*a2x3 + uv * a1x2*a2x3 + uv * a2x2*a1x3; h24 = uu * a1x2*a1x4 + vv * a2x2*a2x4 + uv * a1x2*a2x4 + uv * a2x2*a1x4; h25 = uu * a1x2*a1x5 + vv * a2x2*a2x5 + uv * a1x2*a2x5 + uv * a2x2*a1x5;

		h33 = uu * a1x3*a1x3 + vv * a2x3*a2x3 + uv * a1x3*a2x3 + uv * a2x3*a1x3; h34 = uu * a1x3*a1x4 + vv * a2x3*a2x4 + uv * a1x3*a2x4 + uv * a2x3*a1x4; h35 = uu * a1x3*a1x5 + vv * a2x3*a2x5 + uv * a1x3*a2x5 + uv * a2x3*a1x5;

		h44 = uu * a1x4*a1x4 + vv * a2x4*a2x4 + uv * a1x4*a2x4 + uv * a2x4*a1x4; h45 = uu * a1x4*a1x5 + vv * a2x4*a2x5 + uv * a1x4*a2x5 + uv * a2x4*a1x5;

		h55 = uu * a1x5*a1x5 + vv * a2x5*a2x5 + uv * a1x5*a2x5 + uv * a2x5*a1x5;

		if (walpha >= 0)
		{
			h00 += hwa1*(s1 + s2) + hwa2*a1x0*a1x0; h01 += hwa1*(-s1) + hwa2*a1x0*a1x1; h02 += hwa1*(-s2) + hwa2*a1x0*a1x2; h03 += hwa2*a1x0*a1x3; h04 += hwa1*(h4)+hwa2*a1x0*a1x4; h05 += hwa1*(-h4) + hwa2*a1x0*a1x5;
			h11 += hwa1*(h1)+hwa2*a1x1*a1x1;        h12 += hwa1*(h2)+hwa2*a1x1*a1x2;    h13 += hwa1*(-h4) + hwa2*a1x1*a1x3; h14 += hwa2*a1x1*a1x4; h15 += hwa1*(h4)+hwa2*a1x1*a1x5;
			h22 += hwa1*(h3)+hwa2*a1x2*a1x2;        h23 += hwa1*(h4)+hwa2*a1x2*a1x3;    h24 += hwa1*(-h4) + hwa2*a1x2*a1x4; h25 += hwa2*a1x2*a1x5;
			h33 += hwa1*(s1 + s2) + hwa2*a1x3*a1x3; h34 += hwa1*(-s1) + hwa2*a1x3*a1x4; h35 += hwa1*(-s2) + hwa2*a1x3*a1x5;
			h44 += hwa1*(h1)+hwa2*a1x4*a1x4;        h45 += hwa1*(h2)+hwa2*a1x4*a1x5;
			h55 += hwa1*(h3)+hwa2*a1x5*a1x5;
		}
		h00 += hwb1 * (s1 + s2) + hwb2 * a2x0*a2x0; h01 += hwb1 * (-s1) + hwb2 * a2x0*a2x1; h02 += hwb1 * (-s2) + hwb2 * a2x0*a2x2; h03 += hwb2 * a2x0*a2x3; h04 += hwb1 * (-h4) + hwb2 * a2x0*a2x4; h05 += hwb1 * (h4)+hwb2 * a2x0*a2x5;
		h11 += hwb1 * (h1)+hwb2 * a2x1*a2x1;        h12 += hwb1 * (h2)+hwb2 * a2x1*a2x2;    h13 += hwb1 * (h4)+hwb2 * a2x1*a2x3;    h14 += hwb2 * a2x1*a2x4; h15 += hwb1 * (-h4) + hwb2 * a2x1*a2x5;
		h22 += hwb1 * (h3)+hwb2 * a2x2*a2x2;        h23 += hwb1 * (-h4) + hwb2 * a2x2*a2x3; h24 += hwb1 * (h4)+hwb2 * a2x2*a2x4;    h25 += hwb2 * a2x2*a2x5;
		h33 += hwb1 * (s1 + s2) + hwb2 * a2x3*a2x3; h34 += hwb1 * (-s1) + hwb2 * a2x3*a2x4; h35 += hwb1 * (-s2) + hwb2 * a2x3*a2x5;
		h44 += hwb1 * (h1)+hwb2 * a2x4*a2x4;        h45 += hwb1 * (h2)+hwb2 * a2x4*a2x5;
		h55 += hwb1 * (h3)+hwb2 * a2x5*a2x5;

		u = aa*walpha; v = bb*wbeta;
		pardiso_b[f0] -= (u*a1x0 + v*a2x0);
		pardiso_b[f1] -= (u*a1x1 + v*a2x1);
		pardiso_b[f2] -= (u*a1x2 + v*a2x2);
		pardiso_b[f0 + V_N] -= (u*a1x3 + v*a2x3);
		pardiso_b[f1 + V_N] -= (u*a1x4 + v*a2x4);
		pardiso_b[f2 + V_N] -= (u*a1x5 + v*a2x5);

		if (Plus_area_term)
		{
			double tmp_w = area_term_weight*area_src[i] / Cur_gap / 2.0;

			h00 += tmp_w * ((p00 + p10)*(p00 + p10) + (p01 + p11)*(p01 + p11));
			h01 -= tmp_w * ((p00 + p10)*p00 + (p01 + p11)*p01);
			h02 -= tmp_w * ((p00 + p10)*p10 + (p01 + p11)*p11);
			h03 += tmp_w * 2.0*(p00 + p10)*(p01 + p11);
			h04 -= tmp_w * ((p00 + p10)*p01 + (p01 + p11)*p00);
			h05 -= tmp_w * ((p00 + p10)*p11 + (p01 + p11)*p10);

			h11 += tmp_w * (p00*p00 + p01 * p01);
			h12 += tmp_w * (p00*p10 + p01 * p11);
			h13 -= tmp_w * (p00*(p01 + p11) + p01 * (p00 + p10));
			h14 += tmp_w * 2 * p00*p01;
			h15 += tmp_w * (p00*p11 + p01 * p10);

			h22 += tmp_w * (p10 * p10 + p11 * p11);
			h23 -= tmp_w * (p10 * (p01 + p11) + p11 * (p00 + p10));
			h24 += tmp_w * (p10 * p01 + p11 * p00);
			h25 += tmp_w * 2.0*p10 * p11;

			h33 += tmp_w * ((p01 + p11) * (p01 + p11) + (p00 + p10) * (p00 + p10));
			h34 -= tmp_w * ((p01 + p11) * p01 + (p00 + p10) * p00);
			h35 -= tmp_w * ((p01 + p11) * p11 + (p00 + p10) * p10);

			h44 += tmp_w * (p01 * p01 + p00 * p00);
			h45 += tmp_w * (p01 * p11 + p00 * p10);

			h45 += tmp_w * (p11 * p11 + p10 * p10);
		}

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;

	}


	//scaffold computation
	for (int i = d_.m_T.rows(); i < F_N; i++)
	{
		area_now = area[i];

		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];
		x0 = position[f0];
		y0 = position[f0 + total_num];

		x1 = position[f1];
		y1 = position[f1 + total_num];

		x2 = position[f2];
		y2 = position[f2 + total_num];

		f0 = id2index[f0];
		f1 = id2index[f1];
		f2 = id2index[f2];
		int flag_ = 0;
		if (f0 != -1)flag_++;
		if (f1 != -1)flag_++;
		if (f2 != -1)flag_++;

		if (flag_ == 0)continue;

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;
		//pp
		//p00 = update_p00[i]; p01 = update_p01[i]; p10 = update_p10[i]; p11 = update_p11[i];

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11;  beta_1 = j10 + j01;

		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1 * beta_1);
		if (beta_norm < 1e-30)
		{
			beta_norm = 1e-15;
		}

		s1 = (p00)*(p00 + p10) + (p01)*(p01 + p11);
		s2 = (p10)*(p00 + p10) + (p11)*(p01 + p11);

		double h1 = p00 * p00 + p01 * p01;
		double h2 = p00 * p10 + p01 * p11;
		double h3 = p10 * p10 + p11 * p11;
		double h4 = p00 * p11 - p01 * p10;

		a1x0 = alpha_0 * (-p00 - p10) + alpha_1 * (p01 + p11);  a1x1 = alpha_0 * p00 - alpha_1 * p01; a1x2 = alpha_0 * p10 - alpha_1 * p11;
		a1x3 = alpha_0 * (-p01 - p11) + alpha_1 * (-p00 - p10); a1x4 = alpha_0 * p01 + alpha_1 * p00; a1x5 = alpha_0 * p11 + alpha_1 * p10;

		a2x0 = beta_0 * (-p00 - p10) + beta_1 * (-p01 - p11);   a2x1 = beta_0 * p00 + beta_1 * p01;   a2x2 = beta_0 * p10 + beta_1 * p11;
		a2x3 = beta_0 * (p01 + p11) + beta_1 * (-p00 - p10);    a2x4 = -beta_0 * p01 + beta_1 * p00;  a2x5 = -beta_0 * p11 + beta_1 * p10;

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;
		aa = 0.25 / alpha_norm; bb = 0.25 / beta_norm;

		if (Plus_shrink_term)
		{
			if (inner_scaffold.count(i) > 0)
			{
				area_now *= shrink_weight;
				hi_1 = 2.0 + 6.0 * shrink_constant / (sig1*sig1*sig1*sig1);
				h_v = area_now * (2 * sig1 - 2 * shrink_constant / (sig1*sig1*sig1));
			}
			else
			{
				area_now *= free_scaffold_weight;
				hi_1 = 6.0 / (sig1*sig1*sig1*sig1);
				h_v = -area_now * 2.0 / (sig1*sig1*sig1);
			}
			h_u = 0.0;
			hi_0 = 0.0;
			uu = aa * aa*area_now * hi_1;
			vv = bb * bb*area_now * hi_1;
			uv = -aa * bb*area_now * hi_1;

			walpha = h_v;
			wbeta = -h_v;

			double hwa1 = (walpha * 0.25 / alpha_norm); double hwa2 = -(walpha * 0.25*0.25 / (alpha_norm*alpha_norm*alpha_norm));
			double hwb1 = (wbeta * 0.25 / beta_norm); double hwb2 = -(wbeta *0.25*0.25 / (beta_norm*beta_norm*beta_norm));


			h00 = uu * a1x0*a1x0 + vv * a2x0*a2x0 + uv * a1x0*a2x0 + uv * a2x0*a1x0; h01 = uu * a1x0*a1x1 + vv * a2x0*a2x1 + uv * a1x0*a2x1 + uv * a2x0*a1x1; h02 = uu * a1x0*a1x2 + vv * a2x0*a2x2 + uv * a1x0*a2x2 + uv * a2x0*a1x2; h03 = uu * a1x0*a1x3 + vv * a2x0*a2x3 + uv * a1x0*a2x3 + uv * a2x0*a1x3; h04 = uu * a1x0*a1x4 + vv * a2x0*a2x4 + uv * a1x0*a2x4 + uv * a2x0*a1x4; h05 = uu * a1x0*a1x5 + vv * a2x0*a2x5 + uv * a1x0*a2x5 + uv * a2x0*a1x5;

			h11 = uu * a1x1*a1x1 + vv * a2x1*a2x1 + uv * a1x1*a2x1 + uv * a2x1*a1x1; h12 = uu * a1x1*a1x2 + vv * a2x1*a2x2 + uv * a1x1*a2x2 + uv * a2x1*a1x2; h13 = uu * a1x1*a1x3 + vv * a2x1*a2x3 + uv * a1x1*a2x3 + uv * a2x1*a1x3; h14 = uu * a1x1*a1x4 + vv * a2x1*a2x4 + uv * a1x1*a2x4 + uv * a2x1*a1x4; h15 = uu * a1x1*a1x5 + vv * a2x1*a2x5 + uv * a1x1*a2x5 + uv * a2x1*a1x5;

			h22 = uu * a1x2*a1x2 + vv * a2x2*a2x2 + uv * a1x2*a2x2 + uv * a2x2*a1x2; h23 = uu * a1x2*a1x3 + vv * a2x2*a2x3 + uv * a1x2*a2x3 + uv * a2x2*a1x3; h24 = uu * a1x2*a1x4 + vv * a2x2*a2x4 + uv * a1x2*a2x4 + uv * a2x2*a1x4; h25 = uu * a1x2*a1x5 + vv * a2x2*a2x5 + uv * a1x2*a2x5 + uv * a2x2*a1x5;

			h33 = uu * a1x3*a1x3 + vv * a2x3*a2x3 + uv * a1x3*a2x3 + uv * a2x3*a1x3; h34 = uu * a1x3*a1x4 + vv * a2x3*a2x4 + uv * a1x3*a2x4 + uv * a2x3*a1x4; h35 = uu * a1x3*a1x5 + vv * a2x3*a2x5 + uv * a1x3*a2x5 + uv * a2x3*a1x5;

			h44 = uu * a1x4*a1x4 + vv * a2x4*a2x4 + uv * a1x4*a2x4 + uv * a2x4*a1x4; h45 = uu * a1x4*a1x5 + vv * a2x4*a2x5 + uv * a1x4*a2x5 + uv * a2x4*a1x5;

			h55 = uu * a1x5*a1x5 + vv * a2x5*a2x5 + uv * a1x5*a2x5 + uv * a2x5*a1x5;

			if (walpha >= 0)
			{
				h00 += hwa1 * (s1 + s2) + hwa2 * a1x0*a1x0; h01 += hwa1 * (-s1) + hwa2 * a1x0*a1x1; h02 += hwa1 * (-s2) + hwa2 * a1x0*a1x2; h03 += hwa2 * a1x0*a1x3; h04 += hwa1 * (h4)+hwa2 * a1x0*a1x4; h05 += hwa1 * (-h4) + hwa2 * a1x0*a1x5;
				h11 += hwa1 * (h1)+hwa2 * a1x1*a1x1;        h12 += hwa1 * (h2)+hwa2 * a1x1*a1x2;    h13 += hwa1 * (-h4) + hwa2 * a1x1*a1x3; h14 += hwa2 * a1x1*a1x4; h15 += hwa1 * (h4)+hwa2 * a1x1*a1x5;
				h22 += hwa1 * (h3)+hwa2 * a1x2*a1x2;        h23 += hwa1 * (h4)+hwa2 * a1x2*a1x3;    h24 += hwa1 * (-h4) + hwa2 * a1x2*a1x4; h25 += hwa2 * a1x2*a1x5;
				h33 += hwa1 * (s1 + s2) + hwa2 * a1x3*a1x3; h34 += hwa1 * (-s1) + hwa2 * a1x3*a1x4; h35 += hwa1 * (-s2) + hwa2 * a1x3*a1x5;
				h44 += hwa1 * (h1)+hwa2 * a1x4*a1x4;        h45 += hwa1 * (h2)+hwa2 * a1x4*a1x5;
				h55 += hwa1 * (h3)+hwa2 * a1x5*a1x5;
			}
			else
			{
				h00 += hwb1 * (s1 + s2) + hwb2 * a2x0*a2x0; h01 += hwb1 * (-s1) + hwb2 * a2x0*a2x1; h02 += hwb1 * (-s2) + hwb2 * a2x0*a2x2; h03 += hwb2 * a2x0*a2x3; h04 += hwb1 * (-h4) + hwb2 * a2x0*a2x4; h05 += hwb1 * (h4)+hwb2 * a2x0*a2x5;
				h11 += hwb1 * (h1)+hwb2 * a2x1*a2x1;        h12 += hwb1 * (h2)+hwb2 * a2x1*a2x2;    h13 += hwb1 * (h4)+hwb2 * a2x1*a2x3;    h14 += hwb2 * a2x1*a2x4; h15 += hwb1 * (-h4) + hwb2 * a2x1*a2x5;
				h22 += hwb1 * (h3)+hwb2 * a2x2*a2x2;        h23 += hwb1 * (-h4) + hwb2 * a2x2*a2x3; h24 += hwb1 * (h4)+hwb2 * a2x2*a2x4;    h25 += hwb2 * a2x2*a2x5;
				h33 += hwb1 * (s1 + s2) + hwb2 * a2x3*a2x3; h34 += hwb1 * (-s1) + hwb2 * a2x3*a2x4; h35 += hwb1 * (-s2) + hwb2 * a2x3*a2x5;
				h44 += hwb1 * (h1)+hwb2 * a2x4*a2x4;        h45 += hwb1 * (h2)+hwb2 * a2x4*a2x5;
				h55 += hwb1 * (h3)+hwb2 * a2x5*a2x5;
			}


		}
		else
		{
			hi_0 = 2.0 + 6.0 / (sig0*sig0*sig0*sig0); hi_1 = 2.0 + 6.0 / (sig1*sig1*sig1*sig1);
			h_u = area_now * (2.0 * sig0 - 2.0 / (sig0*sig0*sig0));
			h_v = area_now * (2.0 * sig1 - 2.0 / (sig1*sig1*sig1));

			uu = aa * aa*(area_now*hi_0 + area_now * hi_1);
			vv = bb * bb*(area_now*hi_0 + area_now * hi_1);
			uv = aa * bb*(area_now*hi_0 - area_now * hi_1);

			walpha = h_u + h_v;
			wbeta = h_u - h_v;

			double hwa1 = (walpha * 0.25 / alpha_norm); double hwa2 = -(walpha * 0.25*0.25 / (alpha_norm*alpha_norm*alpha_norm));
			double hwb1 = (wbeta * 0.25 / beta_norm); double hwb2 = -(wbeta *0.25*0.25 / (beta_norm*beta_norm*beta_norm));


			h00 = uu * a1x0*a1x0 + vv * a2x0*a2x0 + uv * a1x0*a2x0 + uv * a2x0*a1x0; h01 = uu * a1x0*a1x1 + vv * a2x0*a2x1 + uv * a1x0*a2x1 + uv * a2x0*a1x1; h02 = uu * a1x0*a1x2 + vv * a2x0*a2x2 + uv * a1x0*a2x2 + uv * a2x0*a1x2; h03 = uu * a1x0*a1x3 + vv * a2x0*a2x3 + uv * a1x0*a2x3 + uv * a2x0*a1x3; h04 = uu * a1x0*a1x4 + vv * a2x0*a2x4 + uv * a1x0*a2x4 + uv * a2x0*a1x4; h05 = uu * a1x0*a1x5 + vv * a2x0*a2x5 + uv * a1x0*a2x5 + uv * a2x0*a1x5;

			h11 = uu * a1x1*a1x1 + vv * a2x1*a2x1 + uv * a1x1*a2x1 + uv * a2x1*a1x1; h12 = uu * a1x1*a1x2 + vv * a2x1*a2x2 + uv * a1x1*a2x2 + uv * a2x1*a1x2; h13 = uu * a1x1*a1x3 + vv * a2x1*a2x3 + uv * a1x1*a2x3 + uv * a2x1*a1x3; h14 = uu * a1x1*a1x4 + vv * a2x1*a2x4 + uv * a1x1*a2x4 + uv * a2x1*a1x4; h15 = uu * a1x1*a1x5 + vv * a2x1*a2x5 + uv * a1x1*a2x5 + uv * a2x1*a1x5;

			h22 = uu * a1x2*a1x2 + vv * a2x2*a2x2 + uv * a1x2*a2x2 + uv * a2x2*a1x2; h23 = uu * a1x2*a1x3 + vv * a2x2*a2x3 + uv * a1x2*a2x3 + uv * a2x2*a1x3; h24 = uu * a1x2*a1x4 + vv * a2x2*a2x4 + uv * a1x2*a2x4 + uv * a2x2*a1x4; h25 = uu * a1x2*a1x5 + vv * a2x2*a2x5 + uv * a1x2*a2x5 + uv * a2x2*a1x5;

			h33 = uu * a1x3*a1x3 + vv * a2x3*a2x3 + uv * a1x3*a2x3 + uv * a2x3*a1x3; h34 = uu * a1x3*a1x4 + vv * a2x3*a2x4 + uv * a1x3*a2x4 + uv * a2x3*a1x4; h35 = uu * a1x3*a1x5 + vv * a2x3*a2x5 + uv * a1x3*a2x5 + uv * a2x3*a1x5;

			h44 = uu * a1x4*a1x4 + vv * a2x4*a2x4 + uv * a1x4*a2x4 + uv * a2x4*a1x4; h45 = uu * a1x4*a1x5 + vv * a2x4*a2x5 + uv * a1x4*a2x5 + uv * a2x4*a1x5;

			h55 = uu * a1x5*a1x5 + vv * a2x5*a2x5 + uv * a1x5*a2x5 + uv * a2x5*a1x5;

			if (walpha >= 0)
			{
				h00 += hwa1 * (s1 + s2) + hwa2 * a1x0*a1x0; h01 += hwa1 * (-s1) + hwa2 * a1x0*a1x1; h02 += hwa1 * (-s2) + hwa2 * a1x0*a1x2; h03 += hwa2 * a1x0*a1x3; h04 += hwa1 * (h4)+hwa2 * a1x0*a1x4; h05 += hwa1 * (-h4) + hwa2 * a1x0*a1x5;
				h11 += hwa1 * (h1)+hwa2 * a1x1*a1x1;        h12 += hwa1 * (h2)+hwa2 * a1x1*a1x2;    h13 += hwa1 * (-h4) + hwa2 * a1x1*a1x3; h14 += hwa2 * a1x1*a1x4; h15 += hwa1 * (h4)+hwa2 * a1x1*a1x5;
				h22 += hwa1 * (h3)+hwa2 * a1x2*a1x2;        h23 += hwa1 * (h4)+hwa2 * a1x2*a1x3;    h24 += hwa1 * (-h4) + hwa2 * a1x2*a1x4; h25 += hwa2 * a1x2*a1x5;
				h33 += hwa1 * (s1 + s2) + hwa2 * a1x3*a1x3; h34 += hwa1 * (-s1) + hwa2 * a1x3*a1x4; h35 += hwa1 * (-s2) + hwa2 * a1x3*a1x5;
				h44 += hwa1 * (h1)+hwa2 * a1x4*a1x4;        h45 += hwa1 * (h2)+hwa2 * a1x4*a1x5;
				h55 += hwa1 * (h3)+hwa2 * a1x5*a1x5;
			}
			h00 += hwb1 * (s1 + s2) + hwb2 * a2x0*a2x0; h01 += hwb1 * (-s1) + hwb2 * a2x0*a2x1; h02 += hwb1 * (-s2) + hwb2 * a2x0*a2x2; h03 += hwb2 * a2x0*a2x3; h04 += hwb1 * (-h4) + hwb2 * a2x0*a2x4; h05 += hwb1 * (h4)+hwb2 * a2x0*a2x5;
			h11 += hwb1 * (h1)+hwb2 * a2x1*a2x1;        h12 += hwb1 * (h2)+hwb2 * a2x1*a2x2;    h13 += hwb1 * (h4)+hwb2 * a2x1*a2x3;    h14 += hwb2 * a2x1*a2x4; h15 += hwb1 * (-h4) + hwb2 * a2x1*a2x5;
			h22 += hwb1 * (h3)+hwb2 * a2x2*a2x2;        h23 += hwb1 * (-h4) + hwb2 * a2x2*a2x3; h24 += hwb1 * (h4)+hwb2 * a2x2*a2x4;    h25 += hwb2 * a2x2*a2x5;
			h33 += hwb1 * (s1 + s2) + hwb2 * a2x3*a2x3; h34 += hwb1 * (-s1) + hwb2 * a2x3*a2x4; h35 += hwb1 * (-s2) + hwb2 * a2x3*a2x5;
			h44 += hwb1 * (h1)+hwb2 * a2x4*a2x4;        h45 += hwb1 * (h2)+hwb2 * a2x4*a2x5;
			h55 += hwb1 * (h3)+hwb2 * a2x5*a2x5;

		}


		u = aa * walpha; v = bb * wbeta;

		if (f0 != -1)
		{
			pardiso_b[f0] -= (u*a1x0 + v * a2x0);
			pardiso_b[f0 + V_N] -= (u*a1x3 + v * a2x3);
			pardiso_a[id_h00[i]] += h00;
			pardiso_a[id_h33[i]] += h33;
			pardiso_a[id_h03[i]] += h03;
		}
		if (f1 != -1)
		{
			pardiso_b[f1] -= (u*a1x1 + v * a2x1);
			pardiso_b[f1 + V_N] -= (u*a1x4 + v * a2x4);

			pardiso_a[id_h11[i]] += h11;
			pardiso_a[id_h44[i]] += h44;
			pardiso_a[id_h14[i]] += h14;
		}
		if (f2 != -1)
		{
			pardiso_b[f2] -= (u*a1x2 + v * a2x2);
			pardiso_b[f2 + V_N] -= (u*a1x5 + v * a2x5);

			pardiso_a[id_h22[i]] += h22;
			pardiso_a[id_h55[i]] += h55;
			pardiso_a[id_h25[i]] += h25;
		}

		if (flag_ == 3)
		{
			pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02;  pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
			pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h15[i]] += h15;
			pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; 
			pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
			pardiso_a[id_h45[i]] += h45;
		}
		else if(flag_==2)
		{
			if (f0 == -1)
			{
				pardiso_a[id_h12[i]] += h12;
				pardiso_a[id_h15[i]] += h15;
				pardiso_a[id_h24[i]] += h24; 
				pardiso_a[id_h45[i]] += h45;
			}
			else if (f1 == -1)
			{
				pardiso_a[id_h02[i]] += h02;
				pardiso_a[id_h05[i]] += h05;
				pardiso_a[id_h23[i]] += h23;
				pardiso_a[id_h35[i]] += h35;
			}
			else
			{
				pardiso_a[id_h01[i]] += h01;
				pardiso_a[id_h04[i]] += h04;
				pardiso_a[id_h13[i]] += h13;
				pardiso_a[id_h34[i]] += h34;
			}
		}		
	}

	//gradient_modified(position_of_mesh, pardiso_b);

	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;

	pardiso->factorize();
	pardiso->pardiso_solver();

	vector<double> result_d = pardiso->result;

	VectorXd negative_grad(2 * total_num), d(2 * total_num);
	negative_grad.setZero();
	d.setZero();

	for (int i = 0; i <V_N; i++)
	{
		negative_grad(var_ids[i]) = pardiso_b[i];
		negative_grad(var_ids[i]+total_num) = pardiso_b[i+V_N];
		d(var_ids[i]) = result_d[i];
		d(var_ids[i]+total_num) = result_d[i+V_N];
	}

	pardiso->free_numerical_factorization_memory();

	//descent_direction_modified(position_of_mesh, d);
	double temp_t;
	max_step(position_of_mesh, d, temp_t);

	double alpha = 0.95 * temp_t;

	backtracking_line_search(position_of_mesh, d, negative_grad, alpha);

	position_of_mesh += alpha * d;
	//std::cout << "cm step length: " << alpha << endl;
	Energysource(alpha);
}

void Parafun::max_step(const VectorXd &xx, const VectorXd &dd, double &step)
{
	double temp_t = numeric_limits<double>::infinity();
	int f0, f1, f2;
	double a, b, c, b1, b2, tt, tt1, tt2;
	double x0, x1, x2, x3, x4, x5, d0, d1, d2, d3, d4, d5;
	const double *x = xx.data();
	const double *d = dd.data();
	for (int i = 0; i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = x[f0]; x1 = x[f1]; x2 = x[f2]; x3 = x[f0 + total_num]; x4 = x[f1 + total_num]; x5 = x[f2 + total_num];
		d0 = d[f0]; d1 = d[f1]; d2 = d[f2]; d3 = d[f0 + total_num]; d4 = d[f1 + total_num]; d5 = d[f2 + total_num];

		a = (d1 - d0) * (d5 - d3) - (d4 - d3) * (d2 - d0);
		b1 = (d1 - d0) * (x5 - x3) + (x1 - x0) * (d5 - d3);
		b2 = (x4 - x3) * (d2 - d0) + (x2 - x0) * (d4 - d3);
		b = b1 - b2;
		c = (x1 - x0) * (x5 - x3) - (x4 - x3) * (x2 - x0);
		tt = get_smallest_pos_quad_zero( a, b, c);

		if (temp_t > tt)
		{
			temp_t = tt;
		}

	}
	if (temp_t == INFINITY)
		temp_t = 100;
	step = temp_t;
}

void Parafun::descent_direction_modified(const VectorXd & pos, VectorXd & d)
{
	printf("descent_direction_modified---------------\n");
	double band = 1e-6;
	vector<int> restrict_top;
	vector<int> restrict_bottom;
	vector<int> restrict_left;
	vector<int> restrict_right;

	//top,bottom,left,right
	double top_line = bbox(1, 1) - band;
	double bottom_line = bbox(1, 0) + band;
	double left_line = bbox(0, 0) + band;
	double right_line = bbox(0, 1) - band;

	double x, y;
	for (auto&loop_b : bnds)
	{
		for (auto&var : loop_b)
		{
			x = pos(var_ids[var]);
			y = pos(var_ids[var] + total_num);
			if (y > top_line)
				restrict_top.push_back(var);
			if (y < bottom_line)
				restrict_bottom.push_back(var);
			if (x < left_line)
				restrict_left.push_back(var);
			if (x > right_line)
				restrict_right.push_back(var);
		}
	}
	printf("top size: %zd; bottom size: %zd; left size: %zd; right size: %zd;\n", restrict_top.size(),restrict_bottom.size(),restrict_left.size(),restrict_right.size());
	double factor = 1e-3;
	for (auto&var : restrict_top)
	{
		double &cur_coord = d(var_ids[var] + total_num);
		if (cur_coord > 0)
			cur_coord *= factor;
	}
	for (auto&var : restrict_bottom)
	{
		double &cur_coord = d(var_ids[var] + total_num);
		if (cur_coord < 0)
			cur_coord *= factor;
	}
	for (auto&var : restrict_left)
	{
		double &cur_coord = d(var_ids[var]);
		if (cur_coord < 0)
			cur_coord *= factor;
	}
	for (auto&var : restrict_right)
	{
		double &cur_coord = d(var_ids[var]);
		if (cur_coord > 0)
			cur_coord *= factor;
	}
}

void Parafun::gradient_modified(const VectorXd & pos, vector<double>& b)
{
	printf("gradient_modified---------------\n");
	double band = 1e-6;
	vector<int> restrict_top;
	vector<int> restrict_bottom;
	vector<int> restrict_left;
	vector<int> restrict_right;

	//top,bottom,left,right
	double top_line = bbox(1, 1) - band;
	double bottom_line = bbox(1, 0) + band;
	double left_line = bbox(0, 0) + band;
	double right_line = bbox(0, 1) - band;

	double x, y;
	for (auto&loop_b : bnds)
	{
		for (auto&var : loop_b)
		{
			x = pos(var_ids[var]);
			y = pos(var_ids[var] + total_num);
			if (y > top_line)
				restrict_top.push_back(var);
			if (y < bottom_line)
				restrict_bottom.push_back(var);
			if (x < left_line)
				restrict_left.push_back(var);
			if (x > right_line)
				restrict_right.push_back(var);
		}
	}
	printf("top size: %zd; bottom size: %zd; left size: %zd; right size: %zd;\n", restrict_top.size(), restrict_bottom.size(), restrict_left.size(), restrict_right.size());
	double factor = 1e-3;
	for (auto&var : restrict_top)
	{
		double &cur_coord = b[var + V_N];
		if (cur_coord > 0)
			cur_coord *= factor;
	}
	for (auto&var : restrict_bottom)
	{
		double &cur_coord = b[var + V_N];
		if (cur_coord < 0)
			cur_coord *= factor;
	}
	for (auto&var : restrict_left)
	{
		double &cur_coord = b[var];
		if (cur_coord < 0)
			cur_coord *= factor;
	}
	for (auto&var : restrict_right)
	{
		double &cur_coord = b[var];
		if (cur_coord > 0)
			cur_coord *= factor;
	}
}


void Parafun::calc_area_H()
{

}

double Parafun::calc_efficiency_area(const VectorXd & x)
{
	double area_sum = 0;
	for (auto&loop_b : bnds)
	{
		double area_part = 0;
		double area_;
		Eigen::Vector2d v0;
		Eigen::Vector2d v1;
		for (size_t i = 1; i < loop_b.size(); i++)
		{
			v0 << x(var_ids[loop_b[i - 1]]), x(var_ids[loop_b[i - 1]] + total_num);
			v1 << x(var_ids[loop_b[i]]), x(var_ids[loop_b[i]] + total_num);
			area_ = v0(0)*v1(1) - v0(1)*v1(0);
			area_part += area_;
		}
		v0 << x(var_ids[loop_b.back()]), x(var_ids[loop_b.back()] + total_num);
		v1 << x(var_ids[loop_b[0]]), x(var_ids[loop_b[0]] + total_num);
		area_ = v0(0)*v1(1) - v0(1)*v1(0);
		area_part += area_;

		area_sum += area_part;
	}
	area_sum /= 2.0;
	return area_sum;
}

double Parafun::calc_Curgap(const VectorXd &x)
{
	area_bbox = calc_cur_bbox(x);
	double area_sum = calc_efficiency_area(x);
	double gap = area_sum - area_bbox * packing_efficiency;
	return gap;
}

void Parafun::calc_areaTerm_hessian_gradient()
{
	Cur_gap = calc_Curgap(position_of_mesh);
	double gap_tmp = -0.5*area_term_weight / Cur_gap;
	const double *x = position_of_mesh.data();

	vector<vector<double>> g_local;
	g_local.resize(bnds.size());

	//calc gradient
	for (size_t j = 0; j < bnds.size(); j++)
//	for (auto&loop_b : bnds)
	{
		auto& loop_b = bnds[j];
		int loop_b_size = loop_b.size();
		g_local[j].resize(2 * loop_b_size);
		double tmp1, tmp2;
		for (size_t i = 0; i < loop_b.size(); i++)
		{
			int v_pre = loop_b[(i - 1 + loop_b_size) % loop_b_size];
			int v_cur = loop_b[i];
			int v_nex = loop_b[(i + 1) % loop_b_size];

			tmp1 = gap_tmp * (x[var_ids[v_nex] + total_num] - x[var_ids[v_pre] + total_num]);
			tmp2 = gap_tmp * (x[var_ids[v_pre]] - x[var_ids[v_nex]]);
			pardiso_b[v_cur] -= tmp1;
			pardiso_b[v_cur + V_N] -= tmp2;

			g_local[j][i] = tmp1;
			g_local[j][i+ loop_b_size] = tmp2;

		}
	}

	//calc hessian

	for (size_t i = 0; i < bnds.size(); i++)
	{
		int bnds_i_size2 = g_local[i].size();
		int index = 0;
		for (size_t j = 0; j < bnds_i_size2; j++)
		{
			for (size_t k = j; k < bnds_i_size2; k++)
			{
				pardiso_a[id_bnd_h[i][index]] += g_local[i][j] * g_local[i][k] / area_term_weight;
				index++;
			}
		}
	}



}

void Parafun::areaCondition(const VectorXd & xx, const VectorXd & dd, Vector3d & quad_equ)
{
	quad_equ.setZero();
	double x0, x1, x2, x3, d0, d1, d2, d3;
	for (auto&loop_b : bnds)
	{
		for (size_t i = 1; i < loop_b.size(); i++)
		{
			x0 = xx(var_ids[loop_b[i - 1]]);
			x2 = xx(var_ids[loop_b[i - 1]] + total_num);
			x1 = xx(var_ids[loop_b[i]]);
			x3=xx(var_ids[loop_b[i]] + total_num);

			d0 = dd(var_ids[loop_b[i - 1]]);
			d2 = dd(var_ids[loop_b[i - 1]] + total_num);
			d1 = dd(var_ids[loop_b[i]]);
			d3 = dd(var_ids[loop_b[i]] + total_num);

			quad_equ(0) += (d0 * d3 - d1 * d2);
			quad_equ(1) += (x0*d3 + x3 * d0 - x2 * d1 - x1 * d2);
			quad_equ(2) += x0 * x3 - x2 * x1;
		}

		x0 = xx(var_ids[loop_b.back()]);
		x2 = xx(var_ids[loop_b.back()] + total_num);
		x1 = xx(var_ids[loop_b[0]]);
		x3 = xx(var_ids[loop_b[0]] + total_num);

		d0 = dd(var_ids[loop_b.back()]);
		d2 = dd(var_ids[loop_b.back()] + total_num);
		d1 = dd(var_ids[loop_b[0]]);
		d3 = dd(var_ids[loop_b[0]] + total_num);

		quad_equ(0) += (d0 * d3 - d1 * d2);
		quad_equ(1) += (x0*d3 + x3 * d0 - x2 * d1 - x1 * d2);
		quad_equ(2) += x0 * x3 - x2 * x1;
	}
	quad_equ /= 2.0;
	quad_equation = quad_equ;
}

double Parafun::calc_cur_bbox(const VectorXd & x)
{
	double xmin = numeric_limits<double>::infinity();
	double ymin = numeric_limits<double>::infinity();
	double xmax = -numeric_limits<double>::infinity();
	double ymax = -numeric_limits<double>::infinity();

	//int mv_n = d_.mv_num;

	for (size_t i = 0; i < d_.m_T.rows(); i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			int vid = d_.m_T(i, j);
			if (xmin > x(vid))
				xmin = x(vid);
			if (xmax < x(vid))
				xmax = x(vid);
			if (ymin > x(vid + total_num))
				ymin = x(vid + total_num);
			if (ymax < x(vid + total_num))
				ymax = x(vid + total_num);

		}
	}



	bbox(0, 0) = xmin;
	bbox(0, 1) = xmax;
	bbox(1, 0) = ymin;
	bbox(1, 1) = ymax;


	double area_bb = (xmax - xmin)*(ymax - ymin);
	return area_bb;
}

void Parafun::separate_scaffold()
{
	inner_scaffold.clear();
	bool isinner = true;
	double x, y;
	int vid;
	int mf_n = d_.m_T.rows();
	for (size_t i = mf_n; i < F_N; i++)
	{
		isinner = true;
		for (size_t j = 0; j < 3; j++)
		{
			vid = d_.s_T(i- mf_n, j);
			x = position_of_mesh(vid);
			y = position_of_mesh(vid + total_num);
			if (x > bbox(0, 1) || x<bbox(0, 0) || y>bbox(1, 1) || y < bbox(1, 0))
			{
				isinner = false;
				break;
			}
		}
		//isinner = true;
		if (isinner)
		{
			inner_scaffold.insert(i);
		}
	}
}

double get_smallest_pos_quad_zero(double a, double b, double c)
{
	using namespace std;
	double t1, t2;
	if (std::abs(a) <= 1.0e-10)
	{
		a *= 1e6;
		b *= 1e6;
		c *= 1e6;
	}
	if (std::abs(a) > 1.0e-10)
	{
		double delta_in = pow(b, 2) - 4 * a * c;
		if (delta_in <= 0)
		{
			return INFINITY;
		}

		double delta = sqrt(delta_in); // delta >= 0
		if (b >= 0) // avoid subtracting two similar numbers
		{
			double bd = -b - delta;
			t1 = 2 * c / bd;
			t2 = bd / (2 * a);
		}
		else
		{
			double bd = -b + delta;
			t1 = bd / (2 * a);
			t2 = (2 * c) / bd;
		}

		assert(std::isfinite(t1));
		assert(std::isfinite(t2));

		if (a < 0) std::swap(t1, t2); // make t1 > t2
		// return the smaller positive root if it exists, otherwise return infinity
		if (t1 > 0)
		{
			return t2 > 0 ? t2 : t1;
		}
		else
		{
			return INFINITY;
		}
	}
	else
	{
		if (b == 0) return INFINITY; // just to avoid divide-by-zero
		t1 = -c / b;
		return t1 > 0 ? t1 : INFINITY;
	}
}


void Parafun::backtracking_line_search(const VectorXd &x, const VectorXd &d, const VectorXd &negetive_grad, double &alpha)
{
	double h = 0.5;
	double tt = -(negetive_grad.transpose()*d)(0, 0);
	double c = 0.2; 
	double ex;
	double e;


	if (Plus_area_term_search)
	{
		Vector3d quad_equ;
		areaCondition(x, d, quad_equ);
		auto judgeArea = [&](double t,double area_bb)->bool {return (quad_equ(0)*t*t + quad_equ(1)*t + quad_equ(2)-area_bb*packing_efficiency) > 0; };

		double cur_bb = calc_cur_bbox(x);
		Vector4d quad_(quad_equ(0), quad_equ(1), quad_equ(2), cur_bb);
		Energy(x, ex, 0.0, quad_);
		VectorXd x_new = x + alpha * d;
		bool area_isok = false;
		cur_bb = calc_cur_bbox(x_new);
		if (judgeArea(alpha, cur_bb))
			area_isok = true;
		if (area_isok)
		{
			quad_(3)=cur_bb;
			Energy(x_new, e, alpha, quad_);
		}

		while (!area_isok||e > ex + alpha * c * tt)
		{
			alpha = h * alpha;
			x_new = x + alpha * d;
			cur_bb = calc_cur_bbox(x_new);
			area_isok = judgeArea(alpha, cur_bb);
			if (area_isok)
			{
				quad_(3) = cur_bb;
				Energy(x_new, e, alpha, quad_);
			}
		}

	}
	else
	{
		Energy(x, ex);
		VectorXd x_new = x + alpha * d;
		Energy(x_new, e);
		while (e > ex + alpha * c * tt)
		{
			alpha = h * alpha;
			x_new = x + alpha * d;
			Energy(x_new, e);
		}
	}


}


void Parafun::Energy(const VectorXd &position, double &energyupdate, double step, const Vector4d& quad)
{
	double energy = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_d;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;
	const double *pos = position.data();
	for (int i = 0; i < d_.m_T.rows(); ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		//pp
		//p00 = update_p00[i]; p01 = update_p01[i]; p10 = update_p10[i]; p11 = update_p11[i];
		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;


		det = j00*j11 - j01*j10;
		E_d = (1 + 1 / (det*det)) * (j00*j00 + j01*j01 + j10*j10 + j11*j11);

		energy += area[i] * E_d;
	}
	double alpha_0, alpha_1, beta_0, beta_1;

	double alpha_norm, beta_norm,sigma2_2;
	for (int i = d_.m_T.rows(); i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		//pp
		//p00 = update_p00[i]; p01 = update_p01[i]; p10 = update_p10[i]; p11 = update_p11[i];
		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;


		det = j00 * j11 - j01 * j10;
		if (Plus_shrink_term)
		{
			alpha_0 = j00 + j11; alpha_1 = j10 - j01;
			beta_0 = j00 - j11;  beta_1 = j10 + j01;
			alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
			beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1 * beta_1);
			sigma2_2 = (alpha_norm - beta_norm)*(alpha_norm - beta_norm);

			if (inner_scaffold.count(i) > 0)
			{
				E_d = shrink_weight * (sigma2_2 + shrink_constant / sigma2_2);
			}
			else
			{
				E_d = free_scaffold_weight / sigma2_2;
			}

		}
		else
		{
			E_d = (1.0 + 1.0 / (det*det)) * (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);

		}
		energy += area[i] * E_d;
	}

	if (Plus_area_term)
	{
		double gap = quad(0)*step*step + quad(1)*step + quad(2) - quad(3)*packing_efficiency;
		energy -= area_term_weight*std::log(gap);
	}

	energyupdate = energy;
}

void Parafun::Energysource(double step)
{
	double end_e_one_temp = 0, end_e_area = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	const double *pos = position_of_mesh.data();
	for (int i = 0; i < d_.m_T.rows(); ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;

		det = j00*j11 - j01*j10;

		E_1 = (j00*j00 + j01*j01 + j10*j10 + j11*j11);
		E_2 = 1.0 / (det*det)* E_1;

		end_e_one_temp += E_1;
		end_e_one_temp += E_2;
		end_e_area += ((E_1 + E_2)*area_src[i]);
	}
	energy_uniform = end_e_one_temp /F_N;

	//printf("energy without AreaTerm: %f;\n", end_e_area);
	if (Plus_area_term)
	{
		area_bbox = calc_cur_bbox(position_of_mesh);
		double gap = quad_equation(0)*step*step + quad_equation(1)*step + quad_equation(2)-area_bbox*packing_efficiency;
		double e_area_term = -area_term_weight * std::log(gap);
		end_e_area += e_area_term;
		printf("e_area_term: %.16f;gap: %.16f;", e_area_term,gap);
	}

	energy_area = end_e_area;
}

double Parafun::compute_energy(const Eigen::MatrixXd & x, bool whole)
{
	double end_e_one_temp = 0, end_e_area = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	const double *pos = x.data();
	int src_t_num = d_.m_T.rows();

	for (int i = 0; i < src_t_num; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		det = j00 * j11 - j01 * j10;

		E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
		E_2 = E_1 / (det*det);

		end_e_one_temp += E_1;
		end_e_one_temp += E_2;
		end_e_area += ((E_1 + E_2)*area_src[i]);
	}
	//cout << "compute energy " << end_e_area << endl;
	if (whole)
	{
		for (int i = src_t_num; i < F_N; ++i)
		{
			f0 = F0[i];
			f1 = F1[i];
			f2 = F2[i];

			x0 = pos[f0];
			y0 = pos[f0 + total_num];

			x1 = pos[f1];
			y1 = pos[f1 + total_num];

			x2 = pos[f2];
			y2 = pos[f2 + total_num];

			q00 = x1 - x0; q01 = x2 - x0;
			q10 = y1 - y0; q11 = y2 - y0;

			p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

			j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

			det = j00 * j11 - j01 * j10;

			E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
			E_2 = 1.0 / (det*det)* E_1;

			end_e_one_temp += E_1;
			end_e_one_temp += E_2;
			end_e_area += ((E_1 + E_2)*area_scaf[i- src_t_num]);
		}
	}


	return end_e_area;
}

void Parafun::local_coordinate_inverse(int i, double &p00, double &p01, double &p10, double &p11)
{
	int f0 = F0[i];
	int f1 = F1[i];
	int f2 = F2[i];
	
	Vector3d x_(d_.m_V(f1,0)- d_.m_V(f0,0), d_.m_V(f1, 1) - d_.m_V(f0, 1), d_.m_V(f1, 2) - d_.m_V(f0, 2));
	double x1_0 = x_.norm();
	x_/=x1_0;
	Vector3d l_(d_.m_V(f2, 0) - d_.m_V(f0, 0), d_.m_V(f2, 1) - d_.m_V(f0, 1), d_.m_V(f2, 2) - d_.m_V(f0, 2));

	Vector3d n_ = x_.cross(l_);
	n_.normalize();
	Vector3d y_ = n_.cross(x_);
	double x2_0 = l_.dot(x_);
	double y2_0 = l_.dot(y_);

	p00 = 1 / x1_0;
	p01 = -x2_0 / (x1_0*y2_0);
	p10 = 0;
	p11 = 1 / y2_0;
}

void Parafun::local_coordinate_inverse_scaf(int i, double & p00, double & p01, double & p10, double & p11)
{
	int f0 = F0[i];
	int f1 = F1[i];
	int f2 = F2[i];

	Vector2d x_(d_.w_uv(f1, 0) - d_.w_uv(f0, 0), d_.w_uv(f1, 1) - d_.w_uv(f0, 1));
	Vector2d l_(d_.w_uv(f2, 0) - d_.w_uv(f0, 0), d_.w_uv(f2, 1) - d_.w_uv(f0, 1));

	double area_tri = abs(x_(0)*l_(1) - x_(1)*l_(0));
	double x1_0, x2_0, y2_0;
	if (area_tri > area_threshold)
	{
		x1_0 = x_.norm();
		x_ /= x1_0;
		Vector2d y_(-x_(1), x_(0));
		x2_0 = l_.dot(x_);
		y2_0 = l_.dot(y_);
	}
	else
	{
		//cout << "area too small!!!!!!!!!!!!! " << endl;
		double h = sqrt((2*area_threshold) / sqrt(3.0));
		x1_0 = h;
		x2_0 = h / 2.0;
		y2_0 = sqrt(3.0)*h / 2.0;
	}
	p00 = 1 / x1_0;
	p01 = -x2_0 / (x1_0*y2_0);
	p10 = 0;
	p11 = 1 / y2_0;
}

void Parafun::local_coordinate_inverse_gap(int i, double &p00, double &p01, double &p10, double &p11)
{
	int f0 = d_.gz.g_F(i,0);
	int f1 = d_.gz.g_F(i, 1);
	int f2 = d_.gz.g_F(i, 2);
	auto &g_V = d_.gz.g_V;

	Vector3d x_(g_V(f1, 0) - g_V(f0, 0), g_V(f1, 1) - g_V(f0, 1), g_V(f1, 2) - g_V(f0, 2));
	double x1_0 = x_.norm();
	x_ /= x1_0;
	Vector3d l_(g_V(f2, 0) - g_V(f0, 0), g_V(f2, 1) - g_V(f0, 1), g_V(f2, 2) - g_V(f0, 2));

	Vector3d n_ = x_.cross(l_);
	n_.normalize();
	Vector3d y_ = n_.cross(x_);
	double x2_0 = l_.dot(x_);
	double y2_0 = l_.dot(y_);

	p00 = 1 / x1_0;
	p01 = -x2_0 / (x1_0*y2_0);
	p10 = 0;
	p11 = 1 / y2_0;
}

double Parafun::newton_equation(const double & a, const double & b, const double & K)
{
	double tt = 1;
	double E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	while (abs(E_d) > 1e-5)
	{
		tt = tt - 1 / (2 * log(a)*pow(a, 2 * tt) + 2 * log(b)* pow(b, 2 * tt) + 2 * log(1 / a)* pow(1 / a, 2 * tt) + 2 * log(1 / b)*pow(1 / b, 2 * tt))*(pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K);
		E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	}
	return tt;
}

void Parafun::adjust_scaf_weight(double new_weight)
{
	d_.scaffold_factor = new_weight;
	d_.update_scaffold();

	inter_weight_1=new_weight*100.0*d_.sf_num;

	//init_area();
}

void Parafun::handle_mintri()
{
	double min_bnd_edge_len = numeric_limits<double>::max();
	int acc_bnd = 0;
	for (int i = 0; i < d_.bnd_sizes.size(); i++)
	{
		int current_size = d_.bnd_sizes[i];

		for (int e = acc_bnd; e < acc_bnd + current_size - 1; e++)
		{
			min_bnd_edge_len = (std::min)(min_bnd_edge_len,
				(d_.w_uv.row(d_.internal_bnd(e)) -
					d_.w_uv.row(d_.internal_bnd(e + 1)))
				.squaredNorm());
		}
		min_bnd_edge_len = (std::min)(min_bnd_edge_len,
			(d_.w_uv.row(d_.internal_bnd(acc_bnd)) -
				d_.w_uv.row(d_.internal_bnd(acc_bnd + current_size -
					1)))
			.squaredNorm());
		acc_bnd += current_size;
	}

	area_threshold = min_bnd_edge_len / 4.0;
}


void Parafun::out_dis(const string & file_str)
{
	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	const double *pos = position_of_mesh.data();
	vector<double> dis_out;
	dis_out.resize(d_.m_T.rows());
	for (int i = 0; i < d_.m_T.rows(); ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		det = j00 * j11 - j01 * j10;

		E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
		E_2 = E_1 / (det*det);

		dis_out[i] = (E_1 + E_2)*0.25;
	}
	ofstream ofs_(file_str, ios::trunc);
	for (size_t i = 0; i < dis_out.size(); i++)
	{
		ofs_ << dis_out[i] << endl;
	}
	ofs_.close();
}
