#include "MyParafun.h"
#include "OMPHelper.h"

#define min_value 1e-8
//Parafun::Parafun(string srcname, string initname, int dim) :n_constraints(0)
//{
//	OpenMesh::IO::read_mesh(source_mesh, srcname);
//	OpenMesh::IO::read_mesh(init_mesh, initname);
//	source_dim = dim;
//	lb = -1;
//	ub = -1;
//	iter_max = 1000;
//	tol_err = 1e-4;
//	tol_energy = 1e-10;
//	use_weighted_metric = false;
//	para_1 = 2.0 / 3.0;
//	para_2 = 2;
//
//	is_optimal = false;
//	iter_equ = 0;
//	iter_project = 0;
//}

MyParafun::MyParafun(Mesh &mesh_src_, Mesh &mesh_tutte_, int dim) :n_constraints(0)
{
	source_mesh = mesh_src_;
	init_mesh = mesh_tutte_;
	source_dim = dim;
	lb = -1;
	ub = -1;
	iter_max = 1000;
	tol_err = 1e-4;
	tol_energy = 1e-10;
	use_weighted_metric = false;
	para_1 = 2.0 / 3.0;
	para_2 = 2;

	is_optimal = false;
	iter_equ = 0;
	iter_project = 0;
}

MyParafun::~MyParafun()
{
	if (pardiso != NULL)
	{
		delete pardiso;
		pardiso = NULL;
	}
}

bool MyParafun::writeobj(Mesh& out_mesh, const string out_filename)
{
	for (auto it1 = out_mesh.vertices_begin(); it1 != out_mesh.vertices_end(); ++it1)
	{
		int id = it1->idx();
		OpenMesh::Vec3d pos(cur_position(id), cur_position(id + V_N), 0);
		out_mesh.set_point(*it1, pos);
	}
	
	return 1;
	//return OpenMesh::IO::write_mesh(out_mesh, out_filename);
}

void MyParafun::init()
{
	V_N = source_mesh.n_vertices();
	F_N = source_mesh.n_faces();

	F.resize(F_N, 3);
	Mesh::FaceIter f_it = source_mesh.faces_begin();
	Mesh::FaceIter f_end = source_mesh.faces_end();
	Mesh::FaceVertexIter fv_it;

	source_position.resize(V_N * source_dim);
	Mesh::VertexIter v_it = source_mesh.vertices_begin();
	Mesh::VertexIter v_end = source_mesh.vertices_end();
	Mesh::Point pos;

	global_position.resize(V_N * 2);
	Mesh::VertexIter iv_it = init_mesh.vertices_begin();
	Mesh::Point ipos;

	Tx.resize(F_N * 2 * 2);
	pTx.resize(F_N * 2 * 2);
	distortions.resize(F_N);
	flips.resize(F_N);
	minsv.resize(F_N);
	maxsv.resize(F_N);
	K.resize(F_N);

	n_constraints = 0;
	int i = 0, k = 0;

	//ÃÓ≥‰Tæÿ’Û
	int T_rows = F_N * 2 * 2;
	int T_cols = V_N * 2;
	T.resize(T_cols, T_rows);
	vector<Tri> triplist_T;
	triplist_T.reserve(T_rows * 3);
	vector<Tri> triplist_A;
	triplist_A.reserve(2 * T_cols);

	MatrixXd currV(3, source_dim), currVFlat(3, 2), currT(source_dim, 3);
	MatrixXd B = MatrixXd::Identity(3, 3);
	B = B.array() - (1.0 / 3);
	Mesh::HalfedgeHandle he_h;
	double x31, x23, x12, y32, y13, y21, area;
	int curr_row = 0;

	if (3 == source_dim)
	{
		OMP_PARALLEL
		{
			OMP_FOR
			for (; f_it != f_end; ++f_it, ++i)
			{
				fv_it = source_mesh.fv_begin(*f_it);
				for (int j = 0; j < 3; ++j, ++fv_it)
				{
					F(i, j) = fv_it.handle().idx();
				}
			}

		OMP_FOR
			for (; v_it != v_end; ++v_it, ++k, ++iv_it)
			{
				pos = source_mesh.point(*v_it);
				ipos = init_mesh.point(*iv_it);
				for (int j = 0; j < source_dim; ++j)
				{
					source_position[k + V_N*j] = pos[j];
				}
				for (int j = 0; j < 2; ++j)
				{
					global_position[k + V_N*j] = ipos[j];
				}
			}

		OMP_FOR
			for (auto it1 = init_mesh.vertices_begin(); it1 != init_mesh.vertices_end(); ++it1)
			{
				if (init_mesh.is_boundary(*it1))
				{
					triplist_A.push_back(Tri(2 * n_constraints, it1->idx(), 1));
					triplist_A.push_back(Tri(2 * n_constraints + 1, it1->idx() + V_N, 1));
					triplist_A.push_back(Tri(2 * n_constraints, T_cols + 2 * n_constraints, 0.0));
					triplist_A.push_back(Tri(2 * n_constraints + 1, T_cols + 2 * n_constraints + 1, 0.0));
					n_constraints++;
				}
			}

		OMP_FOR
			for (int ii = 0; ii < F_N; ii++)
			{
				for (int jj = 0; jj < 3; jj++)
				{
					currV(jj, 0) = source_position(F(ii, jj));
					currV(jj, 1) = source_position(F(ii, jj) + V_N);
					currV(jj, 2) = source_position(F(ii, jj) + 2 * V_N);
				}
				embedTriangle(currV, currVFlat);
				currVFlat = B*currVFlat; // center
				currT = currVFlat.fullPivLu().solve(B); // solver

				triplist_T.push_back(Tri(F(ii, 0), curr_row, currT(0, 0)));
				triplist_T.push_back(Tri(F(ii, 1), curr_row, currT(0, 1)));
				triplist_T.push_back(Tri(F(ii, 2), curr_row, currT(0, 2))); curr_row = curr_row + 1;
				triplist_T.push_back(Tri(F(ii, 0) + V_N, curr_row, currT(0, 0)));
				triplist_T.push_back(Tri(F(ii, 1) + V_N, curr_row, currT(0, 1)));
				triplist_T.push_back(Tri(F(ii, 2) + V_N, curr_row, currT(0, 2))); curr_row = curr_row + 1;

				triplist_T.push_back(Tri(F(ii, 0), curr_row, currT(1, 0)));
				triplist_T.push_back(Tri(F(ii, 1), curr_row, currT(1, 1)));
				triplist_T.push_back(Tri(F(ii, 2), curr_row, currT(1, 2))); curr_row = curr_row + 1;
				triplist_T.push_back(Tri(F(ii, 0) + V_N, curr_row, currT(1, 0)));
				triplist_T.push_back(Tri(F(ii, 1) + V_N, curr_row, currT(1, 1)));
				triplist_T.push_back(Tri(F(ii, 2) + V_N, curr_row, currT(1, 2))); curr_row = curr_row + 1;
			}
		}
	}
	else
	{
		OMP_PARALLEL
		{
			OMP_FOR
			for (; f_it != f_end; ++f_it, ++i)
			{
				fv_it = source_mesh.fv_begin(*f_it);
				for (int j = 0; j < 3; ++j, ++fv_it)
				{
					F(i, j) = fv_it.handle().idx();
				}
			}

		OMP_FOR
			for (; v_it != v_end; ++v_it, ++k, ++iv_it)
			{
				pos = source_mesh.point(*v_it);
				ipos = init_mesh.point(*iv_it);
				for (int j = 0; j < source_dim; ++j)
				{
					source_position[k + V_N*j] = pos[j];
				}
				for (int j = 0; j < 2; ++j)
				{
					global_position[k + V_N*j] = ipos[j];
				}
			}

		OMP_FOR
			for (auto it1 = init_mesh.vertices_begin(); it1 != init_mesh.vertices_end(); ++it1)
			{
				if (init_mesh.is_boundary(*it1))
				{
					triplist_A.push_back(Tri(2 * n_constraints, it1->idx(), 1));
					triplist_A.push_back(Tri(2 * n_constraints + 1, it1->idx() + V_N, 1));
					triplist_A.push_back(Tri(2 * n_constraints, T_cols + 2 * n_constraints, 0.0));
					triplist_A.push_back(Tri(2 * n_constraints + 1, T_cols + 2 * n_constraints + 1, 0.0));
					n_constraints++;
				}
			}

		OMP_FOR
			for (int ii = 0; ii < F_N; ii++)
			{
				for (int jj = 0; jj < 3; jj++)
				{
					currV(jj, 0) = source_position(F(ii, jj));
					currV(jj, 1) = source_position(F(ii, jj) + V_N);
				}
				currV = B*currV;
				currT = currV.fullPivLu().solve(B); // solver

				triplist_T.push_back(Tri(F(ii, 0), curr_row, currT(0, 0)));
				triplist_T.push_back(Tri(F(ii, 1), curr_row, currT(0, 1)));
				triplist_T.push_back(Tri(F(ii, 2), curr_row, currT(0, 2))); curr_row = curr_row + 1;
				triplist_T.push_back(Tri(F(ii, 0) + V_N, curr_row, currT(0, 0)));
				triplist_T.push_back(Tri(F(ii, 1) + V_N, curr_row, currT(0, 1)));
				triplist_T.push_back(Tri(F(ii, 2) + V_N, curr_row, currT(0, 2))); curr_row = curr_row + 1;

				triplist_T.push_back(Tri(F(ii, 0), curr_row, currT(1, 0)));
				triplist_T.push_back(Tri(F(ii, 1), curr_row, currT(1, 1)));
				triplist_T.push_back(Tri(F(ii, 2), curr_row, currT(1, 2))); curr_row = curr_row + 1;
				triplist_T.push_back(Tri(F(ii, 0) + V_N, curr_row, currT(1, 0)));
				triplist_T.push_back(Tri(F(ii, 1) + V_N, curr_row, currT(1, 1)));
				triplist_T.push_back(Tri(F(ii, 2) + V_N, curr_row, currT(1, 2))); curr_row = curr_row + 1;
			}
		}
	}

	SparseMatrix<double> A_extension(2 * n_constraints, T_cols + 2 * n_constraints);
	A_extension.setFromTriplets(triplist_A.begin(), triplist_A.end());
	b.resize(2 * n_constraints);
	for (int i = 0; i < n_constraints; ++i)
	{
		b(2 * i) = global_position(triplist_A[4 * i].col());
		b(2 * i + 1) = global_position(triplist_A[4 * i + 1].col());
	}

	T.setFromTriplets(triplist_T.begin(), triplist_T.end());
	T = T.transpose();
	SparseMatrix<double> Coef;
	SparseMatrix<double> T_extension(T_cols + 2 * n_constraints, T_rows);
	T_extension.setFromTriplets(triplist_T.begin(), triplist_T.end());
	Coef = T_extension * T_extension.transpose();
	Coef.rightCols(2 * n_constraints) = A_extension.transpose();
	Coef = Coef.triangularView<Upper>();
	Coef = Coef.transpose();

	pardiso = new PardisoSolver();

	int * ia = Coef.outerIndexPtr();
	vector<int> pardiso_ia(ia, ia + 2 * V_N + 1 + 2 * n_constraints);
	pardiso->ia = pardiso_ia;

	int * ja = Coef.innerIndexPtr();
	int nonzero = Coef.nonZeros();
	vector<int> pardiso_ja(ja, ja + nonzero);
	pardiso->ja = pardiso_ja;

	double * a = Coef.valuePtr();
	vector<double> pardiso_a(a, a + nonzero);
	pardiso->a = pardiso_a;

	pardiso->nnz = pardiso_ja.size();
	pardiso->num = V_N * 2 + 2 * n_constraints;

	pardiso->pardiso_init();

	pardiso->factorize();
}
void MyParafun::embedTriangle(const MatrixXd& V, MatrixXd& flatV)
{
	VectorXd v1 = V.row(1) - V.row(0);
	VectorXd v2 = V.row(2) - V.row(0);

	double norm_v1 = v1.norm();
	double norm_v2 = v2.norm();
	double cos_theta = v1.dot(v2) / (norm_v1*norm_v2);
	double sin_theta = sqrt(1 - cos_theta*cos_theta);

	flatV << 0, 0,
		norm_v1, 0,
		norm_v2*cos_theta, norm_v2*sin_theta;
}

void MyParafun::pck(Mesh &mesh_)
{
	long time_beg, time_end;
	time_beg = clock();

	init();
	int iter = 0;
	cur_position = global_position;		//Q_0
	Tx = T*cur_position;
	Distortion();
	//while (flips.sum() && iter<20)
	while (flips.sum() && iter<20)
	{
		std::cout << "--------------------------- " << to_string(iter) << " times projection" << " ---------------------------" << std::endl;
		std::cout << "max distortions" << distortions.maxCoeff() << std::endl;
		std::cout << "flip number" << flips.sum() << std::endl;
		//computeK();
		for (int i = 0; i < F_N; ++i)
		{
			K(i) = 4 * pow(2.0, iter);
		}
		std::cout << "Max  K value is: " << K.maxCoeff() << std::endl;
		solve();
		++iter;
	}
	/*
	while (iter < 20)
	{
		is_optimal = true;
		std::cout << "******************** " << to_string(iter) << " times projection for fairness" << " ********************" << std::endl;
		std::cout << "max distortions" << distortions.maxCoeff() << std::endl;
		std::cout << "flip number" << flips.sum() << std::endl;
		computeK();
		std::cout << "Max  K value is: " << K.maxCoeff() << std::endl;
		solve();
		++iter;
	}
	*/
	time_end = clock();
	double time_consumption = (time_end - time_beg) / 1000.0;
	iter_project = iter;
	double avg_distortion;
	for (int i = 0; i < F_N; ++i)
	{
		avg_distortion += distortions(i);
	}
	avg_distortion = avg_distortion / F_N;
	std::cout << "=========================== " << "optimal result data" << " ===========================" << std::endl;
	std::cout << "max distortions" << distortions.maxCoeff() << std::endl;
	std::cout << "average distortions" << avg_distortion << std::endl;
	std::cout << "flip number" << flips.sum() << std::endl;
	std::cout << "time_consumption: " << time_consumption << std::endl;
	std::cout << "iterator number of equation" << iter_equ << std::endl;
	std::cout << "iterator number of projective " << iter_project << std::endl;

	writeobj(init_mesh, "result.obj");
	mesh_ = init_mesh;
}

void MyParafun::computeK()
{
	double d_i, temp_di;
	for (int i = 0; i < F_N; ++i)
	{
		if (!flips[i])
		{
			d_i = (distortions(i) - 1) / (distortions(i) + 1);
			if (d_i > para_1)
			{
				temp_di = ((distortions(i) / para_2) - 1) / (1 + (distortions(i) / para_2));
				d_i = max(para_1, temp_di);
			}
			K[i] = d_i;
		}
	}

	Mesh::FaceHandle f_h;
	double sumK;
	int neiK;
	while (flips.sum())
	{
		for (int ii = 0; ii < F_N; ++ii)
		{
			if (flips(ii))
			{
				f_h = Mesh::FaceHandle(ii);
				sumK = 0;
				neiK = 0;
				Mesh::FaceFaceIter ff_it = source_mesh.ff_iter(f_h);
				for (; ff_it; ++ff_it)
				{
					if (!flips(ff_it.handle().idx()))
					{
						neiK++;
						sumK += K(ff_it.handle().idx());
					}
				}
				if (neiK)
				{
					flips(ii) = 0;
					K(ii) = sumK / neiK;
				}
			}
		}
	}

	for (int j = 0; j < 10; ++j)
	{
		for (int i = 0; i < F_N; ++i)
		{
			f_h = Mesh::FaceHandle(i);
			sumK = 0;
			neiK = 0;
			Mesh::FaceFaceIter ff_it = source_mesh.ff_iter(f_h);
			for (; ff_it; ++ff_it)
			{
				neiK++;
				sumK += K(ff_it.handle().idx());
			}
			sumK = sumK / neiK;
			//K[i] = 0.4 * K[i] + 0.6 * sumK;
			K[i] = sumK;
		}
	}

	for (int i = 0; i < F_N; ++i)
	{
		d_i = K(i);
		K[i] = (1 + d_i) / (1 - d_i);
	}
}

void MyParafun::solve()
{
	Anderson accelerator;
	int Anderson_m = 5;
	accelerator.init(Anderson_m, 2, cur_position);
	pre_energy = std::numeric_limits<double>::max();
	double is_convergence;

	for (int iter = 0; iter < iter_max; ++iter)
	{
		Tx = T*cur_position;
		local();

		if (Energy()>pre_energy)
		{
			cur_position = global_position;
			Tx = T*cur_position;
			local();

			std::cout << "alternative solution" << std::endl;
			accelerator.replace(cur_position);

			is_optimal = true;
		}

		cur_energy = Energy();
		std::cout << iter << "  " << cur_energy << std::endl;

		is_convergence = abs(pre_energy - cur_energy) / pre_energy;
		//if (is_convergence < tol_err || (!is_optimal && (0 == flips.sum())))
		if (is_convergence < tol_err || cur_energy < tol_energy)
		{
			break;
		}
		pre_energy = cur_energy;
		//Anderson acceleration
		global();				//G_k
		cur_position = accelerator.compute(global_position);
	}
}

double MyParafun::Energy()
{
	VectorXd energy = Tx - pTx;
	return energy.norm();
}

void MyParafun::local()
{
	pTx = Tx;
	int block_size = 2 * 2;
	int num_blocks = pTx.size() / block_size;
	double K2plus1;
	JacobiSVD<Matrix2d> svdA(2, 2, (ComputeFullU | ComputeFullV));
	Vector2d s, snew, sinv;
	double dist, t;
	bool flipped;
	Vector2d V1, V2;
	Matrix2d U, V;
	Map<Matrix2d> currA(pTx.data());
	double a, b, c, d;
	double n_ab, n_cd;
	bool wasProjected;

	// project
	for (int ii = 0; ii < num_blocks; ii++)
	{
		K2plus1 = K[ii] * K[ii] + 1;
		new (&currA) Map<Matrix2d>(pTx.data() + ii*block_size);

		if (abs(currA.determinant()) < min_value)
		{
			flipped = true;
			flips(ii) = flipped;
			distortions[ii] = std::numeric_limits<double>::infinity();
			currA.setIdentity();
			currA(0, 0) = K[ii] * min_value;
			currA(1, 1) = min_value;
			continue;
		}


		// 2x2 sv's and sign of det
		a = (currA(0, 0) + currA(1, 1)) / 2;
		b = (currA(1, 0) - currA(0, 1)) / 2;
		c = (currA(0, 0) - currA(1, 1)) / 2;
		d = (currA(0, 1) + currA(1, 0)) / 2;
		n_ab = sqrt(a*a + b*b);
		n_cd = sqrt(c*c + d*d);
		flipped = (n_cd > n_ab);
		s(0) = n_ab + n_cd;
		sinv(0) = 1 / s(0);
		s(1) = (flipped ? (n_cd - n_ab) : (n_ab - n_cd));
		sinv(1) = 1 / s(1);
		// calc distortion
		dist = s(0) / s(1);
		// update sign of last singular value
		if (flipped)
			s(1) = -s(1);

		// -- projections --
		wasProjected = false;
		snew = s;
		// project BD
		if ((K[ii] > 0) & (flipped | (dist > K[ii])))
		{
			t = (K[ii] * s(0) + s(1)) / K2plus1;
			snew << t*K[ii], t;
			wasProjected = true;
		}
		// project LB
		if ((lb > 0) & (s(1) < lb))
		{
			snew = snew.array().max(lb);
			wasProjected = true;
		}
		// project LB
		if ((ub > 0) & (s(0) > ub))
		{
			snew = snew.array().min(ub);
			wasProjected = true;
		}
		// actually project the matrix
		if (wasProjected)
		{
			V1(0) = 1;
			V1(1) = (s(0)*s(0) - currA(0, 0)*currA(0, 0) - currA(1, 0)*currA(1, 0)) / (currA(0, 0)*currA(0, 1) + currA(1, 0)*currA(1, 1));
			V1.normalize();
		//	V2(0) = 1;
		//	V2(1) = (s(1)*s(1) - currA(0, 0)*currA(0, 0) - currA(1, 0)*currA(1, 0)) / (currA(0, 0)*currA(0, 1) + currA(1, 0)*currA(1, 1));
		//	V2.normalize();
			V2(0) = V1(1);
			V2(1) = -V1(0);
			V.col(0) = V1;
			V.col(1) = V2;
			U = currA * V * (sinv.asDiagonal());
		/*	Matrix2d A = U*s.asDiagonal()*V.transpose();
			// full svd
			svdA.compute(currA);
			U = svdA.matrixU();
			V = svdA.matrixV();
			Matrix2d B = U * s.asDiagonal() * V.transpose();
			// ssvd*/
			if (flipped)
				U.col(1) = -U.col(1);
			currA = U*snew.asDiagonal()*V.transpose();
		}

		distortions(ii) = dist;
		flips(ii) = flipped;
		maxsv(ii) = s(0);
		minsv(ii) = s(1);
	}
}

void MyParafun::global()
{
	VectorXd right(2 * V_N + 2 * n_constraints);
	right.head(2 * V_N) = T.transpose()*pTx;
	right.bottomRows(2 * n_constraints) = b;
	vector<double> pardiso_u(right.data(), right.data() + 2 * V_N + 2 * n_constraints);
	pardiso->rhs = pardiso_u;
	pardiso->pardiso_solver();
	++iter_equ;

	for (size_t i = 0; i < 2 * V_N; i++)
	{
		global_position(i) = (pardiso->result)[i];
	}
}

void MyParafun::Distortion()
{
	int block_size = 2 * 2;
	double a, b, c, d;
	double n_ab, n_cd;
	int num_blocks = Tx.size() / block_size;
	Vector2d s;
	double dist, t;
	bool flipped;
	Map<Matrix2d> currA(Tx.data());

	// project
	for (int ii = 0; ii < num_blocks; ii++)
	{
		new (&currA) Map<Matrix2d>(Tx.data() + ii*block_size);
		
		// 2x2 sv's and sign of det
		a = (currA(0, 0) + currA(1, 1)) / 2;
		b = (currA(0, 1) - currA(1, 0)) / 2;
		c = (currA(0, 0) - currA(1, 1)) / 2;
		d = (currA(0, 1) + currA(1, 0)) / 2;
		n_ab = sqrt(a*a + b*b);
		n_cd = sqrt(c*c + d*d);

		if (abs(currA.determinant()) < min_value)
		{
			flipped = true;
			flips(ii) = flipped;
			distortions[ii] = std::numeric_limits<double>::infinity();
		}
		else
		{
			flipped = (n_cd > n_ab);
			flips(ii) = flipped;
			s(0) = n_ab + n_cd;
			s(1) = (flipped ? (n_cd - n_ab) : (n_ab - n_cd));
			// calc distortion
			dist = s(0) / s(1);
			distortions[ii] = dist;
		}
	}
}