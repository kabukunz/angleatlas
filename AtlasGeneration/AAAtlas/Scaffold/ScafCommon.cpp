#include "ScafCommon.h"

void adjacency_matrix(const Eigen::MatrixXi& F, Eigen::SparseMatrix<int>& A)
{
	typedef Eigen::Triplet<int> IJV;
	std::vector<IJV > ijv;
	ijv.reserve(F.size() * 2);
	// Loop over **simplex** (i.e., **not quad**)
	for (int i = 0; i < F.rows(); i++)
	{
		// Loop over this **simplex**
		for (int j = 0; j < F.cols(); j++)
			for (int k = j + 1; k < F.cols(); k++)
			{
				// Get indices of edge: s --> d
				int s = F(i, j);
				int d = F(i, k);
				ijv.push_back(IJV(s, d, 1));
				ijv.push_back(IJV(d, s, 1));
			}
	}

	const int n = F.maxCoeff() + 1;
	A.resize(n, n);
	switch (F.cols())
	{
	case 3:
		A.reserve(6 * (F.maxCoeff() + 1));
		break;
	case 4:
		A.reserve(26 * (F.maxCoeff() + 1));
		break;
	}
	A.setFromTriplets(ijv.begin(), ijv.end());

	// Force all non-zeros to be one

	// Iterate over outside
	for (int k = 0; k < A.outerSize(); ++k)
	{
		// Iterate over inside
		for (Eigen::SparseMatrix<int>::InnerIterator it(A, k); it; ++it)
		{
			assert(it.value() != 0);
			A.coeffRef(it.row(), it.col()) = 1;
		}
	}
}


void components(const Eigen::SparseMatrix<int>& A,Eigen::MatrixXi& C, Eigen::MatrixXi& counts)
{
	using namespace Eigen;
	using namespace std;
	const size_t n = A.rows();
	Array<bool, Dynamic, 1> seen = Array<bool, Dynamic, 1>::Zero(n, 1);
	C.resize(n, 1);
	int id = 0;
	vector<int> vcounts;
	// breadth first search
	for (int k = 0; k < A.outerSize(); ++k)
	{
		if (seen(k))
		{
			continue;
		}
		queue<int> Q;
		Q.push(k);
		vcounts.push_back(0);
		while (!Q.empty())
		{
			const int f = Q.front();
			Q.pop();
			if (seen(f))
			{
				continue;
			}
			seen(f) = true;
			C(f, 0) = id;
			vcounts[id]++;
			// Iterate over inside
			for (SparseMatrix<int>::InnerIterator it(A, f); it; ++it)
			{
				const int g = it.index();
				if (!seen(g) && it.value())
				{
					Q.push(g);
				}
			}
		}
		id++;
	}
	assert((size_t)id == vcounts.size());
	const size_t ncc = vcounts.size();
	assert((size_t)C.maxCoeff() + 1 == ncc);
	counts.resize(ncc, 1);
	for (size_t i = 0; i < ncc; i++)
	{
		counts(i) = vcounts[i];
	}
}

void doublearea(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::VectorXd& dblA)
{

	dblA.resize(F.rows());

	if (V.cols() == 3)
	{
		Eigen::Vector3d e10, e20, a;
		Eigen::Vector3d v0, v1, v2;
		for (size_t i = 0; i < F.rows(); i++)
		{
			v0 = V.row(F(i, 0));
			v1 = V.row(F(i, 1));
			v2 = V.row(F(i, 2));
			e10 = v1 - v0;
			e20 = v2 - v0;
			a = e10.cross(e20);
			dblA(i) = a.norm();
		}
	}
	else
	{
		Eigen::Vector2d e10, e20;
		Eigen::Vector2d v0, v1, v2;
		double a;
		for (size_t i = 0; i < F.rows(); i++)
		{
			v0 = V.row(F(i, 0));
			v1 = V.row(F(i, 1));
			v2 = V.row(F(i, 2));
			e10 = v1 - v0;
			e20 = v2 - v0;
			a = e10(0)*e20(1)-e10(1)*e20(0);
			dblA(i) = abs(a);
		}
	}
}

void remove_unreferenced(const Eigen::MatrixXd & V, const Eigen::MatrixXi & Fsrc, Eigen::MatrixXd& PV, Eigen::MatrixXi& PF, Eigen::VectorXi& VI)
{
	Eigen::VectorXi mark;
	mark.setZero(V.rows());
	for (size_t i = 0; i < Fsrc.rows(); i++)
	{
		for (size_t j = 0; j < Fsrc.cols(); j++)
		{
			mark(Fsrc(i, j)) = 1;
		}
	}
	int newsize = mark.sum();

	VI.resize(V.rows());
	Eigen::VectorXi J;
	J.resize(newsize);

	int counts = 0;

	for (size_t i = 0; i < mark.size(); i++)
	{
		if (mark(i) == 1)
		{
			VI(i) = counts;
			J(counts) = i;
			counts++;
		}
		else
		{
			VI(i) = -1;
		}
	}

	PF.resize(Fsrc.rows(), Fsrc.cols());
	for (size_t i = 0; i < Fsrc.rows(); i++)
	{
		for (size_t j = 0; j < Fsrc.cols(); j++)
		{
			PF(i, j) = VI(Fsrc(i, j));
		}
	}
	
	PV.resize(newsize, 2);

	for (size_t i = 0; i < newsize; i++)
	{
		PV.row(i) = V.row(J(i));
		//PV(i, 0) = V(J(i), 0);
		//PV(i, 1) = V(J(i), 1);
	}
}


void boundary_loop(const Eigen::MatrixXi &F_ref, std::vector<std::vector<int>>& boundaryloop)
{
	using namespace std;
	std::vector<std::vector<int>> boundaryEdges;
	std::vector<std::vector<int>> edges;
	int n_fvs = F_ref.cols();

	for (int it = 0; it < F_ref.rows(); it++)
	{
		for (int i = 0; i < n_fvs; i++)
		{
			int var = F_ref(it, i);
			int var_n = F_ref(it, (i + 1) % n_fvs);
			if (var > var_n) std::swap(var, var_n);
			std::vector<int> edge(4);
			edge[0] = var;
			edge[1] = var_n;
			edge[2] = it;
			edge[3] = i;
			edges.emplace_back(edge);
		}
	}
	std::sort(edges.begin(), edges.end());
	int i = 1;
	for (; i < edges.size();)
	{
		auto& r1 = edges[i - 1];
		auto& r2 = edges[i];
		if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
		{
			i += 2;
		}
		else
		{
			boundaryEdges.emplace_back(edges[i - 1]);
			i++;
		}
	}
	if (i == edges.size())
		boundaryEdges.emplace_back(edges.back());

	for (auto&var : boundaryEdges)
	{
		var[0] = F_ref(var[2], var[3]);
		var[1] = F_ref(var[2], (var[3] + 1) % n_fvs);
	}
	int ev0 = boundaryEdges.front()[0];
	int ev1 = boundaryEdges.front()[1];

	vector<int> visited;
	visited.resize(boundaryEdges.size(), 0);
	visited[0] = 1;
	vector<int> loop0;
	loop0.push_back(ev1);
	while (ev1 != ev0)
	{
		for (int i = 1; i < boundaryEdges.size(); i++)
		{
			if (visited[i] == 1)
				continue;
			if (boundaryEdges[i][0] == ev1)
			{
				visited[i] = 1;
				ev1 = boundaryEdges[i][1];
				loop0.push_back(ev1);
				break;
			}
		}
	}
	boundaryloop.emplace_back(loop0);
}

void readobj(const string & filename, Eigen::MatrixXd & Vpos, Eigen::MatrixXd & UV_V, Eigen::MatrixXi & Fids, Eigen::MatrixXi & UV_F)
{
	vector<vector<double>> V, TC, N;
	vector<vector<int>>F, FTC, FN;

	FILE * obj_file=fopen(filename.c_str(), "r");

	// File open was successful so clear outputs
	V.clear();
	TC.clear();
	N.clear();
	F.clear();
	FTC.clear();
	FN.clear();

	// variables and constants to assist parsing the .obj file
	// Constant strings to compare against
	std::string v("v");
	std::string vn("vn");
	std::string vt("vt");
	std::string f("f");
	std::string tic_tac_toe("#");
#ifndef IGL_LINE_MAX
#  define IGL_LINE_MAX 2048
#endif

	char line[IGL_LINE_MAX];
	int line_no = 1;
	while (fgets(line, IGL_LINE_MAX, obj_file) != NULL)
	{
		char type[IGL_LINE_MAX];
		// Read first word containing type
		if (sscanf(line, "%s", type) == 1)
		{
			// Get pointer to rest of line right after type
			char * l = &line[strlen(type)];
			if (type == v)
			{
				std::istringstream ls(&line[1]);
				std::vector<double > vertex{ std::istream_iterator<double >(ls), std::istream_iterator<double >() };

				if (vertex.size() < 3)
				{
					fprintf(stderr,
						"Error: readOBJ() vertex on line %d should have at least 3 coordinates",
						line_no);
					fclose(obj_file);
					return;
				}

				V.push_back(vertex);
			}
			else if (type == vn)
			{
				double x[3];
				int count =
					sscanf(l, "%lf %lf %lf\n", &x[0], &x[1], &x[2]);
				if (count != 3)
				{
					fprintf(stderr,
						"Error: readOBJ() normal on line %d should have 3 coordinates",
						line_no);
					fclose(obj_file);
					return;
				}
				std::vector<double > normal(count);
				for (int i = 0; i < count; i++)
				{
					normal[i] = x[i];
				}
				N.push_back(normal);
			}
			else if (type == vt)
			{
				double x[3];
				int count =
					sscanf(l, "%lf %lf %lf\n", &x[0], &x[1], &x[2]);
				if (count != 2 && count != 3)
				{
					fprintf(stderr,
						"Error: readOBJ() texture coords on line %d should have 2 "
						"or 3 coordinates (%d)",
						line_no, count);
					fclose(obj_file);
					return;
				}
				std::vector<double > tex(count);
				for (int i = 0; i < count; i++)
				{
					tex[i] = x[i];
				}
				TC.push_back(tex);
			}
			else if (type == f)
			{
				const auto & shift = [&V](const int i)->int
				{
					return i < 0 ? i + V.size() : i - 1;
				};
				const auto & shift_t = [&TC](const int i)->int
				{
					return i < 0 ? i + TC.size() : i - 1;
				};
				const auto & shift_n = [&N](const int i)->int
				{
					return i < 0 ? i + N.size() : i - 1;
				};
				std::vector<int > f;
				std::vector<int > ftc;
				std::vector<int > fn;
				// Read each "word" after type
				char word[IGL_LINE_MAX];
				int offset;
				while (sscanf(l, "%s%n", word, &offset) == 1)
				{
					// adjust offset
					l += offset;
					// Process word
					long int i, it, in;
					if (sscanf(word, "%ld/%ld/%ld", &i, &it, &in) == 3)
					{
						f.push_back(shift(i));
						ftc.push_back(shift_t(it));
						fn.push_back(shift_n(in));
					}
					else if (sscanf(word, "%ld/%ld", &i, &it) == 2)
					{
						f.push_back(shift(i));
						ftc.push_back(shift_t(it));
					}
					else if (sscanf(word, "%ld//%ld", &i, &in) == 2)
					{
						f.push_back(shift(i));
						fn.push_back(shift_n(in));
					}
					else if (sscanf(word, "%ld", &i) == 1)
					{
						f.push_back(shift(i));
					}
					else
					{
						fprintf(stderr,
							"Error: readOBJ() face on line %d has invalid element format\n",
							line_no);
						fclose(obj_file);
						return;
					}
				}
				if (
					(f.size() > 0 && fn.size() == 0 && ftc.size() == 0) ||
					(f.size() > 0 && fn.size() == f.size() && ftc.size() == 0) ||
					(f.size() > 0 && fn.size() == 0 && ftc.size() == f.size()) ||
					(f.size() > 0 && fn.size() == f.size() && ftc.size() == f.size()))
				{
					// No matter what add each type to lists so that lists are the
					// correct lengths
					F.push_back(f);
					FTC.push_back(ftc);
					FN.push_back(fn);
				}
				else
				{
					fprintf(stderr,
						"Error: readOBJ() face on line %d has invalid format\n", line_no);
					fclose(obj_file);
					return;
				}
			}
			else if (strlen(type) >= 1 && (type[0] == '#' ||
				type[0] == 'g' ||
				type[0] == 's' ||
				strcmp("usemtl", type) == 0 ||
				strcmp("mtllib", type) == 0))
			{
				//ignore comments or other shit
			}
			else
			{
				//ignore any other lines
				fprintf(stderr,
					"Warning: readOBJ() ignored non-comment line %d:\n  %s",
					line_no,
					line);
			}
		}
		else
		{
			// ignore empty line
		}
		line_no++;
	}
	fclose(obj_file);

	assert(F.size() == FN.size());
	assert(F.size() == FTC.size());



	Vpos.resize(V.size(), 3);
	UV_V.resize(TC.size(), 2);
	Fids.resize(F.size(), 3);
	UV_F.resize(FTC.size(), 3);
	for (size_t i = 0; i < V.size(); i++)
	{
		auto& ver = V[i];
		for (size_t j = 0; j < ver.size(); j++)
		{
			Vpos(i, j) = ver[j];
		}
	}

	for (size_t i = 0; i < TC.size(); i++)
	{
		auto& ver = TC[i];
		for (size_t j = 0; j < ver.size(); j++)
		{
			UV_V(i, j) = ver[j];
		}
	}


	for (size_t i = 0; i < F.size(); i++)
	{
		auto& ver = F[i];
		for (size_t j = 0; j < ver.size(); j++)
		{
			Fids(i, j) = ver[j];
		}
	}


	for (size_t i = 0; i < FTC.size(); i++)
	{
		auto& ver = FTC[i];
		for (size_t j = 0; j < ver.size(); j++)
		{
			UV_F(i, j) = ver[j];
		}
	}

}
