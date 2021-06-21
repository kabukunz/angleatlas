#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include<math.h>
#include<vector>
#include<iostream>
#include<fstream>
#include"Eigen\src\Eigenvalues\EigenSolver.h"
#include"Eigen/Dense"
#include"Eigen/Sparse"
#include <Eigen\SparseCholesky>
#include "Eigen/SVD"
#include "MeshViewer\MeshDefinition.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "Scaffold\PardisoSolver.h"
#include "MyAnderson.h"

using namespace Eigen;
using namespace std;

typedef Triplet<double> Tri;

class MyParafun
{
public:
	//Parafun(string srcname, string initname, int dim = 2);
	MyParafun(Mesh &mesh_src_, Mesh &mesh_tutte, int dim = 2);
	~MyParafun();

	bool writeobj(Mesh& out_mesh, const string out_filename);
	
	void init();				//存储数据
	void embedTriangle(const MatrixXd& V, MatrixXd& flatV);

	void local();				//把当前的mapping投影到有界区域
	double Energy();

	void pck(Mesh &mesh_);
	void solve();
	void global();

	void Distortion();
	void computeK();

protected:
	double	lb;						// lower bound on SVs(-1 = disabled)
	double	ub;						// upper bound on SVs(-1 = disabled)
	double	tol_err;				// tolerance for stopping BD projection iterations
	double	tol_energy;				// tolerance for stopping BD projection iterations
	bool	use_weighted_metric;	// use a weighted metric
	int		iter_max;				// maximal number of BD projection iterations

	int		n_constraints;			// the number of constraints
	SparseMatrix<double> A;			// coefficient of constraints
	VectorXd			b;			// right side value of constraints

	VectorXd Tx;				//lifted variable differentials
	VectorXd tanNormal;			//normal
	VectorXd pTx;				//projected 
	VectorXd distortions;		//distortions of differentials i.e. condition number
	VectorXd flips;
	VectorXd minsv;
	VectorXd maxsv;
	SimplicialCholesky<SparseMatrix<double>> TCholesky;						//precompile
	SparseMatrix<double> T;		//lifts init_position(d*n) into a higher dimension space of differentials(m*d*d)

	int			F_N;
	int			V_N;
	MatrixXi	F;
	VectorXd	source_position;
	VectorXd	global_position;
	Mesh		source_mesh;
	Mesh		init_mesh;
	int			source_dim;

	VectorXd	cur_position;
	VectorXd	K;					// conformal distortion bound//vector<double>	energy_everyiter;
	double	para_1;
	double	para_2;
	double	cur_energy;
	double  pre_energy;
	PardisoSolver* pardiso;

	bool is_optimal;
	int iter_equ;
	int iter_project;
};

