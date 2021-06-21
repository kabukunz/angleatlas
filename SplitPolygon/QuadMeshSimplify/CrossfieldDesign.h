#pragma once

#include "MeshDefinition.h"

#include <Eigen\Eigen>
#include <Eigen\Sparse>
#include <Eigen/OrderingMethods>
//#include <Eigen\SPQRSupport>
#include <complex>
#include <OpenMesh/Core/IO/MeshIO.hh>

class CrossfieldDesign
{
public:
	CrossfieldDesign(const MyMesh& _mesh, const std::vector<int>& _f_type) :mesh(_mesh), f_type(_f_type){};
	~CrossfieldDesign() = default;

	void compute_mesh();
	bool compute_boundary_constraint(std::vector<OpenMesh::Vec3d>& field);

private:
	bool ready = false;
	const MyMesh& mesh;
	const std::vector<int>& f_type;//0 omitted, 1 constraint, 2 solving

	using face_frame = std::pair<OpenMesh::Vec3d, OpenMesh::Vec3d>;
	std::vector<face_frame> frames;
	std::vector<std::complex<double>> field_complex4;

	void construct_frames();
	void construct_complex4(const std::vector<OpenMesh::Vec3d>& field);
	void recover_field(std::vector<OpenMesh::Vec3d>& field);

	std::complex<double> get_complex(const OpenMesh::Vec3d& vec, OpenMesh::FaceHandle f_h)
	{
		return std::complex<double>(OpenMesh::dot(vec, frames[f_h.idx()].first), OpenMesh::dot(vec, frames[f_h.idx()].second));
	}

	OpenMesh::Vec3d get_vector(const std::complex<double>& c, OpenMesh::FaceHandle f_h)
	{
		return c.real() * frames[f_h.idx()].first + c.imag() * frames[f_h.idx()].second;
	}

	std::complex<double> normal_pow(const std::complex<double>& c, double x)
	{
		return std::polar(1.0, std::arg(c) * x);
	}

	using complex_mat = Eigen::SparseMatrix<std::complex<double>>;
	//Eigen::SparseLU<complex_mat> complex_solver;
	//Eigen::ConjugateGradient<complex_mat> complex_solver;
	//Eigen::SparseQR<complex_mat, Eigen::COLAMDOrdering<int>> complex_solver;
	Eigen::SimplicialCholesky<complex_mat> complex_solver;
	complex_mat mat_B;

	std::vector<int> fid_vec;
	std::vector<int> bvec_fid, ivec_fid;

	void construct_solver();
	void construct_new_solver();
	void solve_complex4();
};

void ComputeExtractField(std::string fileName, std::vector<OpenMesh::Vec3d> &fields, const std::vector<int>& fType);