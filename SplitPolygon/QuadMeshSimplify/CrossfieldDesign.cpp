#include "CrossfieldDesign.h"

#include <iostream>

void CrossfieldDesign::construct_frames()
{
	frames.clear();
	frames.resize(mesh.n_faces());

	for (auto f_h : mesh.faces())
	{
		if (f_type[f_h.idx()] == 0) continue;

		auto h0 = *mesh.cfh_begin(f_h); //const face-halfedge circulator

		frames[f_h.idx()].first = mesh.calc_edge_vector(h0).normalize();
		OpenMesh::Vec3d nor = mesh.normal(f_h);
		OpenMesh::Vec3d cro = OpenMesh::cross(nor, frames[f_h.idx()].first);
		frames[f_h.idx()].second = OpenMesh::cross(mesh.normal(f_h), frames[f_h.idx()].first);
	}
}

void CrossfieldDesign::construct_complex4(const std::vector<OpenMesh::Vec3d>& field)
{
	field_complex4.clear();
	field_complex4.resize(mesh.n_faces());

	for (auto f_h : mesh.faces())
	{
		if (f_type[f_h.idx()] != 1) continue;
		
		field_complex4[f_h.idx()] = normal_pow(get_complex(field[f_h.idx()], f_h), 4.0);
	}
}

void CrossfieldDesign::recover_field(std::vector<OpenMesh::Vec3d>& field)
{
	for (auto f_h : mesh.faces())
	{
		if (f_type[f_h.idx()] != 2) continue;

		field[f_h.idx()] = get_vector(normal_pow(field_complex4[f_h.idx()], 0.25), f_h);
	}
}

void CrossfieldDesign::construct_solver()
{
	bvec_fid.clear();
	ivec_fid.clear();
	fid_vec.clear();
	fid_vec.resize(mesh.n_faces());

	for (int i = 0; i < mesh.n_faces(); i++)
	{
		if (f_type[i] == 0) continue;

		auto& vec_fid = (f_type[i] == 1) ? bvec_fid : ivec_fid;

		fid_vec[i] = vec_fid.size();
		vec_fid.push_back(i);
	}
		
	int n_constraints = bvec_fid.size();
	int n_vars = ivec_fid.size();

	complex_mat mat_L(n_vars, n_vars);
	mat_B.resize(n_vars, n_constraints);
	mat_B.setZero();
	
	std::vector<Eigen::Triplet<std::complex<double>>> mat_L_triplet, mat_B_triplet;
	mat_L_triplet.reserve(6 * n_vars);
	mat_B_triplet.reserve(2 * n_constraints);

	double max_diff = 0.0;

	for (auto f0 : mesh.faces())
	{
		if (f_type[f0.idx()] != 2) continue;
		int vec_id = fid_vec[f0.idx()];

		for (auto fh_h : mesh.fh_range(f0))
		{
			auto f1 = mesh.opposite_face_handle(fh_h);
			if (!f1.is_valid() || f_type[f1.idx()] == 0) continue;

			auto vec_h = mesh.calc_edge_vector(fh_h).normalize();
			vec_h = ((fh_h.idx() & 1) ? -vec_h : vec_h);	//orientation

			double a0 = mesh.calc_sector_angle(mesh.next_halfedge_handle(fh_h));
			double a1 = mesh.calc_sector_angle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(fh_h)));
			//double weight = (1.0 / std::tan(a0) + 1.0 / std::tan(a1)) / 2.0;
			double weight = 1.0;

			std::complex<double> e0 = normal_pow(get_complex(vec_h, f0), -4.0);
			std::complex<double> e1 = normal_pow(get_complex(vec_h, f1), -4.0);

			mat_L_triplet.emplace_back(vec_id, vec_id, weight);
			(f_type[f1.idx()] == 1 ? mat_B_triplet : mat_L_triplet).emplace_back(vec_id, fid_vec[f1.idx()], -weight * e1 / e0);

		}
	}

	mat_L.setFromTriplets(mat_L_triplet.begin(), mat_L_triplet.end());
	mat_B.setFromTriplets(mat_B_triplet.begin(), mat_B_triplet.end());

	complex_solver.compute(mat_L);
	

	/*std::ofstream ofs1("C:\\Users\\ChiZhang\\Desktop\\matL_real.txt");
	for (int i = 0; i < mat_L.rows(); ++i)
	{
		for (int j = 0; j < mat_L.cols(); ++j)
		{
			ofs1 << std::real(mat_L.coeff(i, j)) << " ";
		}
		ofs1 << std::endl;
	}

	std::ofstream ofs2("C:\\Users\\ChiZhang\\Desktop\\matL_image.txt");
	for (int i = 0; i < mat_L.rows(); ++i)
	{
		for (int j = 0; j < mat_L.cols(); ++j)
		{
			ofs2 << std::imag(mat_L.coeff(i, j)) << " ";
		}
		ofs2 << std::endl;
	}*/

	//std::cout/* << "Det " << complex_solver.determinant()*/ << " Success : " << (complex_solver.info() == 0 ? "true" : "false") << std::endl;
	//std::cout << complex_solver.rows() << std::endl;
}

void CrossfieldDesign::construct_new_solver()
{
	
}

void CrossfieldDesign::solve_complex4()
{
	int n_vars = ivec_fid.size();
	int n_constraints = bvec_fid.size();

	Eigen::VectorXcd x1(n_vars), x2(n_constraints), b(n_vars);

	for (int i = 0; i < n_constraints; i++)
	{
		x2[i] = field_complex4[bvec_fid[i]];
	}

	b = -mat_B * x2;
	//std::cout << "b.rows(): " << b.rows() << std::endl;

	x1 = complex_solver.solve(b);

	for (int i = 0; i < n_vars; i++)
	{
		field_complex4[ivec_fid[i]] = x1[i];
	}

	double max_norm = 0;
	double min_norm = std::numeric_limits<double>::max();
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		if (f_type[i] == 0) continue;

		max_norm = std::max(std::norm(field_complex4[i]), max_norm);
		min_norm = std::min(std::norm(field_complex4[i]), min_norm);
	}


	//std::ofstream ofs1("C:\\Users\\ChiZhang\\Desktop\\breal.txt");
	//for (int i = 0; i < b.size(); ++i)
	//{
	//	ofs1 << std::real(b[i]) << std::endl;
	//}

	//std::ofstream ofs2("C:\\Users\\ChiZhang\\Desktop\\bimage.txt");
	//for (int i = 0; i < b.size(); ++i)
	//{
	//	ofs2 << std::imag(b[i]) << std::endl;
	//}
/*
	std::ofstream ofs("C:\\Users\\ChiZhang\\Desktop\\test.txt");
	ofs << field_complex4.size() << std::endl;

	for (int i = 0; i < field_complex4.size(); ++i)
	{
		ofs << field_complex4[i] << std::endl;
	}*/
	std::cout << "Max Norm : " << max_norm << "\nMin Norm : " << min_norm << std::endl;
}

void CrossfieldDesign::compute_mesh()
{
	construct_frames();
	construct_solver();

	ready = true;
}

bool CrossfieldDesign::compute_boundary_constraint(std::vector<OpenMesh::Vec3d>& field)
{
	if (!ready) return false;

	construct_complex4(field);
	solve_complex4();
	recover_field(field);

	return true;
}

void ComputeExtractField(std::string fileName, std::vector<OpenMesh::Vec3d> &fields, const std::vector<int>& fType)
{
	MyMesh mesh;
	OpenMesh::IO::read_mesh(mesh, fileName.c_str());
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_face_normals();

	CrossfieldDesign cfd(mesh, fType);
	cfd.compute_mesh();

	std::vector<OpenMesh::Vec3d> oFields;
	for (int i = 0; i < fields.size(); ++i)
	{
		OpenMesh::Vec3d currValue(fields[i][0], fields[i][1], fields[i][2]);
		oFields.push_back(currValue);
	}
	cfd.compute_boundary_constraint(oFields);

	for (int i = 0; i < oFields.size(); ++i)
	{
		fields[i][0] = oFields[i][0];
		fields[i][1] = oFields[i][1];
		fields[i][2] = oFields[i][2];
	}
}