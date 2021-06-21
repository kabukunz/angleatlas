#include <OpenMesh/Core/IO/MeshIO.hh>

#include "Mesh_doubleIO.h"

bool Mesh_doubleIO::load_mesh(Mesh& _mesh, const char* _filename, bool load_texture)
{
	if (!OpenMesh::IO::read_mesh(_mesh, _filename))
	{
		return false;
	}

	switch (get_file_type(_filename))
	{
	case file_type::obj:
		return load_obj(_mesh, _filename, load_texture);
	case file_type::off:
		return load_off(_mesh, _filename);
	default:
		return true;
	}
}

bool Mesh_doubleIO::save_mesh(const Mesh& _mesh, const char* _filename, bool save_texture)
{
	switch (get_file_type(_filename))
	{
	case file_type::obj:
		return save_obj(_mesh, _filename, save_texture);
	case file_type::off:
		return save_off(_mesh, _filename);
	default:
		return OpenMesh::IO::write_mesh(_mesh, _filename, OpenMesh::IO::Options::Default, std::numeric_limits<Mesh::Scalar>::max_digits10);
	}
}

bool Mesh_doubleIO::save_uv_mesh(const Mesh& _mesh, const char* _filename)
{
	OpenMesh::MPropHandleT<std::string> mstr_tfile;
	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
	OpenMesh::HPropHandleT<int> hvt_index;

	if (!_mesh.get_property_handle(mvt_list, "mvt_list") || !_mesh.get_property_handle(hvt_index, "hvt_index"))
	{
		std::cout << "Texture data is invalid." << std::endl;
		return false;
	}
	if (_mesh.property(mvt_list).empty())
	{
		std::cout << "Texture data is invalid." << std::endl;
		return false;
	}

	Mesh uv_mesh;

	int n_uv = _mesh.property(mvt_list).size();
	for (int i = 0; i < n_uv; i++)
	{
		auto uv0 = _mesh.property(mvt_list)[i];
		uv_mesh.add_vertex(Mesh::Point(uv0[0], uv0[1], 0.0));
	}

	for (int i = 0; i < _mesh.n_faces(); i++)
	{
		auto f_h = _mesh.face_handle(i);

		std::vector<OpenMesh::VertexHandle> f_v;
		f_v.reserve(_mesh.valence(f_h));
		for (auto fh_h : _mesh.fh_range(f_h))
		{
			f_v.emplace_back(_mesh.property(hvt_index, fh_h));
		}

		uv_mesh.add_face(f_v);
	}

	save_mesh(uv_mesh, _filename);
	return true;
}

Mesh_doubleIO::file_type Mesh_doubleIO::get_file_type(const char* _filename)
{
	std::string filetype(_filename);
	filetype = filetype.substr(filetype.length() - 4, 4);

	for (int i = 1; i < filetype.length(); i++)
	{
		filetype[i] = std::tolower(filetype[i]);
	}

	if (filetype.compare(".obj") == 0)
	{
		return file_type::obj;
	}
	else if (filetype.compare(".off") == 0)
	{
		return file_type::off;
	}
	else
	{
		return file_type::others;
	}
}

bool Mesh_doubleIO::load_obj(Mesh& _mesh, const char* _filename, bool load_texture)
{
	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
	OpenMesh::HPropHandleT<int> hvt_index;

	std::ifstream obj_file(_filename);
	std::vector<Mesh::Point> vec_mesh(_mesh.n_vertices());

	if (!obj_file.is_open())
	{
		return false;
	}

	if (load_texture)
	{
		_mesh.add_property(mvt_list, "mvt_list");
		_mesh.add_property(hvt_index, "hvt_index");

		_mesh.property(mvt_list).reserve(2 * _mesh.n_vertices());
	}
	bool properties_added = load_texture;

	std::string line, prefix;
	line.reserve(1024);
	prefix.reserve(16);

	int count_v = 0;
	int count_f = 0;
	while (!obj_file.eof())
	{
		std::getline(obj_file, line);
		if (line == "") continue;
		std::istringstream iss(line);
		iss >> prefix;

		if (prefix == "v")
		{
			iss >> vec_mesh[count_v][0] >> vec_mesh[count_v][1] >> vec_mesh[count_v][2];
			count_v++;
		}
		else if (load_texture && prefix == "vt")
		{
			Mesh::TexCoord2D tex;
			iss >> tex[0] >> tex[1];

			_mesh.property(mvt_list).push_back(tex);
		}
		else if (load_texture && prefix == "f")
		{
			if (_mesh.property(mvt_list).empty())
			{
				load_texture = false;
				continue;
			}

			int v_id, vt_id, vn_id;
			std::map<int, int> vid2vtid;

			while (!iss.eof())
			{
				iss >> v_id;
				if (iss.get() != '/') load_texture = false;
				if (!load_texture) break;
				iss >> vt_id;
				if (iss.get() == '/') iss >> vn_id;
				vid2vtid[v_id] = vt_id;
			}
			if (!load_texture) continue;

			for (auto fh_h : _mesh.fh_range(_mesh.face_handle(count_f)))
			{
				_mesh.property(hvt_index, fh_h) = vid2vtid[_mesh.to_vertex_handle(fh_h).idx() + 1] - 1;
			}
			count_f++;
		}
	}

	for (auto v_h : _mesh.vertices())
	{
		_mesh.point(v_h) = vec_mesh[v_h.idx()];
	}

	obj_file.close();

	if (!load_texture && properties_added)
	{
		_mesh.remove_property(mvt_list);
		_mesh.remove_property(hvt_index);
	}

	return true;
}

bool Mesh_doubleIO::save_obj(const Mesh& _mesh, const char* _filename, bool save_texture)
{
	OpenMesh::MPropHandleT<std::string> mstr_tfile;
	OpenMesh::MPropHandleT<std::vector<Mesh::TexCoord2D>> mvt_list;
	OpenMesh::HPropHandleT<int> hvt_index;

	if (save_texture && (!_mesh.get_property_handle(mvt_list, "mvt_list") || !_mesh.get_property_handle(hvt_index, "hvt_index")))
	{
		std::cout << "Texture data is invalid." << std::endl;
		return false;
	}
	if (save_texture && _mesh.property(mvt_list).empty())
	{
		std::cout << "Texture data is invalid." << std::endl;
		return false;
	}

	std::ofstream obj_file(_filename);

	if (!obj_file.is_open())
	{
		return false;
	}

	std::string str_filename(_filename);
	obj_file << "# " << _mesh.n_vertices() << " vertices, ";
	obj_file << _mesh.n_faces() << " faces\n";

// 	if (save_texture)
// 	{
// 		obj_file << "mtllib ./" << str_filename.substr(str_filename.find_last_of("/\\") + 1) << ".mtl\n";
// 	}
	
	obj_file << std::setprecision(std::numeric_limits<Mesh::Scalar>::max_digits10);

	for (int i = 0; i < _mesh.n_vertices(); i++)
	{
		auto&& p0 = _mesh.point(_mesh.vertex_handle(i));
		obj_file << "v " << p0[0] << " " << p0[1] << " " << p0[2] << "\n";
	}
	if (save_texture)
	{
		int n_uv = _mesh.property(mvt_list).size();
		for (int i = 0; i < n_uv; i++)
		{
			auto uv0 = _mesh.property(mvt_list)[i];
			obj_file << "vt " << uv0[0] << " " << uv0[1] << "\n";
		}
	}
	for (int i = 0; i < _mesh.n_faces(); i++)
	{
		auto f_h = _mesh.face_handle(i);

		obj_file << "f";
		for (auto fh_h : _mesh.fh_range(f_h))
		{
			if (save_texture)
			{
				obj_file << " " << _mesh.to_vertex_handle(fh_h).idx() + 1 << "/" << _mesh.property(hvt_index, fh_h) + 1;
			}
			else
			{
				obj_file << " " << _mesh.to_vertex_handle(fh_h).idx() + 1;
			}
		}
		obj_file << "\n";
	}

	obj_file.close();

// 	if (save_texture)
// 	{
// 		std::ofstream mtl_file(str_filename + ".mtl");
// 		mtl_file << "newmtl material_0\n";
// 		mtl_file << "Ka 1.000000 1.000000 1.000000\n";
// 		mtl_file << "Kd 1.000000 1.000000 1.000000\n";
// 		mtl_file << "Ks 0.000000 0.000000 0.000000\n";
// 		mtl_file << "Tr 1.000000\nillum 2\nNs 0.000000\n\n";
// 		mtl_file << "map_Ka " << _mesh.property(mstr_tfile) << "\n";
// 		mtl_file << "map_Kd " << _mesh.property(mstr_tfile) << "\n";
// 
// 		mtl_file.close();
// 	}

	return true;
}

bool Mesh_doubleIO::load_off(Mesh& _mesh, const char* _filename)
{
	std::ifstream off_file(_filename);
	std::vector<Mesh::Point> vec_mesh(_mesh.n_vertices());

	if (!off_file.is_open())
	{
		return false;
	}

	std::string line, prefix;
	line.reserve(1024);
	prefix.reserve(16);

	std::getline(off_file, line);
	std::getline(off_file, line);
	for (int i = 0; i < _mesh.n_vertices(); i++)
	{
		off_file >> vec_mesh[i][0] >> vec_mesh[i][1] >> vec_mesh[i][2];
	}

	for (auto v_h : _mesh.vertices())
	{
		_mesh.point(v_h) = vec_mesh[v_h.idx()];
	}

	off_file.close();

	return true;
}

bool Mesh_doubleIO::save_off(const Mesh& _mesh, const char* _filename)
{
	int nv = _mesh.n_vertices();
	int nf = _mesh.n_faces();

	std::ofstream off_file(_filename);

	if (!off_file.is_open())
	{
		return false;
	}

	off_file << "OFF\n" << nv << " " << nf << " 0\n";
	off_file << std::setprecision(std::numeric_limits<Mesh::Scalar>::max_digits10);

	for (int i = 0; i < nv; i++)
	{
		auto&& p0 = _mesh.point(_mesh.vertex_handle(i));
		off_file << p0[0] << " " << p0[1] << " " << p0[2] << "\n";
	}

	for (int i = 0; i < nf; i++)
	{
		auto f_h = _mesh.face_handle(i);
		off_file << _mesh.valence(f_h);

		for (auto fh_h : _mesh.fh_range(_mesh.face_handle(i)))
		{
			off_file << " " << _mesh.to_vertex_handle(fh_h).idx();
		}

		off_file << "\n";
	}

	off_file.close();

	return true;
}

void Mesh_doubleIO::copy_mesh(const Mesh& src, Mesh& dst)
{
	dst.clear();

	for (int i = 0; i < src.n_vertices(); i++)
	{
		dst.add_vertex(src.point(src.vertex_handle(i)));
	}
	for (int i = 0; i < src.n_faces(); i++)
	{
		auto f_h = src.face_handle(i);
		std::vector<OpenMesh::VertexHandle> face_v;
		face_v.reserve(src.valence(f_h));

		for (auto fh_iter = src.cfh_begin(f_h); fh_iter.is_valid(); fh_iter++)
		{
			face_v.push_back(dst.vertex_handle(src.to_vertex_handle(*fh_iter).idx()));
		}

		dst.add_face(face_v);
	}
}
