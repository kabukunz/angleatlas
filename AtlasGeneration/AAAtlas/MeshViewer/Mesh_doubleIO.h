#pragma once

#include "MeshDefinition.h"

class Mesh_doubleIO
{
public:
	static bool load_mesh(Mesh& _mesh, const char* _filename, bool load_texture = false);
	static bool save_mesh(const Mesh& _mesh, const char* _filename, bool save_texture = false);

	static bool save_uv_mesh(const Mesh& _mesh, const char* _filename);

	enum class file_type
	{
		others, obj, off
	};

	static file_type get_file_type(const char* _filename);

	//edge indices may change
	static void copy_mesh(const Mesh& src, Mesh& dst);

private:
	static bool load_obj(Mesh& _mesh, const char* _filename, bool load_texture);
	static bool load_off(Mesh& _mesh, const char* _filename);

	static bool save_obj(const Mesh& _mesh, const char* _filename, bool save_texture);
	static bool save_off(const Mesh& _mesh, const char* _filename);
};