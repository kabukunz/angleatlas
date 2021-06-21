#include "MeshDefinition.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <Eigen/Dense>

bool is_flip_ok_openmesh(MyMesh::EdgeHandle& eh, MyMesh& mesh_)
{
	// boundary edges cannot be flipped
	if ( mesh_.is_boundary(eh) ) return false;

	MyMesh::HalfedgeHandle hh = mesh_.halfedge_handle(eh, 0);
	MyMesh::HalfedgeHandle oh = mesh_.halfedge_handle(eh, 1);

	// check if the flipped edge is already present
	// in the mesh

	MyMesh::VertexHandle ah = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(hh));
	MyMesh::VertexHandle bh = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(oh));

	if(ah == bh)   // this is generally a bad sign !!!
		return false;

	for (MyMesh::ConstVertexVertexIter vvi = mesh_.vv_iter(ah); vvi.is_valid(); ++vvi)
		if( *vvi == bh )
			return false;

	return true;
}

bool flip_openmesh(MyMesh::EdgeHandle& eh, MyMesh& mesh_)
{
	// CAUTION : Flipping a halfedge may result in
	// a non-manifold mesh, hence check for yourself
	// whether this operation is allowed or not!
	if( !is_flip_ok_openmesh(eh, mesh_) )
		return false;//let's make it sure it is actually checked
	//assert( is_flip_ok_openmesh(eh, mesh_ ) );

	MyMesh::HalfedgeHandle a0 = mesh_.halfedge_handle(eh, 0);
	MyMesh::HalfedgeHandle b0 = mesh_.halfedge_handle(eh, 1);

	MyMesh::HalfedgeHandle a1 = mesh_.next_halfedge_handle(a0);
	MyMesh::HalfedgeHandle a2 = mesh_.next_halfedge_handle(a1);

	MyMesh::HalfedgeHandle b1 = mesh_.next_halfedge_handle(b0);
	MyMesh::HalfedgeHandle b2 = mesh_.next_halfedge_handle(b1);

	MyMesh::VertexHandle   va0 = mesh_.to_vertex_handle(a0);
	MyMesh::VertexHandle   va1 = mesh_.to_vertex_handle(a1);

	MyMesh::VertexHandle   vb0 = mesh_.to_vertex_handle(b0);
	MyMesh::VertexHandle   vb1 = mesh_.to_vertex_handle(b1);

	MyMesh::FaceHandle     fa  = mesh_.face_handle(a0);
	MyMesh::FaceHandle     fb  = mesh_.face_handle(b0);

	mesh_.set_vertex_handle(a0, va1);
	mesh_.set_vertex_handle(b0, vb1);

	mesh_.set_next_halfedge_handle(a0, a2);
	mesh_.set_next_halfedge_handle(a2, b1);
	mesh_.set_next_halfedge_handle(b1, a0);

	mesh_.set_next_halfedge_handle(b0, b2);
	mesh_.set_next_halfedge_handle(b2, a1);
	mesh_.set_next_halfedge_handle(a1, b0);

	mesh_.set_face_handle(a1, fb);
	mesh_.set_face_handle(b1, fa);

	mesh_.set_halfedge_handle(fa, a0);
	mesh_.set_halfedge_handle(fb, b0);

	if(mesh_.halfedge_handle(va0) == b0)
		mesh_.set_halfedge_handle(va0, a1);
	if(mesh_.halfedge_handle(vb0) == a0)
		mesh_.set_halfedge_handle(vb0, b1);

	return true;
}

bool check_in_triangle_face(const std::vector<OpenMesh::Vec3d>& tri, const OpenMesh::Vec3d& p)
{
	OpenMesh::Vec3d v1 = tri[1] - tri[0]; OpenMesh::Vec3d v2 = tri[2] - tri[0];
	OpenMesh::Vec3d n = OpenMesh::cross(v1, v2);
	double face_area = n.norm(); n.normalize(); double all_area = 0;
	for(unsigned i=0; i < tri.size(); ++i)
	{
		unsigned next_i = ( i+1 )%tri.size(); unsigned prev_i = ( i + tri.size() - 1 )%tri.size();
		v1 = tri[next_i] - p; v2 = tri[prev_i] - p;
		double area = OpenMesh::dot(OpenMesh::cross(v1, v2), n); all_area += area;
		if(area < 0)
		{
			return false;
		}
	}
	if(std::abs(all_area - face_area) < 1e-8) {return true;}
	else {return false;}
}
