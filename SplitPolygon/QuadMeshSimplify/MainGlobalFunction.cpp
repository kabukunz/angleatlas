#include "MainGlobalFunction.h"

double ANGLE_ADD_NEW = PI * 175.0 / 180.0;
std::string  getCurrentTimeStr()
{
	time_t t = time(NULL);
	char ch[64] = { 0 };
	strftime(ch, sizeof(ch) - 1, "%Y_%m_%d_%H_%M_%S", localtime(&t));     //年-月-日 时-分-秒
	return ch;
}

void TestKDTree()
{
	ANNpointArray apa = annAllocPts(4, 3);
	apa[0][0] = 0; apa[0][1] = 0; apa[0][2] = 0;
	apa[1][0] = 1; apa[1][1] = 0; apa[1][2] = 0;
	apa[2][0] = 1; apa[2][1] = 1; apa[2][2] = 0;
	apa[3][0] = 0; apa[3][1] = 1; apa[3][2] = 0;

	ANNkd_tree *kdTree = new ANNkd_tree(apa, 4, 3);
	for (int i = 0; i < 1000000000; ++i)
	{
		ANNpoint ap = annAllocPt(3);
		ap[0] = 0.9; ap[1] = 0.9; ap[2] = 0.9;
		ANNidxArray nnIdx = new ANNidx[4]; ANNdistArray dists = new ANNdist[4];
		kdTree->annkSearch(ap, 4, nnIdx, dists);
		annDeallocPt(ap);
		delete[] nnIdx;
		delete[] dists;
	}
}

void TransVec3d(std::vector<OpenMesh::Vec3d> &in, std::vector<Eigen::Vector3d> &out)
{
	out.clear();
	out.reserve(in.size());
	for (int i = 0; i < in.size(); ++i)
	{
		out.emplace_back(Eigen::Vector3d(in[i][0], in[i][1], in[i][2]));
	}
}

void BuildTriangleFiled(MyMesh &mesh, std::vector<OpenMesh::Vec3d> &fields, std::vector<OpenMesh::Vec3d> &midPoints)
{
	midPoints.clear();
	midPoints.reserve(mesh.n_faces());
	for (MyMesh::FaceIter fi = mesh.faces_begin(); fi != mesh.faces_end(); ++fi)
	{
		midPoints.emplace_back(mesh.calc_face_centroid(*fi));
	}

	std::vector<int> fTypes(mesh.n_faces(), 2);
	//std::vector<OpenMesh::Vec3d> fields(mesh.n_faces());
	fields.resize(mesh.n_faces());
	for (MyMesh::EdgeIter ei = mesh.edges_begin(); ei != mesh.edges_end(); ++ei)
	{
		if (mesh.is_boundary(*ei))
		{
			MyMesh::HalfedgeHandle notBHeh;
			if (mesh.is_boundary(mesh.halfedge_handle(*ei, 0)))
			{
				notBHeh = mesh.halfedge_handle(*ei, 1);
			}
			else
			{
				notBHeh = mesh.halfedge_handle(*ei, 0);
			}

			const OpenMesh::Vec3d &p0 = mesh.point(mesh.from_vertex_handle(notBHeh));
			const OpenMesh::Vec3d &p1 = mesh.point(mesh.to_vertex_handle(notBHeh));
			int fId = mesh.face_handle(notBHeh).idx();
			fields[fId] = (p1 - p0).normalized();
			fTypes[fId] = 1;
		}
	}

	CrossfieldDesign cfd(mesh, fTypes);
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
		fields[i].normalize();
	}
}


void BuildTriangleFiled(BaseDataStructure::QuadMesh &qm, std::vector<OpenMesh::Vec3d> &fields, std::vector<OpenMesh::Vec3d> &midPoints)
{
	Pipeline::RefineQuadToTri(qm, Pipeline::triMesh_);
	midPoints.clear();
	midPoints.reserve(Pipeline::triMesh_.n_faces());
	for (MyMesh::FaceIter fi = Pipeline::triMesh_.faces_begin(); fi != Pipeline::triMesh_.faces_end(); ++fi)
	{
		midPoints.emplace_back(Pipeline::triMesh_.calc_face_centroid(*fi));
	}

#if OUTPUT_TRI_AND_FIELDS
	OpenMesh::IO::write_mesh(Pipeline::triMesh_, "C:\\Users\\ChiZhang\\Desktop\\outTri.obj");
#endif

	std::vector<int> fTypes(Pipeline::triMesh_.n_faces(), 2);
	//std::vector<OpenMesh::Vec3d> fields(Pipeline::triMesh_.n_faces());
	fields.resize(Pipeline::triMesh_.n_faces());
	for (MyMesh::EdgeIter ei = Pipeline::triMesh_.edges_begin(); ei != Pipeline::triMesh_.edges_end(); ++ei)
	{
		if (Pipeline::triMesh_.is_boundary(*ei))
		{
			MyMesh::HalfedgeHandle notBHeh;
			if (Pipeline::triMesh_.is_boundary(Pipeline::triMesh_.halfedge_handle(*ei, 0)))
			{
				notBHeh = Pipeline::triMesh_.halfedge_handle(*ei, 1);
			}
			else
			{
				notBHeh = Pipeline::triMesh_.halfedge_handle(*ei, 0);
			}

			const OpenMesh::Vec3d &p0 = Pipeline::triMesh_.point(Pipeline::triMesh_.from_vertex_handle(notBHeh));
			const OpenMesh::Vec3d &p1 = Pipeline::triMesh_.point(Pipeline::triMesh_.to_vertex_handle(notBHeh));
			int fId = Pipeline::triMesh_.face_handle(notBHeh).idx();
			fields[fId] = (p1 - p0).normalized();
			fTypes[fId] = 1;
		}
	}

	CrossfieldDesign cfd(Pipeline::triMesh_, fTypes);
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
		fields[i].normalize();
	}

#if OUTPUT_TRI_AND_FIELDS
	std::ofstream ofs("C:\\Users\\ChiZhang\\Desktop\\fields.txt");
	ofs << fields.size() << std::endl;
	for (int i = 0; i < fields.size(); ++i)
	{
		const OpenMesh::Vec3d &midP = Pipeline::triMesh_.calc_face_centroid(Pipeline::triMesh_.face_handle(i));
		ofs << midP[0] << " " << midP[1] << " " << fields[i][0] << " " << fields[i][1] << std::endl;
	}

	exit(0);
#endif
}

void MeshRequest(MyMesh &mesh)
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_face_normals();
}

void RefineTri(MyMesh &inMesh, MyMesh &outMesh, int targetFaces)
{
	outMesh = inMesh;
	if (outMesh.n_faces() >= targetFaces)
		return;

	MyMesh tempMesh;
	do
	{
		tempMesh = outMesh;
		outMesh.clear();
		std::vector<MyMesh::VertexHandle> vHandles;
		vHandles.reserve(tempMesh.n_vertices() * 2);
		std::vector<MyMesh::VertexHandle> fHandles;
		fHandles.reserve(3);

		for (int i = 0; i < tempMesh.n_vertices(); ++i)
		{
			vHandles.emplace_back(outMesh.add_vertex(tempMesh.point(tempMesh.vertex_handle(i))));
		}

		int vh[3];
		for (int i = 0; i < tempMesh.n_faces(); ++i)
		{
			OpenMesh::Vec3d avePos(0, 0, 0);
			MyMesh::FaceHandle fh = tempMesh.face_handle(i);
			for (auto it : tempMesh.fv_range(fh))
			{
				avePos += tempMesh.point(it);
			}
			avePos /= 3;

			MyMesh::VertexHandle tempV = outMesh.add_vertex(MyMesh::Point(avePos[0], avePos[1], avePos[2]));

			MyMesh::FaceVertexCCWIter fvcci = tempMesh.fv_ccwiter(fh);
			vh[0] = fvcci->idx();
			++fvcci; vh[1] = fvcci->idx();
			++fvcci; vh[2] = fvcci->idx();

			for (int j = 0; j < 3; ++j)
			{
				fHandles.clear();
				fHandles.emplace_back(tempV);
				fHandles.emplace_back(vHandles[vh[j]]);
				fHandles.emplace_back(vHandles[vh[(j + 1) % 3]]);
				outMesh.add_face(fHandles);
			}
		}
	} while (outMesh.n_faces() < targetFaces);

	MeshRequest(outMesh);
}

void ReadMeshBuildField(const std::string &fileName, std::vector<Eigen::Vector3d> &fields, std::vector<OpenMesh::Vec3d> &midPoints)
{
	//MyMesh sourceMesh;
	Pipeline::triMesh_.clear();
	OpenMesh::IO::read_mesh(Pipeline::triMesh_, fileName.c_str());
	MeshRequest(Pipeline::triMesh_);
	MyMesh tempMesh;
	RefineTri(Pipeline::triMesh_, tempMesh, 15000);

	//OpenMesh::IO::write_mesh(tempMesh, "C:\\Users\\ChiZhang\\Desktop\\outTri.obj");

	std::vector<OpenMesh::Vec3d> tempFields;
	BuildTriangleFiled(tempMesh, tempFields, midPoints);

	TransVec3d(tempFields, fields);

	/*std::ofstream ofs("C:\\Users\\ChiZhang\\Desktop\\fields.txt");
	ofs << fields.size() << std::endl;
	for (int i = 0; i < fields.size(); ++i)
	{
		const OpenMesh::Vec3d &midP = tempMesh.calc_face_centroid(tempMesh.face_handle(i));
		ofs << midP[0] << " " << midP[1] << " " << fields[i][0] << " " << fields[i][1] << std::endl;
	}*/
}

void OutputTriAndField(std::vector<Eigen::Vector3d> &fields, std::vector<OpenMesh::Vec3d> &midPoints)
{
	std::ofstream ofs("C:\\Users\\ChiZhang\\Desktop\\fields.txt");
	ofs << fields.size() << std::endl;
	for (int i = 0; i < fields.size(); ++i)
	{
		const OpenMesh::Vec3d &midP = Pipeline::triMesh_.calc_face_centroid(Pipeline::triMesh_.face_handle(i));
		ofs << midP[0] << " " << midP[1] << " " << fields[i][0] << " " << fields[i][1] << std::endl;
	}
}

Eigen::Vector2d OpenVec3ToEigenVec2(const OpenMesh::Vec3d &vec)
{
	return Eigen::Vector2d(vec[0], vec[1]);
}

void ComputePolyAddNew(std::vector<Eigen::Vector2d> &polys, std::vector<bool> &polysAddNew, double threAngleRad)
{
	polysAddNew.clear();
	double threCosValue = std::cos(threAngleRad);
	for (int i = 0; i < polys.size(); ++i)
	{
		Eigen::Vector2d vPos0 = polys[PolySplit::mmod((int)(i - 1), (int)polys.size())], vPos1 = polys[i], vPos2 = polys[(i + 1) % polys.size()];
		Eigen::Vector2d vec0 = (vPos0 - vPos1).normalized(), vec1 = (vPos2 - vPos1).normalized();
		double currCosValue = vec0.dot(vec1);
		if (currCosValue <= threCosValue)
			polysAddNew.emplace_back(true);
		else
			polysAddNew.emplace_back(false);
	}
}

int GetNextCCWVertex(MyMesh &mesh, MyMesh::VertexHandle &vh)
{
	for (MyMesh::VertexIHalfedgeCCWIter vihcci = mesh.vih_ccwiter(vh); vihcci.is_valid(); ++vihcci)
	{
		if (mesh.is_boundary(*vihcci))
		{
			return mesh.from_vertex_handle(*vihcci).idx();
		}
	}
}

void GetTriMeshPoly(MyMesh &mesh, std::vector<Eigen::Vector2d> &poly, std::vector<bool> &polyAddNew, std::vector<uint32_t> &polyToQuad)
{
	/*std::vector<uint32_t> polyVIds;
	polyVIds.reserve(mesh.n_vertices());*/
	poly.clear();
	poly.reserve(mesh.n_vertices());
	polyToQuad.clear();
	polyToQuad.reserve(mesh.n_vertices());
	polyAddNew.clear();
	polyAddNew.reserve(mesh.n_vertices());

	uint32_t startVId = (uint32_t)-1;
	for (int i = 0; i < mesh.n_vertices(); ++i)
	{
		MyMesh::VertexHandle vh = mesh.vertex_handle(i);
		if (mesh.is_boundary(vh))
		{
			startVId = vh.idx();
			break;
		}
	}

	MyMesh::VertexHandle startVh = mesh.vertex_handle(startVId);
	poly.emplace_back(OpenVec3ToEigenVec2(mesh.point(startVh)));
	polyToQuad.emplace_back(startVh.idx());
	int nextV = GetNextCCWVertex(mesh, startVh);
	MyMesh::VertexHandle nextVh;
	while (nextV != startVId)
	{
		nextVh = mesh.vertex_handle(nextV);
		poly.emplace_back(OpenVec3ToEigenVec2(mesh.point(nextVh)));
		polyToQuad.emplace_back(nextVh.idx());
		nextV = GetNextCCWVertex(mesh, nextVh);
	}

	ComputePolyAddNew(poly, polyAddNew, ANGLE_ADD_NEW);
}

double Cross2(const Eigen::Vector2d &vec0, const Eigen::Vector2d &vec1)
{
	return vec0[0] * vec1[1] - vec0[1] * vec1[0];
}

bool JudgeIfEdgeIntersection(Eigen::Vector2d &pos0, Eigen::Vector2d &pos1, Eigen::Vector2d &pos2, Eigen::Vector2d &pos3)
{
	double d0 = Cross2(pos1 - pos0, pos2 - pos0), d1 = Cross2(pos1 - pos0, pos3 - pos0), d2 = Cross2(pos3 - pos2, pos0 - pos2), d3 = Cross2(pos3 - pos2, pos1 - pos2);

	if (d0*d1 < 0 && d2*d3 < 0)
		return true;
	else
		return false;
}

void GetTriMeshPolyVersion2(MyMesh &mesh, std::vector<Eigen::Vector2d> &poly, std::vector<bool> &polyAddNew)
{
	poly.clear();
	poly.reserve(mesh.n_vertices());
	polyAddNew.clear();
	polyAddNew.reserve(mesh.n_vertices());

	uint32_t startVId = (uint32_t)-1;
	for (int i = 0; i < mesh.n_vertices(); ++i)
	{
		MyMesh::VertexHandle vh = mesh.vertex_handle(i);
		if (mesh.is_boundary(vh))
		{
			startVId = vh.idx();
			break;
		}
	}

	MyMesh::VertexHandle startVh = mesh.vertex_handle(startVId);
	poly.emplace_back(OpenVec3ToEigenVec2(mesh.point(startVh)));
	int nextV = GetNextCCWVertex(mesh, startVh);
	MyMesh::VertexHandle nextVh;
	while (nextV != startVId)
	{
		nextVh = mesh.vertex_handle(nextV);
		poly.emplace_back(OpenVec3ToEigenVec2(mesh.point(nextVh)));
		nextV = GetNextCCWVertex(mesh, nextVh);
	}

	std::vector<int> vsFlag(poly.size(), 0);
	//ComputePolyAddNew(poly, polyAddNew, ANGLE_ADD_NEW);
	std::vector<uint32_t> oldVIds;
	oldVIds.reserve(poly.size());
	double threCosValue = std::cos(PI * 170.0 / 180.0);
	for (int i = 0; i < poly.size(); ++i)
	{
		Eigen::Vector2d &vPos0 = poly[PolySplit::mmod((int)(i - 1), (int)poly.size())], &vPos1 = poly[i], &vPos2 = poly[(i + 1) % poly.size()];
		Eigen::Vector2d vec0 = (vPos0 - vPos1).normalized(), vec1 = (vPos2 - vPos1).normalized();
		double currCosValue = vec0.dot(vec1);
		if (currCosValue > threCosValue)
		{
			oldVIds.emplace_back(i);
			vsFlag[i] = 1;
		}
	}
	if (oldVIds.empty())
		return;

	//std::vector<bool> backupPvsAddNew = pvsAddNew;
	for (int j = 0; j < oldVIds.size(); ++j)
	{
		for (int k = j + 2; k < oldVIds.size(); ++k)
		{
			if (j == k || j == ((k + 1) % oldVIds.size()) || ((j + 1) % oldVIds.size()) == k || ((j + 1) % oldVIds.size()) == ((k + 1) % oldVIds.size()))
				continue;
			uint32_t idJ = oldVIds[j], idJ1 = oldVIds[(j + 1) % oldVIds.size()], idK = oldVIds[k], idK1 = oldVIds[(k + 1) % oldVIds.size()];
			if (JudgeIfEdgeIntersection(poly[idJ], poly[idJ1], poly[idK], poly[idK1]))
			{
				int currId = idJ;
				while (currId != idJ1)
					vsFlag[currId] = 2;

				currId = idK;
				while (currId != idK1)
					vsFlag[currId] = 2;
			}
		}
	}
	for (int i = 0; i < oldVIds.size(); ++i)
	{
		vsFlag[oldVIds[i]] = 1;
	}

	std::vector<Eigen::Vector2d> tempPoly;
	tempPoly.reserve(poly.size());
	for (int i = 0; i < vsFlag.size(); ++i)
	{
		if (vsFlag[i] == 0)
			continue;

		tempPoly.emplace_back(poly[i]);
		if (vsFlag[i] == 1)
			polyAddNew.emplace_back(false);
		else
			polyAddNew.emplace_back(true);
	}
	poly.swap(tempPoly);
}