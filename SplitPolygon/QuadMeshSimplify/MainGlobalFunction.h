#pragma once
#include <iostream>
#include <Eigen/Dense>
#include "ExtraFunc.h"
#include "QuadMesh.h"
#include "QuadMeshIO.h"
#include "BaseComplex.h"
#include "Pipeline.h"
#include "AssistFunc.h"
#include "CrossfieldDesign.h"
#include "GeoPolyReadWrite.h"
#include "PolySplit.h"
#include <ctime>
#include <string>
#include <fstream>
#include <windows.h>
#include <OpenMesh/Core/IO/MeshIO.hh>

using namespace BaseDataStructure;
using namespace DataOperation;


#define OUTPUT_TRI_AND_FIELDS 0
#define TEST_SPLIT 1

static std::string  getCurrentTimeStr();

void TestKDTree();

void TransVec3d(std::vector<OpenMesh::Vec3d> &in, std::vector<Eigen::Vector3d> &out);

void BuildTriangleFiled(BaseDataStructure::QuadMesh &qm, std::vector<OpenMesh::Vec3d> &fields, std::vector<OpenMesh::Vec3d> &midPoints);

void MeshRequest(MyMesh &mesh);

void RefineTri(MyMesh &inMesh, MyMesh &outMesh, int targetFaces);

void BuildTriangleFiled(MyMesh &mesh, std::vector<OpenMesh::Vec3d> &fields, std::vector<OpenMesh::Vec3d> &midPoints);

void ReadMeshBuildField(const std::string &fileName, std::vector<Eigen::Vector3d> &fields, std::vector<OpenMesh::Vec3d> &midPoints);

void OutputTriAndField(std::vector<Eigen::Vector3d> &fields, std::vector<OpenMesh::Vec3d> &midPoints);

void GetTriMeshPoly(MyMesh &mesh, std::vector<Eigen::Vector2d> &poly, std::vector<bool> &polyAddNew, std::vector<uint32_t> &polyToQuad);

void GetTriMeshPolyVersion2(MyMesh &mesh, std::vector<Eigen::Vector2d> &poly, std::vector<bool> &polyAddNew);

Eigen::Vector2d OpenVec3ToEigenVec2(const OpenMesh::Vec3d &vec);

void ComputePolyAddNew(std::vector<Eigen::Vector2d> &polys, std::vector<bool> &polysAddNew, double threAngleRad);

int GetNextCCWVertex(MyMesh &mesh, MyMesh::VertexHandle &vh);

bool JudgeIfEdgeIntersection(Eigen::Vector2d &pos0, Eigen::Vector2d &pos1, Eigen::Vector2d &pos2, Eigen::Vector2d &pos3);

double Cross2(const Eigen::Vector2d &vec0, const Eigen::Vector2d &vec1);