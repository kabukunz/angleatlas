#pragma once
#include <vector>
#include <Eigen/Dense>

namespace BaseDataStructure
{
	class QuadMesh;

	typedef struct FrameVertex 
	{
		uint32_t id, qId; //id is the base complex id, wile qId is the id in source quad mesh.	
		std::vector<uint32_t> b_neighbor_vs;
		std::vector<uint32_t> b_neighbor_es;
		std::vector<uint32_t> b_neighbor_fs;
		bool boundary = false;
		bool singular =true;	
	}FrameVertex;

	typedef struct FrameEdge
	{
		uint32_t id;
		std::vector<uint32_t> b_vs;
		std::vector<uint32_t> b_neighbor_es;
		std::vector<uint32_t> b_neighbor_fs;
		std::vector<uint32_t> vs_link;	//source quad vertices contained.
		std::vector<uint32_t> es_link;	//source quad edges contained.
		bool boundary = false;
	}FrameEdge;

	typedef struct FrameFace
	{
		uint32_t id;
		std::vector<uint32_t> b_vs;
		std::vector<uint32_t> b_es;
		std::vector<uint32_t> b_neighbor_fs;
		
		std::vector<uint32_t> vs_net; //vertices contained in this frame face.
		std::vector<uint32_t> fs_net; //faces contained in this frame face.
		bool boundary = false;
	}FrameFace;

	enum SheetType
	{
		OPEN,
		CLOSED,
		INTERSECTION
	};
	//enum SheetBoundaryType	//sheet 的两侧的边界情况
	//{
	//	NONE,
	//	ONE_SIDE,
	//	TWO_SIDES
	//};
	typedef struct Sheet
	{
		uint32_t id;
		SheetType st;
		//SheetBoundaryType sbt;
		std::vector<uint32_t> b_vs;
		std::vector<uint32_t> b_es;
		std::vector<uint32_t> b_fs;
		std::vector<uint32_t> b_middle_es, b_middle_boundary_es;//, left_es, right_es;
		std::vector<uint32_t> left_vs, right_vs;	//和b_middle_es对应，因此可能有重复。
		std::vector<uint32_t> qf_delete;

		std::vector<std::vector<uint32_t>> vs_pairs, vs_links, vs_group;
		std::vector<uint32_t> target_vs;
		std::vector<Eigen::Vector2d> target_coords;

		double weight;
		bool valenceFilter;
		bool haveAddValenceWidget = false;
	}Sheet; 

	typedef struct Chord
	{
		uint32_t id;
		uint32_t fid; //this is the fid of chord, the same as base complex face.
		std::vector<uint32_t> qf_delete;
		std::vector<std::vector<uint32_t>> vs_group;
		std::vector<uint32_t> target_vs;
		std::vector<Eigen::Vector2d> target_coords;

		uint32_t keepedV0, keepedV1;	//base complex id
		uint32_t removedV0, removedV1;	//base complex id

		uint32_t keepedQV0, keepedQV1;

		double weight;
		bool valenceFilter;
		bool haveAddValenceWidget = false;
	}Chord;

	class BaseComplex
	{
	public:
		BaseComplex();
		~BaseComplex();

		void ExtractBaseComplex(QuadMesh *quadMesh);
		void ExtractBaseComplex(QuadMesh *quadMesh, double angleThre);
		void ExtractBaseComplex(QuadMesh *quadMesh, std::vector<uint32_t> &corners);
		void ExtractBCVerticesEdges(QuadMesh *quadMesh);
		void ExtractBCVerticesEdges(QuadMesh *quadMesh, double angleThre);
		void ExtractBCVerticesEdges(QuadMesh *quadMesh, std::vector<uint32_t> &corners);
		void ExtractBCFaces(QuadMesh *quadMesh);
		void ClearBaseComplex();

		void OutputBCEdges(QuadMesh *quadMesh, std::string ofileName);

	public:
		std::vector<FrameVertex> Bv_;
		std::vector<FrameEdge> Be_;
		std::vector<FrameFace> Bf_;

		bool isError_ = false;

		uint32_t tempVs_[2];
	};

}