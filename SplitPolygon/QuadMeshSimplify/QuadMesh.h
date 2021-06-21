#pragma once
#include <Eigen/Dense>
#include <vector>
#include <stack>
namespace BaseDataStructure
{
	const uint32_t quadTriTable[4][3] = { {0,1,2}, {1,2,3}, {2,3,0}, {3,0,1} };
	typedef struct QuadVertex
	{
		uint32_t id, bId;//bId is the id in base complex, =(uint32_t)-1 means it is not a base complex node.	//ok
		std::vector<uint32_t> neighbor_vs;
		std::vector<uint32_t> neighbor_es;
		std::vector<uint32_t> neighbor_fs;
		bool boundary = false;	

		//bool keeped = true;	//如果collapse之后为内部点且度为2，或非角点度为1，则把这个点干掉。
		
	}QuadVertex;

	typedef struct QuadEdge
	{
		uint32_t id, bId;
		std::vector<uint32_t> vs;
		std::vector<uint32_t> neighbor_es;
		std::vector<uint32_t> neighbor_fs;
		bool boundary = false;

	}QuadEdge;

	typedef struct QuadFace
	{
		uint32_t id, bId;
		std::vector<uint32_t> vs;
		std::vector<uint32_t> es;
		std::vector<uint32_t> neighbor_fs;	//faces which have the same edge
		bool boundary = false;
	}QuadFace;

	typedef struct SingularNode
	{
		uint32_t currV = (uint32_t)-1, vValence = (uint32_t)-1, currE = (uint32_t)-1;
		std::vector<uint32_t> orderEs;
		uint32_t currEId;	//当前的所使用的orderEs的id
	}SingularNode;

	typedef struct SplitVertexInfo
	{
		uint32_t oldV = (uint32_t)-1, newV0 = (uint32_t)-1, newV1 = (uint32_t)-1;
		std::vector<uint32_t> fsNP0, fsNP1;	//在两条edges一侧的面的id， 在另一侧的id
		std::vector<uint32_t> esNP0, esNP1;	//除掉两条edges之外的，和oldV关联的所有边。
	}SplitVertexInfo;

	typedef struct OrderEFs
	{
		std::vector<uint32_t> oEs;
		std::vector<uint32_t> oFs;
	}OrderEFs;

	class QuadMesh
	{
	public:
		QuadMesh();
		~QuadMesh();


		void ComputeQuadMeshBBox(Eigen::Vector2d &minP, Eigen::Vector2d &maxP);
		void ScaleQuadMesh(double targetBBox);
		void ScaleQuadMesh(double targetLength, double &radio, Eigen::Vector2d &minP);
		void ScaleQuadMeshBack(double radio, Eigen::Vector2d &minP);
		bool BuildConnectivity();	//After reading obj file, build the topology information.
		void OritateFace(); //Make the normal of faces to be (0,0,1);
		bool JudgeIfTopologyProblem() const;
		void ClearMesh();
		void ComputeMeshNormal();
		double ComputeAverageEdgeLength();

		void OrderElements();

		bool JudgeQuadMeshJacobi();
		bool JudgeQuadMeshJacobiWithoutBoundary();

		bool SplitBadLine(uint32_t badV, uint32_t badF, uint32_t startV, uint32_t startE, std::vector<uint32_t> &outFs);
		bool FindSplitLine(uint32_t startV, uint32_t startE, std::vector<uint32_t> &es, std::vector<uint32_t> &vs, std::vector<OrderEFs> &oefs);	//oefs边界点没有元素
		void FindOrderEdges(uint32_t startE, uint32_t currV, std::vector<uint32_t> &orderEs, std::vector<uint32_t> &orderFs);	//只能找非边界点
		void FindOrderEdges(uint32_t startE, uint32_t startF, uint32_t currV, std::vector<uint32_t> &orderEs, std::vector<uint32_t> &orderFs);
		//函数返回false表示碰到了已经找到的点，返回true表示找到了下一个奇异点或者边界点，其中，isBoundary为true表示到了边界，否就是到了下一个奇异点
		bool FindSubLine(uint32_t startE, uint32_t startV, std::vector<uint32_t> &es, std::vector<uint32_t> &vs, bool *vFlag, bool &isBoundary, std::vector<OrderEFs> &oefs);
		void ReorderEs(std::vector<uint32_t> &es);	//给边重新排序，中间的靠前面。
		bool FindNewSubStart(std::stack<SingularNode> &snStack, uint32_t &newSE, uint32_t &newSV);	//如果所有可能的情况都被找过了的话，return false, 否则return true.

		int FindVInFId(uint32_t vId, uint32_t fId);
		template <typename T>
		static bool JudgeSign(T tValue0, T tValue1)
		{
			if ( (tValue0 >= 0&& tValue1>=0) || (tValue0<=0&&tValue1<=0) )
				return true;
			else
				return false;
		}

	public:
		//Eigen::MatrixXd V;
		std::vector<Eigen::Vector2d> V_;
		std::vector<QuadVertex> Vs_;
		std::vector<QuadEdge> Es_;
		std::vector<QuadFace> Fs_;

		double minJacobi_ = 0;

		Eigen::Vector3d meshNormal_;

		std::stack<SingularNode> snStack_;

		uint32_t tempEF, otherV;	//不想用临时变量而设立的
		std::vector<uint32_t> tempEs, tempVs; 
		std::vector<OrderEFs> tempOefs;

		double radio_ = 0;
		Eigen::Vector2d minP_;
	};
}

