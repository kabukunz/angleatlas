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

		//bool keeped = true;	//���collapse֮��Ϊ�ڲ����Ҷ�Ϊ2����ǽǵ��Ϊ1����������ɵ���
		
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
		uint32_t currEId;	//��ǰ����ʹ�õ�orderEs��id
	}SingularNode;

	typedef struct SplitVertexInfo
	{
		uint32_t oldV = (uint32_t)-1, newV0 = (uint32_t)-1, newV1 = (uint32_t)-1;
		std::vector<uint32_t> fsNP0, fsNP1;	//������edgesһ������id�� ����һ���id
		std::vector<uint32_t> esNP0, esNP1;	//��������edges֮��ģ���oldV���������бߡ�
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
		bool FindSplitLine(uint32_t startV, uint32_t startE, std::vector<uint32_t> &es, std::vector<uint32_t> &vs, std::vector<OrderEFs> &oefs);	//oefs�߽��û��Ԫ��
		void FindOrderEdges(uint32_t startE, uint32_t currV, std::vector<uint32_t> &orderEs, std::vector<uint32_t> &orderFs);	//ֻ���ҷǱ߽��
		void FindOrderEdges(uint32_t startE, uint32_t startF, uint32_t currV, std::vector<uint32_t> &orderEs, std::vector<uint32_t> &orderFs);
		//��������false��ʾ�������Ѿ��ҵ��ĵ㣬����true��ʾ�ҵ�����һ���������߽߱�㣬���У�isBoundaryΪtrue��ʾ���˱߽磬����ǵ�����һ�������
		bool FindSubLine(uint32_t startE, uint32_t startV, std::vector<uint32_t> &es, std::vector<uint32_t> &vs, bool *vFlag, bool &isBoundary, std::vector<OrderEFs> &oefs);
		void ReorderEs(std::vector<uint32_t> &es);	//�������������м�Ŀ�ǰ�档
		bool FindNewSubStart(std::stack<SingularNode> &snStack, uint32_t &newSE, uint32_t &newSV);	//������п��ܵ���������ҹ��˵Ļ���return false, ����return true.

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

		uint32_t tempEF, otherV;	//��������ʱ������������
		std::vector<uint32_t> tempEs, tempVs; 
		std::vector<OrderEFs> tempOefs;

		double radio_ = 0;
		Eigen::Vector2d minP_;
	};
}

