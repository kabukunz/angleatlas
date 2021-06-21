#pragma once

#include <Eigen/Sparse>
#include <vector>
#include <tuple>
#include <ANN/ANN.h>
#include <list>
#include "CGAL_AABBTree.h"
class PolySplit
{
public:
	PolySplit();
	~PolySplit();

	typedef std::tuple<double, uint32_t, uint32_t, uint32_t, uint32_t, bool> TupleRankingLine;	//从前到后5个元素的含义分别是：与场的平均角度，v0，v1，v0在当前poly中的下标，v1在当前poly中的下标，是否两个点都是凹点。
	typedef struct PolyVertex
	{
		uint32_t id = (uint32_t)-1;
		Eigen::Vector2d pos;

		std::vector<uint32_t> neighbor_pes;
		std::vector<uint32_t> neighbor_pfs;

		bool isBoundary = false;

	}PolyVertex;

	typedef struct PolyEdge
	{
		uint32_t id = (uint32_t)-1;
		uint32_t v0 = (uint32_t)-1, v1 = (uint32_t)-1;
		uint32_t neighbor_f0 = (uint32_t)-1, neighbor_f1 = (uint32_t)-1;

		bool isBoundary = true;
	}PolyEdge;

	typedef struct PolyFace
	{
		uint32_t id = (uint32_t)-1;
		/*std::vector<uint32_t> pvs_new_add;
		std::vector<uint32_t> pvs;*/
		std::list<uint32_t> pvs;
		std::list<bool> pvs_new_add;
		std::vector<uint32_t> pes;
		std::vector<uint32_t> neighbor_pfs;

		bool isBoundary = false;
	}PolyFace;

	void InitializePoly(std::vector<Eigen::Vector2d> &inPoly);
	void InitializePoly(std::vector<Eigen::Vector2d> &inPoly, std::vector<bool> &inPolyAddNew);
	void SplitPoly(std::vector<Eigen::Vector2d> &inPoly);
	void SplitPoly(std::vector<Eigen::Vector2d> &inPoly, std::vector<bool> &inPolyAddNew);
	void SplitPoly(std::vector<Eigen::Vector2d> &inPoly, std::vector<bool> &inPolyAddNew, int polyNumThre);
	void SplitPoly(std::vector<Eigen::Vector2d> &inPoly, std::vector<bool> &inPolyAddNew, int polyNumThre, int resultNum);
	void OutputPolysWithoutSelfIntersection();
	void OutputPolysWithoutSelfIntersection(std::vector<Eigen::Vector2d> &inPoly, std::vector<bool> &inPolyAddNew);

	void SetKdTree(ANNkd_tree *kdTree)
	{
		fieldKdtree_ = kdTree;
	}
	void SetFields(std::vector<Eigen::Vector3d> &fields)
	{
		fields_ = &fields;
	}

	bool JudgePolyCCW(std::vector<uint32_t> &polyVs);
	bool JudgePolyCCW(std::vector<Eigen::Vector2d> &polyVPos);
	void GetConcaveVertices(std::vector<uint32_t> &polyVs, std::vector<uint32_t> &ccVsIds);
	void GetConcaveVerticesLarger(std::vector<uint32_t> &polyVs, std::vector<uint32_t> &ccVsIds);
	void GetConcaveVerticesLarger(std::vector<uint32_t> &polyVs, std::vector<uint32_t> &ccVsIds, std::vector<bool> &polyAddNew);
	double ComputePolyArea(std::vector<uint32_t> &polyVs);
	double ComputePolyArea(std::vector<Eigen::Vector2d> &polyPos);

	bool CheckPointInsidePoly(std::vector<uint32_t> &poly, Eigen::Vector2d &midPos);
	int CheckPointPolyState(std::vector<uint32_t> &poly, Eigen::Vector2d &pos);	//0:out; 1:on boundary; 2: in
	bool CheckPointOnBoundary(Eigen::Vector2d &pos);
	bool JudgeIfEdgeIntersectionOrOutside(std::vector<uint32_t> &polyVs, uint32_t v0, uint32_t v1);
	bool JudgeIfEdgeIntersection(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3);	//v0,v1是一条线段的两个端点，v2,v3是另一条。
	bool JudgeIfEdgeIntersectionLarger(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3);
	bool JudgeIfEdgeIntersection(Eigen::Vector2d &pos0, Eigen::Vector2d &pos1, Eigen::Vector2d &pos2, Eigen::Vector2d &pos3);

	bool JudgeIfOnSameLine(std::vector<bool> &polyAddNew, uint32_t v0Id, uint32_t v1Id);
	bool JudgeIfOnSameLine(std::vector<uint32_t> &poly, uint32_t currId, Eigen::Vector2d &interPos);

	static double Cross2(const Eigen::Vector2d &vec0, const Eigen::Vector2d &vec1)
	{
		return vec0[0] * vec1[1] - vec0[1] * vec1[0];
	}

	double LineFieldMatchAngle(uint32_t v0, uint32_t v1);
	void LineSample(uint32_t v0, uint32_t v1, std::vector<Eigen::Vector2d> &outSamples);
	double FieldMatchAngle(uint32_t v0, uint32_t v1, const Eigen::Vector2d &field);

	void SplitPolyIntoTwo(std::vector<uint32_t> &polys, uint32_t v0Id, uint32_t v1Id, std::vector<uint32_t> &newPoly0, std::vector<uint32_t> &newPoly1);
	void SplitPolyIntoTwo(std::vector<uint32_t> &polys, std::vector<bool> &polyAddNew, uint32_t v0Id, uint32_t v1Id, std::vector<uint32_t> &newPoly0, std::vector<uint32_t> &newPoly1, std::vector<bool> &newPAN0, std::vector<bool> &newPAN1);

	double ComputePolyAngle(std::vector<uint32_t> &polys, int vId);
	double ComputePolyAngle(std::vector<Eigen::Vector2d> &polys, int vId);

	void BuildEFsInfo();
	void BuildTopologyInfoFromVsFs();

	void SetScale(double radio, Eigen::Vector2d &minP);

	void TranslateToEigenPoly(std::vector<std::vector<Eigen::Vector2d>> &polys);
	void TranslateToEigenPolyFromFs(std::vector<std::vector<Eigen::Vector2d>> &polys);

	int GetFalseNum(std::vector<bool> &vecB)
	{
		int resu = 0;
		for (int i = 0; i < vecB.size(); ++i)
		{
			if (!vecB[i])
				++resu;
		}
		return resu;
	}

	template<typename T>
	static T mmod(T n0, T n1)
	{
		return std::fmod(std::fmod(n0, n1) + n1, n1);
	}

	void FindLinePolyNearIntersection(Eigen::Vector2d &vPos0, Eigen::Vector2d &vPos1, std::vector<Eigen::Vector2d> &poly, double &length0, double &length1, uint32_t &interId0, uint32_t &interId1, Eigen::Vector2d &interPos);
	bool CalcLinesIntersection(Eigen::Vector2d &vPos0, Eigen::Vector2d &vPos1, Eigen::Vector2d &vPos2, Eigen::Vector2d &vPos3, double &length0, double &length1,  Eigen::Vector2d &interPos); //return false means parallel
	bool CalcLinesIntersection2(Eigen::Vector2d &vPos0, Eigen::Vector2d &vPos1, Eigen::Vector2d &vPos2, Eigen::Vector2d &vec1, double &length0, double &length1, Eigen::Vector2d &interPos);

	bool GetPolyRayIntersection(std::vector<uint32_t> &poly, uint32_t v0Id, uint32_t v1Id, uint32_t &outV2Id, uint32_t &outV3Id, Eigen::Vector2d &interPos, std::vector<Eigen::Vector2d> &outPoly0, std::vector<Eigen::Vector2d> &outPoly1);	//以v0为顶点，沿着v1减v0的负方向往外延伸的射线与poly的交点
	bool GetPolyRayIntersection(std::vector<uint32_t> &poly, uint32_t v0Id, Eigen::Vector2d &vecDir, uint32_t &outV2Id, uint32_t &outV3Id, Eigen::Vector2d &interPos, std::vector<Eigen::Vector2d> &outPoly0, std::vector<Eigen::Vector2d> &outPoly1);
	void SplitPolyIntoTwo(std::vector<uint32_t> &polys, uint32_t startVId, Eigen::Vector2d &interPos, uint32_t anV0Id, uint32_t anV1Id, std::vector<Eigen::Vector2d> &outPoly0, std::vector<Eigen::Vector2d> &outPoly1);

	double ComputePolyDiagLength(std::vector<Eigen::Vector2d> &poly);

	void ComputePolyEdgesLength(std::vector<Eigen::Vector2d> &poly, std::vector<double> &edgesLength);
	void ComputeShortestAndLongestEdge(std::vector<Eigen::Vector2d> &poly, double &minL, double &maxL);
	void ComputeShortestAndLongestEdge(std::vector<Eigen::Vector2d> &poly, std::vector<bool> &polyAddNew, double &minL, double &maxL);

	bool ComputeVerticalIntersection(Eigen::Vector2d &v0, Eigen::Vector2d &v10, Eigen::Vector2d &v11, double &alpha, Eigen::Vector2d &interPos);
	bool JudgeIfVerticalIntersection(std::vector<uint32_t> &poly, uint32_t startVId, Eigen::Vector2d &interPos, uint32_t interV0Id, uint32_t interV1Id);

	void FindListInsertPos(std::list<uint32_t> &pvs, uint32_t v0, uint32_t v1, std::list<uint32_t>::iterator &it, int &diss);

	int FindLargerOldVertexId(int currId, std::vector<bool> &polyAddNew);
	int FindSmallerOldVertexId(int currId, std::vector<bool> &polyAddNew);

	int FindRealLargerOldVertexId(int currId, std::vector<bool> &polyAddNew);
	int FindRealSmallerOldVertexId(int currId, std::vector<bool> &polyAddNew);

	void OutputInnerEdgesInfo(std::string &fileName);
	void OutputInnerEdgesInfo(std::string &fileName, std::vector<uint32_t> &polyToTri, std::vector<Eigen::Vector2d> &sourcePoly);

	bool JudgeVBound(uint32_t id)
	{
		std::vector<uint32_t> &ves = Vs_[id].neighbor_pes;
		for (int i = 0; i < ves.size(); ++i)
		{
			if (Es_[ves[i]].isBoundary)
				return true;
		}
		return false;
	}

	void GetSubAddNew(int currId, int outV2Id, int outV3Id, std::vector<bool> sourcePolyAddNew, std::vector<bool> &poly0AddNew, std::vector<bool> &poly1AddNew)
	{
		poly0AddNew.clear();
		poly1AddNew.clear();
		poly0AddNew.reserve(sourcePolyAddNew.size());
		poly1AddNew.reserve(sourcePolyAddNew.size());
		int i = currId;
		while (i != outV2Id&&i!=outV3Id)
		{
			if (i == currId)
				poly0AddNew.emplace_back(false);
			else
				poly0AddNew.emplace_back(sourcePolyAddNew[i]);
			i = (i + 1) % sourcePolyAddNew.size();
		}
		poly0AddNew.emplace_back(sourcePolyAddNew[i]);
		poly0AddNew.emplace_back(false);

		poly1AddNew.emplace_back(false);
		if (i == outV2Id)
			i = outV3Id;
		else
			i = outV2Id;

		while (i != currId)
		{
			poly1AddNew.emplace_back(sourcePolyAddNew[i]);
			i = (i + 1) % sourcePolyAddNew.size();
		}
		poly1AddNew.emplace_back(false);
	}

	void GetSubAddNew(int currId, int targetId, std::vector<bool> sourcePolyAddNew, std::vector<bool> &poly0AddNew, std::vector<bool> &poly1AddNew);

	void BuildAABBTree(std::vector<Eigen::Vector2d> &inPoly);

	void ComputeTargetSubQuadsNum(std::vector<std::vector<Eigen::Vector2d>> &polys, uint32_t totalQuadsNum, std::vector<int> &subNum);

	void OutputPolysInfo(std::string &fileName);	//.pti poly topology info

	void OutputSplitEdgesInfo(std::string &fileName);
public:
	template<typename T>
	void SetVecToList(std::vector<T> &vec, std::list<T> &lis)
	{
		lis.clear();
		for (int i = 0; i < vec.size(); ++i)
		{
			lis.emplace_back(vec[i]);
		}
	}

	template<typename T>
	void SetListToVec(std::vector<T> &vec, std::list<T> &lis)
	{
		vec.clear();
		vec.reserve(lis.size());
		for (auto it = lis.begin(); it != lis.end(); ++it)
		{
			vec.emplace_back(*it);
		}
	}

public:
	std::vector<PolyVertex> Vs_;
	std::vector<PolyEdge> Es_;
	std::vector<PolyFace> Fs_;

	ANNkd_tree *fieldKdtree_ = NULL;
	std::vector<Eigen::Vector3d> *fields_ = NULL;

	//std::vector<std::vector<uint32_t> > polysVec_;
	std::vector<uint32_t> currFaceId_;
	double sourcePolyArea_ = 0;
	double diagLength_ = 0;

	std::vector<uint32_t> initialPoly_;
	std::vector<uint32_t> poly_;
	uint32_t outV2Id_; 
	uint32_t outV3Id_; 
	Eigen::Vector2d interPos_;
	std::vector<Eigen::Vector2d> outPoly0_;
	std::vector<Eigen::Vector2d> outPoly1_;

	SegmentTree* aabbTree_ = NULL;
	std::vector<Segment> segList_;

	std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d>> splitEdges_;
};

