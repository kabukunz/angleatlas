#pragma once
#include "QuadMesh.h"
#include "BaseComplex.h"
#include "CGAL_AABBTree.h"
#include "MeshDefinition.h"
#include <tuple>
#include <ANN/ANN.h>
#define CORNER_BC 1

namespace DataOperation
{
	const double PI = 3.141592653;
	const double FEATURE_THRESHOLD = 140.0*PI / 180.0;
	typedef struct MeshFeatures
	{
		std::vector<Eigen::Vector2d> bPoints;	//�߽�������
		std::vector<uint32_t> sourcePIndices;	//Ӧ����bPoints��Ӧ��ԭʼ����id
		
		//corners and curves with the respective corner. the curve id is corresponded with the blew.
		std::vector<uint32_t> corners;
		std::vector<std::tuple<int, int>> cornerCurves;

		//curves �����Ķ��� index, �м���corner���м���curve������index��bPoints��index����ԭ����ġ���Ӧ�ð����������˵�
		std::vector<std::vector<uint32_t>> curveVs;
		std::vector<double> curveLength;

		bool haveFeature;

		std::vector<SegmentTree*> aabbTrees;
		std::vector<std::vector<Segment>> segList;

		ANNkd_tree *kdTree = NULL;	//corner kdTree.
		//std::
	}MeshFeatures;

	typedef struct FeatureConstraints
	{
		std::vector<uint32_t> V_ids;	//���²���ô�ԣ������ã�����Ҳûʲô�ã�
		std::vector<int> V_types;	//-2 means this vertex is not a boundary vertex, -1 means this vertex is a corner vertex, >=0 means the curve index the vertex on.

		//corners
		std::vector<uint32_t> ids_C;
		std::vector<Eigen::Vector2d> C;

		//line
		std::vector<uint32_t> ids_L;	
		std::vector<uint32_t> on_which_L;	//��ǰ�ĵ�������feature line�ϣ� ��curveVs��Ӧ��
		
		std::vector<Eigen::Vector2d> axa_L;
		std::vector<Eigen::Vector2d> origin_L;

		
	}FeatureConstraints;

	typedef struct CollapseInfo
	{
		std::vector<std::vector<uint32_t>> vsGroup;
		std::vector<uint32_t> targetVs;
		std::vector<Eigen::Vector2d> targetCoords;
		
		std::vector<uint32_t> hsToBeKilled;

		std::vector<uint32_t> vsNotMove;
		std::vector<uint32_t> hsSeveralring;

		//uint32_t keepedQVId0 = (uint32_t)-1, keepedQVId1 = (uint32_t)-1;
	}CollapseInfo;

	typedef struct OptimizeInfo
	{
		std::vector<uint32_t> globalToLocalVVec;
		std::vector<uint32_t> localToGlobalVVec;

		std::vector<Eigen::Vector2d> localVPos;
		std::vector<Eigen::Vector3i> localTriIds;
		std::vector<std::vector<Eigen::Vector2d>> localStanTriPos;
		std::vector<std::vector<uint32_t>> localVsGroups;
		std::vector<Eigen::Vector2d> localVsCoords;
		std::vector<bool> vsGroupHasBoundary;
		FeatureConstraints localFc;

		std::vector<std::tuple<double, uint32_t, uint32_t, uint32_t>> idsLDis;	//localFc�е�idsL�����ڵıߵ��׸�����ľ��룬Ҫ����

		std::vector<uint32_t> localVsNotMove;
		std::vector<Eigen::Vector2d> localVsNotMovePos;

	}OptimizeInfo;

	enum ElementType
	{
		SHEET,
		CHORD
	};

	typedef std::tuple<double, ElementType, uint32_t> RankingTuple3;

	class Pipeline
	{
	public:
		Pipeline(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<Eigen::Vector3d> &fields);
		Pipeline(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<Eigen::Vector3d> &fields, bool temp);
		~Pipeline();

		void Initialize();
		void DoPipeline();

		bool TopologyCheck(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc);
		bool TopologyCheck(BaseDataStructure::QuadMesh &qm);
		bool TopologyCheckWithoutValence(BaseDataStructure::QuadMesh &qm);
		bool TopologyCheck(BaseDataStructure::BaseComplex &bc);

		void ClearMeshFeatures(MeshFeatures &mff);
		void ClearFeatureConstraints(FeatureConstraints &fcc);
		void InitializeMeshFeatures();
		void FindSomeMeshFeatures(MeshFeatures &mf, BaseDataStructure::QuadVertex &startQv);
		void InitializeFeatureConstraint();

		void ExtractSheetsAndChords(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords);
		void RankingSheetsAndChords(std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, std::vector<RankingTuple3> &candidates);
		void RankingSheetsAndChordsByFields(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, std::vector<RankingTuple3> &candidates);
		void SubdivideSheets(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, FeatureConstraints &fcc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<RankingTuple3> &candidates);
		void SubdivideSheetsFixV2(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, FeatureConstraints &fcc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<RankingTuple3> &candidates);
		bool SubdivideSheet(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, FeatureConstraints &fcc, BaseDataStructure::Sheet &sheet, std::vector<RankingTuple3> &candidates);
		bool SubdivideSheet(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, FeatureConstraints &fcc, uint32_t startE, uint32_t growNum);
		bool SubdivideSheetFixV2(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, FeatureConstraints &fcc, uint32_t startE, uint32_t growNum);
		void SubdivideSheetByEdgeLength(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, FeatureConstraints &fcc, double aveLength);	//Only handle the open sheet at current time.
		void FindSubdivideSheet(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &startEVec, std::vector<uint32_t> &, double aveLength);
		void FindLongestOpenSheet(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &startEVec, std::vector<uint32_t> &fNumVec);
		void FindNearestOpenSheet(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &startEVec, std::vector<uint32_t> &fNumVec);
		//���������ĺ���ֱ�Ϊ����ʼ�ģ��߽磩��Id�����Ĵ�����һ��sheet�������Ŀ������ʼʱ������µ���ԭʼ�����е�id�� �µĵ�λ�ã��µĵ�id�� �����������id�� һ�ߵĵ�id����һ�ߵĵ�id, ԭʼ����
		//ע�⣬�����������������ģ����ı��εĶ����йأ�������㻻��
		void GrowSplitSheet(uint32_t startE, uint32_t growNum, uint32_t &startVId, std::vector<Eigen::Vector2d> &newVPoss, std::vector<uint32_t> &newVIds, std::vector<uint32_t> &fIds, std::vector<uint32_t> &leftVIds, std::vector<uint32_t> &rightVIds, BaseDataStructure::QuadMesh &qm);
		void FindSourceCorner(Eigen::Vector2d &pos, uint32_t &outCorner);

		void BuildSheet(BaseDataStructure::QuadMesh& qm, BaseDataStructure::BaseComplex &bc, uint32_t startBEId, BaseDataStructure::Sheet &sheet, std::vector<bool> &besFlag);
		void BuildChord(BaseDataStructure::QuadMesh& qm, BaseDataStructure::BaseComplex &bc, BaseDataStructure::Chord &chord);

		bool RemoveSheetsAndChords(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, std::vector<RankingTuple3> &candidates);
		bool RemoveSheetsAndChordsFixV2(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, std::vector<RankingTuple3> &candidates);
		bool FilterTopologyInfo(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, RankingTuple3 &candidate);
		double ComputeQuadFieldWeight(BaseDataStructure::QuadMesh &qm, uint32_t quadId);
		double ComputeEdgeFieldAngle(const Eigen::Vector3d &field, const Eigen::Vector2d &ev0, const Eigen::Vector2d &ev1);
		template <typename T>
		bool JudgeIfRepeatEle(const std::vector<T> &sourceVec)
		{
			std::vector<T> copyVec = sourceVec;
			std::sort(copyVec.begin(), copyVec.end());
			if (std::unique(copyVec.begin(), copyVec.end()) != copyVec.end())
				return true;
			else
				return false;
		}

		template <typename T>
		void SafeDeletePtr(T *ptr)
		{
			if (ptr != NULL)
			{
				delete ptr;
				ptr = NULL;
			}
		}

		void FindSheetToBeCollapsedPair(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, BaseDataStructure::Sheet &sheet);
		void FindChordToBeCollapsedPair(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, BaseDataStructure::Chord &chord);

		void FindFaceOppsiteEdges(BaseDataStructure::BaseComplex &bc, uint32_t bfId, uint32_t *pair0, uint32_t *pair1);
		void FindRegularEdgePoints(BaseDataStructure::QuadMesh &qm, uint32_t preE, uint32_t preV, uint32_t vNum, std::vector<uint32_t> &finalVVec);	//finalVVecһ��ʼ�ǿյģ�vNum��ʾ�ڰ�preVҲpush_back������£�finalVVec���Ҫ�����������㡣

		bool JudgeCollapseCriterion(BaseDataStructure::QuadMesh &qm, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, RankingTuple3 &rt3, FeatureConstraints &fcc);	//�����target_vs��target_coords�㶨��
		void BuildCollapseInfo(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, RankingTuple3 &rt3, CollapseInfo *ci);	//�����collapse_info�㶨�ˡ�

		void GrowRegion(BaseDataStructure::QuadMesh &qm, CollapseInfo *ci);
		bool BuildNewTopologyMesh(BaseDataStructure::QuadMesh &inputQm, BaseDataStructure::QuadMesh *outputQm, BaseDataStructure::BaseComplex *outputBc, CollapseInfo *ci, std::vector<uint32_t> &vNewToOld, std::vector<uint32_t> &vOldToNew, FeatureConstraints &fcc);		//�����µ�mesh�����ж��Ƿ����á�
		bool DeleteUselessElements(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &vNewToNewNew, std::vector<uint32_t> &vNewNewToNew);
		bool DeleteUselessElements(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &vNewToNewNew, std::vector<uint32_t> &vNewNewToNew, std::vector<uint32_t> &fNewToNewNew, std::vector<uint32_t> &fNewNewToNew);
		void UpdateOldNewAfterDelete(std::vector<uint32_t> &vOldToNew, std::vector<uint32_t> &vNewToOld, std::vector<uint32_t> &vNewToNewNew, std::vector<uint32_t> &vNewNewToNew);

		void ComputeInitialGenus(BaseDataStructure::QuadMesh &qm)
		{
			uint32_t vNum = qm.Vs_.size();
			uint32_t eNum = qm.Es_.size();
			uint32_t fNum = qm.Fs_.size();

			initialGenus_ = vNum + fNum - eNum;
		}

		void BuildLocalOptimizeInfo(BaseDataStructure::QuadMesh &qm, CollapseInfo *ci, OptimizeInfo *oi);
		void BuildSubdivideOptimizeInfoV2(BaseDataStructure::QuadMesh &qm, FeatureConstraints &fc, const std::vector<uint32_t> &oldFs, const std::vector<uint32_t> &newFs, OptimizeInfo *oi, std::vector<uint32_t> &locaLToGlobalIdsL);
		void BuildSubdivideOptimizeInfo(BaseDataStructure::QuadMesh &qm, FeatureConstraints &fc, const std::vector<uint32_t> &oldFs, const std::vector<uint32_t> &newFs, OptimizeInfo *oi, std::vector<uint32_t> &locaLToGlobalIdsL);
		void ComputeStandardTriangles(BaseDataStructure::QuadMesh &qm, CollapseInfo *ci, OptimizeInfo *oi);
		void ComputeStandardTriangles(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &fs, OptimizeInfo *oi);
		void ComputeStandardTriangles(BaseDataStructure::QuadMesh &qm, uint32_t faceId, double targetArea, std::vector<std::vector<Eigen::Vector2d>> &stanTris);
		double ComputeFaceArea(BaseDataStructure::QuadMesh &qm, uint32_t faceId);
		double ComputeFaceAreasSum(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &faceIds);
		bool DoCollapse(BaseDataStructure::QuadMesh &sourceQm, BaseDataStructure::QuadMesh *newQm, OptimizeInfo *oi);

		void ProjectToSouceLine(OptimizeInfo *oi, MeshFeatures &mf);
		void UpdateFeatures(OptimizeInfo *oldOi, FeatureConstraints &oldFc, FeatureConstraints &newFc, std::vector<uint32_t> &vNewToOld, std::vector<uint32_t> &vNewNewToNew, std::vector<uint32_t> &vOldToNew, std::vector<uint32_t> &vNewToNewNew);
		void UpdateFeatures(OptimizeInfo *newOi, FeatureConstraints &oldFc, FeatureConstraints &newFc);
		void UpdateFeaturesFromOi(FeatureConstraints &fcc, OptimizeInfo *oi);
		void UpdateFeaturesAfterDelete(FeatureConstraints &fc, std::vector<uint32_t> &vNewToNewNew, std::vector<uint32_t> &vNewNewToNew);
		void BuildSegmentAABBTree(MeshFeatures &mf);
		void TestBuildAABBTrees(Point targetP);
		void BuildNewMeshOptimizeInfo(OptimizeInfo *oldOi, OptimizeInfo *newOi, BaseDataStructure::QuadMesh *newQm);
		void ComputeVNewToOld(std::vector<uint32_t> &vOldToNew, std::vector<uint32_t> &vNewToOld, FeatureConstraints &fcc);
		void SortEdgeFromLongToShort(std::vector<uint32_t> &es, BaseDataStructure::QuadMesh &qm);

		void FinalOptimize(BaseDataStructure::QuadMesh &qm);
		void FinalOptimizeWithCornerChange(BaseDataStructure::QuadMesh &qm);
		void FinalOptimizeFixV2(BaseDataStructure::QuadMesh &qm);
		void FinalOptimizeInner(BaseDataStructure::QuadMesh &qm);
		void FinalOptimizeWithCornerChangeFixV2(BaseDataStructure::QuadMesh &qm);
		void BuildOptimizeInfoHasBoundary(BaseDataStructure::QuadMesh &qm, OptimizeInfo *oi);
		void ComputeIdsLDis(OptimizeInfo *oi);
		void AdjustBoundaryPos(BaseDataStructure::QuadMesh &newQm, OptimizeInfo *newOi);
		void AdjustBoundaryPosInOptimization(BaseDataStructure::QuadMesh &newQm, OptimizeInfo *newOi);

		static void RefineQuadToTri(BaseDataStructure::QuadMesh &qm, MyMesh &triMesh);
		void BuildKDTree(std::vector<OpenMesh::Vec3d> &triMids);

		//���Ըı�ǵ�
		uint32_t FindCurveEdgeNum(BaseDataStructure::QuadMesh &qm, FeatureConstraints &fc, uint32_t startV);
		uint32_t FindCurveEdgeNum(BaseDataStructure::QuadMesh &newQm, OptimizeInfo *newOi, uint32_t startV);
		void ChangeTheCorners(BaseDataStructure::QuadMesh &qm, BaseDataStructure::QuadMesh &newQm, std::vector<uint32_t> &vOldToNew, OptimizeInfo *oi, CollapseInfo *ci, FeatureConstraints &fc);
		void ChangeTheCorners(BaseDataStructure::QuadMesh &newQm, OptimizeInfo *newOi);
		void ChangeTheCornersInTwoCurves(BaseDataStructure::QuadMesh &newQm, OptimizeInfo *newOi);
		void SetQuadMesh(BaseDataStructure::QuadMesh *input, BaseDataStructure::QuadMesh *output);
		void SetOI(OptimizeInfo *input, OptimizeInfo *output);

		//�����Ǹ�ֻ����collapse֮����Ż��иı��˽ǵ㣬�������ϣ��������collapse��ʱ����Կ�ǵ���С�isCollapseSheet����˼�塣
		bool CollapseCorners(BaseDataStructure::QuadMesh &qm, CollapseInfo *ci, OptimizeInfo *oi, FeatureConstraints &currFc, std::vector<uint32_t> &vOldToNew, std::vector<uint32_t> &vNewToOld, bool isCollapseSheet);
		bool CollapseCornersVersion2(BaseDataStructure::QuadMesh &qm, CollapseInfo *ci, OptimizeInfo *oi, FeatureConstraints &currFc, std::vector<uint32_t> &vOldToNew, std::vector<uint32_t> &vNewToOld, bool isCollapseSheet);
		
		//void FindOrder
		void DeleteValence2Vertices(BaseDataStructure::QuadMesh &qm, FeatureConstraints &fcc, std::vector<uint32_t> &outFs);

		//Ϊ�˸ɵ��߽��Ϊ2���õ㣬һϵ�еĺ�������Ҫ��д
		bool DoCollapseVersion2(BaseDataStructure::QuadMesh &sourceQm, BaseDataStructure::QuadMesh *newQm, OptimizeInfo *oi);
		void FindFaceSomeRings(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &inFaces, std::vector<uint32_t> &outFaces, int rings);
		void BuildNewOptimizeInfo(BaseDataStructure::QuadMesh &qm, FeatureConstraints &fcc, std::vector<uint32_t> &faces, OptimizeInfo *newOi);

		void FindMovedFaces(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &faces, FeatureConstraints &fcc);
		void ComputeV2BoundaryPos(BaseDataStructure::QuadMesh &qm, uint32_t badV, double threshold);
		void OutputFcInfo(FeatureConstraints &fcc);

	public:
		BaseDataStructure::QuadMesh backupQm_;
		Eigen::Vector2d minP_, maxP_;
		BaseDataStructure::QuadMesh &qm_;
		BaseDataStructure::BaseComplex &bc_;

		MeshFeatures mf;
		FeatureConstraints fc;		//��ǰ�����feature constraints

		std::vector<BaseDataStructure::Sheet> sheets_;
		std::vector<BaseDataStructure::Chord> chords_;
		std::vector<RankingTuple3> candidates_;

		BaseDataStructure::QuadMesh *newQm_ = NULL;
		BaseDataStructure::BaseComplex *newBc_ = NULL;
		//std::vector<uint32_t> sourceToNewVVec_;
		//std::vector<uint32_t> newToSourceVVec_;
		CollapseInfo *ci_ = NULL;
		OptimizeInfo *oi_ = NULL;

		std::vector<uint32_t> vNewToOld_, vOldToNew_;	//��mesh����mesh֮���ӳ��
		std::vector<uint32_t> newQmOptimizedFs_;	//�µ��������Ҫ���Ż��Ĳ��֡�
		Eigen::Vector2d globalVec2d_;
		Eigen::Vector3d vec0_, vec1_, normalVec_;

		int initialGenus_ = 0;
		int initialQuadNum_ = 0;
		double subdivideThre_ = 0.7;
		double stopSubThre_ = 1.0;
		double initialArea_ = 0;

		double hausdorffRadio, hausdorffThre = 0.01;
		double averageEdgeLength_ = 0;

		int lastCandidatePos = 0;
		bool valenceReturn = false;

		bool noInverse = false;	//���ֱ��collapse���ᷴת�Ļ�����ô�Ͳ���Ҫ��һ�����Ż�
		int subdivideRings_ = 2;
		int minusValue = 0;

		SegmentTree::Point_and_primitive_id pAndP; 

public:
		static MyMesh triMesh_;
		ANNkd_tree *kdTree_ = NULL;
		std::vector<Eigen::Vector3d> &fields_;

		bool collapseCorner_ = false;
		bool isCollapseSheet_ = false;

		std::vector<uint32_t> corners_;
		bool isCornerBCV = false;

		double radio_;
		Eigen::Vector2d minPP_;
	};
}

