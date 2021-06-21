#include "Pipeline.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <queue>
#include <numeric>
#include "MetroDis.h"
#include "ExtraFieldEnergyAlone.h"

#include "QuadMeshIO.h"
#include "AESolver.h"
#include "AESolverSquare.h"
#include "AssistFunc.h"
#define USE_AVERAGE_AREA
#define USE_SQUARE_AREA
#define CHANGE_QUAD_RINGS 1
#define USE_VALENCE_WEIGHT 0
#define USE_NO_INVERSE 0
#define USE_FIELD 0
#define USE_3To0 0
#define SPLIT_LONG_EDGES 0
#define CHANGE_THE_CORNERS 1
#define COLLAPSE_CORNER 0
#define OPTIMIZE_FIELDS 1
//#define OUTPUT_MID_MESH
//#define USE_SQUARE_SOLVER
#define OUTPUT_DEBUG_INFO 0
#define NEW_SUBDIVIDE 0
#define TWO_LOOPS 1
using namespace BaseDataStructure;

const double MIN_ALLOWED_JACOBI = 0.1;
int QUAD_RINGS = 4;
const int AES_ITER = 6;
const double VALENCE_WEIGHT = 1000;
double CORNER_ERROR_THRE = 1E-4;
const double WEITH_WIDGET = 0.01;
namespace DataOperation
{
	MyMesh Pipeline::triMesh_;
	Pipeline::Pipeline(QuadMesh &qm, BaseComplex &bc, std::vector<Eigen::Vector3d> &fields) : qm_(qm), bc_(bc), backupQm_(qm), initialQuadNum_(qm.Fs_.size()), fields_(fields)
	{
		backupQm_.ComputeQuadMeshBBox(minP_, maxP_);

		std::vector<uint32_t> faceIds(qm.Fs_.size());
		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			faceIds[i] = i;
		}
		initialArea_ = ComputeFaceAreasSum(qm, faceIds);

		qm.ComputeMeshNormal();
		averageEdgeLength_ = qm.ComputeAverageEdgeLength();
		std::cout << "AverageEdgeLength: " << averageEdgeLength_ << std::endl;
		if (CORNER_ERROR_THRE > averageEdgeLength_ / 10.0)
			CORNER_ERROR_THRE = averageEdgeLength_ / 10.0;
	}

	Pipeline::Pipeline(QuadMesh &qm, BaseComplex &bc, std::vector<Eigen::Vector3d> &fields, bool temp) : qm_(qm), bc_(bc), backupQm_(qm), initialQuadNum_(qm.Fs_.size()), fields_(fields)
	{
	}


	Pipeline::~Pipeline()
	{
		SafeDeletePtr(newQm_);
		SafeDeletePtr(newBc_);
		SafeDeletePtr(ci_);
		SafeDeletePtr(kdTree_);
		SafeDeletePtr(mf.kdTree);
		mf.kdTree = NULL;
		kdTree_ = NULL;
		newQm_ = NULL;
		newBc_ = NULL;
		ci_ = NULL;

		for (int i = 0; i < mf.aabbTrees.size(); ++i)
		{
			if (mf.aabbTrees[i] != NULL)
			{
				delete mf.aabbTrees[i];
				mf.aabbTrees[i] = NULL;
			}
		}
	}

	void Pipeline::Initialize()
	{
		if (!TopologyCheck(qm_, bc_))
		{
			std::cout << "Didn't pass topology check! " << std::endl;
			return;
		}

		InitializeMeshFeatures();
		BuildSegmentAABBTree(mf);
		InitializeFeatureConstraint();

		std::vector<uint32_t> outFs;
		DeleteValence2Vertices(qm_, fc, outFs);

		bc_.ExtractBaseComplex(&qm_, FEATURE_THRESHOLD);
		if (bc_.isError_)
			exit(-23452);

		ExtractSheetsAndChords(qm_, bc_, sheets_, chords_);
#if !USE_FIELD
		RankingSheetsAndChords(sheets_, chords_, candidates_);
#else
		RankingSheetsAndChordsByFields(qm_, bc_, sheets_, chords_, candidates_);
#endif
	}

	void Pipeline::InitializeMeshFeatures()
	{
#pragma region Old_Version
		//ClearMeshFeatures(mf);

		////Build boundary Points
		//QuadVertex startQv;
		//for (QuadVertex &qv : qm_.Vs_)
		//{
		//	if (qv.boundary)
		//	{
		//		startQv = qv;
		//		break;
		//	}
		//}
		//QuadEdge startQe;
		//for (uint32_t qeIndex : startQv.neighbor_es)
		//{
		//	if (qm_.Es_[qeIndex].boundary)
		//	{
		//		startQe = qm_.Es_[qeIndex];
		//		break;
		//	}
		//}
		//mf.bPoints.emplace_back(qm_.V_[startQv.id]);
		//mf.sourcePIndices.emplace_back(startQv.id);

		//uint32_t preV, preE = startQe.id, nextE;
		//if (startQe.vs[0] == startQv.id)
		//{
		//	preV = startQe.vs[1];
		//	mf.bPoints.emplace_back(qm_.V_[startQe.vs[1]]);
		//	mf.sourcePIndices.emplace_back(startQe.vs[1]);
		//}
		//else
		//{
		//	preV = startQe.vs[0];
		//	mf.bPoints.emplace_back(qm_.V_[startQe.vs[0]]);
		//	mf.sourcePIndices.emplace_back(startQe.vs[0]);
		//}


		//std::function<bool(uint32_t, uint32_t, uint32_t&)> findNextPoint = [&](uint32_t preV, uint32_t preE, uint32_t &nextE)->bool
		//{
		//	for (uint32_t vEs : qm_.Vs_[preV].neighbor_es)
		//	{
		//		if (qm_.Es_[vEs].boundary && vEs != preE)
		//		{
		//			nextE = vEs;
		//			return true;
		//		}
		//	}
		//	return false;
		//};

		//while (findNextPoint(preV, preE, nextE))
		//{
		//	uint32_t nextV = 0;
		//	if (qm_.Es_[nextE].vs[0] == preV)
		//		nextV = qm_.Es_[nextE].vs[1];
		//	else
		//		nextV = qm_.Es_[nextE].vs[0];

		//	if (nextV == startQv.id)
		//		break;
		//	else
		//	{
		//		mf.bPoints.emplace_back(qm_.V_[nextV]);
		//		mf.sourcePIndices.emplace_back(nextV);
		//		preV = nextV;
		//		preE = nextE;
		//	}
		//}

		////Find feature points;
		//for (int i = 0; i < mf.bPoints.size(); ++i)
		//{
		//	Eigen::Vector2d vec1 = mf.bPoints[(i - 1 + mf.bPoints.size()) % mf.bPoints.size()] - mf.bPoints[i];
		//	Eigen::Vector2d vec2 = mf.bPoints[(i + 1) % mf.bPoints.size()] - mf.bPoints[i];
		//	double cosValue = vec1.normalized().dot(vec2.normalized());
		//	if (cosValue > std::cos(FEATURE_THRESHOLD))
		//		mf.corners.emplace_back(i);
		//}
		//if (mf.corners.empty())
		//	mf.corners.emplace_back(0);
		//	/*mf.haveFeature = false;
		//else
		//	mf.haveFeature = true;*/

		//for (int i = 0; i < mf.corners.size(); ++i)
		//{
		//	mf.cornerCurves.emplace_back(std::tuple<uint32_t, uint32_t>((i - 1 + mf.corners.size()) % mf.corners.size(), i));
		//}
		//mf.curveVs.resize(mf.corners.size());
		//if (mf.corners.size() == 1)
		//{
		//	mf.curveVs[0].reserve(mf.bPoints.size());
		//	for (int i = 0; i < mf.bPoints.size(); ++i)
		//	{
		//		mf.curveVs[0].emplace_back((mf.corners[0] + i) % mf.bPoints.size());
		//	}
		//}
		//else if (mf.corners.size() > 1)
		//{
		//	for (int i = 0; i < mf.corners.size() - 1; ++i)
		//	{
		//		uint32_t currCor = mf.corners[i];
		//		uint32_t nextCor = mf.corners[i + 1];
		//		uint32_t currP = currCor;
		//		while (currP <= nextCor)
		//		{
		//			mf.curveVs[i].emplace_back(currP);
		//			++currP;
		//		}
		//	}

		//	uint32_t temp0 = mf.corners[mf.corners.size() - 1];
		//	uint32_t temp1 = mf.corners[0] + mf.bPoints.size();
		//	for (int i = temp0; i <= temp1; ++i)
		//	{
		//		mf.curveVs[mf.corners.size() - 1].emplace_back(i%mf.bPoints.size());
		//	}
		//}
#pragma endregion
		std::vector<uint32_t> boundaryV;
		boundaryV.reserve(qm_.Vs_.size());
		for (int i = 0; i < qm_.Vs_.size(); ++i)
		{
			if (qm_.Vs_[i].boundary)
				boundaryV.emplace_back(i);
		}
		std::vector<bool> vFlag(qm_.Vs_.size(), false);
		std::vector<MeshFeatures> mfVec;

		QuadVertex startQv;
		while (true)
		{
			bool hasNotBuildV = false;
			for (int i = 0; i < boundaryV.size(); ++i)
			{
				if (vFlag[boundaryV[i]] == false)
				{
					hasNotBuildV = true;
					startQv = qm_.Vs_[boundaryV[i]];
					break;
				}
			}
			if (hasNotBuildV)
			{
				mfVec.emplace_back(MeshFeatures());
				MeshFeatures &currMf = mfVec[mfVec.size() - 1];
				FindSomeMeshFeatures(currMf, startQv);	//存在某些边界curve上没有corner，但是有curve。
				for (int i = 0; i < currMf.sourcePIndices.size(); ++i)
				{
					vFlag[currMf.sourcePIndices[i]] = true;
				}
			}
			else
				break;
		}

		uint32_t currVCount = 0, currCurveCount = 0;
		for (int i = 0; i < mfVec.size(); ++i)
		{
			mf.haveFeature = true;
			mf.bPoints.insert(mf.bPoints.end(), mfVec[i].bPoints.begin(), mfVec[i].bPoints.end());
			mf.sourcePIndices.insert(mf.sourcePIndices.end(), mfVec[i].sourcePIndices.begin(), mfVec[i].sourcePIndices.end());

			MeshFeatures &currMf = mfVec[i];
			for (int j = 0; j < currMf.corners.size(); ++j)
			{
				mf.corners.emplace_back(currMf.corners[j] + currVCount);
			}
			for (int j = 0; j < currMf.cornerCurves.size(); ++j)
			{
				int f = std::get<0>(currMf.cornerCurves[j]);
				int s = std::get<1>(currMf.cornerCurves[j]);
				mf.cornerCurves.emplace_back(std::tuple<int, int>(f + currCurveCount, s + currCurveCount));
			}
			for (int j = 0; j < currMf.curveVs.size(); ++j)
			{
				std::vector<uint32_t> tempVec;
				for (int k = 0; k < currMf.curveVs[j].size(); ++k)
				{
					tempVec.emplace_back(currMf.curveVs[j][k] + currVCount);
				}
				mf.curveVs.emplace_back(tempVec);
			}
			currVCount += currMf.bPoints.size();
			currCurveCount += currMf.curveVs.size();
		}
			
		//Build curveLength
		for (int i = 0; i < mf.curveVs.size(); ++i)
		{
			double sumLength = 0;
			std::vector<uint32_t> &currCurve = mf.curveVs[i];
			for (int j = 0; j < currCurve.size() - 1; ++j)
			{
				sumLength += (mf.bPoints[currCurve[j + 1]] - mf.bPoints[currCurve[j]]).norm();
			}
			mf.curveLength.emplace_back(sumLength);
		}

		//Build corner kdTree
		std::vector<Eigen::Vector3d> cornerPos;
		cornerPos.reserve(mf.corners.size());
		for (int i = 0; i < mf.corners.size(); ++i)
		{
			cornerPos.emplace_back(Eigen::Vector3d(mf.bPoints[mf.corners[i]][0], mf.bPoints[mf.corners[i]][1], 0));
		}
		if (mf.kdTree != NULL)
		{
			delete mf.kdTree;
			mf.kdTree = NULL;
		}

		ANNpointArray dataPts = annAllocPts(cornerPos.size(), 3);
		for (int i = 0; i < cornerPos.size(); ++i)
		{
			dataPts[i][0] = cornerPos[i][0];
			dataPts[i][1] = cornerPos[i][1];
			dataPts[i][2] = cornerPos[i][2];
		}
		mf.kdTree = new ANNkd_tree(dataPts, cornerPos.size(), 3);
	}

	void Pipeline::FindSourceCorner(Eigen::Vector2d &pos, uint32_t &outCorner)
	{
		ANNpoint ap = annAllocPt(3);
		ap[0] = pos[0]; ap[1] = pos[1]; ap[2] = 0;
		ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
		mf.kdTree->annkSearch(ap, 1, nnIdx, dists);
		outCorner = nnIdx[0];
		
		annDeallocPt(ap);
		delete[] nnIdx; delete[] dists;
	}

	void Pipeline::FindSomeMeshFeatures(MeshFeatures &mf, QuadVertex &startQv)
	{
		ClearMeshFeatures(mf);

		//Build boundary Points
		/*QuadVertex startQv;
		for (QuadVertex &qv : qm_.Vs_)
		{
			if (qv.boundary)
			{
				startQv = qv;
				break;
			}
		}*/
		QuadEdge startQe;
		for (uint32_t qeIndex : startQv.neighbor_es)
		{
			if (qm_.Es_[qeIndex].boundary)
			{
				startQe = qm_.Es_[qeIndex];
				break;
			}
		}
		mf.bPoints.emplace_back(qm_.V_[startQv.id]);
		mf.sourcePIndices.emplace_back(startQv.id);

		uint32_t preV, preE = startQe.id, nextE;
		if (startQe.vs[0] == startQv.id)
		{
			preV = startQe.vs[1];
			mf.bPoints.emplace_back(qm_.V_[startQe.vs[1]]);
			mf.sourcePIndices.emplace_back(startQe.vs[1]);
		}
		else
		{
			preV = startQe.vs[0];
			mf.bPoints.emplace_back(qm_.V_[startQe.vs[0]]);
			mf.sourcePIndices.emplace_back(startQe.vs[0]);
		}


		std::function<bool(uint32_t, uint32_t, uint32_t&)> findNextPoint = [&](uint32_t preV, uint32_t preE, uint32_t &nextE)->bool
		{
			for (uint32_t vEs : qm_.Vs_[preV].neighbor_es)
			{
				if (qm_.Es_[vEs].boundary && vEs != preE)
				{
					nextE = vEs;
					return true;
				}
			}
			return false;
		};

		while (findNextPoint(preV, preE, nextE))
		{
			uint32_t nextV = 0;
			if (qm_.Es_[nextE].vs[0] == preV)
				nextV = qm_.Es_[nextE].vs[1];
			else
				nextV = qm_.Es_[nextE].vs[0];

			if (nextV == startQv.id)
				break;
			else
			{
				mf.bPoints.emplace_back(qm_.V_[nextV]);
				mf.sourcePIndices.emplace_back(nextV);
				preV = nextV;
				preE = nextE;
			}
		}

		//Find feature points;
		for (int i = 0; i < mf.bPoints.size(); ++i)
		{
			Eigen::Vector2d vec1 = mf.bPoints[(i - 1 + mf.bPoints.size()) % mf.bPoints.size()] - mf.bPoints[i];
			Eigen::Vector2d vec2 = mf.bPoints[(i + 1) % mf.bPoints.size()] - mf.bPoints[i];
			double cosValue = vec1.normalized().dot(vec2.normalized());
			if (cosValue > std::cos(FEATURE_THRESHOLD))
				mf.corners.emplace_back(i);
		}
		if (mf.corners.empty())
		{
			mf.corners.emplace_back(0);
			mf.haveFeature = false;
		}
		else
			mf.haveFeature = true;

		for (int i = 0; i < mf.corners.size(); ++i)
		{
			mf.cornerCurves.emplace_back(std::tuple<uint32_t, uint32_t>((i - 1 + mf.corners.size()) % mf.corners.size(), i));
		}
		mf.curveVs.resize(mf.corners.size());
		if (mf.corners.size() == 1)
		{
			mf.curveVs[0].reserve(mf.bPoints.size());
			for (int i = 0; i < mf.bPoints.size(); ++i)
			{
				mf.curveVs[0].emplace_back((mf.corners[0] + i) % mf.bPoints.size());
			}
		}
		else if (mf.corners.size() > 1)
		{
			for (int i = 0; i < mf.corners.size() - 1; ++i)
			{
				uint32_t currCor = mf.corners[i];
				uint32_t nextCor = mf.corners[i + 1];
				uint32_t currP = currCor;
				while (currP <= nextCor)
				{
					mf.curveVs[i].emplace_back(currP);
					++currP;
				}
			}

			uint32_t temp0 = mf.corners[mf.corners.size() - 1];
			uint32_t temp1 = mf.corners[0] + mf.bPoints.size();
			for (int i = temp0; i <= temp1; ++i)
			{
				mf.curveVs[mf.corners.size() - 1].emplace_back(i%mf.bPoints.size());
			}
		}

		if (!mf.haveFeature)
		{
			mf.corners.clear();
			mf.cornerCurves.clear();
		}
	}

	void Pipeline::ClearMeshFeatures(MeshFeatures &mff)
	{
		/*std::vector<Eigen::Vector2d>().swap(mff.bPoints);

		std::vector<uint32_t>().swap(mff.corners);
		std::vector<std::tuple<int, int>>().swap(mff.cornerCurves);

		std::vector<std::vector<uint32_t>>().swap(mff.curveVs);*/
		mff.bPoints.clear();
		mff.corners.clear();
		mff.cornerCurves.clear();
		mff.curveVs.clear();

		bool haveFeature = true;
	}

	bool Pipeline::TopologyCheck(QuadMesh &qm)
	{
		uint32_t vNum = qm.Vs_.size();
		uint32_t eNum = qm.Es_.size();
		uint32_t fNum = qm.Fs_.size();

		if (vNum + fNum - eNum != initialGenus_)
		{
			return false;
		}

		//Base structure check.
		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fes = qm.Fs_[i].es;
			std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
			if (fes.size() != 4 || fvs.size() != 4)
				return false;
		}
		for (int i = 0; i < qm.Es_.size(); ++i)
		{
			std::vector<uint32_t> &efs = qm.Es_[i].neighbor_fs;
			std::vector<uint32_t> &evs = qm.Es_[i].vs;
			if (efs.size() == 0 || evs.size() != 2)
				return false;
		}
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			std::vector<uint32_t> &vfs = qm.Vs_[i].neighbor_fs;
			std::vector<uint32_t> &ves = qm.Vs_[i].neighbor_es;
			if (vfs.empty() || ves.empty())
				return false;
		}

		//Quad mesh 2-manifold check.
		for (int i = 0; i < qm.Es_.size(); ++i)
		{
			uint32_t eFsNum = qm.Es_[i].neighbor_fs.size();
			if (eFsNum != 1 && eFsNum != 2)
				return false;
		}
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			uint32_t vEsNum = qm.Vs_[i].neighbor_es.size();
			uint32_t vFsNum = qm.Vs_[i].neighbor_fs.size();
			if ((qm.Vs_[i].boundary && vEsNum - vFsNum != 1) || (!qm.Vs_[i].boundary && vEsNum - vFsNum != 0))
				return false;
		}

		//Base complex 2-manifold check.
		/*for (int i = 0; i < bc.Be_.size(); ++i)
		{
			uint32_t beFsNum = bc.Be_[i].b_neighbor_fs.size();
			if (beFsNum != 1 && beFsNum != 2)
				return false;
		}
		for (int i = 0; i < bc.Bv_.size(); ++i)
		{
			uint32_t bvEsNum = bc.Bv_[i].b_neighbor_es.size();
			uint32_t bvFsNum = bc.Bv_[i].b_neighbor_fs.size();
			if ((bc.Bv_[i].boundary&&bvEsNum - bvFsNum != 1) || (!bc.Bv_[i].boundary&&bvEsNum - bvFsNum != 0))
				return false;
		}
*/
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if ((!qm.Vs_[i].boundary && qm.Vs_[i].neighbor_es.size() <= 2)
				|| qm.Vs_[i].neighbor_es.size() <= 1)
				return false;
		}
		return true;
	}

	bool Pipeline::TopologyCheck(BaseComplex &bc)
	{
		uint32_t bvNum = bc.Bv_.size();
		uint32_t beNum = bc.Be_.size();
		uint32_t bfNum = bc.Bf_.size();
		if (bvNum + bfNum - beNum != initialGenus_ )
		{
			return false;
		}

		for (int i = 0; i < bc.Bf_.size(); ++i)
		{
			if (bc.Bf_[i].b_vs.size() != 4)
				return false;
		}

		for (int i = 0; i < bc.Be_.size(); ++i)
		{
			std::vector<uint32_t> &bevs = bc.Be_[i].b_vs;
			if (bevs.size() != 2 || bevs[0] == bevs[1])
				return false;
		}

		return true;
	}

	bool Pipeline::TopologyCheck(QuadMesh &qm, BaseComplex &bc)
	{
		uint32_t vNum = qm.Vs_.size();
		uint32_t eNum = qm.Es_.size();
		uint32_t fNum = qm.Fs_.size();

		uint32_t bvNum = bc.Bv_.size();
		uint32_t beNum = bc.Be_.size();
		uint32_t bfNum = bc.Bf_.size();
		if (bvNum + bfNum - beNum != initialGenus_ || vNum + fNum - eNum != initialGenus_)
		{
			return false;
		}
		
		//Base structure check.
		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fes = qm.Fs_[i].es;
			std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
			if (fes.size() != 4 || fvs.size() != 4)
				return false;
		}
		for (int i = 0; i < qm.Es_.size(); ++i)
		{
			std::vector<uint32_t> &efs = qm.Es_[i].neighbor_fs;
			std::vector<uint32_t> &evs = qm.Es_[i].vs;
			if (efs.size() == 0 || evs.size() != 2)
				return false;
		}
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			std::vector<uint32_t> &vfs = qm.Vs_[i].neighbor_fs;
			std::vector<uint32_t> &ves = qm.Vs_[i].neighbor_es;
			if (vfs.empty() || ves.empty())
				return false;
		}

		//Quad mesh 2-manifold check.
		for (int i = 0; i < qm.Es_.size(); ++i)
		{
			uint32_t eFsNum = qm.Es_[i].neighbor_fs.size();
			if (eFsNum != 1 && eFsNum != 2)
				return false;
		}
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			uint32_t vEsNum = qm.Vs_[i].neighbor_es.size();
			uint32_t vFsNum = qm.Vs_[i].neighbor_fs.size();
			if ( (qm.Vs_[i].boundary && vEsNum - vFsNum != 1) || (!qm.Vs_[i].boundary && vEsNum - vFsNum != 0) )
				return false;
		}

		//Base complex 2-manifold check.
		/*for (int i = 0; i < bc.Be_.size(); ++i)
		{
			uint32_t beFsNum = bc.Be_[i].b_neighbor_fs.size();
			if (beFsNum != 1 && beFsNum != 2)
				return false;
		}
		for (int i = 0; i < bc.Bv_.size(); ++i)
		{
			uint32_t bvEsNum = bc.Bv_[i].b_neighbor_es.size();
			uint32_t bvFsNum = bc.Bv_[i].b_neighbor_fs.size();
			if ((bc.Bv_[i].boundary&&bvEsNum - bvFsNum != 1) || (!bc.Bv_[i].boundary&&bvEsNum - bvFsNum != 0))
				return false;
		}
*/
		for (int i = 0; i < bc.Bf_.size(); ++i)
		{
			if (bc.Bf_[i].b_vs.size() != 4)
				return false;
		}
		for (int i = 0; i < bc.Be_.size(); ++i)
		{
			std::vector<uint32_t> &bevs = bc.Be_[i].b_vs;
			if (bevs.size() != 2 || bevs[0] == bevs[1])
				return false;
		}
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if ((!qm.Vs_[i].boundary && qm.Vs_[i].neighbor_es.size() <= 2)
				|| qm.Vs_[i].neighbor_es.size()<=1)
				return false;
		}
		return true;
	}

	void Pipeline::InitializeFeatureConstraint()
	{
		ClearFeatureConstraints(fc);

		fc.V_types.resize(qm_.Vs_.size());
		std::fill(fc.V_types.begin(), fc.V_types.end(), -2);
		for (int i = 0; i < mf.corners.size(); ++i)
		{
			fc.V_types[mf.sourcePIndices[mf.corners[i]]] = -1;
			fc.ids_C.emplace_back(mf.sourcePIndices[mf.corners[i]]);
			fc.C.emplace_back(qm_.V_[mf.sourcePIndices[mf.corners[i]]]);
		}
		for (int i = 0; i < mf.curveVs.size(); ++i)
		{
			for (int j = 1; j < mf.curveVs[i].size()-1; ++j)
			{
				fc.V_types[mf.sourcePIndices[mf.curveVs[i][j]]] = i;
			}
		}

		/*fc.V_ids.resize(qm_.Vs_.size());
		std::fill(fc.V_ids.begin(), fc.V_ids.end(), (uint32_t)-1);
		for (int i = 0; i < mf.sourcePIndices.size(); ++i)
		{
			fc.V_ids[mf.sourcePIndices[i]] = i;
		}*/

		for (int i = 0; i < mf.curveVs.size(); ++i)
		{
			for (int j = 1; j < mf.curveVs[i].size() - 1; ++j)
			{
				Eigen::Vector2d vec1 = mf.bPoints[mf.curveVs[i][j]] - mf.bPoints[mf.curveVs[i][j - 1]];
				Eigen::Vector2d vec2 = mf.bPoints[mf.curveVs[i][j + 1]] - mf.bPoints[mf.curveVs[i][j]];
				fc.axa_L.emplace_back((vec1.normalized() + vec2.normalized()).normalized());
				fc.ids_L.emplace_back(mf.sourcePIndices[mf.curveVs[i][j]]);
				fc.on_which_L.emplace_back(i);
				fc.origin_L.emplace_back(qm_.V_[mf.sourcePIndices[mf.curveVs[i][j]]]);
			}
		}

	}

	void Pipeline::ClearFeatureConstraints(FeatureConstraints &fcc)
	{
		//std::vector<uint32_t>().swap( fcc.V_ids );
		//std::vector<int>().swap( fcc.V_types );

		////corners
		//std::vector<uint32_t>().swap( fcc.ids_C );
		//std::vector<Eigen::Vector2d>().swap( fcc.C );

		////line
		//std::vector<uint32_t>().swap( fcc.ids_L );
		//std::vector<Eigen::Vector2d>().swap( fcc.axa_L);
		//std::vector<Eigen::Vector2d>().swap( fcc.origin_L);
		fcc.V_ids.clear();
		fcc.V_types.clear();
		fcc.ids_C.clear();
		fcc.C.clear();
		fcc.ids_L.clear();
		fcc.axa_L.clear();
		fcc.on_which_L.clear();
		fcc.origin_L.clear();
	}

	void Pipeline::ExtractSheetsAndChords(QuadMesh &qm, BaseComplex &bc, std::vector<Sheet> &sheets, std::vector<Chord> &chords)
	{
		//std::vector<Sheet>().swap(sheets);
		//std::vector<Chord>().swap(chords);
		sheets.clear();
		chords.clear();

		std::vector<bool> besFlag(bc.Be_.size(), false);

		for (int i = 0; i < bc.Be_.size(); ++i)
		{
			if (besFlag[i])
				continue;
			besFlag[i] = true;

			Sheet currSheet;
			currSheet.id = sheets.size();
			currSheet.st = SheetType::OPEN;
			currSheet.haveAddValenceWidget = false;
			BuildSheet(qm, bc, i, currSheet, besFlag);
			sheets.emplace_back(currSheet);
		}

		std::vector<uint32_t> sameVec;
		for (int i = 0; i < bc.Bf_.size(); ++i)
		{
			std::vector<uint32_t> &v0BEs = bc.Bv_[bc.Bf_[i].b_vs[0]].b_neighbor_es;
			std::sort(v0BEs.begin(), v0BEs.end());
			std::vector<uint32_t> couple0(2,0), couple1;
			couple0[0] = bc.Bf_[i].b_vs[0];
			for (int j = 0; j < 4; ++j)
			{
				std::vector<uint32_t> &v1BEs = bc.Bv_[bc.Bf_[i].b_vs[j]].b_neighbor_es;
				std::sort(v1BEs.begin(), v1BEs.end());
				sameVec.clear();
				std::set_intersection(v0BEs.begin(), v0BEs.end(), v1BEs.begin(), v1BEs.end(), std::inserter(sameVec, sameVec.begin()));
				if (sameVec.empty())
				{
					couple0[1] = bc.Bf_[i].b_vs[j];
					break;
				}
			}
			std::sort(bc.Bf_[i].b_vs.begin(), bc.Bf_[i].b_vs.end());
			if (couple0[0] > couple0[1])
				std::swap(couple0[0], couple0[1]);
			std::set_difference(bc.Bf_[i].b_vs.begin(), bc.Bf_[i].b_vs.end(), couple0.begin(), couple0.end(), std::inserter(couple1, couple1.begin()));

			Chord currChord, currChord1;
			currChord.haveAddValenceWidget = false;
			currChord.fid = i;
			currChord.id = chords.size();
			currChord.keepedV0 = couple0[0]; currChord.keepedV1 = couple0[1];
			currChord.removedV0 = couple1[0]; currChord.removedV1 = couple1[1];
			currChord.keepedQV0 = bc.Bv_[couple0[0]].qId;
			currChord.keepedQV1 = bc.Bv_[couple0[1]].qId;
			BuildChord(qm, bc, currChord);
			chords.emplace_back(currChord);

			currChord1.haveAddValenceWidget = false;
			currChord1.fid = i;
			currChord1.id = chords.size();
			currChord1.keepedV0 = couple1[0]; currChord1.keepedV1 = couple1[1];
			currChord1.removedV0 = couple0[0]; currChord1.removedV1 = couple0[1];
			currChord.keepedQV0 = bc.Bv_[couple1[0]].qId;
			currChord.keepedQV1 = bc.Bv_[couple1[1]].qId;
			BuildChord(qm, bc, currChord1);
			chords.emplace_back(currChord1);
		}
	}

	void Pipeline::BuildSheet(QuadMesh& qm, BaseComplex &bc, uint32_t startBEId, BaseDataStructure::Sheet &sheet, std::vector<bool> &besFlag)
	{

		std::vector<bool> midEVFlag(bc.Bv_.size(), false);
		std::vector<bool> bfFlag(bc.Bf_.size(), false);
		std::vector<bool> beFlag(bc.Be_.size(), false);	//beFlag是判断当前边是否加入了sheet, 而besFlag是判断当前边是否是mid_es。

		//Find preF and preE.
		uint32_t startF = bc.Be_[startBEId].b_neighbor_fs[0];

		uint32_t preF = startF, preE = (uint32_t)-1;
		for (uint32_t nextE : bc.Bf_[startF].b_es)
		{
			uint32_t startV0 = bc.Be_[startBEId].b_vs[0], startV1 = bc.Be_[startBEId].b_vs[1];
			uint32_t nextV0 = bc.Be_[nextE].b_vs[0], nextV1 = bc.Be_[nextE].b_vs[1];
			if (startV0 != nextV0 && startV0 != nextV1 && startV1 != nextV0 && startV1 != nextV1)
			{
				preE = nextE;
				break;
			}
		}
		if (preE == (uint32_t)-1)
			exit(-4);

		//emplace back the initial v e f ....and flags.
		if (!midEVFlag[bc.Be_[startBEId].b_vs[0]])
		{
			sheet.b_vs.emplace_back(bc.Be_[startBEId].b_vs[0]);
			midEVFlag[bc.Be_[startBEId].b_vs[0]] = true;
		}
		if (!midEVFlag[bc.Be_[startBEId].b_vs[1]])
		{
			sheet.b_vs.emplace_back(bc.Be_[startBEId].b_vs[1]);
			midEVFlag[bc.Be_[startBEId].b_vs[1]] = true;
		}
		if (!midEVFlag[bc.Be_[preE].b_vs[0]])
		{
			sheet.b_vs.emplace_back(bc.Be_[preE].b_vs[0]);
			midEVFlag[bc.Be_[preE].b_vs[0]] = true;
		}
		if (!midEVFlag[bc.Be_[preE].b_vs[1]])
		{
			sheet.b_vs.emplace_back(bc.Be_[preE].b_vs[1]);
			midEVFlag[bc.Be_[preE].b_vs[1]] = true;
		}
		for (int i = 0; i < 4; ++i)
		{
			if (!beFlag[bc_.Bf_[startF].b_es[i]])
			{
				sheet.b_es.emplace_back(bc_.Bf_[startF].b_es[i]);
				beFlag[bc_.Bf_[startF].b_es[i]] = true;
			}
		}
		sheet.b_fs.emplace_back(startF);
		bfFlag[startF] = true;
		sheet.b_middle_es.emplace_back(startBEId);
		sheet.b_middle_es.emplace_back(preE);
		besFlag[startBEId] = true;
		//besFlag[preE] = true;
		if (bc.Be_[startBEId].boundary)
			sheet.b_middle_boundary_es.emplace_back(startBEId);
		if (bc.Be_[preE].boundary)
			sheet.b_middle_boundary_es.emplace_back(preE);

		//find left_vs and right_vs
		uint32_t leftV0 = bc.Be_[startBEId].b_vs[0];
		uint32_t rightV0 = bc.Be_[startBEId].b_vs[1];
		uint32_t newLV, newRV;
		sheet.left_vs.emplace_back(leftV0);
		sheet.right_vs.emplace_back(rightV0);

		std::function<void(uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t&, uint32_t&)> findNewLvRv = [&](uint32_t bf, uint32_t bePre, uint32_t beNext, uint32_t lastLv, uint32_t lastRv, uint32_t& newLv, uint32_t& newRv)
		{
			for (int i = 0; i < 4; ++i)
			{
				uint32_t currE = bc.Bf_[bf].b_es[i];
				if (currE != bePre && currE != beNext)
				{
					if (bc.Be_[currE].b_vs[0] == lastLv)
					{
						newLv = bc.Be_[currE].b_vs[1];
						continue;
					}
					else if (bc.Be_[currE].b_vs[1] == lastLv)
					{
						newLv = bc.Be_[currE].b_vs[0];
						continue;
					}
					else if (bc.Be_[currE].b_vs[0] == lastRv)
					{
						newRv = bc.Be_[currE].b_vs[1];
						continue;
					}
					else if (bc.Be_[currE].b_vs[1] == lastRv)
					{
						newRv = bc.Be_[currE].b_vs[0];
						continue;
					}
				}
			}
		};
		findNewLvRv(startF, startBEId, preE, leftV0, rightV0, newLV, newRV);
		sheet.left_vs.emplace_back(newLV); sheet.right_vs.emplace_back(newRV);

		//start mid es
		std::function<bool(uint32_t, uint32_t, uint32_t&, uint32_t&)> findNextFE = [&](uint32_t preF, uint32_t preE, uint32_t &nextF, uint32_t &nextE)->bool
		{
			if (bc.Be_[preE].boundary || besFlag[preE])
				return false;

			if (bc.Be_[preE].b_neighbor_fs[0] == preF)
				nextF = bc.Be_[preE].b_neighbor_fs[1];
			else
				nextF = bc.Be_[preE].b_neighbor_fs[0];

			nextE = (uint32_t)-1;
			for (uint32_t tempE : bc.Bf_[nextF].b_es)
			{
				uint32_t startV0 = bc.Be_[preE].b_vs[0], startV1 = bc.Be_[preE].b_vs[1];
				uint32_t nextV0 = bc.Be_[tempE].b_vs[0], nextV1 = bc.Be_[tempE].b_vs[1];
				if (startV0 != nextV0 && startV0 != nextV1 && startV1 != nextV0 && startV1 != nextV1)
				{
					nextE = tempE;
					break;
				}
			}

			if (nextE == (uint32_t)-1)
				exit(-16);
			return true;
		};

		uint32_t nextF, nextE;
		while (findNextFE(preF, preE, nextF, nextE))
		{
			besFlag[preE] = true;
			for (int i = 0; i < 4; ++i)
			{
				uint32_t currE = bc.Bf_[nextF].b_es[i];
				if (currE != preE && currE != nextE && !beFlag[currE])
				{
					sheet.b_es.emplace_back(currE);
					beFlag[currE] = true;
				}
			}
			if (!bfFlag[nextF])
			{
				sheet.b_fs.emplace_back(nextF);
				bfFlag[nextF] = true;
			}

			if (!besFlag[nextE])
			{
				sheet.b_middle_es.emplace_back(nextE);
				uint32_t nextEV0 = bc.Be_[nextE].b_vs[0], nextEV1 = bc.Be_[nextE].b_vs[1];
				if (midEVFlag[nextEV0] || midEVFlag[nextEV1])
					sheet.st = SheetType::INTERSECTION;

				if (!midEVFlag[nextEV0])
				{
					sheet.b_vs.emplace_back(nextEV0);
					midEVFlag[nextEV0] = true;
				}
				if (!midEVFlag[nextEV1])
				{
					sheet.b_vs.emplace_back(nextEV1);
					midEVFlag[nextEV1] = true;
				}

				if (!beFlag[nextE])
				{
					sheet.b_es.emplace_back(nextE);
					beFlag[nextE] = true;
				}

				if (bc.Be_[nextE].boundary)
					sheet.b_middle_boundary_es.emplace_back(nextE);

				findNewLvRv(nextF, preE, nextE, sheet.left_vs[sheet.left_vs.size() - 1], sheet.right_vs[sheet.right_vs.size() - 1], newLV, newRV);
				sheet.left_vs.emplace_back(newLV); sheet.right_vs.emplace_back(newRV);
			}
			else
			{
				sheet.st = SheetType::CLOSED;
			}

			preF = nextF; preE = nextE;
		}
		besFlag[preE] = true;
		if (sheet.b_middle_boundary_es.empty())
			sheet.st = SheetType::CLOSED;

		//reverse some elements.
		std::reverse(sheet.b_fs.begin(), sheet.b_fs.end());
		std::reverse(sheet.b_middle_es.begin(), sheet.b_middle_es.end());
		std::reverse(sheet.left_vs.begin(), sheet.left_vs.end());
		std::reverse(sheet.right_vs.begin(), sheet.right_vs.end());

		

		if (sheet.st != SheetType::CLOSED)
		{
			//如果开始的startBEId不是边界边且没有闭环，那么就往另一个方向找
			if (!bc.Be_[startBEId].boundary)
			{
				preF = bc.Be_[startBEId].b_neighbor_fs[1];
				preE = (uint32_t)-1;
				for (int i = 0; i < 4; ++i)
				{
					uint32_t currBE = bc.Bf_[preF].b_es[i];
					uint32_t currBV0 = bc.Be_[currBE].b_vs[0];
					uint32_t currBV1 = bc.Be_[currBE].b_vs[1];
					uint32_t preBV0 = bc.Be_[startBEId].b_vs[0];
					uint32_t preBV1 = bc.Be_[startBEId].b_vs[1];
					if (currBV0 != preBV0 && currBV0 != preBV1 && currBV1 != preBV0 && currBV1 != preBV1)
					{
						preE = currBE;
					}
				}
				if (preE == (uint32_t)-1)
					exit(-40);

				if (!midEVFlag[bc.Be_[preE].b_vs[0]])
				{
					sheet.b_vs.emplace_back(bc.Be_[preE].b_vs[0]);
					midEVFlag[bc.Be_[preE].b_vs[0]] = true;
				}
				if (!midEVFlag[bc.Be_[preE].b_vs[1]])
				{
					sheet.b_vs.emplace_back(bc.Be_[preE].b_vs[1]);
					midEVFlag[bc.Be_[preE].b_vs[1]] = true;
				}
				for (int i = 0; i < 4; ++i)
				{
					if (bc_.Bf_[preF].b_es[i] == startBEId)
						continue;
					if (!beFlag[bc_.Bf_[preF].b_es[i]])
					{
						sheet.b_es.emplace_back(bc_.Bf_[preF].b_es[i]);
						beFlag[bc_.Bf_[preF].b_es[i]] = true;
					}
				}
				sheet.b_fs.emplace_back(preF);
				bfFlag[preF] = true;
				sheet.b_middle_es.emplace_back(preE);
				//besFlag[preE] = true;
				if (bc.Be_[preE].boundary)
					sheet.b_middle_boundary_es.emplace_back(preE);

				uint32_t leftV0 = bc.Be_[startBEId].b_vs[0];
				uint32_t rightV0 = bc.Be_[startBEId].b_vs[1];
				uint32_t newLV, newRV;
				findNewLvRv(preF, startBEId, preE, leftV0, rightV0, newLV, newRV);
				sheet.left_vs.emplace_back(newLV); sheet.right_vs.emplace_back(newRV);

				while (findNextFE(preF, preE, nextF, nextE))
				{
					besFlag[preE] = true;
					for (int i = 0; i < 4; ++i)
					{
						uint32_t currE = bc.Bf_[nextF].b_es[i];
						if (currE != preE && currE != nextE && !beFlag[currE])
						{
							sheet.b_es.emplace_back(currE);
							beFlag[currE] = true;
						}
					}
					if (!bfFlag[nextF])
					{
						sheet.b_fs.emplace_back(nextF);
						bfFlag[nextF] = true;
					}

					if (!besFlag[nextE])
					{
						sheet.b_middle_es.emplace_back(nextE);
						uint32_t nextEV0 = bc.Be_[nextE].b_vs[0], nextEV1 = bc.Be_[nextE].b_vs[1];
						if (midEVFlag[nextEV0] || midEVFlag[nextEV1])
							sheet.st = SheetType::INTERSECTION;

						if (!midEVFlag[nextEV0])
						{
							sheet.b_vs.emplace_back(nextEV0);
							midEVFlag[nextEV0] = true;
						}
						if (!midEVFlag[nextEV1])
						{
							sheet.b_vs.emplace_back(nextEV1);
							midEVFlag[nextEV1] = true;
						}

						if (!beFlag[nextE])
						{
							sheet.b_es.emplace_back(nextE);
							beFlag[nextE] = true;
						}

						if (bc.Be_[nextE].boundary)
							sheet.b_middle_boundary_es.emplace_back(nextE);

						findNewLvRv(nextF, preE, nextE, sheet.left_vs[sheet.left_vs.size() - 1], sheet.right_vs[sheet.right_vs.size() - 1], newLV, newRV);
						sheet.left_vs.emplace_back(newLV); sheet.right_vs.emplace_back(newRV);
					}
					else
					{
						sheet.st = SheetType::CLOSED;
					}

					preF = nextF; preE = nextE;
				}
				besFlag[preE] = true;
			}
		}

		//Build sheet.qf_delete
		//std::vector<uint32_t>().swap(sheet.qf_delete);
		sheet.qf_delete.clear();
		for (int i = 0; i < sheet.b_fs.size(); ++i)
		{
			sheet.qf_delete.insert(sheet.qf_delete.end(), bc.Bf_[sheet.b_fs[i]].fs_net.begin(), bc.Bf_[sheet.b_fs[i]].fs_net.end());
		}

		//judge if something repeat.
		std::sort(sheet.b_fs.begin(), sheet.b_fs.end());
		sheet.b_fs.erase(std::unique(sheet.b_fs.begin(), sheet.b_fs.end()), sheet.b_fs.end());
		std::sort(sheet.qf_delete.begin(), sheet.qf_delete.end());
		sheet.qf_delete.erase(std::unique(sheet.qf_delete.begin(), sheet.qf_delete.end()), sheet.qf_delete.end());
		if (JudgeIfRepeatEle(sheet.b_es) /*|| JudgeIfRepeatEle(sheet.b_fs)*/ || JudgeIfRepeatEle(sheet.b_middle_es) || JudgeIfRepeatEle(sheet.b_vs) || JudgeIfRepeatEle(sheet.b_middle_boundary_es) /*|| JudgeIfRepeatEle(sheet.qf_delete)*/)
		{
			QuadMeshIO qmi;
			qmi.WriteQuadMesh(&qm, "C:\\Users\\ChiZhang\\Desktop\\repeatMesh.obj");
			ofstream ofs("C:\\Users\\ChiZhang\\Desktop\\repeatMesh.txt");
			ofs << "b_es: " << JudgeIfRepeatEle(sheet.b_es) << std::endl;
			ofs << "b_fs: " << JudgeIfRepeatEle(sheet.b_fs) << std::endl;
			ofs << "b_middle_es: " << JudgeIfRepeatEle(sheet.b_middle_es) << std::endl;
			ofs << "b_vs: " << JudgeIfRepeatEle(sheet.b_vs) << std::endl;
			ofs << "b_middle_boundary_es: " << JudgeIfRepeatEle(sheet.b_middle_boundary_es) << std::endl;
			ofs << "qf_delete: " << JudgeIfRepeatEle(sheet.qf_delete) << std::endl;
			ofs << std::endl;
			for (int i = 0; i < sheet.qf_delete.size(); ++i)
			{
				ofs << sheet.qf_delete[i] << std::endl;
			}
			exit(-400);
		}

		std::function<double(uint32_t)> computeBELength = [&](uint32_t bme)->double
		{
			std::vector<uint32_t> &eslinks = bc.Be_[bme].es_link;
			double sumLength = 0;
			for (int i = 0; i < eslinks.size(); ++i)
			{
				uint32_t v0 = qm.Es_[eslinks[i]].vs[0], v1 = qm.Es_[eslinks[i]].vs[1];
				sumLength += (qm.V_[v0] - qm.V_[v1]).norm();
			}
			return sumLength;
		};

		//build sheet widget
		sheet.weight = 0;
		for (int i = 0; i < sheet.b_middle_es.size(); ++i)
		{
			sheet.weight += computeBELength(sheet.b_middle_es[i]);
		}
		sheet.weight /= sheet.b_middle_es.size();
	}

	void Pipeline::BuildChord(QuadMesh& qm, BaseComplex &bc, BaseDataStructure::Chord &chord)
	{
		uint32_t qv0 = bc.Bv_[chord.removedV0].qId;
		uint32_t qv1 = bc.Bv_[chord.removedV1].qId;

		FrameFace &currFc = bc.Bf_[chord.fid];

		chord.weight = (qm.V_[qv0] - qm.V_[qv1]).norm();
	}

	void Pipeline::RankingSheetsAndChords(std::vector<Sheet> &sheets, std::vector<Chord> &chords, std::vector<RankingTuple3> &candidates)
	{
		//std::vector<RankingTuple3>().swap(candidates);
		candidates.clear();
		candidates.reserve(sheets.size() + chords.size());

		for (int i = 0; i < sheets.size(); ++i)
		{
			candidates.emplace_back(RankingTuple3(sheets[i].weight, ElementType::SHEET, sheets[i].id));
		}
		for (int i = 0; i < chords.size(); ++i)
		{
			candidates.emplace_back(RankingTuple3(chords[i].weight, ElementType::CHORD, chords[i].id));
		}
		std::sort(candidates.begin(), candidates.end());
	}

	void Pipeline::DoPipeline()
	{
		/*FinalOptimizeInner(qm_);
		QuadMeshIO qmio;
		qmio.WriteQuadMesh(&qm_, "C:\\Users\\ChiZhang\\Desktop\\errorG00.obj");
		exit(0);*/

		/*FinalOptimize(qm_);
		QuadMeshIO qmio;
		double finalAveLength0 = qm_.ComputeAverageEdgeLength();
		SubdivideSheetByEdgeLength(qm_, bc_, fc, finalAveLength0);
		qmio.WriteQuadMesh(&qm_, "C:\\Users\\ChiZhang\\Desktop\\errorG00.obj");
		FinalOptimize(qm_);
		qmio.WriteQuadMesh(&qm_, "C:\\Users\\ChiZhang\\Desktop\\errorG01.obj");
		exit(0);*/

		isCornerBCV = false;
		int sucCount = 0;
		while (true)
		{
			/*SubdivideSheetByEdgeLength(qm_, bc_, fc, 0.3);
			exit(0);*/
			/*if (!qm_.JudgeQuadMeshJacobi() || qm_.minJacobi_<MIN_ALLOWED_JACOBI)
			{
				std::cout << "Inverse Jacobi! " << std::endl;
				return;
			}*/

			if (!RemoveSheetsAndChords(qm_, bc_, sheets_, chords_, candidates_))
				break;

			std::cout << "Success Count: " << ++sucCount << std::endl;
			ExtractSheetsAndChords(qm_, bc_, sheets_, chords_);
#if !USE_FIELD
			RankingSheetsAndChords(sheets_, chords_, candidates_);
#else
			RankingSheetsAndChordsByFields(qm_, bc_, sheets_, chords_, candidates_);
#endif

#if OUTPUT_DEBUG_INFO
			std::cout << "Before SubdivideSheets" << std::endl;
#endif
			SubdivideSheets(qm_, bc_, fc, sheets_, candidates_);

#if OUTPUT_DEBUG_INFO
			std::cout << "After SubdivideSheets" << std::endl;
#endif
			/*if (sucCount == 102)
			{
				QuadMeshIO qmio;
				qmio.WriteQuadMesh(&qm_, "C:\\Users\\ChiZhang\\Desktop\\error7.obj");
				exit(0);
			}*/
		}

#if TWO_LOOPS
		isCornerBCV = true;

#define CHANGE_THE_CORNERS 0
		corners_.clear();
		for (int i = 0; i < fc.ids_C.size(); ++i)
		{
			corners_.emplace_back(fc.ids_C[i]);
		}
		bc_.ExtractBaseComplex(&qm_, corners_);
		ExtractSheetsAndChords(qm_, bc_, sheets_, chords_);
		RankingSheetsAndChords(sheets_, chords_, candidates_);
		sucCount = 0;
		while (true)
		{
			/*SubdivideSheetByEdgeLength(qm_, bc_, fc, 0.3);
			exit(0);*/
			/*if (!qm_.JudgeQuadMeshJacobi() || qm_.minJacobi_<MIN_ALLOWED_JACOBI)
			{
				std::cout << "Inverse Jacobi! " << std::endl;
				return;
			}*/

			if (!RemoveSheetsAndChords(qm_, bc_, sheets_, chords_, candidates_))
				break;

			std::cout << "Success Count: " << ++sucCount << std::endl;
			ExtractSheetsAndChords(qm_, bc_, sheets_, chords_);
#if !USE_FIELD
			RankingSheetsAndChords(sheets_, chords_, candidates_);
#else
			RankingSheetsAndChordsByFields(qm_, bc_, sheets_, chords_, candidates_);
#endif

#if OUTPUT_DEBUG_INFO
			std::cout << "Before SubdivideSheets" << std::endl;
#endif
			SubdivideSheets(qm_, bc_, fc, sheets_, candidates_);

#if OUTPUT_DEBUG_INFO
			std::cout << "After SubdivideSheets" << std::endl;
#endif
			/*if (sucCount == 102)
			{
				QuadMeshIO qmio;
				qmio.WriteQuadMesh(&qm_, "C:\\Users\\ChiZhang\\Desktop\\error7.obj");
				exit(0);
			}*/
		}
#endif

		double finalAveLength = qm_.ComputeAverageEdgeLength();
#if SPLIT_LONG_EDGES
		FinalOptimize(qm_);
		SubdivideSheetByEdgeLength(qm_, bc_, fc, finalAveLength);
#endif
		FinalOptimizeInner(qm_);
		FinalOptimize(qm_);
		FinalOptimizeWithCornerChangeFixV2(qm_);

		qm_.ScaleQuadMeshBack(radio_, minPP_);
	}

	void Pipeline::OutputFcInfo(FeatureConstraints &fcc)
	{
		ofstream ofs("C:\\Users\\ChiZhang\\Desktop\\fccInfo.txt");
		ofs << "mf.corner size: " << mf.corners.size() << std::endl;

		ofs << "fcc.corner size: " << fcc.ids_C.size() << std::endl;
		for (int i = 0; i < fcc.ids_C.size(); ++i)
		{
			ofs << fcc.ids_C[i] << " ";
		}
		ofs << std::endl;

		ofs << "fcc.idsL size: " << fcc.ids_L.size() << std::endl;
		for (int i = 0; i < fcc.ids_L.size(); ++i)
		{
			ofs << fcc.ids_L[i] << " ";
		}
		ofs << std::endl;

		for (int i = 0; i < fcc.ids_C.size(); ++i)
		{
			if (fcc.V_types[fcc.ids_C[i]] != -1)
				ofs << "Error in V_types C! " << std::endl;
		}
		for (int i = 0; i < fcc.ids_L.size(); ++i)
		{
			if (fcc.V_types[fcc.ids_L[i]] != fcc.on_which_L[i])
				ofs << "Error in V_types L! " << std::endl;
		}
	}

	bool Pipeline::RemoveSheetsAndChordsFixV2(QuadMesh &qm, BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, std::vector<RankingTuple3> &candidates)
	{
		std::vector<RankingTuple3> tempCandidates;
		tempCandidates.insert(tempCandidates.end(), candidates.begin() + lastCandidatePos, candidates.end());
		tempCandidates.insert(tempCandidates.end(), candidates.begin(), candidates.begin() + lastCandidatePos);
		tempCandidates.swap(candidates);

#ifdef OUTPUT_MID_MESH
		QuadMeshIO qmi;
		qmi.WriteQuadMesh(&qm, "C:\\Users\\ChiZhang\\Desktop\\startMesh.obj");
		OutputFcInfo(fc);
#endif

#if OUTPUT_DEBUG_INFO
		std::cout << std::endl << "Start Candidates Num: " << lastCandidatePos << std::endl;
#endif

		for (int i = 0; i < candidates.size(); ++i)
		{
#if OUTPUT_DEBUG_INFO
			std::cout << std::get<0>(candidates[i]) << " " << std::get<1>(candidates[i]) << " " << std::get<2>(candidates[i]) << endl;
#endif
			//i = 8;
			//std::cout << std::get<1>(candidates[i]) << " " << std::get<2>(candidates[i]) << std::endl;
			if (std::get<1>(candidates[i]) == ElementType::SHEET)
				isCollapseSheet_ = true;
			else
				isCollapseSheet_ = false;

#if OUTPUT_DEBUG_INFO
			std::cout << "Before FilterTopologyInfo" << std::endl;
#endif
			if (!FilterTopologyInfo(qm, bc, sheets, chords, candidates[i]))
			{
#if USE_VALENCE_WEIGHT
				if (valenceReturn)
				{
					std::sort(candidates.begin(), candidates.end());
					tempCandidates.clear();
					tempCandidates.insert(tempCandidates.end(), candidates.begin() + lastCandidatePos, candidates.end());
					tempCandidates.insert(tempCandidates.end(), candidates.begin(), candidates.begin() + lastCandidatePos);
					tempCandidates.swap(candidates);
					--i;
				}
#endif
				continue;
			}

#if OUTPUT_DEBUG_INFO
			std::cout << "Before BuildLocalOptimizeInfo" << std::endl;
#endif
			SafeDeletePtr(oi_);
			oi_ = NULL;
			oi_ = new OptimizeInfo();
			BuildLocalOptimizeInfo(qm, ci_, oi_);

#if OUTPUT_DEBUG_INFO
			std::cout << "Vs Pair: " << ci_->vsGroup[0][0] << " " << ci_->vsGroup[0][1] << std::endl;
#endif
#if OUTPUT_DEBUG_INFO
			std::cout << "Before Collapse" << std::endl;
#endif

			std::cout << "Curr Ite: " << i << std::endl;
			if (!DoCollapseVersion2(qm, newQm_, oi_))
			{
				newQm_ = NULL;
				continue;
			}

			/*for (int i = 0; i < fc.C.size(); ++i)
			{
				std::cout << "Corner: " << i << " " << fc.C[i][0] << ", " << fc.C[i][1] << std::endl;
			}*/

			if (i >= candidates.size() - lastCandidatePos)
				lastCandidatePos = i - candidates.size() + lastCandidatePos;
			else if (i < candidates.size() - lastCandidatePos)
				lastCandidatePos = i + lastCandidatePos;

			if (lastCandidatePos > sheets.size()*0.3 || lastCandidatePos > 40) lastCandidatePos = 0;

			newQm_ = NULL;
			return true;
		}
		return false;
	}

	bool Pipeline::RemoveSheetsAndChords(QuadMesh &qm, BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, std::vector<RankingTuple3> &candidates)
	{
		std::vector<RankingTuple3> tempCandidates;
		tempCandidates.insert(tempCandidates.end(), candidates.begin() + lastCandidatePos, candidates.end());
		tempCandidates.insert(tempCandidates.end(), candidates.begin(), candidates.begin() + lastCandidatePos);
		tempCandidates.swap(candidates);

#ifdef OUTPUT_MID_MESH
		QuadMeshIO qmi;
		qmi.WriteQuadMesh(&qm, "C:\\Users\\ChiZhang\\Desktop\\startMesh.obj");
		OutputFcInfo(fc);
#endif

#if OUTPUT_DEBUG_INFO
		std::cout << std::endl << "Start Candidates Num: " << lastCandidatePos << std::endl;
#endif

		for (int i = 0; i < candidates.size(); ++i)
		{
#if OUTPUT_DEBUG_INFO
			std::cout << std::get<0>(candidates[i]) << " " << std::get<1>(candidates[i]) << " " << std::get<2>(candidates[i]) << endl;
#endif
			//std::cout << std::get<1>(candidates[i]) << " " << std::get<2>(candidates[i]) << std::endl;
			if (std::get<1>(candidates[i]) == ElementType::SHEET)
				isCollapseSheet_ = true;
			else
				isCollapseSheet_ = false;

#if OUTPUT_DEBUG_INFO
			std::cout << "Before FilterTopologyInfo" << std::endl;
#endif
			if (!FilterTopologyInfo(qm, bc, sheets, chords, candidates[i]))
			{
#if USE_VALENCE_WEIGHT
				if (valenceReturn)
				{
					std::sort(candidates.begin(), candidates.end());
					tempCandidates.clear();
					tempCandidates.insert(tempCandidates.end(), candidates.begin() + lastCandidatePos, candidates.end());
					tempCandidates.insert(tempCandidates.end(), candidates.begin(), candidates.begin() + lastCandidatePos);
					tempCandidates.swap(candidates);
					--i;
				}
#endif
				continue;
			}

#if OUTPUT_DEBUG_INFO
			std::cout << "Before BuildLocalOptimizeInfo" << std::endl;
#endif
			SafeDeletePtr(oi_);
			oi_ = NULL;
			oi_ = new OptimizeInfo();
			BuildLocalOptimizeInfo(qm, ci_, oi_);

#if OUTPUT_DEBUG_INFO
			std::cout << "Vs Pair: " << ci_->vsGroup[0][0] << " " << ci_->vsGroup[0][1] << std::endl;
#endif
#if OUTPUT_DEBUG_INFO
			std::cout << "Before Collapse" << std::endl;
#endif

			std::cout << "Curr Ite: " << i << std::endl;
			if (!DoCollapse/*Version2*/(qm, newQm_, oi_))
			{
				newQm_ = NULL;
				continue;
			}

			/*for (int i = 0; i < fc.C.size(); ++i)
			{
				std::cout << "Corner: " << i << " " << fc.C[i][0] << ", " << fc.C[i][1] << std::endl;
			}*/

			if (i >= candidates.size() - lastCandidatePos)
				lastCandidatePos = i - candidates.size() + lastCandidatePos;
			else if (i < candidates.size() - lastCandidatePos)
				lastCandidatePos = i + lastCandidatePos;

			if (lastCandidatePos >sheets.size()*0.3 || lastCandidatePos > 40) lastCandidatePos = 0;

			newQm_ = NULL;
			return true;
		}
		return false;
	}

	bool Pipeline::FilterTopologyInfo(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, RankingTuple3 &candidate)
	{

		ElementType currET = std::get<1>(candidate);
		uint32_t currId = std::get<2>(candidate);
		/*if (currId == 20 && currET == ElementType::CHORD)
			double ddds = 34;*/
		if (currET == ElementType::SHEET)
		{
			FindSheetToBeCollapsedPair(qm, bc, sheets[currId]);
		}
		else
		{
			FindChordToBeCollapsedPair(qm, bc, chords[currId]);
		}

		if (!JudgeCollapseCriterion(qm, sheets, chords, candidate, fc))
			return false;
		SafeDeletePtr(ci_);
		ci_ = NULL;
		ci_ = new CollapseInfo();
		BuildCollapseInfo(qm, bc, sheets, chords, candidate, ci_);

		/*ofstream ofs("C:\\Users\\ChiZhang\\Desktop\\faces3.txt");
		ofs << sheets[4].qf_delete.size() << std::endl;
		for (int j = 0; j < sheets[4].qf_delete.size(); ++j)
		{
			ofs << sheets[4].qf_delete[j] << std::endl;
		}
		exit(0);*/
		
		SafeDeletePtr(newQm_);
		SafeDeletePtr(newBc_);
		newQm_ = NULL;
		newBc_ = NULL;
		newQm_ = new QuadMesh();
		newBc_ = new BaseComplex();
#if USE_NO_INVERSE
		noInverse = false;
#endif
		if (!BuildNewTopologyMesh(qm, newQm_, newBc_, ci_, vNewToOld_, vOldToNew_, fc))	//这里只检查了拓扑
			return false;

#ifdef OUTPUT_MID_MESH
		QuadMeshIO qmio;
		qmio.WriteQuadMesh(&qm, "C:\\Users\\ChiZhang\\Desktop\\error17.obj");
		qmio.WriteQuadMesh(newQm_, "C:\\Users\\ChiZhang\\Desktop\\error18.obj");
#endif

		return true;
	}

	void Pipeline::FindSheetToBeCollapsedPair(QuadMesh &qm, BaseComplex &bc, Sheet &sheet)
	{
		/*std::vector<std::vector<uint32_t>>().swap(sheet.vs_pairs);
		std::vector<std::vector<uint32_t>>().swap(sheet.vs_group);
		std::vector<std::vector<uint32_t>>().swap(sheet.vs_links);
		std::vector<uint32_t>().swap(sheet.target_vs);
		std::vector<Eigen::Vector2d>().swap(sheet.target_coords);*/
		sheet.vs_pairs.clear();
		sheet.vs_group.clear();
		sheet.vs_links.clear();
		sheet.target_vs.clear();
		sheet.target_coords.clear();

		std::vector<uint32_t> vsCounts(qm.Vs_.size(), 0);
		std::vector<std::vector<uint32_t>> pairVsInfo(qm.Vs_.size());
		std::vector<uint32_t> tempVsPair(2);
		std::vector<bool> bMiddleEsFlag(bc.Be_.size(), false);
		//先把b_middle_es搞定，并且碰到的点次数加一。
		for (int i = 0; i < sheet.b_middle_es.size(); ++i)
		{
			uint32_t be = sheet.b_middle_es[i];
			tempVsPair[0] = bc.Bv_[bc.Be_[be].b_vs[0]].qId;
			tempVsPair[1] = bc.Bv_[bc.Be_[be].b_vs[1]].qId;
			pairVsInfo[tempVsPair[0]].emplace_back(tempVsPair[1]);
			pairVsInfo[tempVsPair[1]].emplace_back(tempVsPair[0]);
			sheet.vs_pairs.emplace_back(tempVsPair);
			sheet.vs_links.emplace_back(bc.Be_[be].vs_link);
			for (int j = 0; j < bc.Be_[be].vs_link.size(); ++j)
			{
				++vsCounts[bc.Be_[be].vs_link[j]];
			}
			bMiddleEsFlag[be] = true;
		}

		//找每个base complex内部的，位于source quad中的点。
		uint32_t vsLinkNum = sheet.vs_links[0].size();
		std::vector<uint32_t> currVsLink(vsLinkNum, 0);
		for (int i = 0; i < sheet.b_fs.size(); ++i)
		{
			std::vector<uint32_t> currSideBE;
			for (int j = 0; j < 4; ++j)
			{
				uint32_t tempBE = bc.Bf_[sheet.b_fs[i]].b_es[j];
				if (!bMiddleEsFlag[tempBE])
				{
					currSideBE.emplace_back(tempBE);
					break;
				}
			}
			if (currSideBE.empty())
			{
				currSideBE.emplace_back(bc.Bf_[sheet.b_fs[i]].b_es[0]);
				std::vector<uint32_t> &evs0 = bc.Be_[bc.Bf_[sheet.b_fs[i]].b_es[0]].b_vs;
				for (int j = 1; j < 4; ++j)
				{
					std::vector<uint32_t> &evs1 = bc.Be_[bc.Bf_[sheet.b_fs[i]].b_es[j]].b_vs;
					if (evs0[0] == evs1[0] || evs0[0] == evs1[1] || evs0[1] == evs1[0] || evs0[1] == evs1[1])
					{
						currSideBE.emplace_back(bc.Bf_[sheet.b_fs[i]].b_es[j]);
						break;
					}
				}
			}

			for (int n = 0; n < currSideBE.size(); ++n)
			{
				std::vector<uint32_t> &sideBeVs = bc.Be_[currSideBE[n]].vs_link;
				for (int j = 1; j < sideBeVs.size() - 1; ++j)
				{
					std::vector<uint32_t> &currVEs = qm.Vs_[sideBeVs[j]].neighbor_es;
					uint32_t rightE = (uint32_t)-1;
					for (int k = 0; k < currVEs.size(); ++k)
					{
						std::vector<uint32_t> &tempFs = qm.Es_[currVEs[k]].neighbor_fs;
						if (tempFs.size() == 2 && qm.Fs_[tempFs[0]].bId == sheet.b_fs[i] && qm.Fs_[tempFs[1]].bId == sheet.b_fs[i])
						{
							rightE = currVEs[k];
							break;
						}
					}
					if (rightE == (uint32_t)-1)
						exit(-6);

					tempVsPair[0] = sideBeVs[j];
					currVsLink[0] = sideBeVs[j];
					++vsCounts[sideBeVs[j]];
					uint32_t preV = sideBeVs[j];
					for (int k = 1; k < vsLinkNum; ++k)
					{
						if (qm.Es_[rightE].vs[0] == preV)
							currVsLink[k] = qm.Es_[rightE].vs[1];
						else
							currVsLink[k] = qm.Es_[rightE].vs[0];
						++vsCounts[currVsLink[k]];

						if (k != vsLinkNum - 1)
						{
							preV = currVsLink[k];
							std::vector<uint32_t> &ves = qm.Vs_[preV].neighbor_es;
							for (int m = 0; m < ves.size(); ++m)
							{
								std::vector<uint32_t> &fs0 = qm.Es_[rightE].neighbor_fs;
								std::vector<uint32_t> &fs1 = qm.Es_[ves[m]].neighbor_fs;
								if (fs0[0] != fs1[0] && fs0[0] != fs1[1] && fs0[1] != fs1[0] && fs0[1] != fs1[1])
								{
									rightE = ves[m];
									break;
								}
							}
						}
					}
					tempVsPair[1] = currVsLink[currVsLink.size() - 1];
					pairVsInfo[tempVsPair[0]].emplace_back(tempVsPair[1]);
					pairVsInfo[tempVsPair[1]].emplace_back(tempVsPair[0]);
					sheet.vs_links.emplace_back(currVsLink);
					sheet.vs_pairs.emplace_back(tempVsPair);
				}
			}
		}

		//在考虑到一些Base Complex Quad 自交的情况下，确定vs_group
		//先确定一个quad四条边都是middle es的情况，且要考虑这样的quad相邻的情况。
		std::vector<uint32_t> vInPairFlag(bc.Bv_.size(), (uint32_t)-1);
		std::vector<std::vector<uint32_t>> tempVsGroups;
		std::vector<uint32_t> onePointFaces;
		std::vector<std::vector<uint32_t>> fsInOntPoint;
		for (int i = 0; i < sheet.b_fs.size(); ++i)
		{
			bool isFAllMiddleEs = true;
			for (int j = 0; j < 4; ++j)
			{
				if (!bMiddleEsFlag[bc.Bf_[sheet.b_fs[i]].b_es[j]])
				{
					isFAllMiddleEs = false;
					break;
				}
			}

			if (isFAllMiddleEs)
			{
				onePointFaces.emplace_back(sheet.b_fs[i]);

				std::vector<uint32_t> &currVsNet = bc.Bf_[sheet.b_fs[i]].vs_net;
				for (int j = 0; j < currVsNet.size(); ++j)
				{
					vsCounts[currVsNet[j]] = 0;
				}
			}
		}

		std::vector<uint32_t> resultVec;
		std::vector<uint32_t> fes0, fes1, fvs0, fvs1;
		std::function<bool(uint32_t, uint32_t)> judgeIfAdj = [&](uint32_t f0, uint32_t f1)->bool
		{
			resultVec.clear();
			//即使不共边，单单只有共点的情况也会出现，这时也是邻接的。
			/*fes0 = bc.Bf_[f0].b_es;
			fes1 = bc.Bf_[f1].b_es;
			std::sort(fes0.begin(), fes0.end());
			std::sort(fes1.begin(), fes1.end());
			std::set_intersection(fes0.begin(), fes0.end(), fes1.begin(), fes1.end(), std::inserter(resultVec, resultVec.begin()));*/
			fvs0 = bc.Bf_[f0].b_vs;
			fvs1 = bc.Bf_[f1].b_vs;
			std::sort(fvs0.begin(), fvs0.end());
			std::sort(fvs1.begin(), fvs1.end());
			std::set_intersection(fvs0.begin(), fvs0.end(), fvs1.begin(), fvs1.end(), std::inserter(resultVec, resultVec.begin()));
			
			if (resultVec.empty())
				return false;
			else
				return true;
		};

		if (!onePointFaces.empty())
		{
			//找相邻“全被干掉四边形”找错了，有些没有找到。 下面vs_group的更新也错了，全被干掉四边形也有可能与其他被干掉点相邻。
			/*std::vector<uint32_t> tempF(1, onePointFaces[0]);
			fsInOntPoint.emplace_back(tempF);

			for (int i = 1; i < onePointFaces.size(); ++i)
			{
				bool isAdj = false;
				for (int j = 0; j < fsInOntPoint.size(); ++j)
				{
					for (int k = 0; k < fsInOntPoint[j].size(); ++k)
					{
						if (judgeIfAdj(onePointFaces[i], fsInOntPoint[j][k]))
						{
							fsInOntPoint[j].emplace_back(onePointFaces[i]);
							isAdj = true;
							break;
						}
					}
					if (isAdj)
						break;
				}
				if (!isAdj)
				{
					tempF[0] = onePointFaces[i];
					fsInOntPoint.emplace_back(tempF);
				}
			}*/

			std::vector<std::vector<uint32_t>> adjFsPairs(onePointFaces.size());
			for (int i = 0; i < onePointFaces.size(); ++i)
			{
				for (int j = i + 1; j < onePointFaces.size(); ++j)
				{
					if (judgeIfAdj(onePointFaces[i], onePointFaces[j]))
					{
						adjFsPairs[i].emplace_back(j);
						adjFsPairs[j].emplace_back(i);
					}
				}
			}
			std::vector<bool> adjFsFlag(onePointFaces.size(), false);
			for (int i = 0; i < adjFsPairs.size(); ++i)
			{
				if (adjFsPairs[i].empty())
				{
					fsInOntPoint.emplace_back(std::vector<uint32_t>{onePointFaces[i]});
					adjFsFlag[i] = true;
				}
				else
				{
					if (adjFsFlag[i] == true)
					{
						continue;
					}
					else
					{
						std::queue<uint32_t> fsPool;
						std::vector<uint32_t> currAdjFs;
						currAdjFs.emplace_back(onePointFaces[i]);
						adjFsFlag[i] = true;
						fsPool.push(i);
						while (!fsPool.empty())
						{
							uint32_t currEle = fsPool.front();
							fsPool.pop();
							
							std::vector<uint32_t> &eleVec = adjFsPairs[currEle];
							for (int j = 0; j < eleVec.size(); ++j)
							{
								if (!adjFsFlag[eleVec[j]])
								{
									adjFsFlag[eleVec[j]] = true;
									currAdjFs.emplace_back(onePointFaces[eleVec[j]]);
									fsPool.push(eleVec[j]);
								}
							}
						}
						fsInOntPoint.emplace_back(currAdjFs);
					}
				}
			}
		}

		for (int i = 0; i < fsInOntPoint.size(); ++i)
		{
			std::vector<uint32_t> tempVec;
			for (int j = 0; j < fsInOntPoint[i].size(); ++j)
			{
				std::vector<uint32_t> &fvs = bc.Bf_[fsInOntPoint[i][j]].b_vs;
				for (int k = 0; k < fvs.size(); ++k)
				{
					tempVec.emplace_back(bc.Bv_[fvs[k]].qId);
				}
			}

			std::vector<uint32_t> tempVecc = tempVec;
			for (int j = 0; j < tempVecc.size(); ++j)
			{
				std::queue<uint32_t> vsPool;
				vsPool.push(tempVecc[j]);
				while (!vsPool.empty())
				{
					uint32_t currEle = vsPool.front();
					vsPool.pop();
					std::vector<uint32_t> &adjVs = pairVsInfo[currEle];
					for (int k = 0; k < adjVs.size(); ++k)
					{
						if (vsCounts[adjVs[k]] != 0)
						{
							vsCounts[adjVs[k]] = 0;
							tempVec.emplace_back(adjVs[k]);
							vsPool.push(adjVs[k]);
						}
					}
				}
			}
			std::sort(tempVec.begin(), tempVec.end());
			tempVec.erase(std::unique(tempVec.begin(), tempVec.end()), tempVec.end());
			sheet.vs_group.emplace_back(tempVec);
		}


		for (int j = 0; j < sheet.vs_pairs.size(); ++j)
		{
			std::vector<uint32_t> &currVsPair = sheet.vs_pairs[j];
			std::vector<uint32_t> tempVsGroup;
			if (vsCounts[currVsPair[0]] != 0 && vsCounts[currVsPair[1]] != 0)
			{
				std::queue<uint32_t> groupPool;
				/*tempVsGroup.emplace_back(currVsPair[0]); tempVsGroup.emplace_back(currVsPair[1]);
				vsCounts[currVsPair[0]] = 0; vsCounts[currVsPair[1]] = 0;*/
				groupPool.push(currVsPair[0]); groupPool.push(currVsPair[1]);
				vsCounts[currVsPair[0]] = 0; vsCounts[currVsPair[1]] = 0;
				while (!groupPool.empty())
				{
					uint32_t currV = groupPool.front(); groupPool.pop();
					//if (pairVsInfo[currV].size() == 1)
						tempVsGroup.emplace_back(currV);
					for (int k = 0; k < pairVsInfo[currV].size(); ++k)
					{
						if (vsCounts[pairVsInfo[currV][k]] != 0)
						{
							groupPool.push(pairVsInfo[currV][k]);
							vsCounts[pairVsInfo[currV][k]] = 0;
						}
					}
				}
				std::sort(tempVsGroup.begin(), tempVsGroup.end());
				tempVsGroup.erase(std::unique(tempVsGroup.begin(), tempVsGroup.end()), tempVsGroup.end());
				sheet.vs_group.emplace_back(tempVsGroup);
			}
		}
	}

	void Pipeline::FindFaceOppsiteEdges(BaseComplex &bc, uint32_t bfId, uint32_t *pair0, uint32_t *pair1)
	{
		std::vector<uint32_t> &fes = bc.Bf_[bfId].b_es;
		pair0[0] = fes[0];
		std::vector<uint32_t> &evs0 = bc.Be_[fes[0]].b_vs;
		for (int i = 1; i < fes.size(); ++i)
		{
			std::vector<uint32_t> &evs1 = bc.Be_[fes[i]].b_vs;
			if (evs0[0] != evs1[0] && evs0[0] != evs1[1] && evs0[1] != evs1[0] && evs0[1] != evs1[1])
			{
				pair0[1] = fes[i];
				pair1[0] = fes[1 + (i) % 3];
				pair1[1] = fes[1 + (i + 1) % 3];
				break; 
			}
		}
	}

	void Pipeline::FindRegularEdgePoints(BaseDataStructure::QuadMesh &qm, uint32_t preE, uint32_t preV, uint32_t vNum, std::vector<uint32_t> &finalVVec)
	{
		uint32_t nextE, nextV;
		finalVVec.emplace_back(preV);
		std::vector<uint32_t> resultVec;

		for (int i = 1; i < vNum; ++i)
		{
			std::vector<uint32_t> &efs0 = qm.Es_[preE].neighbor_fs;
			std::vector<uint32_t> &ves = qm.Vs_[preV].neighbor_es;
			for (int j = 0; j < ves.size(); ++j)
			{
				resultVec.clear();
				std::vector<uint32_t> &efs1 = qm.Es_[ves[j]].neighbor_fs;
				std::sort(efs0.begin(), efs0.end());
				std::sort(efs1.begin(), efs1.end());
				std::set_intersection(efs0.begin(), efs0.end(), efs1.begin(), efs1.end(), std::inserter(resultVec, resultVec.begin()));
				if (resultVec.empty())
				{
					nextE = ves[j];
					if (qm.Es_[nextE].vs[0] == preV)
						nextV = qm.Es_[nextE].vs[1];
					else
						nextV = qm.Es_[nextE].vs[0];
					finalVVec.emplace_back(nextV);
					break;
				}
			}
			preE = nextE; preV = nextV;
		}
	}

	void Pipeline::FindChordToBeCollapsedPair(QuadMesh &qm, BaseComplex &bc, Chord &chord)
	{
		/*std::vector<uint32_t>().swap(chord.qf_delete);
		std::vector<std::vector<uint32_t>>().swap( chord.vs_group);
		std::vector<uint32_t>().swap( chord.target_vs);
		std::vector<Eigen::Vector2d>().swap( chord.target_coords);*/
		chord.qf_delete.clear();
		chord.vs_group.clear();
		chord.target_vs.clear();
		chord.target_coords.clear();

		if (bc.Bv_[chord.keepedV1].boundary)
		{
			std::swap(chord.keepedV0, chord.keepedV1);
			std::swap(chord.keepedQV0, chord.keepedQV1);
		}

		uint32_t currId = chord.fid;
		uint32_t edgePair0[2], edgePair1[2];

		std::vector<bool> fFlag(qm.Fs_.size(), false);	//如果面已经被考虑过了，则标记为true;
		for (int i = 0; i < bc.Bf_[currId].fs_net.size(); ++i)
		{
			chord.qf_delete.emplace_back(bc.Bf_[currId].fs_net[i]);
			fFlag[bc.Bf_[currId].fs_net[i]] = true;
		}

		FindFaceOppsiteEdges(bc, currId, edgePair0, edgePair1);
		std::function<void(uint32_t *, uint32_t *)> findKeepedEdges = [&](uint32_t *keepedEdges0, uint32_t *keepedEdges1)
		{
			std::vector<uint32_t> &evs0 = bc.Be_[edgePair0[0]].b_vs;
			if (evs0[0] == chord.keepedV0 || evs0[1] == chord.keepedV0)
			{
				keepedEdges0[0] = edgePair0[0];
				keepedEdges1[0] = edgePair0[1];
			}
			else
			{
				keepedEdges0[0] = edgePair0[1];
				keepedEdges1[0] = edgePair0[0];
			}

			std::vector<uint32_t> &evs1 = bc.Be_[edgePair1[0]].b_vs;
			if (evs1[0] == chord.keepedV0 || evs1[1] == chord.keepedV0)
			{
				keepedEdges0[1] = edgePair1[0];
				keepedEdges1[1] = edgePair1[1];
			}
			else
			{
				keepedEdges0[1] = edgePair1[1];
				keepedEdges1[1] = edgePair1[0];
			}


		};
		uint32_t keepedEdges0[2], keepedEdges1[2];
		findKeepedEdges(keepedEdges0, keepedEdges1);

		std::function<void(uint32_t, uint32_t, uint32_t&, uint32_t&)> findNextEV = [&](uint32_t preE, uint32_t preV, uint32_t &nextE, uint32_t &nextV)	//here the index is in quad mesh (not base complex).
		{
			nextV = (uint32_t)-1;
			nextE = (uint32_t)-1;
			std::vector<uint32_t> &ves = qm.Vs_[preV].neighbor_es;
			std::vector<uint32_t> &efs0 = qm.Es_[preE].neighbor_fs;
			std::vector<uint32_t> resultVec;
			for (int i = 0; i < ves.size(); ++i)
			{
				resultVec.clear();
				std::vector<uint32_t> &efs1 = qm.Es_[ves[i]].neighbor_fs;
				std::sort(efs0.begin(), efs0.end());
				std::sort(efs1.begin(), efs1.end());
				std::set_intersection(efs0.begin(), efs0.end(), efs1.begin(), efs1.end(), std::inserter(resultVec, resultVec.begin()));
				if (resultVec.empty())
				{
					nextE = ves[i];
					if (qm.Es_[nextE].vs[0] == preV)
						nextV = qm.Es_[nextE].vs[1];
					else
						nextV = qm.Es_[nextE].vs[0];
				}
			}
		};

		std::function<void(uint32_t, uint32_t, uint32_t, std::vector<uint32_t>&, std::vector<uint32_t>&)> findPairVsInChord = [&](uint32_t startBv, uint32_t be0, uint32_t vsNum, std::vector<uint32_t>& finalPairs, std::vector<uint32_t> &finalEdges)	//vsNum是一条边上总共的点数
		{
			uint32_t startQv = bc.Bv_[startBv].qId;
			std::vector<uint32_t> &ves = qm.Vs_[startQv].neighbor_es;
			uint32_t qe0 = (uint32_t)-1;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (qm.Es_[ves[i]].bId == be0)
					qe0 = ves[i];
			}

			if (qe0 == (uint32_t)-1)
				exit(-7);

			uint32_t preE = qe0, preV;
			if (qm.Es_[preE].vs[0] == startQv)
				preV = qm.Es_[preE].vs[1];
			else
				preV = qm.Es_[preE].vs[0];

			if (vsNum >= 2)
			{
				finalPairs.emplace_back(preV);
				finalEdges.emplace_back(preE);
			}
			uint32_t nextV = (uint32_t)-1, nextE = (uint32_t)-1;
			for (int i = 1; i < vsNum - 1; ++i)
			{
				findNextEV(preE, preV, nextE, nextV);
				if (nextE == (uint32_t)-1 || nextV == (uint32_t)-1)
					exit(-9);
				finalPairs.emplace_back(nextV);
				finalEdges.emplace_back(nextE);
				preE = nextE; preV = nextV;
			}
		};

		if (bc.Be_[keepedEdges0[0]].vs_link.size() == bc.Be_[keepedEdges0[1]].vs_link.size())	//如果两条相交的边顶点数一致
		{
			std::vector<uint32_t> tempLV0, tempLV1, tempLE0, tempLE1;
			uint32_t vsNum = bc.Be_[edgePair0[0]].vs_link.size();
			findPairVsInChord(chord.keepedV0, keepedEdges0[0], vsNum, tempLV0, tempLE0);
			findPairVsInChord(chord.keepedV0, keepedEdges0[1], vsNum, tempLV1, tempLE1);
			findPairVsInChord(chord.keepedV1, keepedEdges1[1], vsNum - 1, tempLV0, tempLE0);
			findPairVsInChord(chord.keepedV1, keepedEdges1[0], vsNum - 1, tempLV1, tempLE1);
			
			if (tempLV0.size() != tempLV1.size())
				exit(-10);

			chord.vs_group.resize(tempLV0.size());
			for (int i = 0; i < tempLV0.size(); ++i)
			{
				chord.vs_group[i].resize(2);
				chord.vs_group[i][0] = tempLV0[i];
				chord.vs_group[i][1] = tempLV1[i];
			}
		}
		else
		{
			//由于开始的时候已经确定了如果chord.keepedV1是边界就换一下V0和V1，所以只要搞V1就行了。
			uint32_t keeped0Shorter, keeped0Longer, keeped1Shorter, keeped1Longer;
			uint32_t shortVsLinkNum, longVsLinkNum;
			if (bc.Be_[keepedEdges0[0]].vs_link.size() < bc.Be_[keepedEdges0[1]].vs_link.size())
			{
				shortVsLinkNum = bc.Be_[keepedEdges0[0]].vs_link.size();
				longVsLinkNum = bc.Be_[keepedEdges0[1]].vs_link.size();
				keeped0Shorter = keepedEdges0[0]; keeped0Longer = keepedEdges0[1];
				keeped1Shorter = keepedEdges1[0]; keeped1Longer = keepedEdges1[1];
			}
			else
			{
				shortVsLinkNum = bc.Be_[keepedEdges0[1]].vs_link.size();
				longVsLinkNum = bc.Be_[keepedEdges0[0]].vs_link.size();
				keeped0Shorter = keepedEdges0[1]; keeped0Longer = keepedEdges0[0];
				keeped1Shorter = keepedEdges1[1]; keeped1Longer = keepedEdges1[0];
			}

			std::vector<uint32_t> shortK0Vec, longK0Vec, shortK1Vec, longK1Vec;
			std::vector<uint32_t> shortK0EVec, longK0EVec, shortK1EVec, longK1EVec;
			findPairVsInChord(chord.keepedV0, keeped0Shorter, shortVsLinkNum, shortK0Vec, shortK0EVec);
			findPairVsInChord(chord.keepedV0, keeped0Longer, longVsLinkNum, longK0Vec, longK0EVec);
			findPairVsInChord(chord.keepedV1, keeped1Shorter, shortVsLinkNum, shortK1Vec, shortK1EVec);
			findPairVsInChord(chord.keepedV1, keeped1Longer, longVsLinkNum, longK1Vec, longK1EVec);

			uint32_t newStartP = longK1Vec[longVsLinkNum - shortVsLinkNum - 1];
			uint32_t newStartE = (uint32_t)-1;
			std::vector<uint32_t> &ves = qm.Vs_[newStartP].neighbor_es;
			for (int i = 0; i < ves.size(); ++i)
			{
				std::vector<uint32_t> &efs = qm.Es_[ves[i]].neighbor_fs;
				if (efs.size() != 2)
					continue;
				if (qm.Fs_[efs[0]].bId == chord.fid && qm.Fs_[efs[1]].bId == chord.fid)
				{
					newStartE = ves[i];
					break;
				}
			}
			if (newStartE == (uint32_t)-1)
				exit(-10);

			uint32_t preE = newStartE, preV = (uint32_t)-1, nextE = (uint32_t)-1, nextV = (uint32_t)-1;
			if (qm.Es_[newStartE].vs[0] == newStartP)
				preV = qm.Es_[newStartE].vs[1];
			else
				preV = qm.Es_[newStartE].vs[0];

			std::vector<uint32_t> newVs;
			FindRegularEdgePoints(qm, preE, preV, shortVsLinkNum - 1, newVs);

			//记录每个顶点的配对信息
			std::vector<uint32_t> adjVNum(qm.Vs_.size(), 0);
			std::vector<std::vector<uint32_t>> adjVPair(qm.Vs_.size(), std::vector<uint32_t>());
			std::vector<uint32_t> shouldGroupVs;
			for (int i = 0; i < shortVsLinkNum-2; ++i)
			{
				++adjVNum[shortK0Vec[i]]; ++adjVNum[longK0Vec[i]];
				adjVPair[shortK0Vec[i]].emplace_back(longK0Vec[i]);
				adjVPair[longK0Vec[i]].emplace_back(shortK0Vec[i]);
				shouldGroupVs.emplace_back(shortK0Vec[i]); shouldGroupVs.emplace_back(longK0Vec[i]);
			}
			for (int i = 0; i < shortVsLinkNum - 1; ++i)
			{
				uint32_t vLeft = shortK1Vec[i], vMid = newVs[i], vRight = longK1Vec[longVsLinkNum - shortVsLinkNum + i];
				++adjVNum[vLeft]; adjVNum[vMid] += 2; ++adjVNum[vRight];
				adjVPair[vLeft].emplace_back(vMid);
				adjVPair[vMid].emplace_back(vLeft);
				adjVPair[vRight].emplace_back(vMid);
				adjVPair[vMid].emplace_back(vRight);
				shouldGroupVs.emplace_back(vLeft);
				shouldGroupVs.emplace_back(vMid);
				shouldGroupVs.emplace_back(vRight);
			}
			uint32_t kq1id = bc.Bv_[chord.keepedV1].qId;
			++adjVNum[kq1id]; ++adjVNum[newStartP];
			adjVPair[kq1id].emplace_back(newStartP);
			adjVPair[newStartP].emplace_back(kq1id);
			shouldGroupVs.emplace_back(kq1id);
			shouldGroupVs.emplace_back(newStartP);

			std::vector<uint32_t> startVVec0, startEVec0, startVVec1, startEVec1;
			for (int i = shortVsLinkNum - 2; i < longVsLinkNum - 1; ++i)
			{
				startVVec0.emplace_back(longK0Vec[i]);
			}
			for (int i = shortVsLinkNum - 1; i < longVsLinkNum - 1; ++i)
			{
				startEVec0.emplace_back(longK0EVec[i]);
			}
			for (int i = longVsLinkNum - shortVsLinkNum-1; i >= 0; --i)
			{
				startVVec1.emplace_back(longK1Vec[i]);
				startEVec1.emplace_back(longK1EVec[i]);
			}
			startVVec1.emplace_back(bc.Bv_[chord.keepedV1].qId);

			std::function<bool(std::vector<uint32_t>&, std::vector<uint32_t>&)> findAllPairsInSubsheet = [&](std::vector<uint32_t>& startVVec, std::vector<uint32_t> &startEVec)->bool
			{
				if (startVVec.empty() || startEVec.empty())
					exit(-11);
				uint32_t firstE = startEVec[0], lastE = startEVec[startEVec.size() - 1];

				for (int i=0; i < startEVec.size(); ++i)
				{
					if (qm.Es_[startEVec[i]].boundary)
						return false;
				}

				std::vector<uint32_t> nextFs(startEVec.size()), nextEs(startEVec.size()), nextVs(startVVec.size());
				for (int i = 0; i < startEVec.size(); ++i)
				{
					std::vector<uint32_t> &efs = qm.Es_[startEVec[i]].neighbor_fs;
					if (fFlag[efs[0]] && !fFlag[efs[1]])
					{
						nextFs[i] = efs[1];
						fFlag[efs[1]] = true;
						chord.qf_delete.emplace_back(efs[1]);
					}
					else if (fFlag[efs[1]] && !fFlag[efs[0]])
					{
						nextFs[i] = efs[0];
						fFlag[efs[0]] = true;
						chord.qf_delete.emplace_back(efs[0]);
					}
					else
					{
						//exit(-12);
						return false;
					}
				}
				for (int i = 0; i < startEVec.size(); ++i)
				{
					uint32_t currF = nextFs[i];
					std::vector<uint32_t> &fes = qm.Fs_[currF].es;
					std::vector<uint32_t> &evs = qm.Es_[startEVec[i]].vs;
					for (int j = 0; j < 4; ++j)
					{
						std::vector<uint32_t> &currEvs = qm.Es_[fes[j]].vs;
						if (currEvs[0] != evs[0] && currEvs[0] != evs[1] && currEvs[1] != evs[0] && currEvs[1] != evs[1])
						{
							nextEs[i] = fes[j];
							break;
						}
					}
				}
				for (int i = 0; i < startEVec.size(); ++i)
				{
					std::vector<uint32_t> &efs = qm.Es_[nextEs[i]].neighbor_fs;
					if (efs.size() != 2)
						continue;
					if (fFlag[efs[0]] && fFlag[efs[1]])
						return false;
				}
				std::vector<uint32_t> tempVec;
				std::sort(qm.Vs_[startVVec[0]].neighbor_vs.begin(), qm.Vs_[startVVec[0]].neighbor_vs.end());
				std::sort(qm.Es_[nextEs[0]].vs.begin(), qm.Es_[nextEs[0]].vs.end());
				std::set_intersection(qm.Vs_[startVVec[0]].neighbor_vs.begin(), qm.Vs_[startVVec[0]].neighbor_vs.end(), qm.Es_[nextEs[0]].vs.begin(), qm.Es_[nextEs[0]].vs.end(), std::inserter(tempVec, tempVec.begin()));
				if (tempVec.empty())
					exit(-13);
				nextVs[0] = tempVec[0];
				for (int i = 0; i < startEVec.size(); ++i)
				{
					std::vector<uint32_t> &evs = qm.Es_[nextEs[i]].vs;
					if (evs[0] == nextVs[i])
						nextVs[i + 1] = evs[1];
					else
						nextVs[i+1] = evs[0];
				}

				++adjVNum[nextVs[0]]; ++adjVNum[nextVs[nextVs.size()-1]];
				adjVPair[nextVs[0]].emplace_back(nextVs[nextVs.size() - 1]);
				adjVPair[nextVs[nextVs.size() - 1]].emplace_back(nextVs[0]);
				shouldGroupVs.emplace_back(nextVs[0]);
				shouldGroupVs.emplace_back(nextVs[nextVs.size() - 1]);

				startEVec = nextEs;
				startVVec = nextVs;
				return true;
			};

			while (findAllPairsInSubsheet(startVVec0, startEVec0))
				;
			while (findAllPairsInSubsheet(startVVec1, startEVec1))
				;

			//所有配对信息都找到之后，把他们配对
			shouldGroupVs.erase(std::unique(shouldGroupVs.begin(), shouldGroupVs.end()), shouldGroupVs.end());
			for (int i = 0; i < shouldGroupVs.size(); ++i)
			{
				if (adjVNum[shouldGroupVs[i]] == 0)
					continue;

				std::vector<uint32_t> tempGroup;
				adjVNum[shouldGroupVs[i]] = 0;
				std::queue<uint32_t> groupPool;
				groupPool.push(shouldGroupVs[i]);
				while (!groupPool.empty())
				{
					uint32_t currEle = groupPool.front();
					groupPool.pop();
					//if (adjVPair[currEle].size() == 1)
						tempGroup.emplace_back(currEle);

					for (int j = 0; j < adjVPair[currEle].size(); ++j)
					{
						if (adjVNum[adjVPair[currEle][j]] != 0)
						{
							groupPool.push(adjVPair[currEle][j]);
							adjVNum[adjVPair[currEle][j]] = 0;
						}
					}
				}
				std::sort(tempGroup.begin(), tempGroup.end());
				tempGroup.erase(std::unique(tempGroup.begin(), tempGroup.end()), tempGroup.end());
				chord.vs_group.emplace_back(tempGroup);
			}
		}
	}

	bool Pipeline::JudgeCollapseCriterion(BaseDataStructure::QuadMesh &qm, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, RankingTuple3 &rt3, FeatureConstraints &fcc)
	{
		ElementType eleType = std::get<1>(rt3);
		uint32_t currId = std::get<2>(rt3);
		std::vector<bool> fFlag(qm.Fs_.size(), false);

		collapseCorner_ = false;
		if (eleType == ElementType::SHEET)
		{
			Sheet &currObj = sheets[currId];
			for (int i = 0; i < currObj.qf_delete.size(); ++i)
			{
				fFlag[currObj.qf_delete[i]] = true;
			}

			//build valenceFilter
			currObj.valenceFilter = false;
			valenceReturn = false;
			for (int i = 0; i < currObj.vs_group.size(); ++i)
			{
				std::vector<uint32_t> &currVsGroup = currObj.vs_group[i];

				uint32_t minValence = 1000;
				std::vector<uint32_t> fsAfterCo;
				for (int j = 0; j < currVsGroup.size(); ++j)
				{
					uint32_t currV = currVsGroup[j];
					if (minValence > qm.Vs_[currV].neighbor_fs.size())
						minValence = qm.Vs_[currV].neighbor_fs.size();

					std::vector<uint32_t> &vFs = qm.Vs_[currV].neighbor_fs;
					for (int k = 0; k < vFs.size(); ++k)
					{
						if (!fFlag[vFs[k]])
							fsAfterCo.emplace_back(vFs[k]);
					}
				}
				std::sort(fsAfterCo.begin(), fsAfterCo.end());
				fsAfterCo.erase(std::unique(fsAfterCo.begin(), fsAfterCo.end()), fsAfterCo.end());
				if (fsAfterCo.size() < minValence)
				{
					currObj.valenceFilter = true;
//					return false;
				}
				
#if USE_VALENCE_WEIGHT
				if (!currObj.haveAddValenceWidget)
				{
					//currObj.haveAddValenceWidget = true;
					//uint32_t currTarget = currObj.target_vs[i];
					bool hasBoundary = false;
					for (int j = 0; j < currVsGroup.size(); ++j)
					{
						if (qm.Vs_[currVsGroup[j]].boundary)
						{
							hasBoundary = true;
							break;
						}
					}
					if (hasBoundary && fsAfterCo.size() > 2)
					{
						currObj.weight += VALENCE_WEIGHT * (fsAfterCo.size() - 2) * averageEdgeLength_;
						std::get<0>(rt3) += VALENCE_WEIGHT * (fsAfterCo.size() - 2) * averageEdgeLength_;
						valenceReturn = true;
					}
					else if (!hasBoundary && fsAfterCo.size() > 4)
					{
						currObj.weight += VALENCE_WEIGHT * (fsAfterCo.size() - 4) * averageEdgeLength_;
						std::get<0>(rt3) += VALENCE_WEIGHT * (fsAfterCo.size() - 4) * averageEdgeLength_;
						valenceReturn = true;
					}
				}
#endif
			}
#if USE_VALENCE_WEIGHT
			currObj.haveAddValenceWidget = true;
			if (valenceReturn)
				return false;
#endif

			//target_vs and target_coords
			//std::vector<uint32_t>().swap(currObj.target_vs);
			//std::vector<Eigen::Vector2d>().swap(currObj.target_coords);
			currObj.target_coords.clear();
			currObj.target_vs.clear();

			currObj.target_vs.resize(currObj.vs_group.size());
			currObj.target_coords.resize(currObj.vs_group.size());
			for (int i = 0; i < currObj.vs_group.size(); ++i)
			{
				std::vector<uint32_t> &currVsGroup = currObj.vs_group[i];
				std::vector<uint32_t> curvesVs;
				std::vector<uint32_t> cornerVs;
				std::vector<uint32_t> curvesVsIds;
				//如果存在两个都是coner或者两个位于不同的边上，则不搞。
				for (int j = 0; j < currVsGroup.size(); ++j)
				{
					if (fcc.V_types[currVsGroup[j]] == -1)
						cornerVs.emplace_back(currVsGroup[j]);
					else if (fcc.V_types[currVsGroup[j]] >= 0)
					{
						curvesVs.emplace_back(fcc.V_types[currVsGroup[j]]);
						curvesVsIds.emplace_back(currVsGroup[j]);
					}
				}
				std::sort(curvesVs.begin(), curvesVs.end());
				curvesVs.erase(std::unique(curvesVs.begin(), curvesVs.end()), curvesVs.end());
				std::sort(cornerVs.begin(), cornerVs.end());
				cornerVs.erase(std::unique(cornerVs.begin(), cornerVs.end()), cornerVs.end());

#if !COLLAPSE_CORNER
				if (cornerVs.size() > 1 || curvesVs.size() > 1)
					return false;

				if (cornerVs.size()!=0)	//如果角点数量不为0（其实就是1）
				{
					currObj.target_vs[i] = cornerVs[0];
					currObj.target_coords[i] = qm.V_[cornerVs[0]];
				}
				else if (curvesVs.size()!=0)
				{
					currObj.target_vs[i] = curvesVsIds[0];
					Eigen::Vector2d averagePos(0, 0);
					for (int j = 0; j < curvesVsIds.size(); ++j)
					{
						averagePos += qm.V_[curvesVsIds[j]];
					}
					averagePos /= curvesVsIds.size();
					currObj.target_coords[i] = averagePos;
				}
				else
				{
					currObj.target_vs[i] = currVsGroup[0];
					Eigen::Vector2d currCoord(0, 0);
					for (int j = 0; j < currVsGroup.size(); ++j)
					{
						currCoord += qm.V_[currVsGroup[j]];
					}
					currCoord /= currVsGroup.size();
					currObj.target_coords[i] = currCoord;
				}
#else
				if (cornerVs.size() > 1 || curvesVs.size() > 1)
					collapseCorner_ = true;
				if (!cornerVs.empty())
					currObj.target_vs[i] = cornerVs[0];
				else if (!curvesVs.empty())
					currObj.target_vs[i] = curvesVsIds[0];
				else
					currObj.target_vs[i] = currVsGroup[0];

				Eigen::Vector2d currCoord(0, 0);
				uint32_t countt = 0;
				if (cornerVs.size() > 0 || curvesVs.size() > 0)
				{
					for (int j = 0; j < currVsGroup.size(); ++j)
					{
						if (!qm.Vs_[currVsGroup[j]].boundary)
							continue;
						currCoord += qm.V_[currVsGroup[j]];
						++countt;
					}
				}
				else
				{
					for (int j = 0; j < currVsGroup.size(); ++j)
					{
						currCoord += qm.V_[currVsGroup[j]];
						++countt;
					}
				}
				currCoord /= countt;
				currObj.target_coords[i] = currCoord;
#endif
			}
		}
		else if (eleType == ElementType::CHORD)
		{
			Chord &currObj = chords[currId];
			for (int i = 0; i < currObj.qf_delete.size(); ++i)
			{
				fFlag[currObj.qf_delete[i]] = true;
			}

			//build valenceFilter
			currObj.valenceFilter = false;
			valenceReturn = false;
			for (int i = 0; i < currObj.vs_group.size(); ++i)
			{
				std::vector<uint32_t> &currVsGroup = currObj.vs_group[i];

				uint32_t minValence = 1000;
				std::vector<uint32_t> fsAfterCo;
				for (int j = 0; j < currVsGroup.size(); ++j)
				{
					uint32_t currV = currVsGroup[j];
					if (minValence > qm.Vs_[currV].neighbor_fs.size())
						minValence = qm.Vs_[currV].neighbor_fs.size();

					std::vector<uint32_t> &vFs = qm.Vs_[currV].neighbor_fs;
					for (int k = 0; k < vFs.size(); ++k)
					{
						if (!fFlag[vFs[k]])
							fsAfterCo.emplace_back(vFs[k]);
					}
				}
				std::sort(fsAfterCo.begin(), fsAfterCo.end());
				fsAfterCo.erase(std::unique(fsAfterCo.begin(), fsAfterCo.end()), fsAfterCo.end());
				uint32_t keeped1QId = currObj.keepedQV1;
				if (fsAfterCo.size() < minValence && std::find(currVsGroup.begin(), currVsGroup.end(), keeped1QId)==currVsGroup.end())
				{
					currObj.valenceFilter = true;
//					return false;
				}
				/*std::uint32_t currTargetV = currObj.target_vs[i];
				if (!qm.Vs_[currTargetV].boundary && fsAfterCo.size() == 2)
					qm.Vs_[currTargetV].keeped = false;
				else if (qm.Vs_[currTargetV].boundary && !fcc.V_types[currTargetV] == -1 && fsAfterCo.size() == 0)
					qm.Vs_[currTargetV].keeped = false;*/

#if USE_VALENCE_WEIGHT
				if (!currObj.haveAddValenceWidget)
				{
					//currObj.haveAddValenceWidget = true;
					bool hasBoundary = false;
					for (int j = 0; j < currVsGroup.size(); ++j)
					{
						if (qm.Vs_[currVsGroup[j]].boundary)
						{
							hasBoundary = true;
							break;
						}
					}
					if (hasBoundary && fsAfterCo.size() > 2)
					{
						currObj.weight += VALENCE_WEIGHT * (fsAfterCo.size() - 2) * averageEdgeLength_;
						std::get<0>(rt3) += VALENCE_WEIGHT * (fsAfterCo.size() - 2) * averageEdgeLength_;
						valenceReturn = true;
					}
					else if (!hasBoundary && fsAfterCo.size() > 4)
					{
						currObj.weight += VALENCE_WEIGHT * (fsAfterCo.size() - 4) * averageEdgeLength_;
						std::get<0>(rt3) += VALENCE_WEIGHT * (fsAfterCo.size() - 4) * averageEdgeLength_;
						valenceReturn = true;
					}
				}
#endif
			}
#if USE_VALENCE_WEIGHT
			currObj.haveAddValenceWidget = true;
			if (valenceReturn)
				return false;
#endif
			/*uint32_t keeped0QId = currObj.keepedQV0;
			if (!qm.Vs_[keeped0QId].boundary && qm.Vs_[keeped0QId].neighbor_es.size() == 3)
				qm.Vs_[keeped0QId].keeped = false;
			else if (qm.Vs_[keeped0QId].boundary && fcc.V_types[keeped0QId] != -1 && qm.Vs_[keeped0QId].neighbor_fs.size() == 1)
				qm.Vs_[keeped0QId].keeped = false;*/

			//target_vs and target_coords
			//std::vector<uint32_t>().swap(currObj.target_vs);
			//std::vector<Eigen::Vector2d>().swap(currObj.target_coords);
			currObj.target_vs.clear();
			currObj.target_coords.clear();

			currObj.target_vs.resize(currObj.vs_group.size());
			currObj.target_coords.resize(currObj.vs_group.size());
			for (int i = 0; i < currObj.vs_group.size(); ++i)
			{
				std::vector<uint32_t> &currVsGroup = currObj.vs_group[i];
				std::vector<uint32_t> curvesVs;
				std::vector<uint32_t> cornerVs;
				std::vector<uint32_t> curvesVsIds;
				//如果存在两个都是coner或者两个位于不同的边上，则不搞。
				for (int j = 0; j < currVsGroup.size(); ++j)
				{
					if (fcc.V_types[currVsGroup[j]] == -1)
						cornerVs.emplace_back(currVsGroup[j]);
					else if (fcc.V_types[currVsGroup[j]] >= 0)
					{
						curvesVs.emplace_back(fcc.V_types[currVsGroup[j]]);
						curvesVsIds.emplace_back(currVsGroup[j]);
					}
				}
				std::sort(curvesVs.begin(), curvesVs.end());
				curvesVs.erase(std::unique(curvesVs.begin(), curvesVs.end()), curvesVs.end());
				std::sort(cornerVs.begin(), cornerVs.end());
				cornerVs.erase(std::unique(cornerVs.begin(), cornerVs.end()), cornerVs.end());
				if (cornerVs.size() > 1 || curvesVs.size() > 1)
					return false;

				if (cornerVs.size())	//如果角点数量不为0（其实就是1）
				{
					currObj.target_vs[i] = cornerVs[0];
					currObj.target_coords[i] = qm.V_[cornerVs[0]];
				}
				else if (curvesVs.size())
				{
					currObj.target_vs[i] = curvesVsIds[0];
					Eigen::Vector2d averagePos(0, 0);
					for (int j = 0; j < curvesVsIds.size(); ++j)
					{
						averagePos += qm.V_[curvesVsIds[j]];
					}
					averagePos /= curvesVsIds.size();
					currObj.target_coords[i] = averagePos;
				}
				else
				{
					currObj.target_vs[i] = currVsGroup[0];
					Eigen::Vector2d currCoord(0, 0);
					for (int j = 0; j < currVsGroup.size(); ++j)
					{
						currCoord += qm.V_[currVsGroup[j]];
					}
					currCoord /= currVsGroup.size();
					currObj.target_coords[i] = currCoord;
				}
			}
		}
		
		return true;
	}

	void Pipeline::GrowRegion(QuadMesh &qm, CollapseInfo *ci)
	{
		std::vector<int> fFlag(qm.Fs_.size(), -2);	//-2 表示初始情况，即没有被标记；-1表示该四边形是要被干掉的四边形；0表示该四边形是后面长出来的四边形。
		for (int i = 0; i < ci->hsToBeKilled.size(); ++i)
		{
			fFlag[ci->hsToBeKilled[i]] = -1;
		}

		std::vector<uint32_t> currQVec, nextQVec;
		for (int i = 0; i < ci->hsToBeKilled.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsToBeKilled[i]].vs;
			for (int j = 0; j < 4; ++j)
			{
				std::vector<uint32_t> &vfs = qm.Vs_[fvs[j]].neighbor_fs;
				for (int k = 0; k < vfs.size(); ++k)
				{
					if (fFlag[vfs[k]] == -2)
					{
						fFlag[vfs[k]] = 0;
						ci->hsSeveralring.emplace_back(vfs[k]);
						currQVec.emplace_back(vfs[k]);
					}
				}
			}
		}

		uint32_t iterateNum = 0;
		while (iterateNum < QUAD_RINGS - 1)
		{
			for (int i = 0; i < currQVec.size(); ++i)
			{
				std::vector<uint32_t> &fvs = qm.Fs_[currQVec[i]].vs;
				for (int j = 0; j < 4; ++j)
				{
					std::vector<uint32_t> &vfs = qm.Vs_[fvs[j]].neighbor_fs;
					for (int k = 0; k < vfs.size(); ++k)
					{
						if (fFlag[vfs[k]] == -2)
						{
							fFlag[vfs[k]] = 0;
							ci->hsSeveralring.emplace_back(vfs[k]);
							nextQVec.emplace_back(vfs[k]);
						}
					}
				}
			}
			currQVec.swap(nextQVec);
			nextQVec.clear();
			++iterateNum;
		}

		std::sort(ci->hsSeveralring.begin(), ci->hsSeveralring.end());
		ci->hsSeveralring.erase(std::unique(ci->hsSeveralring.begin(), ci->hsSeveralring.end()), ci->hsSeveralring.end());

		//找固定点。
		for (int i = 0; i < ci->hsSeveralring.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsSeveralring[i]].vs;
			for (int j = 0; j < 4; ++j)
			{
				std::vector<uint32_t> &vfs = qm.Vs_[fvs[j]].neighbor_fs;
				int vState = 0;
				for (int k = 0; k < vfs.size(); ++k)
				{
					if (fFlag[vfs[k]] == -2)
					{
						vState = 1;
						break;
					}
				}

				if (vState == 1)
					ci->vsNotMove.emplace_back(fvs[j]);
			}
		}
		std::sort(ci->vsNotMove.begin(), ci->vsNotMove.end());
		ci->vsNotMove.erase(std::unique(ci->vsNotMove.begin(), ci->vsNotMove.end()), ci->vsNotMove.end());
	}

	void Pipeline::BuildCollapseInfo(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, RankingTuple3 &rt3, CollapseInfo *ci)
	{
#if CHANGE_QUAD_RINGS
		std::function<void(Sheet &)> getSheetNewRing = [&](Sheet &currSheet)
		{
			uint32_t midBE = currSheet.b_middle_es[0];
			int midBELength = bc.Be_[midBE].es_link.size();
			QUAD_RINGS = 4 * midBELength;
		};
		std::function<void(Chord &)> getChordNewRing = [&](Chord &currChord)
		{
			std::vector<uint32_t> &fes = bc.Bf_[currChord.fid].b_es;
			uint32_t be0 = fes[0], be1;
			std::vector<uint32_t> &evs0 = bc.Be_[be0].b_vs;
			for (int i = 1; i < fes.size(); ++i)
			{
				std::vector<uint32_t> &evs1 = bc.Be_[fes[i]].b_vs;
				if (evs0[0] == evs1[0] || evs0[0] == evs1[1] || evs0[1] == evs1[0] || evs0[1] == evs1[1])
				{
					be1 = fes[i];
					break;
				}
			}
			int e0Length = bc.Be_[be0].es_link.size();
			int e1Length = bc.Be_[be1].es_link.size();
			QUAD_RINGS = 4 * std::max(e0Length, e1Length);
		};
#endif

		ElementType eleType = std::get<1>(rt3);
		uint32_t currId = std::get<2>(rt3);
		if (eleType == ElementType::SHEET)
		{
			Sheet &currObj = sheets[currId];
#if CHANGE_QUAD_RINGS
			getSheetNewRing(currObj);
#endif
			ci->vsGroup = currObj.vs_group;
			ci->targetCoords = currObj.target_coords;
			ci->targetVs = currObj.target_vs;
			ci->hsToBeKilled = currObj.qf_delete;
		}
		else if (eleType == ElementType::CHORD)
		{
			Chord &currObj = chords[currId];
#if CHANGE_QUAD_RINGS
			getChordNewRing(currObj);
#endif
			ci->vsGroup = currObj.vs_group;
			ci->targetCoords = currObj.target_coords;
			ci->targetVs = currObj.target_vs;
			ci->hsToBeKilled = currObj.qf_delete;

		}
		GrowRegion(qm, ci);
	}

	bool Pipeline::BuildNewTopologyMesh(QuadMesh &inputQm, QuadMesh *outputQm, BaseComplex *outputBc, CollapseInfo *ci, std::vector<uint32_t> &vNewToOld, std::vector<uint32_t> &vOldToNew, FeatureConstraints &fcc)
	{
		//标记面是否保留
		std::vector<bool> isFRemain(inputQm.Fs_.size(), true);
		for (int i = 0; i < ci->hsToBeKilled.size(); ++i)
		{
			isFRemain[ci->hsToBeKilled[i]] = false;
		}

		//标记点是否保留
		std::vector<bool> isVRemain(inputQm.Vs_.size(), true);
		for (int i = 0; i < ci->hsToBeKilled.size(); ++i)
		{
			std::vector<uint32_t> &fvs = inputQm.Fs_[ci->hsToBeKilled[i]].vs;
			for (int j = 0; j < 4; ++j)
			{
				std::vector<uint32_t> &vfs = inputQm.Vs_[fvs[j]].neighbor_fs;
				bool isAllAdjFaceKilled = true;
				for (int k = 0; k < vfs.size(); ++k)
				{
					if (isFRemain[vfs[k]])
					{
						isAllAdjFaceKilled = false;
						break;
					}
				}
				if (isAllAdjFaceKilled)
					isVRemain[fvs[j]] = false;
			}
		}
		for (int i = 0; i < ci->vsGroup.size(); ++i)
		{
			for (int j = 0; j < ci->vsGroup[i].size(); ++j)
			{
				isVRemain[ci->vsGroup[i][j]] = false;
			}
		}
		for (int i = 0; i < ci->targetVs.size(); ++i)
		{
			isVRemain[ci->targetVs[i]] = true;
		}

		//std::vector<uint32_t>().swap(vNewToOld);
		//std::vector<uint32_t>().swap(vOldToNew);
		vNewToOld.clear();
		vOldToNew.clear();
		vOldToNew.resize(inputQm.Vs_.size(), (uint32_t)-1);

		uint32_t currVIndex = 0;
		for (int i = 0; i < inputQm.Vs_.size(); ++i)
		{
			if (isVRemain[i] == true)
			{
				vOldToNew[i] = currVIndex++;
				vNewToOld.emplace_back(i);
				outputQm->V_.emplace_back(inputQm.V_[i]);

				QuadVertex qv;
				qv.id = currVIndex - 1;
				qv.bId = (uint32_t)-1;
				outputQm->Vs_.emplace_back(qv);
			}
		}
		for (int i = 0; i < ci->targetVs.size(); ++i)
		{
			outputQm->V_[vOldToNew[ci->targetVs[i]]] = ci->targetCoords[i];
		}
		for (int i = 0; i < ci->vsGroup.size(); ++i)
		{
			std::vector<uint32_t> &vsGroup = ci->vsGroup[i];
			for (int j = 0; j < vsGroup.size(); ++j)
			{
				vOldToNew[vsGroup[j]] = vOldToNew[ci->targetVs[i]];
			}
		}

		ComputeVNewToOld(vOldToNew, vNewToOld, fcc);	//为了特征映回去还是特征，特此写了这个函数。

		newQmOptimizedFs_.clear();
		newQmOptimizedFs_.reserve(inputQm.Fs_.size());
		uint32_t *fOldToNew = new uint32_t[inputQm.Fs_.size()];
		std::memset(fOldToNew, (uint32_t)-1, inputQm.Fs_.size() * sizeof(uint32_t));
		uint32_t currFIndex = 0;
		for (int i = 0; i < inputQm.Fs_.size(); ++i)
		{
			if (isFRemain[i] == false)
				continue;

			fOldToNew[i] = currFIndex;
			QuadFace qf;
			qf.id = currFIndex++;
			qf.bId = (uint32_t)-1;
			std::vector<uint32_t> &fvs = inputQm.Fs_[i].vs;
			for (int j = 0; j < 4; ++j)
			{
				uint32_t newVId = vOldToNew[fvs[j]];
				if (newVId == (uint32_t)-1)
				{
					std::cout << "Error in build new mesh! " << std::endl;
					exit(-13);
				}
				qf.vs.emplace_back(newVId);
				outputQm->Vs_[newVId].neighbor_fs.emplace_back(currFIndex - 1);
			}
			outputQm->Fs_.emplace_back(qf);
		}
		for (int i = 0; i < ci_->hsSeveralring.size(); ++i)
		{
			uint32_t currNewFId = fOldToNew[ci_->hsSeveralring[i]];
			if (currNewFId == (uint32_t)-1)
				exit(-908);

			newQmOptimizedFs_.emplace_back(currNewFId);
		}
		delete[] fOldToNew;
		if (!outputQm->BuildConnectivity())
			return false;
		/*outputBc->ExtractBaseComplex(outputQm);
		if (outputBc->isError_)
			return false;*/

		if (!TopologyCheckWithoutValence(*outputQm))
		{
			//std::cout << "Didn't pass topology check! " << std::endl;
			return false;
		}

#if USE_NO_INVERSE
		if (outputQm->JudgeQuadMeshJacobi() && outputQm->minJacobi_ > MIN_ALLOWED_JACOBI)
		{
			noInverse = true;
		}
#endif
#ifdef OUTPUT_MID_MESH
		QuadMeshIO qmio;
		qmio.WriteQuadMesh(outputQm, "C:\\Users\\ChiZhang\\Desktop\\collapsedMesh.obj");
#endif
		return true;
	}

	bool Pipeline::TopologyCheckWithoutValence(QuadMesh &qm)
	{
		uint32_t vNum = qm.Vs_.size();
		uint32_t eNum = qm.Es_.size();
		uint32_t fNum = qm.Fs_.size();

		if (vNum + fNum - eNum != initialGenus_)
		{
			return false;
		}

		//Base structure check.
		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fes = qm.Fs_[i].es;
			std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
			if (fes.size() != 4 || fvs.size() != 4)
				return false;
		}
		for (int i = 0; i < qm.Es_.size(); ++i)
		{
			std::vector<uint32_t> &efs = qm.Es_[i].neighbor_fs;
			std::vector<uint32_t> &evs = qm.Es_[i].vs;
			if (efs.size() == 0 || evs.size() != 2)
				return false;
		}
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			std::vector<uint32_t> &vfs = qm.Vs_[i].neighbor_fs;
			std::vector<uint32_t> &ves = qm.Vs_[i].neighbor_es;
			if (vfs.empty() || ves.empty())
				return false;
		}

		//Quad mesh 2-manifold check.
		for (int i = 0; i < qm.Es_.size(); ++i)
		{
			uint32_t eFsNum = qm.Es_[i].neighbor_fs.size();
			if (eFsNum != 1 && eFsNum != 2)
				return false;
		}
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			uint32_t vEsNum = qm.Vs_[i].neighbor_es.size();
			uint32_t vFsNum = qm.Vs_[i].neighbor_fs.size();
			if ((qm.Vs_[i].boundary && vEsNum - vFsNum != 1) || (!qm.Vs_[i].boundary && vEsNum - vFsNum != 0))
				return false;
		}

		//Base complex 2-manifold check.
		/*for (int i = 0; i < bc.Be_.size(); ++i)
		{
			uint32_t beFsNum = bc.Be_[i].b_neighbor_fs.size();
			if (beFsNum != 1 && beFsNum != 2)
				return false;
		}
		for (int i = 0; i < bc.Bv_.size(); ++i)
		{
			uint32_t bvEsNum = bc.Bv_[i].b_neighbor_es.size();
			uint32_t bvFsNum = bc.Bv_[i].b_neighbor_fs.size();
			if ((bc.Bv_[i].boundary&&bvEsNum - bvFsNum != 1) || (!bc.Bv_[i].boundary&&bvEsNum - bvFsNum != 0))
				return false;
		}
*/
		/*for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if ((!qm.Vs_[i].boundary && qm.Vs_[i].neighbor_es.size() <= 2)
				|| qm.Vs_[i].neighbor_es.size() <= 1)
				return false;
		}*/
		return true;
	}

	bool Pipeline::DeleteUselessElements(QuadMesh &qm, std::vector<uint32_t> &vNewToNewNew, std::vector<uint32_t> &vNewNewToNew, std::vector<uint32_t> &fNewToNewNew, std::vector<uint32_t> &fNewNewToNew)
	{
		fNewToNewNew.resize(qm.Fs_.size(), (uint32_t)-1);
		fNewNewToNew.reserve(qm.Fs_.size());
		std::vector<int> vFlag(qm.Vs_.size(), -2); //-2表示初始状态， -1表示这是一个内部度为2的点
		std::vector<uint32_t> valence2P;
		bool haveUselessEle = false;
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if (!qm.Vs_[i].boundary && qm.Vs_[i].neighbor_fs.size() == 2)
			{
				vFlag[i] = -1;
				valence2P.emplace_back(i);
				haveUselessEle = true;
			}
		}
		if (!haveUselessEle)
			return true;

		vNewToNewNew.resize(qm.Vs_.size(), (uint32_t)-1);
		uint32_t startCount = 0;
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if (vFlag[i] == -2)
			{
				vNewToNewNew[i] = startCount++;
				vNewNewToNew.emplace_back(i);
			}
		}

		std::vector<bool> fFlag(qm.Fs_.size(), false); //判断是否已经找过了。
		std::queue<uint32_t> fV2Pool;
		std::vector<std::vector<uint32_t>> fV2SetVec;
		std::vector<int> fFindFlag(qm.Fs_.size(), -2);//-2表示这个面不是一个包含度2点的面；>=0记录了一个包含度2顶点的面所在的fV2SetVec的下标，-1表示虽然这个面包含度2顶点，但是已经被找过了。
		for (int i = 0; i < valence2P.size(); ++i)
		{
			uint32_t f0 = qm.Vs_[valence2P[i]].neighbor_fs[0];
			if (fFlag[f0])
				continue;

			std::vector<uint32_t> currAdjF;
			fV2Pool.emplace(f0);
			fFlag[f0] = true;
			while (!fV2Pool.empty())
			{
				uint32_t currF = fV2Pool.front();
				fV2Pool.pop();
				//fFlag[currF] = true;
				currAdjF.emplace_back(currF);
				fFindFlag[currF] = fV2SetVec.size();

				std::vector<uint32_t> &fvs = qm.Fs_[currF].vs;
				for (int j = 0; j < 4; ++j)
				{
					std::vector<uint32_t> &vfs = qm.Vs_[fvs[j]].neighbor_fs;
					if (vfs.size() == 2 && !qm.Vs_[fvs[j]].boundary)
					{
						uint32_t newF0 = vfs[0];
						if (!fFlag[newF0])
						{
							fV2Pool.emplace(newF0);
							fFlag[newF0] = true;
						}
						uint32_t newF1 = vfs[1];
						if (!fFlag[newF1])
						{
							fFlag[newF1] = true;
							fV2Pool.emplace(newF1);
						}
					}
				}
			}
			fV2SetVec.emplace_back(currAdjF);
		}

		std::vector<Eigen::Vector4i> fOldVs;
		Eigen::Vector4i currFvs;
		for (int i = 0; i < fFindFlag.size(); ++i)
		{
			currFvs = Eigen::Vector4i(-1, -1, -1, -1);
			if (fFindFlag[i] == -2)
			{
				std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
				fNewToNewNew[i] = fOldVs.size();
				fNewNewToNew.emplace_back(i);
				for (int j = 0; j < 4; ++j)
				{
					currFvs[j] = fvs[j];
				}
			}
			else if (fFindFlag[i] >= 0)
			{
				std::vector<uint32_t> &adjFs = fV2SetVec[fFindFlag[i]];
				for (int j = 0; j < adjFs.size(); ++j)
				{
					fNewToNewNew[adjFs[j]] = fOldVs.size();
				}
				fNewNewToNew.emplace_back(adjFs[0]);
				uint32_t sideFs[2];
				int currSF = 0;
				for (int j = 0; j < adjFs.size(); ++j)
				{
					fFindFlag[adjFs[j]] = -1;

					std::vector<uint32_t> &fvs = qm.Fs_[adjFs[j]].vs;
					uint32_t vMoreThan2 = 0;
					for (int n = 0; n < 4; ++n)
					{
						if (!(!qm.Vs_[fvs[n]].boundary &&qm.Vs_[fvs[n]].neighbor_fs.size() <= 2))
							++vMoreThan2;
					}
					if (vMoreThan2 >= 3)
					{
						sideFs[currSF++] = adjFs[j];
					}
				}

				uint32_t anotherV = 0;
				for (int j = 0; j < 4; ++j)
				{
					uint32_t currVv = qm.Fs_[sideFs[0]].vs[j];
					if (!(!qm.Vs_[currVv].boundary &&qm.Vs_[currVv].neighbor_fs.size() <= 2))
					{
						currFvs[j] = qm.Fs_[sideFs[0]].vs[j];
					}
					else
					{
						anotherV = j;
					}
				}
				for (int j = 0; j < 4; ++j)
				{
					uint32_t currV = qm.Fs_[sideFs[1]].vs[j];
					if (currV != currFvs[0] && currV != currFvs[1] && currV != currFvs[2] && currV != currFvs[3] && !(!qm.Vs_[currV].boundary &&qm.Vs_[currV].neighbor_fs.size() <= 2))
					{
						currFvs[anotherV] = currV;
						break;
					}
				}
			}
			else
				continue;
			fOldVs.emplace_back(currFvs);
		}


		//生成新的拓扑关系和网格。

		QuadMesh newnewQm;
		QuadVertex qv;
		for (int i = 0; i < vNewNewToNew.size(); ++i)
		{
			newnewQm.V_.emplace_back(qm.V_[vNewNewToNew[i]]);
			qv.id = i;
			qv.bId = (uint32_t)-1;
			newnewQm.Vs_.emplace_back(qv);
		}
		for (int i = 0; i < fOldVs.size(); ++i)
		{
			QuadFace qf;
			qf.id = i;
			qf.bId = (uint32_t)-1;
			for (int j = 0; j < 4; ++j)
			{
				qf.vs.emplace_back(vNewToNewNew[fOldVs[i][j]]);
				newnewQm.Vs_[vNewToNewNew[fOldVs[i][j]]].neighbor_fs.emplace_back(i);
			}
			newnewQm.Fs_.emplace_back(qf);
		}
		if (!newnewQm.BuildConnectivity())
			return false;
		if (!TopologyCheck(newnewQm))
			return false;

		qm = newnewQm;
		return true;
	}

	bool Pipeline::DeleteUselessElements(QuadMesh &qm, std::vector<uint32_t> &vNewToNewNew, std::vector<uint32_t> &vNewNewToNew)
	{
		std::vector<int> vFlag(qm.Vs_.size(), -2); //-2表示初始状态， -1表示这是一个内部度为2的点
		std::vector<uint32_t> valence2P;
		bool haveUselessEle = false;
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if (!qm.Vs_[i].boundary && qm.Vs_[i].neighbor_fs.size() == 2)
			{
				//也许要用到。
//				std::vector<uint32_t> vvs = qm.Vs_[i].neighbor_vs;
//				if (!qm.Vs_[vvs[0]].boundary && !qm.Vs_[vvs[1]].boundary && qm.Vs_[vvs[0]].neighbor_es.size() == 4 && qm.Vs_[vvs[1]].neighbor_es.size() == 4)
//					return false;
				vFlag[i] = -1;
				valence2P.emplace_back(i);
				haveUselessEle = true;
			}
			/*else if (fcc.V_types[vNewToOld[i]] != -1 && qm.Vs_[i].neighbor_fs.size() == 0)
			{
				vFlag[i] = 0;
				haveUselessEle = true;
			}*/
		}
		if (!haveUselessEle)
			return true;

		vNewToNewNew.resize(qm.Vs_.size(), (uint32_t)-1);
		uint32_t startCount = 0;
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if (vFlag[i] == -2)
			{
				vNewToNewNew[i] = startCount++;
				vNewNewToNew.emplace_back(i);
			}
		}

		std::vector<bool> fFlag(qm.Fs_.size(), false); //判断是否已经找过了。
		std::queue<uint32_t> fV2Pool;
		std::vector<std::vector<uint32_t>> fV2SetVec;
		std::vector<int> fFindFlag(qm.Fs_.size(), -2);//-2表示这个面不是一个包含度2点的面；>=0记录了一个包含度2顶点的面所在的fV2SetVec的下标，-1表示虽然这个面包含度2顶点，但是已经被找过了。
		for (int i = 0; i < valence2P.size(); ++i)
		{
			uint32_t f0 = qm.Vs_[valence2P[i]].neighbor_fs[0];
			if (fFlag[f0])
				continue;

			std::vector<uint32_t> currAdjF;
			fV2Pool.emplace(f0);
			fFlag[f0] = true;
			while (!fV2Pool.empty())
			{
				uint32_t currF = fV2Pool.front();
				fV2Pool.pop();
				//fFlag[currF] = true;
				currAdjF.emplace_back(currF);
				fFindFlag[currF] = fV2SetVec.size();

				std::vector<uint32_t> &fvs = qm.Fs_[currF].vs;
				for (int j = 0; j < 4; ++j)
				{
					std::vector<uint32_t> &vfs = qm.Vs_[fvs[j]].neighbor_fs;
					if (vfs.size() == 2 && !qm.Vs_[fvs[j]].boundary)
					{
						uint32_t newF0 = vfs[0];
						if (!fFlag[newF0])
						{
							fV2Pool.emplace(newF0);
							fFlag[newF0] = true;
						}
						uint32_t newF1 = vfs[1];
						if (!fFlag[newF1])
						{
							fFlag[newF1] = true;
							fV2Pool.emplace(newF1);
						}
					}
				}
			}
			fV2SetVec.emplace_back(currAdjF);
		}

		std::vector<Eigen::Vector4i> fOldVs;
		Eigen::Vector4i currFvs;
		for (int i = 0; i < fFindFlag.size(); ++i)
		{
			currFvs = Eigen::Vector4i(-1, -1, -1, -1);
			if (fFindFlag[i] == -2)
			{
				std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
				for (int j = 0; j < 4; ++j)
				{
					currFvs[j] = fvs[j];
				}
			}
			else if (fFindFlag[i] >= 0)
			{
				std::vector<uint32_t> &adjFs = fV2SetVec[fFindFlag[i]];
				uint32_t sideFs[2];
				int currSF = 0;
				for (int j = 0; j < adjFs.size(); ++j)
				{
					fFindFlag[adjFs[j]] = -1;

					std::vector<uint32_t> &fvs = qm.Fs_[adjFs[j]].vs;
					uint32_t vMoreThan2 = 0;
					for (int n = 0; n < 4; ++n)
					{
						if (!(!qm.Vs_[fvs[n]].boundary &&qm.Vs_[fvs[n]].neighbor_fs.size() <= 2))
							++vMoreThan2;
					}
					if (vMoreThan2 >= 3)
					{
						sideFs[currSF++] = adjFs[j];
					}
				}

				uint32_t anotherV = 0;
				for (int j = 0; j < 4; ++j)
				{
					uint32_t currVv = qm.Fs_[sideFs[0]].vs[j];
					if (!(!qm.Vs_[currVv].boundary &&qm.Vs_[currVv].neighbor_fs.size() <= 2))
					{
						currFvs[j] = qm.Fs_[sideFs[0]].vs[j];
					}
					else
					{
						anotherV = j;
					}
				}
				for (int j = 0; j < 4; ++j)
				{
					uint32_t currV = qm.Fs_[sideFs[1]].vs[j];
					if (currV != currFvs[0] && currV != currFvs[1] && currV != currFvs[2] && currV != currFvs[3] && !(!qm.Vs_[currV].boundary &&qm.Vs_[currV].neighbor_fs.size() <= 2))
					{
						currFvs[anotherV] = currV;
						break;
					}
				}
			}
			else
				continue;
			fOldVs.emplace_back(currFvs);
		}


		//生成新的拓扑关系和网格。

		QuadMesh newnewQm;
		QuadVertex qv;
		for (int i = 0; i < vNewNewToNew.size(); ++i)
		{
			newnewQm.V_.emplace_back(qm.V_[vNewNewToNew[i]]);
			qv.id = i;
			qv.bId = (uint32_t)-1;
			newnewQm.Vs_.emplace_back(qv);
		}
		for (int i = 0; i < fOldVs.size(); ++i)
		{
			QuadFace qf;
			qf.id = i;
			qf.bId = (uint32_t) - 1;
			for (int j = 0; j < 4; ++j)
			{
				qf.vs.emplace_back(vNewToNewNew[fOldVs[i][j]]);
				newnewQm.Vs_[vNewToNewNew[fOldVs[i][j]]].neighbor_fs.emplace_back(i);
			}
			newnewQm.Fs_.emplace_back(qf);
		}
		if (!newnewQm.BuildConnectivity())
			return false;
		if (!TopologyCheck(newnewQm))
			return false;

		/*std::vector<uint32_t> tempOldToNew(vOldToNew.size(), (uint32_t)-1), tempNewToOld;
		for (int i = 0; i < vOldToNew.size(); ++i)
		{
			uint32_t midV = vOldToNew[i];
			if (midV != (uint32_t)-1)
				tempOldToNew[i] = vNewToNewNew[midV];
		}
		for (int i = 0; i < vNewNewToNew.size(); ++i)
		{
			tempNewToOld.emplace_back(vNewToOld[vNewNewToNew[i]]);
		}*/
		qm = newnewQm;
		//vNewToOld = tempNewToOld;
		//vOldToNew = tempOldToNew;
		return true;
	}

	void Pipeline::BuildLocalOptimizeInfo(QuadMesh &qm, CollapseInfo *ci, OptimizeInfo *oi)
	{
		std::vector<bool> vFlag(qm.Vs_.size(), false);
		for (int i = 0; i < ci->hsSeveralring.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsSeveralring[i]].vs;
			for (int j = 0; j < fvs.size(); ++j)
			{
				vFlag[fvs[j]] = true;
			}
		}
		for (int i = 0; i < ci->vsGroup.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = ci->vsGroup[i];
			for (int j = 0; j < currVsGroup.size(); ++j)
			{
				if (fc.V_types[currVsGroup[j]] != -2)
					vFlag[currVsGroup[j]] = true;
			}
		}

		oi->globalToLocalVVec.resize(qm.Vs_.size(), (uint32_t)-1);
		uint32_t currNewVId = 0;
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if (vFlag[i])
			{
				oi->globalToLocalVVec[i] = currNewVId++;
				oi->localToGlobalVVec.emplace_back(i);
				oi->localVPos.emplace_back(qm.V_[i]);
			}
		}

		Eigen::Vector3i tempVi;
		for (int i = 0; i < ci->hsSeveralring.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsSeveralring[i]].vs;
			for (int j = 0; j < 4; ++j)
			{
				tempVi[0] = oi->globalToLocalVVec[fvs[quadTriTable[j][0]]];
				tempVi[1] = oi->globalToLocalVVec[fvs[quadTriTable[j][1]]];
				tempVi[2] = oi->globalToLocalVVec[fvs[quadTriTable[j][2]]];
				oi->localTriIds.emplace_back(tempVi);
			}
		}

		ComputeStandardTriangles(qm, ci, oi);

		oi->localVsGroups.resize(ci->vsGroup.size());
		oi->localVsCoords.resize(ci->targetCoords.size());
		if (ci->vsGroup.size() != ci->targetCoords.size())
			exit(-30);

		for (int i = 0; i < ci->vsGroup.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = ci->vsGroup[i];
			for (int j = 0; j < currVsGroup.size(); ++j)
			{
				uint32_t localVsGroup = oi->globalToLocalVVec[currVsGroup[j]];
				if (localVsGroup == (uint32_t)-1)
					continue;

				oi->localVsGroups[i].emplace_back(localVsGroup);
			}
			oi->localVsCoords[i] = ci->targetCoords[i];
		}

		for (int i = 0; i < ci->vsNotMove.size(); ++i)
		{
			if (oi->globalToLocalVVec[ci->vsNotMove[i]] == (uint32_t)-1)
				continue;

			oi->localVsNotMove.emplace_back(oi->globalToLocalVVec[ci->vsNotMove[i]]);
			oi->localVsNotMovePos.emplace_back(qm.V_[ci->vsNotMove[i]]);
		}

		//build new fc
		FeatureConstraints &newFc = oi->localFc;
		newFc.V_types.reserve(oi->localToGlobalVVec.size());
		for (int i = 0; i < oi->localToGlobalVVec.size(); ++i)
		{
			uint32_t refeVId = oi->localToGlobalVVec[i];
			newFc.V_types.emplace_back(fc.V_types[refeVId]);
		}

		for (int i = 0; i < fc.ids_C.size(); ++i)
		{
			uint32_t newV = oi->globalToLocalVVec[fc.ids_C[i]];
			if (newV == (uint32_t)-1)
				continue;

			newFc.ids_C.emplace_back(newV);
			newFc.C.emplace_back(fc.C[i]);
		}

		for (int i = 0; i < fc.ids_L.size(); ++i)
		{
			uint32_t newV = oi->globalToLocalVVec[fc.ids_L[i]];
			if (newV == (uint32_t)-1)
				continue;

			newFc.ids_L.emplace_back(newV);
			newFc.on_which_L.emplace_back(fc.on_which_L[i]);
			newFc.axa_L.emplace_back(fc.axa_L[i]);
			newFc.origin_L.emplace_back(fc.origin_L[i]);
		}
	}

	void Pipeline::BuildOptimizeInfoHasBoundary(QuadMesh &qm, OptimizeInfo *oi)
	{
		std::vector<std::vector<uint32_t>> &localVsGroups = oi->localVsGroups;
		//std::vector<bool>().swap(oi->vsGroupHasBoundary);
		oi->vsGroupHasBoundary.clear();
		oi->vsGroupHasBoundary.reserve(localVsGroups.size());
		for (int i = 0; i < localVsGroups.size(); ++i)
		{
			bool hasBoundary = false;
			for (int j = 0; j < localVsGroups[i].size(); ++j)
			{
				bool tempBS = qm.Vs_[oi->localToGlobalVVec[localVsGroups[i][j]]].boundary;
				if (tempBS)
				{
					hasBoundary = true;
					break;
				}
			}
			oi->vsGroupHasBoundary.emplace_back(hasBoundary);
		}
	}

	void Pipeline::ComputeStandardTriangles(QuadMesh &qm, CollapseInfo *ci, OptimizeInfo *oi)
	{
		std::function<double(Eigen::Vector2d&, Eigen::Vector2d&)> crossLength = [&](Eigen::Vector2d& vec1, Eigen::Vector2d& vec2)->double
		{
			return std::abs(vec1[0] * vec2[1] - vec1[1] * vec2[0]);
		};
		std::function<double(uint32_t)> computeFaceArea = [&](uint32_t faceId)->double
		{
			std::vector<uint32_t> &fvs = qm.Fs_[faceId].vs;
			Eigen::Vector2d e0 = qm.V_[fvs[1]] - qm.V_[fvs[0]];
			Eigen::Vector2d e1 = qm.V_[fvs[1]] - qm.V_[fvs[2]];
			Eigen::Vector2d e2 = qm.V_[fvs[3]] - qm.V_[fvs[2]];
			Eigen::Vector2d e3 = qm.V_[fvs[3]] - qm.V_[fvs[0]];

			double sumArea = 0;
			sumArea += crossLength(e0, e3) / 2;
			sumArea += crossLength(e1, e2) / 2;
			/*sumArea += e0.cross(e3).norm() / 2;
			sumArea += e1.cross(e2).norm() / 2;*/
			return sumArea;
		};

		std::vector<Eigen::Vector2d> standVVec(4);
		standVVec[0][0] = 0; standVVec[0][1] = 0;
		for (int i = 0; i < ci->hsSeveralring.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsSeveralring[i]].vs;
			double e0 = 0, e1 = 0;
#ifndef USE_SQUARE_AREA
			e0 += (qm.V_[fvs[1]] - qm.V_[fvs[0]]).norm();
			e0 += (qm.V_[fvs[3]] - qm.V_[fvs[2]]).norm();
			e0 /= 2;
			e1 += (qm.V_[fvs[2]] - qm.V_[fvs[1]]).norm();
			e1 += (qm.V_[fvs[3]] - qm.V_[fvs[0]]).norm();
			e1 /= 2;
#else
			e0 = 1;
			e1 = 1;
#endif

#ifndef USE_AVERAGE_AREA
			double referenceArea = computeFaceArea(ci->hsSeveralring[i]);
#else
			double referenceArea = initialArea_/qm.Fs_.size();
#endif
			double radio = std::sqrt(referenceArea / (e0*e1));
			e0 *= radio; e1 *= radio;

			standVVec[1][0] = e0; standVVec[1][1] = 0;
			standVVec[2][0] = e0; standVVec[2][1] = e1;
			standVVec[3][0] = 0; standVVec[3][1] = e1;

			std::vector<Eigen::Vector2d> tempVVec(3);
			for (int j = 0; j < 4; ++j)
			{
				tempVVec[0] = standVVec[quadTriTable[j][0]];
				tempVVec[1] = standVVec[quadTriTable[j][1]];
				tempVVec[2] = standVVec[quadTriTable[j][2]];
				oi->localStanTriPos.emplace_back(tempVVec);
			}
		}
	}

	void Pipeline::AdjustBoundaryPos(BaseDataStructure::QuadMesh &newQm, OptimizeInfo *newOi)
	{
		newQm.ComputeMeshNormal();
		uint32_t bfvs[4];
		for (int i = 0; i < newQm.Vs_.size(); ++i)
		{
			if (!newQm.Vs_[i].boundary || newQm.Vs_[i].neighbor_fs.size()!=1)
				continue;

			uint32_t currF = newQm.Vs_[i].neighbor_fs[0];
			uint32_t startV = (uint32_t)-1;
			const std::vector<uint32_t> &fvs = newQm.Fs_[currF].vs;
			for (int j = 0; j < 4; ++j)
			{
				if (fvs[j] == i)
				{
					startV = j;
					break;
				}
			}
			for (int j = 0; j < 4; ++j)
			{
				bfvs[j] = fvs[(startV + j) % 4];	//bfvs中保存的是以度二边界边为起始的，逆时针的顶点集合。
			}

			vec0_[0] = newQm.V_[bfvs[1]][0] - newQm.V_[bfvs[0]][0]; vec0_[1] = newQm.V_[bfvs[1]][1] - newQm.V_[bfvs[0]][1]; vec0_[2] = 0;
			vec1_[0] = newQm.V_[bfvs[3]][0] - newQm.V_[bfvs[0]][0]; vec1_[1] = newQm.V_[bfvs[3]][1] - newQm.V_[bfvs[0]][1]; vec1_[2] = 0;
			vec0_.normalize(); vec1_.normalize();
			double cosValue = vec0_[0] * vec1_[0] + vec0_[1] * vec1_[1];
			normalVec_ = vec0_.cross(vec1_);
			if (std::abs(normalVec_[2]) >= MIN_ALLOWED_JACOBI || cosValue>0)
				continue;

			double sinTheta = MIN_ALLOWED_JACOBI * 1.1;
			double theta = std::asin(sinTheta);
			if (std::cos(theta) >= 0)
				theta = PI - theta;

			double res0, res1;
			double x0 = newQm.V_[bfvs[3]][0], y0 = newQm.V_[bfvs[3]][1], x1 = newQm.V_[bfvs[1]][0], y1 = newQm.V_[bfvs[1]][1];
			AssistFunctions::AssistFunc::ComputeMinVPos(x0, y0, x1, y1, theta, res0, res1);
			if (res0 > 1E20)
				continue;

			Eigen::Vector2d testP((x0 + x1)*0.5 + res0 * (y0 - y1), (y0 + y1)*0.5 + res0 * (x1 - x0));
			vec0_[0] = newQm.V_[bfvs[1]][0] - testP[0]; vec0_[1] = newQm.V_[bfvs[1]][1] - testP[1]; vec0_[2] = 0;
			vec1_[0] = newQm.V_[bfvs[3]][0] - testP[0]; vec1_[1] = newQm.V_[bfvs[3]][1] - testP[1]; vec1_[2] = 0;
			double testV = vec0_[0] * vec1_[1] - vec0_[1] * vec1_[0];
			vec0_[0] = newQm.V_[bfvs[3]][0] - newQm.V_[bfvs[2]][0]; vec0_[1] = newQm.V_[bfvs[3]][1] - newQm.V_[bfvs[2]][1]; vec0_[2] = 0;
			vec1_[0] = newQm.V_[bfvs[1]][0] - newQm.V_[bfvs[2]][0]; vec1_[1] = newQm.V_[bfvs[1]][1] - newQm.V_[bfvs[2]][1]; vec1_[2] = 0;
			double testV2 = vec0_[0] * vec1_[1] - vec0_[1] * vec1_[0];
			if ((testV > 0 && testV2 < 0) || (testV < 0 && testV2>0))
			{
				testP[0] = (x0 + x1)*0.5 + res1 * (y0 - y1);
				testP[1] = (y0 + y1)*0.5 + res1 * (x1 - x0);
			}
			newQm.V_[i] = testP;

			/*uint32_t localId = newOi->globalToLocalVVec[i];
			if (localId == (uint32_t)-1)
				continue;
			if (newOi->localFc.V_types[localId] == -2)
				continue;
			else if (newOi->localFc.V_types[localId] == -1)
			{
				auto it = std::find(newOi->localFc.ids_C.begin(), newOi->localFc.ids_C.end(), localId);
				if (it == newOi->localFc.ids_C.end())
					continue;

				newQm.V_[i] = testP;
				newOi->localVPos[localId] = testP;
				int dis = std::distance(newOi->localFc.ids_C.begin(), it);
				newOi->localFc.C[dis] = testP;
			}
			else if (newOi->localFc.V_types[localId] >= 0)
			{
				auto it = std::find(newOi->localFc.ids_L.begin(), newOi->localFc.ids_L.end(), localId);
				if (it == newOi->localFc.ids_L.end())
					continue;

				newQm.V_[i] = testP;
				newOi->localVPos[localId] = testP;
				int dis = std::distance(newOi->localFc.ids_L.begin(), it);
				newOi->localFc.axa_L[dis] = Eigen::Vector2d(newQm.V_[bfvs[1]][0] - newQm.V_[bfvs[3]][0], newQm.V_[bfvs[1]][1] - newQm.V_[bfvs[3]][1]).normalized();
				newOi->localFc.origin_L[dis] = testP;
			}*/
		}
	}

	void Pipeline::ComputeV2BoundaryPos(QuadMesh &qm, uint32_t badV, double threshold)
	{
		qm.ComputeMeshNormal();
		uint32_t bfvs[4];
		if (!qm.Vs_[badV].boundary || qm.Vs_[badV].neighbor_fs.size() != 1)
			return;

		uint32_t currF = qm.Vs_[badV].neighbor_fs[0];
		uint32_t startV = (uint32_t)-1;
		std::vector<uint32_t> &fvs = qm.Fs_[currF].vs;
		for (int j = 0; j < 4; ++j)
		{
			if (fvs[j] == badV)
			{
				startV = j;
				break;
			}
		}
		for (int j = 0; j < 4; ++j)
		{
			bfvs[j] = fvs[(startV + j) % 4];	//bfvs中保存的是以度二边界边为起始的，逆时针的顶点集合。
		}

		vec0_[0] = qm.V_[bfvs[1]][0] - qm.V_[bfvs[0]][0]; vec0_[1] = qm.V_[bfvs[1]][1] - qm.V_[bfvs[0]][1]; vec0_[2] = 0;
		vec1_[0] = qm.V_[bfvs[3]][0] - qm.V_[bfvs[0]][0]; vec1_[1] = qm.V_[bfvs[3]][1] - qm.V_[bfvs[0]][1]; vec1_[2] = 0;
		vec0_.normalize(); vec1_.normalize();
		double cosValue = vec0_[0] * vec1_[0] + vec0_[1] * vec1_[1];
		normalVec_ = vec0_.cross(vec1_);
		if (std::abs(normalVec_[2]) >= threshold || cosValue > 0)
			return;

		double sinTheta = threshold * 1.1;
		double theta = std::asin(sinTheta);
		if (std::cos(theta) >= 0)
			theta = PI - theta;

		double res0, res1;
		double x0 = qm.V_[bfvs[3]][0], y0 = qm.V_[bfvs[3]][1], x1 = qm.V_[bfvs[1]][0], y1 = qm.V_[bfvs[1]][1];
		AssistFunctions::AssistFunc::ComputeMinVPos(x0, y0, x1, y1, theta, res0, res1);
		if (res0 > 1E20)
			return;

		Eigen::Vector2d testP((x0 + x1)*0.5 + res0 * (y0 - y1), (y0 + y1)*0.5 + res0 * (x1 - x0));
		vec0_[0] = qm.V_[bfvs[1]][0] - testP[0]; vec0_[1] = qm.V_[bfvs[1]][1] - testP[1]; vec0_[2] = 0;
		vec1_[0] = qm.V_[bfvs[3]][0] - testP[0]; vec1_[1] = qm.V_[bfvs[3]][1] - testP[1]; vec1_[2] = 0;
		double testV = vec0_[0] * vec1_[1] - vec0_[1] * vec1_[0];
		vec0_[0] = qm.V_[bfvs[3]][0] - qm.V_[bfvs[2]][0]; vec0_[1] = qm.V_[bfvs[3]][1] - qm.V_[bfvs[2]][1]; vec0_[2] = 0;
		vec1_[0] = qm.V_[bfvs[1]][0] - qm.V_[bfvs[2]][0]; vec1_[1] = qm.V_[bfvs[1]][1] - qm.V_[bfvs[2]][1]; vec1_[2] = 0;
		double testV2 = vec0_[0] * vec1_[1] - vec0_[1] * vec1_[0];
		if ((testV > 0 && testV2 < 0) || (testV < 0 && testV2>0))
		{
			testP[0] = (x0 + x1)*0.5 + res1 * (y0 - y1);
			testP[1] = (y0 + y1)*0.5 + res1 * (x1 - x0);
		}
		qm.V_[badV] = testP;
	}

	void Pipeline::AdjustBoundaryPosInOptimization(BaseDataStructure::QuadMesh &newQm, OptimizeInfo *newOi)
	{
		for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
		{
			uint32_t globalId = newOi->localToGlobalVVec[i];
			newQm.V_[globalId] = newOi->localVPos[i];
		}

		newQm.ComputeMeshNormal();
		uint32_t bfvs[4];
		for (int i = 0; i < newQm.Vs_.size(); ++i)
		{
			if (!newQm.Vs_[i].boundary || newQm.Vs_[i].neighbor_fs.size() != 1)
				continue;

			uint32_t currF = newQm.Vs_[i].neighbor_fs[0];
			uint32_t startV = (uint32_t)-1;
			const std::vector<uint32_t> &fvs = newQm.Fs_[currF].vs;
			for (int j = 0; j < 4; ++j)
			{
				if (fvs[j] == i)
				{
					startV = j;
					break;
				}
			}
			for (int j = 0; j < 4; ++j)
			{
				bfvs[j] = fvs[(startV + j) % 4];	//bfvs中保存的是以度二边界边为起始的，逆时针的顶点集合。
			}

			vec0_[0] = newQm.V_[bfvs[1]][0] - newQm.V_[bfvs[0]][0]; vec0_[1] = newQm.V_[bfvs[1]][1] - newQm.V_[bfvs[0]][1]; vec0_[2] = 0;
			vec1_[0] = newQm.V_[bfvs[3]][0] - newQm.V_[bfvs[0]][0]; vec1_[1] = newQm.V_[bfvs[3]][1] - newQm.V_[bfvs[0]][1]; vec1_[2] = 0;
			vec0_.normalize(); vec1_.normalize();
			double cosValue = vec0_[0] * vec1_[0] + vec0_[1] * vec1_[1];
			normalVec_ = vec0_.cross(vec1_);
			if (std::abs(normalVec_[2]) >= MIN_ALLOWED_JACOBI || cosValue>0)
				continue;

			double sinTheta = MIN_ALLOWED_JACOBI * 1.1;
			double theta = std::asin(sinTheta);
			if (std::cos(theta) >= 0)
				theta = PI - theta;

			double res0, res1;
			double x0 = newQm.V_[bfvs[3]][0], y0 = newQm.V_[bfvs[3]][1], x1 = newQm.V_[bfvs[1]][0], y1 = newQm.V_[bfvs[1]][1];
			AssistFunctions::AssistFunc::ComputeMinVPos(x0, y0, x1, y1, theta, res0, res1);
			if (res0 > 1E20)
				continue;

			Eigen::Vector2d testP((x0 + x1)*0.5 + res0 * (y0 - y1), (y0 + y1)*0.5 + res0 * (x1 - x0));
			vec0_[0] = newQm.V_[bfvs[1]][0] - testP[0]; vec0_[1] = newQm.V_[bfvs[1]][1] - testP[1]; vec0_[2] = 0;
			vec1_[0] = newQm.V_[bfvs[3]][0] - testP[0]; vec1_[1] = newQm.V_[bfvs[3]][1] - testP[1]; vec1_[2] = 0;
			double testV = vec0_[0] * vec1_[1] - vec0_[1] * vec1_[0];
			vec0_[0] = newQm.V_[bfvs[3]][0] - newQm.V_[bfvs[2]][0]; vec0_[1] = newQm.V_[bfvs[3]][1] - newQm.V_[bfvs[2]][1]; vec0_[2] = 0;
			vec1_[0] = newQm.V_[bfvs[1]][0] - newQm.V_[bfvs[2]][0]; vec1_[1] = newQm.V_[bfvs[1]][1] - newQm.V_[bfvs[2]][1]; vec1_[2] = 0;
			double testV2 = vec0_[0] * vec1_[1] - vec0_[1] * vec1_[0];
			if ((testV > 0 && testV2 < 0) || (testV < 0 && testV2>0))
			{
				testP[0] = (x0 + x1)*0.5 + res1 * (y0 - y1);
				testP[1] = (y0 + y1)*0.5 + res1 * (x1 - x0);
			}
			newQm.V_[i] = testP;

			uint32_t localId = newOi->globalToLocalVVec[i];
			if (localId == (uint32_t)-1)
				continue;
			if (newOi->localFc.V_types[localId] == -2)
				continue;
			else
			{
				newOi->localVPos[localId] = testP;
			}
		}
	}

	void Pipeline::ChangeTheCorners(QuadMesh &qm, QuadMesh &newQm, std::vector<uint32_t> &vOldToNew, OptimizeInfo *oi, CollapseInfo *ci, FeatureConstraints &fc)
	{
		bool *ciVFlag = new bool[qm.Vs_.size()];	//标记所有ci中的边界点
		std::memset(ciVFlag, 0, qm.Vs_.size() * sizeof(bool));
		//bool *removeVFlag = new bool[qm.Vs_.size()];	//标记hsToBeKilled中的点
		//std::memset(removeVFlag, 0, qm.Vs_.size() * sizeof(bool));
		std::vector<uint32_t> bCiVs;	//保存所有ci中的边界点
		bCiVs.reserve(qm.Vs_.size());
		bool *findVFlag = new bool[qm.Vs_.size()];	//标记一个边界点有没有被找过。
		std::memset(findVFlag, 0, qm.Vs_.size() * sizeof(bool));
		bool *deleteVFlag = new bool[qm.Vs_.size()];
		std::memset(deleteVFlag, 0, qm.Vs_.size() * sizeof(bool));

		for (int i = 0; i < ci->hsToBeKilled.size(); ++i)
		{
			if (!qm.Fs_[ci->hsToBeKilled[i]].boundary)
				continue;

			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsToBeKilled[i]].vs;
			for (int j = 0; j < 4; ++j)
			{
				if (qm.Vs_[fvs[j]].boundary)
				{
					ciVFlag[fvs[j]] = true;
					bCiVs.emplace_back(fvs[j]);
					deleteVFlag[fvs[j]] = true;
				}
			}
		}
		for (int i = 0; i < ci->hsSeveralring.size(); ++i)
		{
			if (!qm.Fs_[ci->hsSeveralring[i]].boundary)
				continue;

			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsSeveralring[i]].vs;
			for (int j = 0; j < 4; ++j)
			{
				if (qm.Vs_[fvs[j]].boundary)
				{
					ciVFlag[fvs[j]] = true;
					bCiVs.emplace_back(fvs[j]);
				}
			}
		}

		bool closeLoop = false;
		std::function<bool(uint32_t, uint32_t, uint32_t&, uint32_t&)> findNextVE = [&](uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE) -> bool
		{
			

			std::vector<uint32_t> &ves = qm.Vs_[preV].neighbor_es;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (qm.Es_[ves[i]].boundary && ves[i] != preE)
				{
					nextE = ves[i];
					break;
				}
			}

			if (preV == qm.Es_[nextE].vs[0])
				nextV = qm.Es_[nextE].vs[1];
			else
				nextV = qm.Es_[nextE].vs[0];

			if (findVFlag[nextV])
			{
				closeLoop = true;
				return false;
			}
			if (!ciVFlag[nextV])
				return false;
			
			return true;
		};

		std::vector<std::vector<uint32_t>> vVecVec, eVecVec;
		std::vector<bool> isCloseVec;
		for (int i = 0; i < bCiVs.size(); ++i)
		{
			if (findVFlag[bCiVs[i]] || deleteVFlag[bCiVs[i]])
				continue;

			std::vector<uint32_t> vVec, eVec;
			closeLoop = false;
			vVec.emplace_back(bCiVs[i]);	//如果是loop的话，第一个点没有重复
			findVFlag[bCiVs[i]] = true;

			uint32_t startV = bCiVs[i];
			uint32_t preE = (uint32_t)-1, preV = (uint32_t)-1;
			uint32_t preE2[2];
			std::vector<uint32_t> &ves = qm.Vs_[startV].neighbor_es;
			int currCount = 0;
			for (int j = 0; j < ves.size(); ++j)
			{
				if (qm.Es_[ves[j]].boundary)
				{
					preE2[currCount++] = ves[j];
				}
			}

			preE = preE2[0];
			std::vector<uint32_t> &evs = qm.Es_[preE].vs;
			if (evs[0] == startV)
				preV = evs[1];
			else
				preV = evs[0];

			bool isBegin = false;
			if (!ciVFlag[preV])
			{
				isBegin = true;
				preE = preE2[1];
				std::vector<uint32_t> &newEvs = qm.Es_[preE].vs;
				if (newEvs[0] == startV)
					preV = newEvs[1];
				else
					preV = newEvs[0];
			}
			if (!ciVFlag[preV])
				continue;

			vVec.emplace_back(preV);
			findVFlag[preV] = true;
			eVec.emplace_back(preE);
			uint32_t nextV = (uint32_t)-1, nextE = (uint32_t)-1;
			while (findNextVE(preV, preE, nextV, nextE))
			{
				vVec.emplace_back(nextV);
				findVFlag[nextV] = true;
				eVec.emplace_back(nextE);

				preV = nextV;
				preE = nextE;
			}
			if (closeLoop)
				eVec.emplace_back(nextE);

			if (!closeLoop && !isBegin)
			{
				std::reverse(vVec.begin(), vVec.end());
				std::reverse(eVec.begin(), eVec.end());
				preE = preE2[1];
				std::vector<uint32_t> &newEvs = qm.Es_[preE].vs;
				if (newEvs[0] == startV)
					preV = newEvs[1];
				else
					preV = newEvs[0];
				if (!ciVFlag[preV])
				{
					isCloseVec.emplace_back(closeLoop);
					vVecVec.emplace_back(vVec);
					eVecVec.emplace_back(eVec);
					continue;
				}

				vVec.emplace_back(preV);
				findVFlag[preV] = true;
				eVec.emplace_back(preE);
				//uint32_t nextV = (uint32_t)-1, nextE = (uint32_t)-1;
				while (findNextVE(preV, preE, nextV, nextE))
				{
					vVec.emplace_back(nextV);
					findVFlag[nextV] = true;
					eVec.emplace_back(nextE);

					preV = nextV;
					preE = nextE;
				}
			}
			isCloseVec.emplace_back(closeLoop);
			vVecVec.emplace_back(vVec);
			eVecVec.emplace_back(eVec);
		}

		//找到了许多条curve之后，可以在每条curve上找新的corner点。
		for (int i = 0; i < vVecVec.size(); ++i)
		{
			std::vector<uint32_t> &currVVec = vVecVec[i];
			std::vector<uint32_t> &currEVec = eVecVec[i];
			bool isLoop = isCloseVec[i];
			
			if (!isLoop)
			{
				std::vector<uint32_t> cornerPos;
				std::vector<uint32_t> onWhichC;
				std::vector<uint32_t> toBeDeleteV, dVOnWhichC;	//toBeDeleteV保存了要删掉的那些V，数量和要被删掉的边数相当。dVOnWhichC保存了要被删掉的点所在的边。
				bool findTDSV = false;
				uint32_t toBeDeleteStartV = (uint32_t)-1;
				for (int j = 0; j < currVVec.size(); ++j)
				{
					onWhichC.emplace_back(fc.V_types[currVVec[j]]);
					if (fc.V_types[currVVec[j]] == -1)
						cornerPos.emplace_back(j);

					uint32_t newVId = vOldToNew[currVVec[j]];
					if (newVId == (uint32_t)-1)
					{
						if (findTDSV == false)
						{
							toBeDeleteStartV = j;
							findTDSV = true;
						}
						toBeDeleteV.emplace_back(currVVec[j]);
						dVOnWhichC.emplace_back(fc.V_types[currVVec[j]]);
					}
				}

				std::sort(dVOnWhichC.begin(), dVOnWhichC.end());
				dVOnWhichC.erase(std::unique(dVOnWhichC.begin(), dVOnWhichC.end()), dVOnWhichC.end());
				if (cornerPos.empty() || dVOnWhichC.size() > 2 || dVOnWhichC.empty() || toBeDeleteStartV==0)
					continue;
				else if (dVOnWhichC.size() == 2)
				{
					if (dVOnWhichC[0] != -1 && dVOnWhichC[1] != -1)
						continue;
				}

				//uint32_t dVC = dVOnWhichC[0];	//被删掉的那些点所在的curve
				//if (dVOnWhichC.size() == 2)
				//{
				//	if (dVC == -1)
				//		dVC = dVOnWhichC[1];
				//}

				std::vector<Eigen::Vector2d> ori_L, axa_L;
				std::vector<Eigen::Vector2d> cornersP;
				ori_L.reserve(cornerPos.size() + 1);
				axa_L.reserve(cornerPos.size() + 1);
				cornersP.reserve(cornerPos.size() + 1);
				std::vector<uint32_t> realCurveNums;	//按顺序保存所有碰到的边界边实际上的段数。
				std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> onWhichCInfo; //找到所有位于不同边上的点。三个元分别代表左起点(currVVec的下标），在优化中的边的段数，所在的curve编号
				std::tuple<uint32_t, uint32_t, uint32_t> currEle;
				if (cornerPos[0] != 0)
				{
					std::get<0>(currEle) = 0;
					std::get<1>(currEle) = cornerPos[0];
					std::get<2>(currEle) = onWhichC[0];
					realCurveNums.emplace_back(FindCurveEdgeNum(qm, fc, currVVec[0]));
					onWhichCInfo.emplace_back(currEle);

					ori_L.emplace_back((qm.V_[currVVec[1]] + qm.V_[currVVec[0]])*0.5);
					axa_L.emplace_back((qm.V_[currVVec[1]] - qm.V_[currVVec[0]]).normalized());

					//cornersP.emplace_back(qm.V_[currVVec[cornerPos[0]]]);
				}
				bool isFindRightCorner = true;
				for (int j = 0; j < cornerPos.size()-1; ++j)
				{
					ori_L.emplace_back((qm.V_[currVVec[cornerPos[j]+1]] + qm.V_[currVVec[cornerPos[j]]])*0.5);
					axa_L.emplace_back((qm.V_[currVVec[cornerPos[j]+1]] - qm.V_[currVVec[cornerPos[j]]]).normalized());
					cornersP.emplace_back(qm.V_[currVVec[cornerPos[j]]]);
					std::get<0>(currEle) = cornerPos[j];
					std::get<1>(currEle) = cornerPos[j + 1] - cornerPos[j];
					realCurveNums.emplace_back(cornerPos[j + 1] - cornerPos[j]);
					std::get<2>(currEle) = (uint32_t)-1;
					if (cornerPos[j + 1] - cornerPos[j] != 1)
					{
						std::get<2>(currEle) = onWhichC[cornerPos[j] + 1];
					}
					else
					{
						Eigen::Vector2d &cPos0 = qm.V_[currVVec[cornerPos[j]]], &cPos1 = qm.V_[currVVec[cornerPos[j + 1]]];
						uint32_t cor0 = (uint32_t)-1, cor1 = (uint32_t)-1;
						FindSourceCorner(cPos0, cor0);
						FindSourceCorner(cPos1, cor1);
						std::tuple<int, int> &corCur0 = mf.cornerCurves[cor0], &corCur1 = mf.cornerCurves[cor1];
						if (std::get<0>(corCur0) == std::get<0>(corCur1) || std::get<0>(corCur0) == std::get<1>(corCur1))
							std::get<2>(currEle) = std::get<0>(corCur0);
						else if (std::get<1>(corCur0) == std::get<0>(corCur1) || std::get<1>(corCur0) == std::get<1>(corCur1))
							std::get<2>(currEle) = std::get<1>(corCur0);
					}
					if (std::get<2>(currEle) == (uint32_t)-1)
					{
						isFindRightCorner = false;
						break;
					}
					onWhichCInfo.emplace_back(currEle);
				}
				if (!isFindRightCorner)
					continue;

				if (cornerPos[cornerPos.size() - 1] != currVVec.size() - 1)
				{
					ori_L.emplace_back((qm.V_[currVVec[currVVec.size() - 1]] + qm.V_[currVVec[currVVec.size() - 2]])*0.5);
					axa_L.emplace_back((qm.V_[currVVec[currVVec.size() - 1]] - qm.V_[currVVec[currVVec.size() - 2]]).normalized());
					cornersP.emplace_back(qm.V_[currVVec[cornerPos[cornerPos.size() - 1]]]);
					std::get<0>(currEle) = cornerPos[cornerPos.size() - 1];
					std::get<1>(currEle) = currVVec.size() - 1 - cornerPos[cornerPos.size() - 1];
					std::get<2>(currEle) = onWhichC[currVVec.size()-1];
					realCurveNums.emplace_back(FindCurveEdgeNum(qm, fc, currVVec[currVVec.size() - 1]));
					onWhichCInfo.emplace_back(currEle);
				}

				if (std::get<2>(onWhichCInfo[0]) == std::get<2>(onWhichCInfo[onWhichCInfo.size() - 1]))
				{
					isLoop = true;
					--i;
					continue;
				}

				//计算每个curve应该保留多少
				int finalSumEdges = 0;
				std::vector<double> curvesWeight;
				curvesWeight.reserve(realCurveNums.size());
				for (int j = 0; j < realCurveNums.size(); ++j)
				{
					finalSumEdges += realCurveNums[j];
					curvesWeight.emplace_back(mf.curveLength[std::get<2>(onWhichCInfo[j])]);
				}
				finalSumEdges -= toBeDeleteV.size();
				if (finalSumEdges <= 0)
					continue;
				
				double sumWeight = std::accumulate(curvesWeight.begin(), curvesWeight.end(), 0);
				std::vector<uint32_t> targetENum;
				targetENum.reserve(curvesWeight.size());

				int remainEN = finalSumEdges;
				for (int j = 0; j < curvesWeight.size()-1; ++j)
				{
					int currENum = std::floor(finalSumEdges * curvesWeight[j] / sumWeight);
					if (currENum == 0)
						currENum = 1;
					remainEN = remainEN - currENum;
					targetENum.emplace_back(currENum);
				}
				if (remainEN <= 0)
					continue;
				targetENum.emplace_back(remainEN);

				//找新的corners和curve vertices.
				std::vector<uint32_t> newCorners; //保存的是currVVec的下标，表示currVVec的第几个元素是新顶点。
				int newId = 0;
				if (cornerPos[0] != 0)
				{
					int remainNotMove = realCurveNums[0] - std::get<1>(onWhichCInfo[0]);
					newId = targetENum[0] - remainNotMove;
					if (newId <= 0)
						newId = 1;
					newCorners.emplace_back(newId);
				}
				int offset = cornerPos[0] != 0 ? 0 : 1;
				for (int j = 1; j < cornerPos.size() - 1; ++j)
				{
					newId += targetENum[j - offset];
					newCorners.emplace_back(newId);
				}
				if (cornerPos[cornerPos.size() - 1] != currVVec.size() - 1)
				{
					newId += targetENum[cornerPos.size() - 1 - offset];
					newCorners.emplace_back(newId);
				}
				for (int j = 0; j < newCorners.size(); ++j)
				{
					if (newCorners[j] >= toBeDeleteStartV)
						newCorners[j] += toBeDeleteV.size();
				}
				if (newId > currVVec.size()-1)
					continue;

				//重建oi->localFc
				bool *vFFlag = new bool[oi->localVPos.size()];
				std::memset(vFFlag, 0, oi->localVPos.size() * sizeof(bool));
				FeatureConstraints newLocalFc;
				newLocalFc.V_ids.resize(oi->localFc.V_ids.size());
				newLocalFc.V_types.resize(oi->localFc.V_types.size());
				uint32_t currCurveCount = 0, currCorner = 0;
				for (int j = 1; j < currVVec.size() - 1; ++j)
				{
					uint32_t localVId = oi->globalToLocalVVec[currVVec[j]];
					if (localVId == (uint32_t)-1)
						continue;

					vFFlag[localVId] = true;
					if (j == newCorners[currCorner])
					{
						newLocalFc.ids_C.emplace_back(localVId);
						newLocalFc.C.emplace_back(cornersP[currCorner]);
						newLocalFc.V_types[localVId] = -1;
						++currCorner;
						++currCurveCount;
					}
					else
					{
						newLocalFc.ids_L.emplace_back(localVId);
						newLocalFc.on_which_L.emplace_back(std::get<2>(onWhichCInfo[currCurveCount]));
						newLocalFc.axa_L.emplace_back(axa_L[currCurveCount]);
						newLocalFc.origin_L.emplace_back(ori_L[currCurveCount]);
						newLocalFc.V_types[localVId] = std::get<2>(onWhichCInfo[currCurveCount]);
					}
				}
				for (int j = 0; j < oi->localFc.V_types.size(); ++j)
				{
					if (vFFlag[j])
						continue;
					
					newLocalFc.V_types[j] = oi->localFc.V_types[j];
				}
				for (int j = 0; j < oi->localFc.ids_C.size(); ++j)
				{
					if (vFFlag[j])
						continue;

					newLocalFc.C.emplace_back(oi->localFc.C[j]);
					newLocalFc.ids_C.emplace_back(oi->localFc.ids_C[j]);
				}
				for (int j = 0; j < oi->localFc.ids_L.size(); ++j)
				{
					if (vFFlag[j])
						continue;
					newLocalFc.ids_L.emplace_back(oi->localFc.ids_L[j]);
					newLocalFc.origin_L.emplace_back(oi->localFc.origin_L[j]);
					newLocalFc.axa_L.emplace_back(oi->localFc.axa_L[j]);
				}

				oi->localFc = newLocalFc;
				delete[] vFFlag;
			}
			else
			{
				std::vector<uint32_t> cornerPos;
				std::vector<uint32_t> onWhichC;
				std::vector<uint32_t> toBeDeleteV, dVOnWhichC;	//toBeDeleteV保存了要删掉的那些V，数量和要被删掉的边数相当。dVOnWhichC保存了要被删掉的点所在的边。
				bool findTDSV = false;
				uint32_t toBeDeleteStartV = (uint32_t)-1;
				cornerPos.emplace_back(0);
				for (int j = 0; j < currVVec.size(); ++j)
				{
					onWhichC.emplace_back(fc.V_types[currVVec[j]]);
					if (fc.V_types[currVVec[j]] == -1)
						cornerPos.emplace_back(j);

					uint32_t newVId = vOldToNew[currVVec[j]];
					if (newVId == (uint32_t)-1)
					{
						if (findTDSV == false)
						{
							toBeDeleteStartV = j;
							findTDSV = true;
						}
						toBeDeleteV.emplace_back(currVVec[j]);
						dVOnWhichC.emplace_back(fc.V_types[currVVec[j]]);
					}
				}
				cornerPos.emplace_back(currVVec.size() - 1);

				std::sort(dVOnWhichC.begin(), dVOnWhichC.end());
				dVOnWhichC.erase(std::unique(dVOnWhichC.begin(), dVOnWhichC.end()), dVOnWhichC.end());
				if (cornerPos.empty() || dVOnWhichC.size() > 2 || dVOnWhichC.empty() || toBeDeleteStartV==0)
					continue;
				else if (dVOnWhichC.size() == 2)
				{
					if (dVOnWhichC[0] != -1 && dVOnWhichC[1] != -1)
						continue;
				}

				uint32_t dVC = dVOnWhichC[0];	//被删掉的那些点所在的curve
				if (dVOnWhichC.size() == 2)
				{
					if (dVC == -1)
						dVC = dVOnWhichC[1];
				}

				std::vector<Eigen::Vector2d> ori_L, axa_L;
				std::vector<Eigen::Vector2d> cornersP;
				ori_L.reserve(cornerPos.size() + 1);
				axa_L.reserve(cornerPos.size() + 1);
				cornersP.reserve(cornerPos.size() + 1);
				std::vector<uint32_t> realCurveNums;	//按顺序保存所有碰到的边界边实际上的段数。
				std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> onWhichCInfo; //找到所有位于不同边上的点。三个元分别代表左起点(currVVec的下标），在优化中的边的段数，所在的curve编号
				std::tuple<uint32_t, uint32_t, uint32_t> currEle;
				/*if (cornerPos[0] != 0)
				{
					std::get<0>(currEle) = 0;
					std::get<1>(currEle) = cornerPos[0];
					std::get<2>(currEle) = onWhichC[0];
					realCurveNums.emplace_back(FindCurveEdgeNum(qm, fc, currVVec[0]));
					onWhichCInfo.emplace_back(currEle);

					ori_L.emplace_back((qm.V_[currVVec[1]] + qm.V_[currVVec[0]])*0.5);
					axa_L.emplace_back((qm.V_[currVVec[1]] - qm.V_[currVVec[0]]).normalized());
				}*/
				bool isFindRightCorner = true;
				for (int j = 0; j < cornerPos.size() - 1; ++j)
				{
					ori_L.emplace_back((qm.V_[currVVec[cornerPos[j] + 1]] + qm.V_[currVVec[cornerPos[j]]])*0.5);
					axa_L.emplace_back((qm.V_[currVVec[cornerPos[j] + 1]] - qm.V_[currVVec[cornerPos[j]]]).normalized());
					cornersP.emplace_back(qm.V_[currVVec[cornerPos[j]]]);
					std::get<0>(currEle) = cornerPos[j];
					std::get<1>(currEle) = cornerPos[j + 1] - cornerPos[j];
					realCurveNums.emplace_back(cornerPos[j + 1] - cornerPos[j]);
					std::get<2>(currEle) = (uint32_t)-1;
					if (cornerPos[j + 1] - cornerPos[j] != 1)
					{
						std::get<2>(currEle) = onWhichC[cornerPos[j] + 1];
					}
					else
					{
						Eigen::Vector2d &cPos0 = qm.V_[currVVec[cornerPos[j]]], &cPos1 = qm.V_[currVVec[cornerPos[j + 1]]];
						uint32_t cor0 = (uint32_t)-1, cor1 = (uint32_t)-1;
						FindSourceCorner(cPos0, cor0);
						FindSourceCorner(cPos1, cor1);
						std::tuple<int, int> &corCur0 = mf.cornerCurves[cor0], &corCur1 = mf.cornerCurves[cor1];
						if (std::get<0>(corCur0) == std::get<0>(corCur1) || std::get<0>(corCur0) == std::get<1>(corCur1))
							std::get<2>(currEle) = std::get<0>(corCur0);
						else if (std::get<1>(corCur0) == std::get<0>(corCur1) || std::get<1>(corCur0) == std::get<1>(corCur1))
							std::get<2>(currEle) = std::get<1>(corCur0);
					}
					if (std::get<2>(currEle) == (uint32_t)-1)
					{
						isFindRightCorner = false;
						break;
					}
					onWhichCInfo.emplace_back(currEle);
				}
				if (!isFindRightCorner)
					continue;

				/*if (cornerPos[cornerPos.size() - 1] != currVVec.size() - 1)
				{
					ori_L.emplace_back((qm.V_[currVVec[currVVec.size() - 1]] + qm.V_[currVVec[currVVec.size() - 2]])*0.5);
					axa_L.emplace_back((qm.V_[currVVec[currVVec.size() - 1]] - qm.V_[currVVec[currVVec.size() - 2]]).normalized());
					std::get<0>(currEle) = cornerPos[cornerPos.size() - 1];
					std::get<1>(currEle) = currVVec.size() - 1 - cornerPos[cornerPos.size() - 1];
					std::get<2>(currEle) = onWhichC[currVVec.size() - 1];
					realCurveNums.emplace_back(FindCurveEdgeNum(qm, fc, currVVec[currVVec.size() - 1]));
					onWhichCInfo.emplace_back(currEle);
				}*/

				//计算每个curve应该保留多少
				int finalSumEdges = 0;
				std::vector<double> curvesWeight;
				curvesWeight.reserve(realCurveNums.size());
				//这里的计算和上面不一样，收尾两端不能直接采用mf.curveLength。
				double sumLength = 0;
				for (int j = 0; j < cornerPos[1]; ++j)
				{
					sumLength += (qm.V_[currVVec[j + 1]] - qm.V_[currVVec[j]]).norm();
				}
				finalSumEdges += realCurveNums[0];
				curvesWeight.emplace_back(sumLength);
				for (int j = 1; j < realCurveNums.size()-1; ++j)
				{
					finalSumEdges += realCurveNums[j];
					curvesWeight.emplace_back(mf.curveLength[std::get<2>(onWhichCInfo[j])]);
				}
				sumLength = 0;
				for (int j = cornerPos[cornerPos.size()-2]; j < cornerPos[cornerPos.size() - 1]; ++j)
				{
					sumLength += (qm.V_[currVVec[j + 1]] - qm.V_[currVVec[j]]).norm();
				}
				finalSumEdges += realCurveNums[realCurveNums.size()-1];
				curvesWeight.emplace_back(sumLength);
				finalSumEdges -= toBeDeleteV.size();
				if (finalSumEdges <= 0)
					continue;

				double sumWeight = std::accumulate(curvesWeight.begin(), curvesWeight.end(), 0);
				std::vector<uint32_t> targetENum;
				targetENum.reserve(curvesWeight.size());

				int remainEN = finalSumEdges;
				for (int j = 0; j < curvesWeight.size() - 1; ++j)
				{
					int currENum = std::floor(finalSumEdges * curvesWeight[j] / sumWeight);
					if (currENum == 0)
						currENum = 1;
					remainEN = remainEN - currENum;
					targetENum.emplace_back(currENum);
				}
				if (remainEN <= 0)
					continue;
				targetENum.emplace_back(remainEN);

				//找新的corners和curve vertices.
				std::vector<uint32_t> newCorners; //保存的是currVVec的下标，表示currVVec的第几个元素是新顶点。
				int newId = 0;
				/*if (cornerPos[0] != 0)
				{
					int remainNotMove = realCurveNums[0] - std::get<1>(onWhichCInfo[0]);
					newId = targetENum[0] - remainNotMove;
					if (newId <= 0)
						newId = 1;
					newCorners.emplace_back(newId);
				}*/
				int offset = cornerPos[0] != 0 ? 0 : 1;
				for (int j = 1; j < cornerPos.size() - 1; ++j)
				{
					newId += targetENum[j - offset];
					newCorners.emplace_back(newId);
				}
				/*if (cornerPos[cornerPos.size() - 1] != currVVec.size() - 1)
				{
					newId += targetENum[cornerPos.size() - 1 - offset];
					newCorners.emplace_back(newId);
				}*/
				for (int j = 0; j < newCorners.size(); ++j)
				{
					if (newCorners[j] >= toBeDeleteStartV)
						newCorners[j] += toBeDeleteV.size();
				}
				if (newId > currVVec.size() - 1)
					continue;

				//重建oi->localFc
				bool *vFFlag = new bool[oi->localVPos.size()];
				std::memset(vFFlag, 0, oi->localVPos.size() * sizeof(bool));
				FeatureConstraints newLocalFc;
				newLocalFc.V_ids.resize(oi->localFc.V_ids.size());
				newLocalFc.V_types.resize(oi->localFc.V_types.size());
				uint32_t currCurveCount = 0, currCorner = 0;
				for (int j = 1; j < currVVec.size() - 1; ++j)
				{
					uint32_t localVId = oi->globalToLocalVVec[currVVec[j]];
					if (localVId == (uint32_t)-1)
						continue;

					vFFlag[localVId] = true;
					if (j == newCorners[currCorner])
					{
						newLocalFc.ids_C.emplace_back(localVId);
						newLocalFc.C.emplace_back(cornersP[currCorner]);
						newLocalFc.V_types[localVId] = -1;
						++currCorner;
						++currCurveCount;
					}
					else
					{
						newLocalFc.ids_L.emplace_back(localVId);
						newLocalFc.on_which_L.emplace_back(std::get<2>(onWhichCInfo[currCurveCount]));
						newLocalFc.axa_L.emplace_back(axa_L[currCurveCount]);
						newLocalFc.origin_L.emplace_back(ori_L[currCurveCount]);
						newLocalFc.V_types[localVId] = std::get<2>(onWhichCInfo[currCurveCount]);
					}
				}
				for (int j = 0; j < oi->localFc.V_types.size(); ++j)
				{
					if (vFFlag[j])
						continue;

					newLocalFc.V_types[j] = oi->localFc.V_types[j];
				}
				for (int j = 0; j < oi->localFc.ids_C.size(); ++j)
				{
					if (vFFlag[j])
						continue;

					newLocalFc.C.emplace_back(oi->localFc.C[j]);
					newLocalFc.ids_C.emplace_back(oi->localFc.ids_C[j]);
				}
				for (int j = 0; j < oi->localFc.ids_L.size(); ++j)
				{
					if (vFFlag[j])
						continue;
					newLocalFc.ids_L.emplace_back(oi->localFc.ids_L[j]);
					newLocalFc.origin_L.emplace_back(oi->localFc.origin_L[j]);
					newLocalFc.axa_L.emplace_back(oi->localFc.axa_L[j]);
				}

				oi->localFc = newLocalFc;
				delete[] vFFlag;
			}
		}

		delete[] ciVFlag;
		ciVFlag = NULL;
		//delete[] removeVFlag;
		delete[] findVFlag;
		findVFlag = NULL;
		delete[] deleteVFlag;
		deleteVFlag = NULL;
	}

	uint32_t Pipeline::FindCurveEdgeNum(QuadMesh &qm, FeatureConstraints &fc, uint32_t startV)
	{
		bool *vFlag = new bool[qm.Vs_.size()];
		std::memset(vFlag, 0, qm.Vs_.size() * sizeof(bool));

		bool isLoop = false;
		std::function<bool(uint32_t, uint32_t, uint32_t&, uint32_t&)> findNextVE = [&](uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE) -> bool
		{
			std::vector<uint32_t> &evs = qm.Es_[preE].vs;
			if (evs[0] == preV)
				nextV = evs[1];
			else
				nextV = evs[0];

			std::vector<uint32_t> &ves = qm.Vs_[nextV].neighbor_es;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (qm.Es_[ves[i]].boundary && ves[i] != preE)
				{
					nextE = ves[i];
					break;
				}
			}

			if (vFlag[nextV])
			{
				isLoop = true;
				return false;
			}
			if (fc.V_types[nextV] == -1)
				return false;

			return true;
		};

		uint32_t preE2[2];
		uint32_t preV = startV, preE = (uint32_t)-1;
		vFlag[preV] = true;
		int currCount = 0;
		std::vector<uint32_t> &ves = qm.Vs_[preV].neighbor_es;
		for (int i = 0; i < ves.size(); ++i)
		{
			if (qm.Es_[ves[i]].boundary)
				preE2[currCount++] = ves[i];
		}

		uint32_t finalENum = 1;
		preE = preE2[0];
		uint32_t nextV, nextE;
		while (findNextVE(preV, preE, nextV, nextE))
		{
			++finalENum;
			preV = nextV;
			preE = nextE;
		}

		if (!isLoop)
		{
			preV = startV;
			preE = preE2[1];

			while (findNextVE(preV, preE, nextV, nextE))
			{
				++finalENum;
				preV = nextV;
				preE = nextE;
			}
		}

		delete[] vFlag;
		return finalENum;
	}

	uint32_t Pipeline::FindCurveEdgeNum(BaseDataStructure::QuadMesh &newQm, OptimizeInfo *newOi, uint32_t startV)
	{
		bool *vFlag = new bool[newQm.Vs_.size()];
		std::memset(vFlag, 0, newQm.Vs_.size() * sizeof(bool));
		//FeatureConstraints fcc = newOi->localFc;

		bool isLoop = false;
		std::function<bool(uint32_t, uint32_t, uint32_t&, uint32_t&)> findNextVE = [&](uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE) -> bool
		{
			std::vector<uint32_t> &evs = newQm.Es_[preE].vs;
			if (evs[0] == preV)
				nextV = evs[1];
			else
				nextV = evs[0];

			std::vector<uint32_t> &ves = newQm.Vs_[nextV].neighbor_es;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (newQm.Es_[ves[i]].boundary && ves[i] != preE)
				{
					nextE = ves[i];
					break;
				}
			}

			if (vFlag[nextV])
			{
				isLoop = true;
				return false;
			}
			if (fc.V_types[vNewToOld_[nextV]] == -1)
				return false;

			return true;
		};

		uint32_t preE2[2];
		uint32_t preV = startV, preE = (uint32_t)-1;
		vFlag[preV] = true;
		int currCount = 0;
		std::vector<uint32_t> &ves = newQm.Vs_[preV].neighbor_es;
		for (int i = 0; i < ves.size(); ++i)
		{
			if (newQm.Es_[ves[i]].boundary)
				preE2[currCount++] = ves[i];
		}

		uint32_t finalENum = 1;
		preE = preE2[0];
		uint32_t nextV, nextE;
		while (findNextVE(preV, preE, nextV, nextE))
		{
			++finalENum;
			vFlag[nextV] = true;
			preV = nextV;
			preE = nextE;
		}

		if (!isLoop)
		{
			preV = startV;
			preE = preE2[1];
			++finalENum;

			while (findNextVE(preV, preE, nextV, nextE))
			{
				++finalENum;
				vFlag[nextV] = true;
				preV = nextV;
				preE = nextE;
			}
		}

		delete[] vFlag;
		return finalENum;
	}
	void Pipeline::ChangeTheCorners(QuadMesh &newQm, OptimizeInfo *newOi)
	{
		bool *oiVFlag = new bool[newQm.Vs_.size()];	//判断一个点是否是待优化的点
		std::memset(oiVFlag, 0, newQm.Vs_.size() * sizeof(bool));
		bool *findVFlag = new bool[newQm.Vs_.size()];	//判断一个点是否被找过
		std::memset(findVFlag, 0, newQm.Vs_.size() * sizeof(bool));

		std::vector<uint32_t> oiVs;
		oiVs.reserve(newOi->localToGlobalVVec.size());
		uint32_t currVId;
		for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
		{
			currVId = newOi->localToGlobalVVec[i];
			if (newQm.Vs_[currVId].boundary)
			{
				oiVs.emplace_back(currVId);
				oiVFlag[currVId] = true;
			}
		}
		std::sort(oiVs.begin(), oiVs.end());
		oiVs.erase(std::unique(oiVs.begin(), oiVs.end()), oiVs.end());

		bool isLoop = false;
		std::function<bool(uint32_t, uint32_t, uint32_t &, uint32_t &)> findNextVE = [&](uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE) -> bool
		{
			std::vector<uint32_t> &evs = newQm.Es_[preE].vs;
			if (evs[0] == preV)
				nextV = evs[1];
			else
				nextV = evs[0];

			std::vector<uint32_t> &ves = newQm.Vs_[nextV].neighbor_es;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (newQm.Es_[ves[i]].boundary && ves[i] != preE)
				{
					nextE = ves[i];
					break;
				}
			}

			if (findVFlag[nextV])
			{
				isLoop = true;
				return false;
			}
			if (!oiVFlag[nextV])
				return false;

			return true;
		};

		std::vector<std::vector<uint32_t> > vVecVec;
		std::vector<bool> isLoopVec;
		for (int i = 0; i < oiVs.size(); ++i)
		{
			currVId = oiVs[i];
			if (findVFlag[oiVs[i]])
				continue;

			isLoop = false;
			std::vector<uint32_t> vVec;
			uint32_t startV = currVId, preV = startV, preE = (uint32_t)-1, nextV = (uint32_t)-1, nextE = (uint32_t)-1;
			vVec.emplace_back(preV);
			findVFlag[preV] = true;

			uint32_t currCount = 0;
			uint32_t preE2[2];
			std::vector<uint32_t> &ves = newQm.Vs_[startV].neighbor_es;
			for (int j = 0; j < ves.size(); ++j)
			{
				if (newQm.Es_[ves[j]].boundary)
					preE2[currCount++] = ves[j];
			}

			preE = preE2[0];
			while (findNextVE(preV, preE, nextV, nextE))
			{
				vVec.emplace_back(nextV);
				findVFlag[nextV] = true;
				preV = nextV;
				preE = nextE;
			}

			if (!isLoop)
			{
				std::reverse(vVec.begin(), vVec.end());
				preV = startV;
				preE = preE2[1];

				while (findNextVE(preV, preE, nextV, nextE))
				{
					vVec.emplace_back(nextV);
					findVFlag[nextV] = true;
					preV = nextV;
					preE = nextE;
				}
			}

			vVecVec.emplace_back(vVec);
			isLoopVec.emplace_back(isLoop);
		}

		FeatureConstraints &fcc = newOi->localFc;
		FeatureConstraints newLocalFc;
		newLocalFc.V_ids.resize(newOi->localFc.V_ids.size());
		newLocalFc.V_types.resize(newOi->localFc.V_types.size());
		bool *vFFlag = new bool[newOi->localVPos.size()];
		std::memset(vFFlag, 0, newOi->localVPos.size() * sizeof(bool));
		bool *curvesFound = new bool[mf.curveLength.size()];
		std::memset(curvesFound, 0, mf.curveLength.size() * sizeof(bool));
		for (int i = 0; i < vVecVec.size(); ++i)
		{
			std::vector<uint32_t> &currVVec = vVecVec[i];
			bool iLoop = isLoopVec[i];
			if (currVVec.size() < 2)
				continue;

			if (false)
			{
				std::vector<uint32_t> cornerPos;
				std::vector<int> onWhichC;
				for (int j = 0; j < currVVec.size(); ++j)
				{
					uint32_t localV = newOi->globalToLocalVVec[currVVec[j]];
					onWhichC.emplace_back(fcc.V_types[localV]);
					if (fcc.V_types[localV] == -1)
						cornerPos.emplace_back(j);
				}

				std::vector<int> tempOWC = onWhichC;
				std::sort(tempOWC.begin(), tempOWC.end());
				tempOWC.erase(std::unique(tempOWC.begin(), tempOWC.end()), tempOWC.end());
				if (tempOWC.size() == 1 || (tempOWC.size() == 2 && tempOWC[0] == -1))
					continue;
				bool haveFound = false;
				for (int j = 0; j < tempOWC.size(); ++j)
				{
					if (tempOWC[j] >= 0 && curvesFound[tempOWC[j]])
					{
						haveFound = true;
						break;
					}
				}
				if (haveFound)
					continue;

				if (cornerPos.empty())
					continue;

				std::vector<Eigen::Vector2d> ori_L, axa_L;
				std::vector<Eigen::Vector2d> cornersP;
				ori_L.reserve(cornerPos.size() + 1);
				axa_L.reserve(cornerPos.size() + 1);
				cornersP.reserve(cornerPos.size() + 1);
				std::vector<uint32_t> realCurveNums;	//按顺序保存所有碰到的边界边实际上的段数。
				std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> onWhichCInfo; //找到所有位于不同边上的点。三个元分别代表左起点(currVVec的下标），在优化中的边的段数，所在的curve编号
				std::tuple<uint32_t, uint32_t, uint32_t> currEle;
				if (cornerPos[0] != 0)
				{
					std::get<0>(currEle) = 0;
					std::get<1>(currEle) = cornerPos[0];
					std::get<2>(currEle) = onWhichC[0];
					realCurveNums.emplace_back(FindCurveEdgeNum(newQm, newOi, currVVec[0]));
					onWhichCInfo.emplace_back(currEle);

					ori_L.emplace_back((newQm.V_[currVVec[1]] + newQm.V_[currVVec[0]])*0.5);
					axa_L.emplace_back((newQm.V_[currVVec[1]] - newQm.V_[currVVec[0]]).normalized());

					//cornersP.emplace_back(qm.V_[currVVec[cornerPos[0]]]);
				}
				bool isFindRightCorner = true;
				for (int j = 0; j < cornerPos.size() - 1; ++j)
				{
					ori_L.emplace_back((newQm.V_[currVVec[cornerPos[j] + 1]] + newQm.V_[currVVec[cornerPos[j]]])*0.5);
					axa_L.emplace_back((newQm.V_[currVVec[cornerPos[j] + 1]] - newQm.V_[currVVec[cornerPos[j]]]).normalized());
					if (cornerPos[j] != 0)
						cornersP.emplace_back(newQm.V_[currVVec[cornerPos[j]]]);
					std::get<0>(currEle) = cornerPos[j];
					std::get<1>(currEle) = cornerPos[j + 1] - cornerPos[j];
					realCurveNums.emplace_back(cornerPos[j + 1] - cornerPos[j]);
					std::get<2>(currEle) = (uint32_t)-1;
					if (cornerPos[j + 1] - cornerPos[j] != 1)
					{
						std::get<2>(currEle) = onWhichC[cornerPos[j] + 1];
					}
					else
					{
						Eigen::Vector2d &cPos0 = newQm.V_[currVVec[cornerPos[j]]], &cPos1 = newQm.V_[currVVec[cornerPos[j + 1]]];
						uint32_t cor0 = (uint32_t)-1, cor1 = (uint32_t)-1;
						FindSourceCorner(cPos0, cor0);
						FindSourceCorner(cPos1, cor1);
						std::tuple<int, int> &corCur0 = mf.cornerCurves[cor0], &corCur1 = mf.cornerCurves[cor1];
						if (std::get<0>(corCur0) == std::get<0>(corCur1) || std::get<0>(corCur0) == std::get<1>(corCur1))
							std::get<2>(currEle) = std::get<0>(corCur0);
						else if (std::get<1>(corCur0) == std::get<0>(corCur1) || std::get<1>(corCur0) == std::get<1>(corCur1))
							std::get<2>(currEle) = std::get<1>(corCur0);
					}
					if (std::get<2>(currEle) == (uint32_t)-1)
					{
						isFindRightCorner = false;
						break;
					}
					onWhichCInfo.emplace_back(currEle);
				}
				if (!isFindRightCorner)
					continue;

				if (cornerPos[cornerPos.size() - 1] != currVVec.size() - 1)
				{
					ori_L.emplace_back((newQm.V_[currVVec[currVVec.size() - 1]] + newQm.V_[currVVec[currVVec.size() - 2]])*0.5);
					axa_L.emplace_back((newQm.V_[currVVec[currVVec.size() - 1]] - newQm.V_[currVVec[currVVec.size() - 2]]).normalized());
					cornersP.emplace_back(newQm.V_[currVVec[cornerPos[cornerPos.size() - 1]]]);
					std::get<0>(currEle) = cornerPos[cornerPos.size() - 1];
					std::get<1>(currEle) = currVVec.size() - 1 - cornerPos[cornerPos.size() - 1];
					std::get<2>(currEle) = onWhichC[currVVec.size() - 1];
					realCurveNums.emplace_back(FindCurveEdgeNum(newQm, newOi, currVVec[currVVec.size() - 1]));
					onWhichCInfo.emplace_back(currEle);
				}

				if (std::get<2>(onWhichCInfo[0]) == std::get<2>(onWhichCInfo[onWhichCInfo.size() - 1]))
				{
					isLoopVec[i] = true;
					--i;
					continue;
				}

				//计算每个curve应该保留多少
				int finalSumEdges = 0;
				std::vector<double> curvesWeight;
				curvesWeight.reserve(realCurveNums.size());
				for (int j = 0; j < realCurveNums.size(); ++j)
				{
					finalSumEdges += realCurveNums[j];
					curvesWeight.emplace_back(mf.curveLength[std::get<2>(onWhichCInfo[j])]);
				}
				if (finalSumEdges == 0)
					continue;

				double sumWeight = std::accumulate(curvesWeight.begin(), curvesWeight.end(), 0.0);
				std::vector<uint32_t> targetENum;
				targetENum.reserve(curvesWeight.size());

				int remainEN = finalSumEdges;
				for (int j = 0; j < curvesWeight.size() - 1; ++j)
				{
					int currENum = std::floor(finalSumEdges * curvesWeight[j] / sumWeight);
					if (currENum == 0)
						currENum = 1;
					remainEN = remainEN - currENum;
					targetENum.emplace_back(currENum);
				}
				if (remainEN <= 0)
					continue;
				targetENum.emplace_back(remainEN);

				//找新的corners和curve vertices.
				std::vector<uint32_t> newCorners; //保存的是currVVec的下标，表示currVVec的第几个元素是新顶点。
				int newId = 0;
				if (cornerPos[0] != 0)
				{
					int remainNotMove = realCurveNums[0] - std::get<1>(onWhichCInfo[0]);
					newId = targetENum[0] - remainNotMove;
					if (newId <= 0)
						newId = 1;
					newCorners.emplace_back(newId);
				}
				int offset = cornerPos[0] != 0 ? 0 : 1;
				for (int j = 1; j < cornerPos.size() - 1; ++j)
				{
					newId += targetENum[j - offset];
					newCorners.emplace_back(newId);
				}
				if (cornersP.size()!=1 && cornerPos[cornerPos.size() - 1] != currVVec.size() - 1)
				{
					newId += targetENum[cornerPos.size() - 1 - offset];
					newCorners.emplace_back(newId);
				}

				if (newId >= currVVec.size() - 1)
				{

				}
				/*if (newId > currVVec.size() - 1)
					continue;*/

				//重建oi->localFc
				uint32_t currCurveCount = 0, currCorner = 0;
				for (int j = 1; j < currVVec.size() - 1; ++j)
				{
					uint32_t localVId = newOi->globalToLocalVVec[currVVec[j]];
					if (localVId == (uint32_t)-1 || vFFlag[localVId])
						continue;

					vFFlag[localVId] = true;
					if (currCorner<newCorners.size() && j == newCorners[currCorner])
					{
						newLocalFc.ids_C.emplace_back(localVId);
						newLocalFc.C.emplace_back(cornersP[currCorner]);
						newLocalFc.V_types[localVId] = -1;
						++currCorner;
						++currCurveCount;
					}
					else
					{
						newLocalFc.ids_L.emplace_back(localVId);
						newLocalFc.on_which_L.emplace_back(std::get<2>(onWhichCInfo[currCurveCount]));
						newLocalFc.axa_L.emplace_back(axa_L[currCurveCount]);
						newLocalFc.origin_L.emplace_back(ori_L[currCurveCount]);
						newLocalFc.V_types[localVId] = std::get<2>(onWhichCInfo[currCurveCount]);
					}
				}
				
				for (int j = 0; j < tempOWC.size(); ++j)
				{
					if (tempOWC[j] >= 0)
						curvesFound[tempOWC[j]] = true;
				}
			}
			else
			{
				std::vector<uint32_t> cornerPos;
				std::vector<int> onWhichC;
				cornerPos.emplace_back(0);
				onWhichC.emplace_back(fcc.V_types[newOi->globalToLocalVVec[currVVec[0]]]);
				for (int j = 1; j < currVVec.size()-1; ++j)
				{
					uint32_t localV = newOi->globalToLocalVVec[currVVec[j]];
					onWhichC.emplace_back(fcc.V_types[localV]);
					if (fcc.V_types[localV] == -1)
						cornerPos.emplace_back(j);
				}
				cornerPos.emplace_back(currVVec.size() - 1);
				onWhichC.emplace_back(fcc.V_types[newOi->globalToLocalVVec[currVVec[currVVec.size()-1]]]);

				bool haveMidM1 = false;
				for (int j = 1; j < onWhichC.size() - 1; ++j)
				{
					if (onWhichC[j] == -1)
					{
						haveMidM1 = true;
						break;
					}
				}
				if (!haveMidM1)
					continue;

				std::vector<int> tempOWC = onWhichC;
				std::sort(tempOWC.begin(), tempOWC.end());
				tempOWC.erase(std::unique(tempOWC.begin(), tempOWC.end()), tempOWC.end());
				/*if (tempOWC.size() == 1 || (tempOWC.size() == 2 && tempOWC[0] == -1))
					continue;*/
				bool haveFound = false;
				for (int j = 0; j < tempOWC.size(); ++j)
				{
					if (tempOWC[j] >= 0 && curvesFound[tempOWC[j]])
					{
						haveFound = true;
						break;
					}
				}
				if (haveFound)
					continue;

				if (cornerPos.empty())
					continue;

				std::vector<Eigen::Vector2d> ori_L, axa_L;
				std::vector<Eigen::Vector2d> cornersP;
				ori_L.reserve(cornerPos.size() + 1);
				axa_L.reserve(cornerPos.size() + 1);
				cornersP.reserve(cornerPos.size() + 1);
				std::vector<uint32_t> realCurveNums;	//按顺序保存所有碰到的边界边实际上的段数。
				std::vector<std::tuple<uint32_t, uint32_t, int>> onWhichCInfo; //找到所有位于不同边上的点。三个元分别代表左起点(currVVec的下标），在优化中的边的段数，所在的curve编号
				std::tuple<uint32_t, uint32_t, int> currEle;
				//if (cornerPos[0] != 0)
				//{
				//	std::get<0>(currEle) = 0;
				//	std::get<1>(currEle) = cornerPos[0];
				//	std::get<2>(currEle) = onWhichC[0];
				//	realCurveNums.emplace_back(FindCurveEdgeNum(newQm, newOi, currVVec[0]));
				//	onWhichCInfo.emplace_back(currEle);

				//	ori_L.emplace_back((newQm.V_[currVVec[1]] + newQm.V_[currVVec[0]])*0.5);
				//	axa_L.emplace_back((newQm.V_[currVVec[1]] - newQm.V_[currVVec[0]]).normalized());

				//	//cornersP.emplace_back(qm.V_[currVVec[cornerPos[0]]]);
				//}
				bool isFindRightCorner = true;
				for (int j = 0; j < cornerPos.size() - 1; ++j)
				{
					ori_L.emplace_back((newQm.V_[currVVec[cornerPos[j] + 1]] + newQm.V_[currVVec[cornerPos[j]]])*0.5);
					axa_L.emplace_back((newQm.V_[currVVec[cornerPos[j] + 1]] - newQm.V_[currVVec[cornerPos[j]]]).normalized());
					if (cornerPos[j]!=0)
						cornersP.emplace_back(newQm.V_[currVVec[cornerPos[j]]]);
					std::get<0>(currEle) = cornerPos[j];
					std::get<1>(currEle) = cornerPos[j + 1] - cornerPos[j];
					realCurveNums.emplace_back(cornerPos[j + 1] - cornerPos[j]);
					std::get<2>(currEle) = -100;
					if (cornerPos[j + 1] - cornerPos[j] != 1)
					{
						std::get<2>(currEle) = onWhichC[cornerPos[j] + 1];
					}
					/*else if (j == 0 && onWhichC[0] != -1 && onWhichC[1] == -1)
					{
						std::get<2>(currEle) = onWhichC[0];
					}
					else if (j == cornerPos.size() - 2 && onWhichC[currVVec.size() - 2] == -1 && onWhichC[currVVec.size() - 1] != -1)
					{
						std::get<2>(currEle) = onWhichC[currVVec.size() - 1];
					}*/
					else
					{
						Eigen::Vector2d &cPos0 = newQm.V_[currVVec[cornerPos[j]]], &cPos1 = newQm.V_[currVVec[cornerPos[j + 1]]];
						uint32_t cor0 = (uint32_t)-1, cor1 = (uint32_t)-1;
						FindSourceCorner(cPos0, cor0);
						FindSourceCorner(cPos1, cor1);
						std::tuple<int, int> &corCur0 = mf.cornerCurves[cor0], &corCur1 = mf.cornerCurves[cor1];
						if (std::get<0>(corCur0) == std::get<0>(corCur1) || std::get<0>(corCur0) == std::get<1>(corCur1))
							std::get<2>(currEle) = std::get<0>(corCur0);
						else if (std::get<1>(corCur0) == std::get<0>(corCur1) || std::get<1>(corCur0) == std::get<1>(corCur1))
							std::get<2>(currEle) = std::get<1>(corCur0);
					}
					if (std::get<2>(currEle) == -1 || std::get<2>(currEle) == -100)
					{
						isFindRightCorner = false;
						break;
					}
					onWhichCInfo.emplace_back(currEle);
				}
				if (!isFindRightCorner)
					continue;

				/*if (cornerPos[cornerPos.size() - 1] != currVVec.size() - 1)
				{
					ori_L.emplace_back((newQm.V_[currVVec[currVVec.size() - 1]] + newQm.V_[currVVec[currVVec.size() - 2]])*0.5);
					axa_L.emplace_back((newQm.V_[currVVec[currVVec.size() - 1]] - newQm.V_[currVVec[currVVec.size() - 2]]).normalized());
					cornersP.emplace_back(newQm.V_[currVVec[cornerPos[cornerPos.size() - 1]]]);
					std::get<0>(currEle) = cornerPos[cornerPos.size() - 1];
					std::get<1>(currEle) = currVVec.size() - 1 - cornerPos[cornerPos.size() - 1];
					std::get<2>(currEle) = onWhichC[currVVec.size() - 1];
					realCurveNums.emplace_back(FindCurveEdgeNum(newQm, newOi, currVVec[currVVec.size() - 1]));
					onWhichCInfo.emplace_back(currEle);
				}*/

				//计算每个curve应该保留多少
				int finalSumEdges = 0;
				std::vector<double> curvesWeight;
				curvesWeight.reserve(realCurveNums.size());
				double sumLength = 0;
				for (int j = 0; j < cornerPos[1]; ++j)
				{
					sumLength += (newQm.V_[currVVec[j + 1]] - newQm.V_[currVVec[j]]).norm();
				}
				finalSumEdges += realCurveNums[0];
				curvesWeight.emplace_back(sumLength);
				for (int j = 1; j < realCurveNums.size() - 1; ++j)
				{
					finalSumEdges += realCurveNums[j];
					curvesWeight.emplace_back(mf.curveLength[std::get<2>(onWhichCInfo[j])]);
				}
				sumLength = 0;
				for (int j = cornerPos[cornerPos.size() - 2]; j < cornerPos[cornerPos.size() - 1]; ++j)
				{
					sumLength += (newQm.V_[currVVec[j + 1]] - newQm.V_[currVVec[j]]).norm();
				}
				finalSumEdges += realCurveNums[realCurveNums.size() - 1];
				curvesWeight.emplace_back(sumLength);
				if (finalSumEdges == 0)
					continue;

				double sumWeight = std::accumulate(curvesWeight.begin(), curvesWeight.end(), 0.0);
				std::vector<uint32_t> targetENum;
				targetENum.reserve(curvesWeight.size());

				int remainEN = finalSumEdges;
				for (int j = 0; j < curvesWeight.size() - 1; ++j)
				{
					int currENum = std::floor(finalSumEdges * curvesWeight[j] / sumWeight);
					if (currENum == 0)
						currENum = 1;
					remainEN = remainEN - currENum;
					targetENum.emplace_back(currENum);
				}
				if (remainEN <= 0)
					continue;
				targetENum.emplace_back(remainEN);

				//找新的corners和curve vertices.
				std::vector<uint32_t> newCorners; //保存的是currVVec的下标，表示currVVec的第几个元素是新顶点。
				int newId = 0;
				/*if (cornerPos[0] != 0)
				{
					int remainNotMove = realCurveNums[0] - std::get<1>(onWhichCInfo[0]);
					newId = targetENum[0] - remainNotMove;
					if (newId <= 0)
						newId = 1;
					newCorners.emplace_back(newId);
				}*/
				int offset = 1;
				for (int j = 1; j < cornerPos.size() - 1; ++j)
				{
					newId += targetENum[j - offset];
					newCorners.emplace_back(newId);
				}
				if (newId >= currVVec.size() - 1)
					exit(-3442);
				/*if (cornerPos[cornerPos.size() - 1] != currVVec.size() - 1)
				{
					newId += targetENum[cornerPos.size() - 1 - offset];
					newCorners.emplace_back(newId);
				}*/
				/*if (newId > currVVec.size() - 1)
					continue;*/

				//重建oi->localFc
				uint32_t currCurveCount = 0, currCorner = 0;
				for (int j = 1; j < currVVec.size() - 1; ++j)
				{
					uint32_t localVId = newOi->globalToLocalVVec[currVVec[j]];
					if (localVId == (uint32_t)-1 || vFFlag[localVId])
						continue;

					vFFlag[localVId] = true;
					if (currCorner < newCorners.size() && j == newCorners[currCorner])
					{
						newLocalFc.ids_C.emplace_back(localVId);
						newLocalFc.C.emplace_back(cornersP[currCorner]);
						newLocalFc.V_types[localVId] = -1;
						++currCorner;
						++currCurveCount;
					}
					else
					{
						newLocalFc.ids_L.emplace_back(localVId);
						newLocalFc.on_which_L.emplace_back(std::get<2>(onWhichCInfo[currCurveCount]));
						newLocalFc.axa_L.emplace_back(axa_L[currCurveCount]);
						newLocalFc.origin_L.emplace_back(ori_L[currCurveCount]);
						newLocalFc.V_types[localVId] = std::get<2>(onWhichCInfo[currCurveCount]);
					}
				}
				for (int j = 0; j < tempOWC.size(); ++j)
				{
					if (tempOWC[j] >= 0)
						curvesFound[tempOWC[j]] = true;
				}
			}
		}

		for (int j = 0; j < newOi->localFc.V_types.size(); ++j)
		{
			if (vFFlag[j])
				continue;

			newLocalFc.V_types[j] = newOi->localFc.V_types[j];
		}
		for (int j = 0; j < newOi->localFc.ids_C.size(); ++j)
		{
			if (vFFlag[newOi->localFc.ids_C[j]])
				continue;

			newLocalFc.C.emplace_back(newOi->localFc.C[j]);
			newLocalFc.ids_C.emplace_back(newOi->localFc.ids_C[j]);
		}
		for (int j = 0; j < newOi->localFc.ids_L.size(); ++j)
		{
			if (vFFlag[newOi->localFc.ids_L[j]])
				continue;
			newLocalFc.ids_L.emplace_back(newOi->localFc.ids_L[j]);
			newLocalFc.on_which_L.emplace_back(newOi->localFc.on_which_L[j]);
			newLocalFc.origin_L.emplace_back(newOi->localFc.origin_L[j]);
			newLocalFc.axa_L.emplace_back(newOi->localFc.axa_L[j]);
		}

		newOi->localFc = newLocalFc;
		delete[] vFFlag;
		vFFlag = NULL;
		delete[] oiVFlag;
		oiVFlag = NULL;
		delete[] findVFlag;
		findVFlag = NULL;
	}

	void Pipeline::ChangeTheCornersInTwoCurves(QuadMesh &newQm, OptimizeInfo *newOi)
	{
		bool *oiVFlag = new bool[newQm.Vs_.size()];	//判断一个点是否是待优化的点
		std::memset(oiVFlag, 0, newQm.Vs_.size() * sizeof(bool));
		bool *findVFlag = new bool[newQm.Vs_.size()];	//判断一个点是否被找过
		std::memset(findVFlag, 0, newQm.Vs_.size() * sizeof(bool));

		std::vector<uint32_t> oiVs;
		oiVs.reserve(newOi->localToGlobalVVec.size());
		uint32_t currVId;
		for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
		{
			currVId = newOi->localToGlobalVVec[i];
			if (newQm.Vs_[currVId].boundary)
			{
				oiVs.emplace_back(currVId);
				oiVFlag[currVId] = true;
			}
		}
		std::sort(oiVs.begin(), oiVs.end());
		oiVs.erase(std::unique(oiVs.begin(), oiVs.end()), oiVs.end());

		bool isLoop = false;
		std::function<bool(uint32_t, uint32_t, uint32_t &, uint32_t &)> findNextVE = [&](uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE) -> bool
		{
			std::vector<uint32_t> &evs = newQm.Es_[preE].vs;
			if (evs[0] == preV)
				nextV = evs[1];
			else
				nextV = evs[0];

			std::vector<uint32_t> &ves = newQm.Vs_[nextV].neighbor_es;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (newQm.Es_[ves[i]].boundary && ves[i] != preE)
				{
					nextE = ves[i];
					break;
				}
			}

			if (findVFlag[nextV])
			{
				isLoop = true;
				return false;
			}
			if (!oiVFlag[nextV])
				return false;

			return true;
		};

		std::vector<std::vector<uint32_t> > vVecVec;
		std::vector<bool> isLoopVec;
		for (int i = 0; i < oiVs.size(); ++i)
		{
			currVId = oiVs[i];
			if (findVFlag[oiVs[i]])
				continue;

			isLoop = false;
			std::vector<uint32_t> vVec;
			uint32_t startV = currVId, preV = startV, preE = (uint32_t)-1, nextV = (uint32_t)-1, nextE = (uint32_t)-1;
			vVec.emplace_back(preV);
			findVFlag[preV] = true;

			uint32_t currCount = 0;
			uint32_t preE2[2];
			std::vector<uint32_t> &ves = newQm.Vs_[startV].neighbor_es;
			for (int j = 0; j < ves.size(); ++j)
			{
				if (newQm.Es_[ves[j]].boundary)
					preE2[currCount++] = ves[j];
			}

			preE = preE2[0];
			while (findNextVE(preV, preE, nextV, nextE))
			{
				vVec.emplace_back(nextV);
				findVFlag[nextV] = true;
				preV = nextV;
				preE = nextE;
			}

			if (!isLoop)
			{
				std::reverse(vVec.begin(), vVec.end());
				preV = startV;
				preE = preE2[1];

				while (findNextVE(preV, preE, nextV, nextE))
				{
					vVec.emplace_back(nextV);
					findVFlag[nextV] = true;
					preV = nextV;
					preE = nextE;
				}
			}

			vVecVec.emplace_back(vVec);
			isLoopVec.emplace_back(isLoop);
		}

		FeatureConstraints &fcc = newOi->localFc;
		FeatureConstraints newLocalFc;
		newLocalFc.V_ids.resize(newOi->localFc.V_ids.size());
		newLocalFc.V_types.resize(newOi->localFc.V_types.size());
		bool *vFFlag = new bool[newOi->localVPos.size()];
		std::memset(vFFlag, 0, newOi->localVPos.size() * sizeof(bool));
		bool *curvesFound = new bool[mf.curveLength.size()];
		std::memset(curvesFound, 0, mf.curveLength.size() * sizeof(bool));
		for (int i = 0; i < vVecVec.size(); ++i)
		{
			std::vector<uint32_t> &currVVec = vVecVec[i];
			bool iLoop = isLoopVec[i];
			if (currVVec.size() < 2)
				continue;

			std::vector<uint32_t> cornerPos;
			std::vector<int> onWhichC;
			cornerPos.emplace_back(0);
			onWhichC.emplace_back(fcc.V_types[newOi->globalToLocalVVec[currVVec[0]]]);
			for (int j = 1; j < currVVec.size() - 1; ++j)
			{
				uint32_t localV = newOi->globalToLocalVVec[currVVec[j]];
				onWhichC.emplace_back(fcc.V_types[localV]);
				if (fcc.V_types[localV] == -1)
					cornerPos.emplace_back(j);
			}
			cornerPos.emplace_back(currVVec.size() - 1);
			onWhichC.emplace_back(fcc.V_types[newOi->globalToLocalVVec[currVVec[currVVec.size() - 1]]]);

			bool haveMidM1 = false;
			for (int j = 1; j < onWhichC.size() - 1; ++j)
			{
				if (onWhichC[j] == -1)
				{
					haveMidM1 = true;
					break;
				}
			}
			if (!haveMidM1)
				continue;

			//std::vector<int> tempOWC = onWhichC;
			//std::sort(tempOWC.begin(), tempOWC.end());
			//tempOWC.erase(std::unique(tempOWC.begin(), tempOWC.end()), tempOWC.end());
			///*if (tempOWC.size() == 1 || (tempOWC.size() == 2 && tempOWC[0] == -1))
			//	continue;*/
			//bool haveFound = false;
			//for (int j = 0; j < tempOWC.size(); ++j)
			//{
			//	if (tempOWC[j] >= 0 && curvesFound[tempOWC[j]])
			//	{
			//		haveFound = true;
			//		break;
			//	}
			//}
			//if (haveFound)
			//	continue;

			if (cornerPos.size()<3)
				continue;

			//首先更新一下各个段所在的curves
			std::vector<int> curveC;
			std::vector<Eigen::Vector2d> ori_L, axa_L;
			ori_L.reserve(cornerPos.size() + 1);
			axa_L.reserve(cornerPos.size() + 1);
			std::vector<Eigen::Vector2d> cornersP;
			cornersP.reserve(cornerPos.size() + 1);
			std::vector<uint32_t> realCurveNums;	//按顺序保存所有碰到的边界边实际上的段数。
			realCurveNums.reserve(cornerPos.size() + 1);
			std::vector<double> realLength;
			realLength.reserve(cornerPos.size() + 1);
			std::vector<bool> cornerV2;
			cornerV2.reserve(cornerPos.size() + 1);
			cornerV2.emplace_back(true);
			bool shouldCon0 = false;
			for (int j = 0; j < cornerPos.size() - 1; ++j)
			{
				//ori_L.emplace_back((newQm.V_[currVVec[cornerPos[j] + 1]] + newQm.V_[currVVec[cornerPos[j]]])*0.5);
				ori_L.emplace_back(newQm.V_[currVVec[cornerPos[j]]]);
				axa_L.emplace_back((newQm.V_[currVVec[cornerPos[j] + 1]] - newQm.V_[currVVec[cornerPos[j]]]).normalized());
				realCurveNums.emplace_back(cornerPos[j + 1] - cornerPos[j]);
				double len = 0;
				for (int k = cornerPos[j]; k < cornerPos[j + 1]; ++k)
				{
					len += (newQm.V_[currVVec[k + 1]] - newQm.V_[currVVec[k]]).norm();
				}
				realLength.emplace_back(len);
				if (j != 0)
				{
					cornersP.emplace_back(newQm.V_[currVVec[cornerPos[j]]]);
					if (newQm.Vs_[currVVec[cornerPos[j]]].neighbor_fs.size() == 1)
						cornerV2.emplace_back(true);
					else
						cornerV2.emplace_back(false);
				}
				int currC = -1;
				if (cornerPos[j + 1] - cornerPos[j] != 1)
				{
					currC = onWhichC[cornerPos[j] + 1];
				}
				else if (j == 0 && onWhichC[0] != -1 && onWhichC[1] == -1)
				{
					currC = onWhichC[0];
				}
				else if (j == cornerPos.size() - 2 && onWhichC[currVVec.size() - 2] == -1 && onWhichC[currVVec.size() - 1] != -1)
				{
					currC = onWhichC[currVVec.size() - 1];
				}
				else
				{
					Eigen::Vector2d &cPos0 = newQm.V_[currVVec[cornerPos[j]]], &cPos1 = newQm.V_[currVVec[cornerPos[j + 1]]];
					uint32_t cor0 = (uint32_t)-1, cor1 = (uint32_t)-1;
					FindSourceCorner(cPos0, cor0);
					FindSourceCorner(cPos1, cor1);
					std::tuple<int, int> &corCur0 = mf.cornerCurves[cor0], &corCur1 = mf.cornerCurves[cor1];
					if (std::get<0>(corCur0) == std::get<0>(corCur1) || std::get<0>(corCur0) == std::get<1>(corCur1))
						currC = std::get<0>(corCur0);
					else if (std::get<1>(corCur0) == std::get<0>(corCur1) || std::get<1>(corCur0) == std::get<1>(corCur1))
						currC = std::get<1>(corCur0);
				}
				
				if (currC < 0)
				{
					//exit(-4396);
					shouldCon0 = true;
					break;
				}
				curveC.emplace_back(currC);
			}
			cornerV2.emplace_back(true);
			if (shouldCon0)
				continue;

			//计算边的中点
			std::vector<bool> curveState; //如果是两个相邻的corner，中间没有curve点，则为true，否则false。
			curveState.reserve(cornerPos.size() - 1);
			std::vector<uint32_t> midVVec;
			midVVec.reserve(cornerPos.size() - 1);
			midVVec.emplace_back(cornerPos[0]);
			if (cornerPos[1] - cornerPos[0]!=1)
				curveState.emplace_back(false);
			else
				curveState.emplace_back(true);
			for (int j = 1; j < cornerPos.size() - 2; ++j)
			{
				if (cornerPos[j + 1] - cornerPos[j] == 1)
				{
					midVVec.emplace_back(cornerPos[j]);
					curveState.emplace_back(true);
				}
				else
				{
					midVVec.emplace_back(cornerPos[j] + (cornerPos[j + 1] - cornerPos[j]) / 2);
					curveState.emplace_back(false);
				}
			}
			midVVec.emplace_back(cornerPos[cornerPos.size()-1]);
			if (cornerPos[cornerPos.size() - 1] - cornerPos[cornerPos.size() - 2] != 1)
				curveState.emplace_back(false);
			else
				curveState.emplace_back(true);

			std::vector<std::tuple<uint32_t, uint32_t>> lrNum;
			bool shCon1 = false;
			for (int j = 0; j < midVVec.size() - 1; ++j)
			{
				uint32_t tempNum = midVVec[j + 1] - midVVec[j];
				if (curveState[j+1]&& j!=midVVec.size()-2)
				{
					++tempNum;
				}
				if (tempNum < 2)
				{
					//exit(-6324);
					shCon1 = true;
					break;
				}

				double sumLength0 = 0;
				if (curveState[j] && j!=0)
				{
					for (int k = midVVec[j]; k < cornerPos[j + 1]; ++k)
					{
						sumLength0 += (newQm.V_[currVVec[k + 1]] - newQm.V_[currVVec[k]]).norm()/2;
					}
				}
				else
				{
					for (int k = midVVec[j]; k < cornerPos[j + 1]; ++k)
					{
						sumLength0 += (newQm.V_[currVVec[k + 1]] - newQm.V_[currVVec[k]]).norm();
					}
				}
				double sumLength1 = 0;
				if (curveState[j+1] && j != midVVec.size() - 2)
				{
					for (int k = cornerPos[j + 1]; k < midVVec[j] + tempNum; ++k)
					{
						sumLength1 += (newQm.V_[currVVec[k + 1]] - newQm.V_[currVVec[k]]).norm()/2;
					}
				}
				else
				{
					for (int k = cornerPos[j + 1]; k < midVVec[j] + tempNum; ++k)
					{
						sumLength1 += (newQm.V_[currVVec[k + 1]] - newQm.V_[currVVec[k]]).norm();
					}
				}

				int leftNum = std::round(tempNum * sumLength0 / (sumLength0 + sumLength1));
				if (leftNum < 1)
					leftNum = 1;
				else if (leftNum >= tempNum)
					leftNum = tempNum - 1;

				lrNum.emplace_back(std::tuple<uint32_t, uint32_t>(leftNum, tempNum - leftNum));
			}
			if (shCon1)
				continue;

			std::vector<uint32_t> targetENum;
			uint32_t tempENum = 0;
			targetENum.reserve(cornerPos.size() - 1);
			targetENum.emplace_back(std::get<0>(lrNum[0]));
			tempENum += targetENum[targetENum.size() - 1];
			for (int j = 0; j < lrNum.size()-1; ++j)
			{
				if (curveState[j + 1])
					targetENum.emplace_back(std::get<1>(lrNum[j]) + std::get<0>(lrNum[j + 1]) - 1);
				else
					targetENum.emplace_back(std::get<1>(lrNum[j]) + std::get<0>(lrNum[j + 1]));
				tempENum += targetENum[targetENum.size() - 1];
			}
			targetENum.emplace_back(currVVec.size() - 1 - tempENum);

			//找新的corners和curve vertices.
			std::vector<uint32_t> newCorners; //保存的是currVVec的下标，表示currVVec的第几个元素是新顶点。
			int newId = 0;
			/*if (cornerPos[0] != 0)
			{
				int remainNotMove = realCurveNums[0] - std::get<1>(onWhichCInfo[0]);
				newId = targetENum[0] - remainNotMove;
				if (newId <= 0)
					newId = 1;
				newCorners.emplace_back(newId);
			}*/
			int offset = 1;
			for (int j = 1; j < cornerPos.size() - 1; ++j)
			{
				newId += targetENum[j - offset];
				newCorners.emplace_back(newId);
			}
			if (newId >= currVVec.size() - 1)
				//exit(-3442);
				continue;
			for (int j = 1; j < cornerPos.size() - 1; ++j)
			{
				if (cornerV2[j])
					newCorners[j - 1] = cornerPos[j];
			}

			std::vector<uint32_t> tempNC = newCorners;
			std::sort(tempNC.begin(), tempNC.end());
			if (std::unique(tempNC.begin(), tempNC.end()) != tempNC.end())
				exit(-2800);
			/*if (cornerPos[cornerPos.size() - 1] != currVVec.size() - 1)
			{
				newId += targetENum[cornerPos.size() - 1 - offset];
				newCorners.emplace_back(newId);
			}*/
			/*if (newId > currVVec.size() - 1)
				continue;*/

				//重建oi->localFc
			uint32_t currInC = 0;
			uint32_t currCurveCount = 0, currCorner = 0;
			for (int j = 1; j < currVVec.size() - 1; ++j)
			{
				uint32_t localVId = newOi->globalToLocalVVec[currVVec[j]];
				if (localVId == (uint32_t)-1 || vFFlag[localVId])
					continue;

				vFFlag[localVId] = true;
				if (currCorner < newCorners.size() && j == newCorners[currCorner])
				{
					currInC = 0;
					newLocalFc.ids_C.emplace_back(localVId);
					newLocalFc.C.emplace_back(cornersP[currCorner]);
					newLocalFc.V_types[localVId] = -1;
					++currCorner;
					++currCurveCount;
				}
				else
				{
					newLocalFc.ids_L.emplace_back(localVId);
					newLocalFc.on_which_L.emplace_back(curveC[currCurveCount]);
					newLocalFc.axa_L.emplace_back(axa_L[currCurveCount]);
					newLocalFc.origin_L.emplace_back(ori_L[currCurveCount] + axa_L[currCurveCount] * (realLength[currCurveCount]/targetENum[currCurveCount])*(currInC+1));
					newLocalFc.V_types[localVId] = curveC[currCurveCount];
					++currInC;
				}
			}
			/*for (int j = 0; j < tempOWC.size(); ++j)
			{
				if (tempOWC[j] >= 0)
					curvesFound[tempOWC[j]] = true;
			}*/
			
		}

		for (int j = 0; j < newOi->localFc.V_types.size(); ++j)
		{
			if (vFFlag[j])
				continue;

			newLocalFc.V_types[j] = newOi->localFc.V_types[j];
		}
		for (int j = 0; j < newOi->localFc.ids_C.size(); ++j)
		{
			if (vFFlag[newOi->localFc.ids_C[j]])
				continue;

			newLocalFc.C.emplace_back(newOi->localFc.C[j]);
			newLocalFc.ids_C.emplace_back(newOi->localFc.ids_C[j]);
		}
		for (int j = 0; j < newOi->localFc.ids_L.size(); ++j)
		{
			if (vFFlag[newOi->localFc.ids_L[j]])
				continue;
			newLocalFc.ids_L.emplace_back(newOi->localFc.ids_L[j]);
			newLocalFc.on_which_L.emplace_back(newOi->localFc.on_which_L[j]);
			newLocalFc.origin_L.emplace_back(newOi->localFc.origin_L[j]);
			newLocalFc.axa_L.emplace_back(newOi->localFc.axa_L[j]);
		}

		newOi->localFc = newLocalFc;
		delete[] vFFlag;
		vFFlag = NULL;
		delete[] oiVFlag;
		oiVFlag = NULL;
		delete[] findVFlag;
		findVFlag = NULL;
	}

	void Pipeline::FindFaceSomeRings(QuadMesh &qm, std::vector<uint32_t> &inFaces, std::vector<uint32_t> &outFaces, int rings)
	{

		bool *fFlag = new bool[qm.Fs_.size()];
		std::memset(fFlag, 0, qm.Fs_.size() * sizeof(bool));
		std::stack<uint32_t> fsPool, tempFsPool;

		for (int i = 0; i < inFaces.size(); ++i)
		{
			fsPool.emplace(inFaces[i]);
			outFaces.emplace_back(inFaces[i]);
			fFlag[inFaces[i]] = true;
		}

		for (int i = 0; i < rings; ++i)
		{
			while (!fsPool.empty())
			{
				uint32_t currF = fsPool.top();
				fsPool.pop();

				std::vector<uint32_t> &ffs = qm.Fs_[currF].neighbor_fs;
				for (int j = 0; j < ffs.size(); ++j)
				{
					if (!fFlag[ffs[j]])
					{
						tempFsPool.emplace(ffs[j]);
						outFaces.emplace_back(ffs[j]);
						fFlag[ffs[j]] = true;
					}
				}
			}
			fsPool.swap(tempFsPool);
			std::stack<uint32_t>().swap(tempFsPool);
		}

		delete[] fFlag;
		fFlag = NULL;
	}

	void Pipeline::BuildNewOptimizeInfo(QuadMesh &qm, FeatureConstraints &fcc, std::vector<uint32_t> &faces, OptimizeInfo *newOi)
	{
		bool *vFlag = new bool[qm.Vs_.size()];
		std::memset(vFlag, 0, qm.Vs_.size() * sizeof(bool));
		bool *fFlag = new bool[qm.Fs_.size()];
		std::memset(fFlag, 0, qm.Fs_.size() * sizeof(bool));
		
		
		for (int i = 0; i < faces.size(); ++i)
		{
			fFlag[faces[i]] = true;
			std::vector<uint32_t> &fvs = qm.Fs_[faces[i]].vs;
			for (int j = 0; j < fvs.size(); ++j)
			{
				vFlag[fvs[j]] = true;
			}
		}

		newOi->globalToLocalVVec.resize(qm.Vs_.size(), (uint32_t)-1);
		newOi->localToGlobalVVec.reserve(qm.Vs_.size());
		uint32_t currNewId = 0;
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if (vFlag[i])
			{
				newOi->globalToLocalVVec[i] = currNewId++;
				newOi->localToGlobalVVec.emplace_back(i);
				newOi->localVPos.emplace_back(qm.V_[i]);
			}
		}

		Eigen::Vector3i tempVi;
		for (int i = 0; i < faces.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[faces[i]].vs;
			for (int j = 0; j < 4; ++j)
			{
				tempVi[0] = newOi->globalToLocalVVec[fvs[quadTriTable[j][0]]];
				tempVi[1] = newOi->globalToLocalVVec[fvs[quadTriTable[j][1]]];
				tempVi[2] = newOi->globalToLocalVVec[fvs[quadTriTable[j][2]]];
				newOi->localTriIds.emplace_back(tempVi);
			}
		}
		ComputeStandardTriangles(qm, faces, newOi);

		newOi->localVsNotMove.reserve(newOi->localToGlobalVVec.size());
		newOi->localVsNotMovePos.reserve(newOi->localToGlobalVVec.size());
		for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
		{
			uint32_t currVId = newOi->localToGlobalVVec[i];
			std::vector<uint32_t> &vfs = qm.Vs_[currVId].neighbor_fs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (!fFlag[vfs[j]])
				{
					newOi->localVsNotMove.emplace_back(i);
					newOi->localVsNotMovePos.emplace_back(qm.V_[currVId]);
				}
			}
		}

		//Compute local Fc
		FeatureConstraints &localFc = newOi->localFc;
		localFc.V_types.reserve(newOi->localToGlobalVVec.size());
		for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
		{
			uint32_t currVId = newOi->localToGlobalVVec[i];
			localFc.V_types.emplace_back(fcc.V_types[currVId]);
		}

		for (int i = 0; i < fcc.ids_C.size(); ++i)
		{
			uint32_t newV = newOi->globalToLocalVVec[fcc.ids_C[i]];
			if (newV == (uint32_t)-1)
				continue;

			localFc.ids_C.emplace_back(newV);
			localFc.C.emplace_back(fcc.C[i]);
		}

		for (int i = 0; i < fcc.ids_L.size(); ++i)
		{
			uint32_t newV = newOi->globalToLocalVVec[fcc.ids_L[i]];
			if (newV == (uint32_t)-1)
				continue;

			localFc.ids_L.emplace_back(newV);
			localFc.on_which_L.emplace_back(fcc.on_which_L[i]);
			localFc.axa_L.emplace_back(fcc.axa_L[i]);
			localFc.origin_L.emplace_back(fcc.origin_L[i]);
		}

		delete[] vFlag;
		vFlag = NULL;
		delete[] fFlag;
		fFlag = NULL;
	}

	void Pipeline::ComputeStandardTriangles(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &fs, OptimizeInfo *oi)
	{
		std::vector<Eigen::Vector2d> standVVec(4);
		standVVec[0][0] = 0; standVVec[0][1] = 0;
		for (int i = 0; i < fs.size(); ++i)
		{
			double e0 = 0, e1 = 0;

			e0 = 1;
			e1 = 1;

			double referenceArea = initialArea_ / qm.Fs_.size();
			double radio = std::sqrt(referenceArea / (e0*e1));
			e0 *= radio; e1 *= radio;

			standVVec[1][0] = e0; standVVec[1][1] = 0;
			standVVec[2][0] = e0; standVVec[2][1] = e1;
			standVVec[3][0] = 0; standVVec[3][1] = e1;

			std::vector<Eigen::Vector2d> tempVVec(3);
			for (int j = 0; j < 4; ++j)
			{
				tempVVec[0] = standVVec[quadTriTable[j][0]];
				tempVVec[1] = standVVec[quadTriTable[j][1]];
				tempVVec[2] = standVVec[quadTriTable[j][2]];
				oi->localStanTriPos.emplace_back(tempVVec);
			}
		}
	}


	bool Pipeline::DoCollapseVersion2(QuadMesh &sourceQm, QuadMesh *newQm, OptimizeInfo *oi)
	{
		AESolver aes(oi, kdTree_, fields_);
		if (aes.isInverse)
		{
			/*SafeDeletePtr(newQm);
			newQm = NULL;
			SafeDeletePtr(newBc_);
			newBc_ = NULL;
			return false;*/
		}
		else
		{
			BuildOptimizeInfoHasBoundary(sourceQm, oi);
			ComputeIdsLDis(oi);
#if OPTIMIZE_FIELDS
			ExtraFieldEnergyAlone efea(sourceQm, oi, kdTree_, fields_, 1, 1, sourceQm.ComputeAverageEdgeLength());
#endif
			for (int i = 0; i < AES_ITER; ++i)
			{
				//				aes.efe_.ExtractFields();
				aes.Optimize();
				aes.SetValueBack(oi);
				ProjectToSouceLine(oi, mf);
#if USE_3To0
				//AdjustBoundaryPosInOptimization(sourceQm, oi);
#endif
				aes.SetValueFront(oi);
				if (aes.isInverse)
				{
					/*for (int i = 0; i < oi->localToGlobalVVec.size(); ++i)
					{
						uint32_t newVId = oi->localToGlobalVVec[i];
						sourceQm.V_[newVId] = oi->localVPos[i];
					}
					QuadMeshIO qmi;
					qmi.WriteQuadMesh(&sourceQm, "C:\\Users\\ChiZhang\\Desktop\\testOut.obj");*/

					/*SafeDeletePtr(newQm);
					newQm = NULL;
					SafeDeletePtr(newBc_);
					newBc_ = NULL;
					return false;*/
					break;
				}

#if OPTIMIZE_FIELDS
				efea.Initialize();
				efea.SetValueFront(oi);
				efea.Optimize();
				efea.SetValueBack(oi);
#endif
			}
		}

		for (int i = 0; i < oi->localToGlobalVVec.size(); ++i)
		{
			uint32_t newVId = vOldToNew_[oi->localToGlobalVVec[i]];
			if (newVId == (uint32_t)-1)
				exit(-934);
			newQm->V_[newVId] = oi->localVPos[i];
		}
		for (int i = 0; i < oi->localFc.ids_C.size(); ++i)		//这里要加这个循环的原因是在上面的代码中在新老映射的过程中角点可能没有保持住
		{
			uint32_t oldVId = oi->localFc.ids_C[i];
			uint32_t newVId = vOldToNew_[oi->localToGlobalVVec[oldVId]];
			if (newVId == (uint32_t)-1)
				exit(-935);
			newQm->V_[newVId] = oi->localVPos[oldVId];
		}

#ifdef OUTPUT_MID_MESH
		QuadMeshIO qmi;
		qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\testOut1.obj");
#endif
		std::vector<uint32_t> vNewNewToNew, vNewToNewNew, fNewNewToNew, fNewToNewNew;
		DeleteUselessElements(*newQm, vNewToNewNew, vNewNewToNew, fNewToNewNew, fNewNewToNew);
		FeatureConstraints newFc;
		UpdateFeatures(oi, fc, newFc, vNewToOld_, vNewNewToNew, vOldToNew_, vNewToNewNew);
		if (!vNewNewToNew.empty())		//表明删掉了度为2的内点，要更新面的信息。
		{
			std::vector<uint32_t> tempNewNewF;
			uint32_t tempF;
			for (int i = 0; i < newQmOptimizedFs_.size(); ++i)
			{
				tempF = fNewToNewNew[newQmOptimizedFs_[i]];
				if (tempF == (uint32_t)-1)
					exit(-4657);

				tempNewNewF.emplace_back(tempF);
			}
			newQmOptimizedFs_.swap(tempNewNewF);
		}

		//std::vector<uint32_t> outFs;
		//DeleteValence2Vertices(*newQm, newFc, outFs);
		//if (!outFs.empty())
		//{
		//	FindFaceSomeRings(*newQm, outFs, newQmOptimizedFs_, 4);
		//}
		std::sort(newQmOptimizedFs_.begin(), newQmOptimizedFs_.end());
		newQmOptimizedFs_.erase(std::unique(newQmOptimizedFs_.begin(), newQmOptimizedFs_.end()), newQmOptimizedFs_.end());
		FindMovedFaces(*newQm, newQmOptimizedFs_, newFc);	//由于计算误差等原因，可能出现结果随机的情况，所以要把允许最小角结果调的大一点。

		OptimizeInfo *newOi = new OptimizeInfo();
		BuildNewOptimizeInfo(*newQm, newFc, newQmOptimizedFs_, newOi);

#ifdef OUTPUT_MID_MESH
		qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\testOut3.obj");
#endif

#if CHANGE_THE_CORNERS
		bool suc = true;
		if (!newOi->localFc.ids_C.empty())
		{
			OptimizeInfo nnewOi = *newOi;
			QuadMesh nnewQm = *newQm;
			ChangeTheCornersInTwoCurves(*newQm, newOi);
#ifndef USE_SQUARE_SOLVER
			AESolver aes2(newOi, kdTree_, fields_);
			if (aes2.isInverse)
			{
				/*SafeDeletePtr(newOi);
				newOi = NULL;
				SafeDeletePtr(newQm);
				newQm = NULL;
				SafeDeletePtr(newBc_);
				newBc_ = NULL;
				return false;*/
			}
			ComputeIdsLDis(newOi);
#if OPTIMIZE_FIELDS
			ExtraFieldEnergyAlone efea2(*newQm, newOi, kdTree_, fields_, 1, 1, newQm->ComputeAverageEdgeLength());
#endif
			for (int i = 0; i < AES_ITER; ++i)
			{
				//			aes2.efe_.ExtractFields();
				aes2.Optimize();
				aes2.SetValueBack(newOi);
				ProjectToSouceLine(newOi, mf);
#if USE_3To0
				AdjustBoundaryPosInOptimization(*newQm, newOi);
#endif
				aes2.SetValueFront(newOi);
				if (aes2.isInverse)
				{
#ifdef OUTPUT_MID_MESH
					for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
					{
						newQm->V_[newOi->localToGlobalVVec[i]] = newOi->localVPos[i];
					}
					qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\error9.obj");
#endif
					suc = false;
					SetQuadMesh(&nnewQm, newQm);
					SetOI(&nnewOi, newOi);
					/*newQm = &nnewQm;
					newOi = &nnewOi;*/
					break;
				}

#if OPTIMIZE_FIELDS
				//QuadMeshIO qmi;
				//qmi.WriteQuadMesh(&sourceQm, "C:\\Users\\ChiZhang\\Desktop\\bfo.obj"); 
				//if (i != AES_ITER - 1)
				//{
				efea2.Initialize();
				efea2.SetValueFront(newOi);
				efea2.Optimize();
				efea2.SetValueBack(newOi);
				//}
				//qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\afo.obj");
#endif
			}
#else
			AESolverSquare aes2(newOi);
			if (aes2.isInverse)
				return false;
			for (int i = 0; i < AES_ITER; ++i)
			{
				aes2.Optimize();
				aes2.SetValueBack(newOi);
				ProjectToSouceLine(newOi, mf);
				aes2.SetValueFront(newOi);
				if (aes2.isInverse)
					return false;
			}
#endif
			if (suc)
			{
				for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
				{
					newQm->V_[newOi->localToGlobalVVec[i]] = newOi->localVPos[i];
				}
#if USE_3To0
				AdjustBoundaryPos(*newQm, newOi);
#endif
				/*vNewToNewNew.clear(); vNewNewToNew.clear();
				DeleteUselessElements(*newQm, vNewToNewNew, vNewNewToNew);*/
				if (!newQm->JudgeQuadMeshJacobi() || newQm->minJacobi_ < MIN_ALLOWED_JACOBI)
				{
					suc = false;
					SetQuadMesh(&nnewQm, newQm);
					SetOI(&nnewOi, newOi);
				}
			}
		}
		else
		{
			suc = false;
		}
#endif 

#if  CHANGE_THE_CORNERS
		if (newOi->localFc.ids_C.empty() || suc == false)
		{
#endif
#ifndef USE_SQUARE_SOLVER
			AESolver aes2(newOi, kdTree_, fields_);
			if (aes2.isInverse)
			{
				SafeDeletePtr(newOi);
				newOi = NULL;
				SafeDeletePtr(newQm);
				newQm = NULL;
				SafeDeletePtr(newBc_);
				newBc_ = NULL;
				return false;
			}
			ComputeIdsLDis(newOi);
			for (int i = 0; i < AES_ITER; ++i)
			{
				//			aes2.efe_.ExtractFields();
				aes2.Optimize();
				aes2.SetValueBack(newOi);
				ProjectToSouceLine(newOi, mf);
#if USE_3To0
				AdjustBoundaryPosInOptimization(*newQm, newOi);
#endif
				aes2.SetValueFront(newOi);
				if (aes2.isInverse)
				{
#ifdef OUTPUT_MID_MESH
					for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
					{
						newQm->V_[newOi->localToGlobalVVec[i]] = newOi->localVPos[i];
					}
					qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\error9.obj");
#endif
					SafeDeletePtr(newOi);
					newOi = NULL;
					SafeDeletePtr(newQm);
					newQm = NULL;
					SafeDeletePtr(newBc_);
					newBc_ = NULL;
					return false;
				}
			}
#else
			AESolverSquare aes2(newOi);
			if (aes2.isInverse)
				return false;
			for (int i = 0; i < AES_ITER; ++i)
			{
				aes2.Optimize();
				aes2.SetValueBack(newOi);
				ProjectToSouceLine(newOi, mf);
				aes2.SetValueFront(newOi);
				if (aes2.isInverse)
					return false;
			}
#endif

			for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
			{
				newQm->V_[newOi->localToGlobalVVec[i]] = newOi->localVPos[i];
			}
#if USE_3To0
			AdjustBoundaryPos(*newQm, newOi);
#endif
			/*vNewNewToNew.clear(); vNewToNewNew.clear();
			DeleteUselessElements(*newQm, vNewToNewNew, vNewNewToNew);*/
			if (!newQm->JudgeQuadMeshJacobi() || newQm->minJacobi_ < MIN_ALLOWED_JACOBI)
			{
#ifdef OUTPUT_MID_MESH
				QuadMeshIO qmi;
				qmi.WriteQuadMesh(newQm, ".\\inverseModel.obj");
#endif
				std::cout << "Inverse Jacobi! " << std::endl;
				SafeDeletePtr(newOi);
				newOi = NULL;
				SafeDeletePtr(newQm);
				newQm = NULL;
				SafeDeletePtr(newBc_);
				newBc_ = NULL;
				return false;
			}

#if  CHANGE_THE_CORNERS
		}
#endif

		if (!ComputeHausdorff(backupQm_, *newQm, hausdorffRadio, hausdorffThre))
		{
			SafeDeletePtr(newOi);
			newOi = NULL;
			SafeDeletePtr(newQm);
			newQm = NULL;
			SafeDeletePtr(newBc_);
			newBc_ = NULL;
			return false;
		}

#ifdef OUTPUT_MID_MESH
		qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\testOut2.obj");
#endif
		//exit(0);
		/*FeatureConstraints newFc;
		UpdateFeatures(newOi, fc, newFc);*/
		UpdateFeaturesFromOi(newFc, newOi);
		SafeDeletePtr(newOi);
		newOi = NULL;

		if (!isCornerBCV)
		{
#if !CORNER_BC
			newBc_->ExtractBaseComplex(newQm);
#else
			newBc_->ExtractBaseComplex(newQm, FEATURE_THRESHOLD);
#endif
		}
		else
		{
			corners_.clear();
			for (int i = 0; i < newFc.ids_C.size(); ++i)
			{
				corners_.emplace_back(newFc.ids_C[i]);
			}
			newBc_->ExtractBaseComplex(newQm, corners_);
		}
		if (newBc_->isError_)
		{
			SafeDeletePtr(newQm);
			newQm = NULL;
			SafeDeletePtr(newBc_);
			newBc_ = NULL;
			return false;
		}

		if (!TopologyCheck(*newBc_))
		{
			std::cout << "Didn't pass topology check! " << std::endl;
			SafeDeletePtr(newQm);
			newQm = NULL;
			SafeDeletePtr(newBc_);
			newBc_ = NULL;
			return false;
		}

		/*if (!vNewNewToNew.empty())
		{
			UpdateFeaturesAfterDelete(newFc, vNewToNewNew, vNewNewToNew);
			UpdateOldNewAfterDelete(vOldToNew_, vNewToOld_, vNewToNewNew, vNewNewToNew);
		}*/

		if (newFc.ids_C.size() != mf.corners.size())
		{
			std::cout << "Coners number incorrect!" << std::endl;
			SafeDeletePtr(newQm);
			newQm = NULL;
			SafeDeletePtr(newBc_);
			newBc_ = NULL;
			return false;
		}

		//如果存在边界度为2的点，那么在下面重新构造，再加一步优化。注意oi要重新构造，并且要更新fc。

		sourceQm = *newQm;
		SafeDeletePtr(newQm);
		newQm = NULL;
		bc_ = *newBc_;
		SafeDeletePtr(newBc_);
		newBc_ = NULL;
		fc = newFc;
		return true;
	}

	void Pipeline::UpdateFeaturesFromOi(FeatureConstraints &fcc, OptimizeInfo *oi)
	{
		bool *vFlag = new bool[oi->globalToLocalVVec.size()];
		std::memset(vFlag, 0, oi->globalToLocalVVec.size() * sizeof(bool));

		FeatureConstraints &localFc = oi->localFc;
		for (int i = 0; i < oi->localToGlobalVVec.size(); ++i)
		{
			uint32_t currVId = oi->localToGlobalVVec[i];
			fcc.V_types[currVId] = localFc.V_types[i];
		}

		std::vector<uint32_t> newIdsC;
		std::vector<Eigen::Vector2d> newC;
		newIdsC.reserve(fcc.ids_C.size());
		newC.reserve(fcc.ids_C.size());

		std::vector<uint32_t> newIdsL, newOnWhichL;
		std::vector<Eigen::Vector2d> newAxaL, newOriginL;
		newIdsL.reserve(fcc.ids_L.size());
		newOnWhichL.reserve(fcc.ids_L.size());
		newAxaL.reserve(fcc.ids_L.size());
		newOriginL.reserve(fcc.ids_L.size());

		for (int i = 0; i < localFc.ids_C.size(); ++i)
		{
			uint32_t currVId = oi->localToGlobalVVec[localFc.ids_C[i]];
			if (vFlag[currVId])
				continue;

			vFlag[currVId] = true;
			newIdsC.emplace_back(currVId);
			newC.emplace_back(localFc.C[i]);
		}
		for (int i = 0; i < localFc.ids_L.size(); ++i)
		{
			uint32_t currVId = oi->localToGlobalVVec[localFc.ids_L[i]];
			if (vFlag[currVId])
				continue;
			vFlag[currVId] = true;
			newIdsL.emplace_back(currVId);
			newOnWhichL.emplace_back(localFc.on_which_L[i]);
			newAxaL.emplace_back(localFc.axa_L[i]);
			newOriginL.emplace_back(localFc.origin_L[i]);
		}
		for (int i = 0; i < fcc.ids_C.size(); ++i)
		{
			if (vFlag[fcc.ids_C[i]])
				continue;
			vFlag[fcc.ids_C[i]] = true;

			newIdsC.emplace_back(fcc.ids_C[i]);
			newC.emplace_back(fcc.C[i]);
		}
		for (int i = 0; i < fcc.ids_L.size(); ++i)
		{
			if (vFlag[fcc.ids_L[i]])
				continue;
			vFlag[fcc.ids_L[i]] = true;

			newIdsL.emplace_back(fcc.ids_L[i]);
			newOnWhichL.emplace_back(fcc.on_which_L[i]);
			newAxaL.emplace_back(fcc.axa_L[i]);
			newOriginL.emplace_back(fcc.origin_L[i]);
		}

		fcc.ids_C = newIdsC;
		fcc.C = newC;
		fcc.ids_L = newIdsL;
		fcc.on_which_L = newOnWhichL;
		fcc.axa_L = newAxaL;
		fcc.origin_L = newOriginL;

		delete[] vFlag;
		vFlag = NULL;
	}

	bool Pipeline::DoCollapse(QuadMesh &sourceQm, QuadMesh *newQm, OptimizeInfo *oi)
	{
#ifdef OUTPUT_MID_MESH
		QuadMeshIO qmi;
#endif

#if USE_NO_INVERSE
		if (!noInverse)
		{
#endif
#if COLLAPSE_CORNER
			if (collapseCorner_)
			{
				/*QuadMeshIO qmi;
				qmi.WriteQuadMesh(&sourceQm, "C:\\Users\\ChiZhang\\Desktop\\testOut.obj");*/
				if (!CollapseCornersVersion2(sourceQm, ci_, oi, fc, vOldToNew_, vNewToOld_, isCollapseSheet_))
					return false;
			}
#endif
			AESolver aes(oi, kdTree_, fields_);
			if (aes.isInverse)
			{
				SafeDeletePtr(newQm);
				newQm = NULL;
				SafeDeletePtr(newBc_);
				newBc_ = NULL;
				return false;
			}
			BuildOptimizeInfoHasBoundary(sourceQm, oi);
			ComputeIdsLDis(oi);
#if OPTIMIZE_FIELDS
			ExtraFieldEnergyAlone efea(sourceQm, oi, kdTree_, fields_, 1, 1, sourceQm.ComputeAverageEdgeLength());
#endif
			for (int i = 0; i < AES_ITER; ++i)
			{
				//				aes.efe_.ExtractFields();
				aes.Optimize();
				aes.SetValueBack(oi);
				ProjectToSouceLine(oi, mf);
#if USE_3To0
				//AdjustBoundaryPosInOptimization(sourceQm, oi);
#endif
				aes.SetValueFront(oi);
				if (aes.isInverse)
				{
					/*for (int i = 0; i < oi->localToGlobalVVec.size(); ++i)
					{
						uint32_t newVId = oi->localToGlobalVVec[i];
						sourceQm.V_[newVId] = oi->localVPos[i];
					}
					QuadMeshIO qmi;
					qmi.WriteQuadMesh(&sourceQm, "C:\\Users\\ChiZhang\\Desktop\\testOut.obj");*/
					SafeDeletePtr(newQm);
					newQm = NULL;
					SafeDeletePtr(newBc_);
					newBc_ = NULL;
					return false;
				}

#if OPTIMIZE_FIELDS
				//if (i != AES_ITER - 1)
				//{
				efea.Initialize();
				efea.SetValueFront(oi);
				efea.Optimize();
				efea.SetValueBack(oi);
				//}
#endif
			}

			//Set new pos to new mesh.
			for (int i = 0; i < oi->localToGlobalVVec.size(); ++i)
			{
				uint32_t newVId = vOldToNew_[oi->localToGlobalVVec[i]];
				if (newVId == (uint32_t)-1)
					exit(-934);
				newQm->V_[newVId] = oi->localVPos[i];
			}
			for (int i = 0; i < oi->localFc.ids_C.size(); ++i)		//这里要加这个循环的原因是在上面的代码中在新老映射的过程中角点可能没有保持住
			{
				uint32_t oldVId = oi->localFc.ids_C[i];
				uint32_t newVId = vOldToNew_[oi->localToGlobalVVec[oldVId]];
				if (newVId == (uint32_t)-1)
					exit(-935);
				newQm->V_[newVId] = oi->localVPos[oldVId];
			}
#if USE_NO_INVERSE
	}
#endif
		/*if (!newQm->JudgeQuadMeshJacobi() || newQm->minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			QuadMeshIO qmi;
			qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\testOut2.obj");
			std::cout << "Inverse Jacobi! " << std::endl;
			return false;
		}*/
#ifdef OUTPUT_MID_MESH
		qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\testOut1.obj");
#endif

		//对新网格进行局部优化，去除E_T这一项。
		OptimizeInfo *newOi = new OptimizeInfo();
		BuildNewMeshOptimizeInfo(oi, newOi, newQm);

		std::vector<uint32_t> vNewNewToNew, vNewToNewNew;
#if CHANGE_THE_CORNERS
		bool suc = true;
		if (!newOi->localFc.ids_C.empty())
		{
			OptimizeInfo nnewOi = *newOi;
			QuadMesh nnewQm = *newQm;
			ChangeTheCornersInTwoCurves(*newQm, newOi);
#ifndef USE_SQUARE_SOLVER
			AESolver aes2(newOi, kdTree_, fields_);
			if (aes2.isInverse)
			{
				/*SafeDeletePtr(newOi);
				newOi = NULL;
				SafeDeletePtr(newQm);
				newQm = NULL;
				SafeDeletePtr(newBc_);
				newBc_ = NULL;
				return false;*/
			}
			ComputeIdsLDis(newOi);
#if OPTIMIZE_FIELDS
			ExtraFieldEnergyAlone efea2(*newQm, newOi, kdTree_, fields_, 1, 1, newQm->ComputeAverageEdgeLength());
#endif
			for (int i = 0; i < AES_ITER; ++i)
			{
				//			aes2.efe_.ExtractFields();
				aes2.Optimize();
				aes2.SetValueBack(newOi);
				ProjectToSouceLine(newOi, mf);
#if USE_3To0
				AdjustBoundaryPosInOptimization(*newQm, newOi);
#endif
				aes2.SetValueFront(newOi);
				if (aes2.isInverse)
				{
#ifdef OUTPUT_MID_MESH
					for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
					{
						newQm->V_[newOi->localToGlobalVVec[i]] = newOi->localVPos[i];
					}
					qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\error9.obj");
#endif
					suc = false;
					SetQuadMesh(&nnewQm, newQm);
					SetOI(&nnewOi, newOi);
					/*newQm = &nnewQm;
					newOi = &nnewOi;*/
					break;
			}

#if OPTIMIZE_FIELDS
				//QuadMeshIO qmi;
				//qmi.WriteQuadMesh(&sourceQm, "C:\\Users\\ChiZhang\\Desktop\\bfo.obj"); 
				//if (i != AES_ITER - 1)
				//{
				efea2.Initialize();
				efea2.SetValueFront(newOi);
				efea2.Optimize();
				efea2.SetValueBack(newOi);
				//}
				//qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\afo.obj");
#endif
		}
#else
			AESolverSquare aes2(newOi);
			if (aes2.isInverse)
				return false;
			for (int i = 0; i < AES_ITER; ++i)
			{
				aes2.Optimize();
				aes2.SetValueBack(newOi);
				ProjectToSouceLine(newOi, mf);
				aes2.SetValueFront(newOi);
				if (aes2.isInverse)
					return false;
			}
#endif
			if (suc)
			{
				for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
				{
					newQm->V_[newOi->localToGlobalVVec[i]] = newOi->localVPos[i];
				}
#if USE_3To0
				AdjustBoundaryPos(*newQm, newOi);
#endif
				vNewToNewNew.clear(); vNewNewToNew.clear();
				DeleteUselessElements(*newQm, vNewToNewNew, vNewNewToNew);
				if (!newQm->JudgeQuadMeshJacobi() || newQm->minJacobi_ < MIN_ALLOWED_JACOBI)
				{
					suc = false;
					SetQuadMesh(&nnewQm, newQm);
					SetOI(&nnewOi, newOi);
				}
			}
		}
		else
		{
			suc = false;
		}
#endif 

#if  CHANGE_THE_CORNERS
		if (suc == false)
		{
#endif
#ifndef USE_SQUARE_SOLVER
			AESolver aes2(newOi, kdTree_, fields_);
			if (aes2.isInverse)
			{
				SafeDeletePtr(newOi);
				newOi = NULL;
				SafeDeletePtr(newQm);
				newQm = NULL;
				SafeDeletePtr(newBc_);
				newBc_ = NULL;
				return false;
			}
			ComputeIdsLDis(newOi);
			for (int i = 0; i < AES_ITER; ++i)
			{
				//			aes2.efe_.ExtractFields();
				aes2.Optimize();
				aes2.SetValueBack(newOi);
				ProjectToSouceLine(newOi, mf);
#if USE_3To0
				AdjustBoundaryPosInOptimization(*newQm, newOi);
#endif
				aes2.SetValueFront(newOi);
				if (aes2.isInverse)
				{
#ifdef OUTPUT_MID_MESH
					for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
					{
						newQm->V_[newOi->localToGlobalVVec[i]] = newOi->localVPos[i];
					}
					qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\error9.obj");
#endif
					SafeDeletePtr(newOi);
					newOi = NULL;
					SafeDeletePtr(newQm);
					newQm = NULL;
					SafeDeletePtr(newBc_);
					newBc_ = NULL;
					return false;
			}
		}
#else
			AESolverSquare aes2(newOi);
			if (aes2.isInverse)
				return false;
			for (int i = 0; i < AES_ITER; ++i)
			{
				aes2.Optimize();
				aes2.SetValueBack(newOi);
				ProjectToSouceLine(newOi, mf);
				aes2.SetValueFront(newOi);
				if (aes2.isInverse)
					return false;
			}
#endif

			for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
			{
				newQm->V_[newOi->localToGlobalVVec[i]] = newOi->localVPos[i];
			}
#if USE_3To0
			AdjustBoundaryPos(*newQm, newOi);
#endif
			vNewNewToNew.clear(); vNewToNewNew.clear();
			DeleteUselessElements(*newQm, vNewToNewNew, vNewNewToNew);
			if (!newQm->JudgeQuadMeshJacobi() || newQm->minJacobi_ < MIN_ALLOWED_JACOBI)
			{
#ifdef OUTPUT_MID_MESH
				QuadMeshIO qmi;
				qmi.WriteQuadMesh(newQm, ".\\inverseModel.obj");
#endif
				std::cout << "Inverse Jacobi! " << std::endl;
				SafeDeletePtr(newOi);
				newOi = NULL;
				SafeDeletePtr(newQm);
				newQm = NULL;
				SafeDeletePtr(newBc_);
				newBc_ = NULL;
				return false;
			}

#if  CHANGE_THE_CORNERS
		}
#endif

		if (!ComputeHausdorff(backupQm_, *newQm, hausdorffRadio, hausdorffThre))
		{
			SafeDeletePtr(newOi);
			newOi = NULL;
			SafeDeletePtr(newQm);
			newQm = NULL;
			SafeDeletePtr(newBc_);
			newBc_ = NULL;
			return false;
		}

#ifdef OUTPUT_MID_MESH
		qmi.WriteQuadMesh(newQm, "C:\\Users\\ChiZhang\\Desktop\\testOut2.obj");
#endif
		//exit(0);
		FeatureConstraints newFc;
		UpdateFeatures(newOi, fc, newFc);
		SafeDeletePtr(newOi);
		newOi = NULL;

		if (!vNewNewToNew.empty())
		{
			UpdateFeaturesAfterDelete(newFc, vNewToNewNew, vNewNewToNew);
			UpdateOldNewAfterDelete(vOldToNew_, vNewToOld_, vNewToNewNew, vNewNewToNew);
		}


		if (!isCornerBCV)
		{
#if !CORNER_BC
			newBc_->ExtractBaseComplex(newQm);
#else
			newBc_->ExtractBaseComplex(newQm, FEATURE_THRESHOLD);
#endif
		}
		else
		{
			corners_.clear();
			for (int i = 0; i < newFc.ids_C.size(); ++i)
			{
				corners_.emplace_back(newFc.ids_C[i]);
			}
			newBc_->ExtractBaseComplex(newQm, corners_);
		}
		if (newBc_->isError_)
		{
			SafeDeletePtr(newOi);
			newOi = NULL;
			SafeDeletePtr(newQm);
			newQm = NULL;
			SafeDeletePtr(newBc_);
			newBc_ = NULL;
			return false;
		}

		if (!TopologyCheck(*newBc_))
		{
			std::cout << "Didn't pass topology check! " << std::endl;
			SafeDeletePtr(newOi);
			newOi = NULL;
			SafeDeletePtr(newQm);
			newQm = NULL;
			SafeDeletePtr(newBc_);
			newBc_ = NULL;
			return false;
		}

		if (newFc.ids_C.size() != mf.corners.size())
		{
			std::cout << "Coners number incorrect!" << std::endl;
			SafeDeletePtr(newQm);
			newQm = NULL;
			SafeDeletePtr(newBc_);
			newBc_ = NULL;
			return false;
		}

		//如果存在边界度为2的点，那么在下面重新构造，再加一步优化。注意oi要重新构造，并且要更新fc。

		sourceQm = *newQm;
		SafeDeletePtr(newQm);
		newQm = NULL;
		bc_ = *newBc_;
		SafeDeletePtr(newBc_);
		newBc_ = NULL;
		fc = newFc;
		std::cout << "Collapse! " << std::endl;
		return true;
	}

	void Pipeline::SetQuadMesh(BaseDataStructure::QuadMesh *input, BaseDataStructure::QuadMesh *output)
	{
		output->V_ = input->V_;
		output->Vs_ = input->Vs_;
		output->Es_ = input->Es_;
		output->Fs_ = input->Fs_;
		output->minJacobi_ = input->minJacobi_;
		output->meshNormal_ = input->meshNormal_;
	}
	void Pipeline::SetOI(OptimizeInfo *input, OptimizeInfo *output)
	{
		output->globalToLocalVVec = input->globalToLocalVVec;
		output->localToGlobalVVec = input->localToGlobalVVec;
		output->localVPos = input->localVPos;
		output->localTriIds = input->localTriIds;
		output->localStanTriPos = input->localStanTriPos;
		output->localVsGroups = input->localVsGroups;
		output->localVsCoords = input->localVsCoords;
		output->vsGroupHasBoundary = input->vsGroupHasBoundary;
		output->localFc = input->localFc;
		output->idsLDis = input->idsLDis;
		output->localVsNotMove = input->localVsNotMove;
		output->localVsNotMovePos = input->localVsNotMovePos;
	}

	bool Pipeline::CollapseCorners(BaseDataStructure::QuadMesh &qm, CollapseInfo *ci, OptimizeInfo *oi, FeatureConstraints &currFc, std::vector<uint32_t> &vOldToNew, std::vector<uint32_t> &vNewToOld, bool isCollapseSheet)
	{
		if (!isCollapseSheet)
			return true;

		bool *ciVFlag = new bool[qm.Vs_.size()];	//该点是否是ci中的点
		std::memset(ciVFlag, 0, qm.Vs_.size() * sizeof(bool));
		bool *killVFlag = new bool[qm.Vs_.size()];	//该点是否是包含于killed quads
		std::memset(killVFlag, 0, qm.Vs_.size() * sizeof(bool));
		bool *findVFlag = new bool[qm.Vs_.size()]; //该点是否是被找过的
		std::memset(findVFlag, 0, qm.Vs_.size() * sizeof(bool));
		bool *leftInNewVs = new bool[qm.Vs_.size()];
		std::memset(leftInNewVs, 0, qm.Vs_.size() * sizeof(bool));
		bool *vFFlag = new bool[oi->localVPos.size()];
		std::memset(vFFlag, 0, oi->localVPos.size() * sizeof(bool));

		std::vector<uint32_t> ciVs;
		ciVs.reserve(qm.Vs_.size());
		for (int i = 0; i < ci->hsToBeKilled.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsToBeKilled[i]].vs;
			for (int j = 0; j < fvs.size(); ++j)
			{
				if (!qm.Vs_[fvs[j]].boundary)
					continue;

				ciVFlag[fvs[j]] = true;
				killVFlag[fvs[j]] = true;
				ciVs.emplace_back(fvs[j]);
			}
		}
		for (int i = 0; i < ci->hsSeveralring.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsSeveralring[i]].vs;
			for (int j = 0; j < fvs.size(); ++j)
			{
				if (!qm.Vs_[fvs[j]].boundary)
					continue;

				ciVFlag[fvs[j]] = true;
				ciVs.emplace_back(fvs[j]);
			}
		}
		std::sort(ciVs.begin(), ciVs.end());
		ciVs.erase(std::unique(ciVs.begin(), ciVs.end()), ciVs.end());

		bool isLoop = false;
		std::function<bool(uint32_t, uint32_t, uint32_t&, uint32_t&)> findNextVE = [&](uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE)->bool
		{
			std::vector<uint32_t> evs = qm.Es_[preE].vs;
			if (evs[0] == preV)
				nextV = evs[1];
			else
				nextV = evs[0];

			std::vector<uint32_t> ves = qm.Vs_[nextV].neighbor_es;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (qm.Es_[ves[i]].boundary && ves[i] != preE)
				{
					nextE = ves[i];
					break;
				}
			}

			if (findVFlag[nextV])
			{
				isLoop = true;
				return false;
			}
			if (!ciVFlag[nextV])
				return false;

			return true;
		};

		std::vector<std::vector<uint32_t>> ciVVecVec;
		for (int i = 0; i < ciVs.size(); ++i)
		{
			uint32_t startV = ciVs[i];
			if (killVFlag[startV] || !ciVFlag[startV] || findVFlag[startV])
				continue;

			isLoop = false;
			std::vector<uint32_t> vVec;
			uint32_t preV = startV, preE = (uint32_t)-1, nextV = (uint32_t)-1, nextE = (uint32_t)-1;
			vVec.emplace_back(preV);
			findVFlag[preV] = true;
			uint32_t preE2[2], currC = 0;
			std::vector<uint32_t> &ves = qm.Vs_[preV].neighbor_es;
			for (int j = 0; j < ves.size(); ++j)
			{
				if (qm.Es_[ves[j]].boundary)
					preE2[currC++] = ves[j];
			}

			preE = preE2[0];
			while (findNextVE(preV, preE, nextV, nextE))
			{
				vVec.emplace_back(nextV);
				findVFlag[nextV] = true;
				preV = nextV;
				preE = nextE;
			}

			if (!isLoop)
			{
				std::reverse(vVec.begin(), vVec.end());
				preV = startV;
				preE = preE2[1];

				while (findNextVE(preV, preE, nextV, nextE))
				{
					vVec.emplace_back(nextV);
					findVFlag[nextV] = true;
					preV = nextV;
					preE = nextE;
				}
				if (vVec.size() < 3)
					continue;

				if (killVFlag[vVec[vVec.size() - 1]])
				{
					if (killVFlag[vVec[1]])
					{
						delete[] ciVFlag;
						delete[] killVFlag;
						delete[] findVFlag;
						delete[] vFFlag;
						return false;
					}
					else
					{
						std::vector<uint32_t> vVec2;
						vVec2.resize(vVec.size());
						std::memcpy(&vVec2[0], &vVec[1], (vVec.size() - 1) * sizeof(uint32_t));
						vVec2.emplace_back(vVec[0]);
						ciVVecVec.emplace_back(vVec2);
						continue;
					}
				}
			}
			ciVVecVec.emplace_back(vVec);
		}


		FeatureConstraints newLocalFc;
		for (int i = 0; i < ciVVecVec.size(); ++i)
		{
			std::vector<uint32_t> &currVVec = ciVVecVec[i];
			std::vector<uint32_t> corners;
			std::vector<int> onWhichC;
			std::vector<uint32_t> killedVVec, newLeftVVec;
			killedVVec.reserve(currVVec.size());
			newLeftVVec.reserve(currVVec.size());
			onWhichC.reserve(currVVec.size());
			corners.reserve(currVVec.size());
			corners.emplace_back(0);
			onWhichC.emplace_back(currFc.V_types[currVVec[0]]);

			if (vOldToNew[currVVec[0]] != (uint32_t)-1)
				newLeftVVec.emplace_back(vOldToNew[currVVec[0]]);
			for (int j = 1; j < currVVec.size() - 1; ++j)
			{
				onWhichC.emplace_back(currFc.V_types[currVVec[j]]);
				if (currFc.V_types[currVVec[j]] == -1)
					corners.emplace_back(j);

				if (killVFlag[currVVec[j]])
					killedVVec.emplace_back(j);

				if (vOldToNew[currVVec[j]] != (uint32_t)-1)
					newLeftVVec.emplace_back(vOldToNew[currVVec[j]]);
			}
			corners.emplace_back(currVVec.size() - 1);
			onWhichC.emplace_back(currFc.V_types[currVVec[currVVec.size() - 1]]);
			if (vOldToNew[currVVec[currVVec.size() - 1]] != (uint32_t)-1)
				newLeftVVec.emplace_back(vOldToNew[currVVec[currVVec.size() - 1]]);

			if (killedVVec.empty())
				continue;

			bool haveMid = false;
			for (int j = 0; j < killedVVec.size(); ++j)
			{
				if (onWhichC[killedVVec[j]] == -1)
				{
					haveMid = true;
					break;
				}
			}
			if (!haveMid)
				continue;

			//std::memset(leftInNewVs, 0, qm.Vs_.size() * sizeof(bool));
			for (int j = 0; j < newLeftVVec.size(); ++j)
			{
				if (vNewToOld[newLeftVVec[j]] != (uint32_t)-1)
					leftInNewVs[vNewToOld[newLeftVVec[j]]] = true;
			}
			std::vector<uint32_t> realLeftVVec;
			realLeftVVec.reserve(currVVec.size());
			realLeftVVec.emplace_back(0);
			for (int j = 1; j < currVVec.size() - 1; ++j)
			{
				if (leftInNewVs[currVVec[j]])
					realLeftVVec.emplace_back(j);
			}
			realLeftVVec.emplace_back(currVVec.size() - 1);


			//更新各个段的信息
			std::vector<int> curves;
			std::vector<Eigen::Vector2d> ori_L, axa_L, cornerP;
			ori_L.reserve(corners.size() + 1);
			axa_L.reserve(corners.size() + 1);
			cornerP.reserve(corners.size() + 1);
			std::vector<double> realLength;
			realLength.reserve(corners.size() + 1);
			bool shouldCon0 = false;
			for (int j = 0; j < corners.size() - 1; ++j)
			{
				ori_L.emplace_back(qm.V_[currVVec[corners[j]]]);
				axa_L.emplace_back((qm.V_[currVVec[corners[j] + 1]] - qm.V_[currVVec[corners[j]]]).normalized());
				double len = 0;
				for (int k = corners[j]; k < corners[j + 1]; ++k)
				{
					len += (qm.V_[currVVec[k + 1]] - qm.V_[currVVec[k]]).norm();
				}
				realLength.emplace_back(len);
				if (j != 0)
				{
					cornerP.emplace_back(qm.V_[currVVec[corners[j]]]);
				}

				int currC = -1;
				if (corners[j + 1] - corners[j] != 1)
				{
					currC = onWhichC[corners[j] + 1];
				}
				else if (j == 0 && onWhichC[0] != -1 && onWhichC[1] == -1)
				{
					currC = onWhichC[0];
				}
				else if (j == corners.size() - 2 && onWhichC[currVVec.size() - 2] == -1 && onWhichC[currVVec.size() - 1] != -1)
				{
					currC = onWhichC[currVVec.size() - 1];
				}
				else
				{
					Eigen::Vector2d &cPos0 = qm.V_[currVVec[corners[j]]], &cPos1 = qm.V_[currVVec[corners[j + 1]]];
					uint32_t cor0 = (uint32_t)-1, cor1 = (uint32_t)-1;
					FindSourceCorner(cPos0, cor0);
					FindSourceCorner(cPos1, cor1);
					std::tuple<int, int> &corCur0 = mf.cornerCurves[cor0], &corCur1 = mf.cornerCurves[cor1];
					if (std::get<0>(corCur0) == std::get<0>(corCur1) || std::get<0>(corCur0) == std::get<1>(corCur1))
						currC = std::get<0>(corCur0);
					else if (std::get<1>(corCur0) == std::get<0>(corCur1) || std::get<1>(corCur0) == std::get<1>(corCur1))
						currC = std::get<1>(corCur0);
				}
				if (currC < 0)
				{
					//exit(-4396);
					shouldCon0 = true;
					break;
				}
				curves.emplace_back(currC);
			}
			if (shouldCon0)
			{
				delete[] ciVFlag;
				delete[] killVFlag;
				delete[] findVFlag;
				delete[] leftInNewVs;
				delete[] vFFlag;
				return false;
			}

			uint32_t allEdgeNum = realLeftVVec.size() - 1;
			if (cornerP.size() > realLeftVVec.size() - 2)
			{
				delete[] ciVFlag;
				delete[] killVFlag;
				delete[] findVFlag;
				delete[] leftInNewVs;
				delete[] vFFlag;
				return false;
			}

			double sumLength = std::accumulate(realLength.begin(), realLength.end(), 0.0);
			std::vector<uint32_t> targetNumVec;
			targetNumVec.reserve(realLength.size() + 1);
			uint32_t currSumENum = 0;
			for (int j = 0; j < realLength.size()-1; ++j)
			{
				int currENum = std::round(allEdgeNum *realLength[j] / sumLength);
				if (currENum == 0)
					currENum = 1;
				targetNumVec.emplace_back(currENum);
				currSumENum += currENum;
			}
			if (allEdgeNum - currSumENum <= 0)
			{
				delete[] ciVFlag;
				delete[] killVFlag;
				delete[] findVFlag;
				delete[] leftInNewVs;
				return false;
			}
			targetNumVec.emplace_back(allEdgeNum - currSumENum);

			std::vector<uint32_t> newCorners;
			newCorners.reserve(cornerP.size());
			uint32_t newIdPos = 0;
			for (int j = 0; j < targetNumVec.size()-1; ++j)
			{
				newIdPos += targetNumVec[j];
				newCorners.emplace_back(realLeftVVec[newIdPos]);
			}
			if (newIdPos >= realLeftVVec.size() - 1)
				exit(-3442);

			//处理度为2的点
			std::vector<bool> isCornerOcc(newCorners.size(), false);
			std::vector<uint32_t> va2V, corIdx;
			va2V.reserve(realLeftVVec.size());
			corIdx.reserve(realLeftVVec.size());
			for (int j = 0; j < realLeftVVec.size(); ++j)
			{
				int currV = currVVec[realLeftVVec[j]];
				if (qm.Vs_[currV].neighbor_fs.size() == 1)
				{
					int minDis = 1E10;
					int currIdx = -1;
					for (int k = 0; k < newCorners.size(); ++k)
					{
						int currDis = std::abs((int)newCorners[k] - (int)realLeftVVec[j]);
						if (currDis < minDis && !isCornerOcc[k])
						{
							minDis = currDis;
							currIdx = k;
						}
					}
					if (currIdx == -1)
					{
						delete[] ciVFlag;
						delete[] killVFlag;
						delete[] findVFlag;
						delete[] leftInNewVs;
						delete[] vFFlag;
						return false;
					}
					va2V.emplace_back(realLeftVVec[j]);
					corIdx.emplace_back(currIdx);
					isCornerOcc[currIdx] = true;
				}
			}
			if (!va2V.empty())
			{
				std::sort(corIdx.begin(), corIdx.end());
				for (int j = 0; j < corIdx.size(); ++j)
				{
					newCorners[corIdx[j]] = va2V[j];
				}
			}

			newLocalFc.V_types.resize(oi->localFc.V_types.size());
			uint32_t currInC = 0;
			uint32_t currCurveCount = 0, currCorner = 0;
			for (int j = 1; j < currVVec.size() - 1; ++j)
			{
				uint32_t localVId = oi->globalToLocalVVec[currVVec[j]];
				if (localVId == (uint32_t)-1 || vFFlag[localVId])
					continue;

				vFFlag[localVId] = true;
				if (currCorner < newCorners.size() && j == newCorners[currCorner])
				{
					currInC = 0;
					newLocalFc.ids_C.emplace_back(localVId);
					newLocalFc.C.emplace_back(cornerP[currCorner]);
					newLocalFc.V_types[localVId] = -1;
					++currCorner;
					++currCurveCount;
				}
				else
				{
					newLocalFc.ids_L.emplace_back(localVId);
					newLocalFc.on_which_L.emplace_back(curves[currCurveCount]);
					newLocalFc.axa_L.emplace_back(axa_L[currCurveCount]);
					newLocalFc.origin_L.emplace_back(ori_L[currCurveCount] + axa_L[currCurveCount] * (realLength[currCurveCount] / targetNumVec[currCurveCount])*(currInC + 1));
					newLocalFc.V_types[localVId] = curves[currCurveCount];
					++currInC;
				}
			}
		}

		for (int j = 0; j < oi->localFc.V_types.size(); ++j)
		{
			if (vFFlag[j])
				continue;

			newLocalFc.V_types[j] = oi->localFc.V_types[j];
		}
		for (int j = 0; j < oi->localFc.ids_C.size(); ++j)
		{
			if (vFFlag[oi->localFc.ids_C[j]])
				continue;

			newLocalFc.C.emplace_back(oi->localFc.C[j]);
			newLocalFc.ids_C.emplace_back(oi->localFc.ids_C[j]);
		}
		for (int j = 0; j < oi->localFc.ids_L.size(); ++j)
		{
			if (vFFlag[oi->localFc.ids_L[j]])
				continue;
			newLocalFc.ids_L.emplace_back(oi->localFc.ids_L[j]);
			newLocalFc.on_which_L.emplace_back(oi->localFc.on_which_L[j]);
			newLocalFc.origin_L.emplace_back(oi->localFc.origin_L[j]);
			newLocalFc.axa_L.emplace_back(oi->localFc.axa_L[j]);
		}

		oi->localFc = newLocalFc;

		delete[] ciVFlag;
		delete[] killVFlag;
		delete[] findVFlag;
		delete[] leftInNewVs;
		delete[] vFFlag;
		return true;
	}

	bool Pipeline::CollapseCornersVersion2(BaseDataStructure::QuadMesh &qm, CollapseInfo *ci, OptimizeInfo *oi, FeatureConstraints &currFc, std::vector<uint32_t> &vOldToNew, std::vector<uint32_t> &vNewToOld, bool isCollapseSheet)
	{
		if (!isCollapseSheet)
			return true;

		bool *ciVFlag = new bool[qm.Vs_.size()];	//该点是否是ci中的点
		std::memset(ciVFlag, 0, qm.Vs_.size() * sizeof(bool));
		bool *killVFlag = new bool[qm.Vs_.size()];	//该点是否是包含于killed quads
		std::memset(killVFlag, 0, qm.Vs_.size() * sizeof(bool));
		bool *findVFlag = new bool[qm.Vs_.size()]; //该点是否是被找过的
		std::memset(findVFlag, 0, qm.Vs_.size() * sizeof(bool));
		bool *leftInNewVs = new bool[qm.Vs_.size()];
		std::memset(leftInNewVs, 0, qm.Vs_.size() * sizeof(bool));
		bool *vFFlag = new bool[oi->localVPos.size()];
		std::memset(vFFlag, 0, oi->localVPos.size() * sizeof(bool));

		std::vector<uint32_t> ciVs;
		ciVs.reserve(qm.Vs_.size());
		for (int i = 0; i < ci->hsToBeKilled.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsToBeKilled[i]].vs;
			for (int j = 0; j < fvs.size(); ++j)
			{
				if (!qm.Vs_[fvs[j]].boundary)
					continue;

				ciVFlag[fvs[j]] = true;
				killVFlag[fvs[j]] = true;
				ciVs.emplace_back(fvs[j]);
			}
		}
		for (int i = 0; i < ci->hsSeveralring.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[ci->hsSeveralring[i]].vs;
			for (int j = 0; j < fvs.size(); ++j)
			{
				if (!qm.Vs_[fvs[j]].boundary)
					continue;

				ciVFlag[fvs[j]] = true;
				ciVs.emplace_back(fvs[j]);
			}
		}
		std::sort(ciVs.begin(), ciVs.end());
		ciVs.erase(std::unique(ciVs.begin(), ciVs.end()), ciVs.end());

		bool isLoop = false;
		std::function<bool(uint32_t, uint32_t, uint32_t&, uint32_t&)> findNextVE = [&](uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE)->bool
		{
			std::vector<uint32_t> evs = qm.Es_[preE].vs;
			if (evs[0] == preV)
				nextV = evs[1];
			else
				nextV = evs[0];

			std::vector<uint32_t> ves = qm.Vs_[nextV].neighbor_es;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (qm.Es_[ves[i]].boundary && ves[i] != preE)
				{
					nextE = ves[i];
					break;
				}
			}

			if (findVFlag[nextV])
			{
				isLoop = true;
				return false;
			}
			if (!ciVFlag[nextV])
				return false;

			return true;
		};

		std::vector<std::vector<uint32_t>> ciVVecVec;
		for (int i = 0; i < ciVs.size(); ++i)
		{
			uint32_t startV = ciVs[i];
			if (killVFlag[startV] || !ciVFlag[startV] || findVFlag[startV])
				continue;

			isLoop = false;
			std::vector<uint32_t> vVec;
			uint32_t preV = startV, preE = (uint32_t)-1, nextV = (uint32_t)-1, nextE = (uint32_t)-1;
			vVec.emplace_back(preV);
			findVFlag[preV] = true;
			uint32_t preE2[2], currC = 0;
			std::vector<uint32_t> &ves = qm.Vs_[preV].neighbor_es;
			for (int j = 0; j < ves.size(); ++j)
			{
				if (qm.Es_[ves[j]].boundary)
					preE2[currC++] = ves[j];
			}

			preE = preE2[0];
			while (findNextVE(preV, preE, nextV, nextE))
			{
				vVec.emplace_back(nextV);
				findVFlag[nextV] = true;
				preV = nextV;
				preE = nextE;
			}

			if (!isLoop)
			{
				std::reverse(vVec.begin(), vVec.end());
				preV = startV;
				preE = preE2[1];

				while (findNextVE(preV, preE, nextV, nextE))
				{
					vVec.emplace_back(nextV);
					findVFlag[nextV] = true;
					preV = nextV;
					preE = nextE;
				}
				if (vVec.size() < 3)
					continue;

				if (killVFlag[vVec[vVec.size() - 1]])
				{
					if (killVFlag[vVec[1]])
					{
						delete[] ciVFlag;
						delete[] killVFlag;
						delete[] findVFlag;
						delete[] vFFlag;
						return false;
					}
					else
					{
						std::vector<uint32_t> vVec2;
						vVec2.resize(vVec.size());
						std::memcpy(&vVec2[0], &vVec[1], (vVec.size() - 1) * sizeof(uint32_t));
						vVec2.emplace_back(vVec[0]);
						ciVVecVec.emplace_back(vVec2);
						continue;
					}
				}
			}
			ciVVecVec.emplace_back(vVec);
		}

		std::vector<std::vector<uint32_t>> tempCiVVecVec;
		for (int i = 0; i < ciVVecVec.size(); ++i)
		{
			if (ciVVecVec[i].size() < 2)
				continue;
			std::vector<uint32_t> va2VIdx;
			va2VIdx.emplace_back(0);
			for (int j = 1; j < ciVVecVec[i].size()-1; ++j)
			{
				if (!killVFlag[ciVVecVec[i][j]] && qm.Vs_[ciVVecVec[i][j]].neighbor_fs.size() == 1)
					va2VIdx.emplace_back(j);
			}
			va2VIdx.emplace_back(ciVVecVec[i].size() - 1);

			for (int j = 0; j < va2VIdx.size() - 1; ++j)
			{
				tempCiVVecVec.emplace_back(std::vector<uint32_t>());
				uint32_t currId = tempCiVVecVec.size() - 1;
				for (int k = va2VIdx[j]; k <= va2VIdx[j + 1]; ++k)
				{
					tempCiVVecVec[currId].emplace_back(ciVVecVec[i][k]);
				}
			}
		}
		ciVVecVec.swap(tempCiVVecVec);

		FeatureConstraints newLocalFc;
		for (int i = 0; i < ciVVecVec.size(); ++i)
		{
			std::vector<uint32_t> &currVVec = ciVVecVec[i];
			std::vector<uint32_t> corners;
			std::vector<int> onWhichC;
			std::vector<uint32_t> killedVVec, newLeftVVec;
			killedVVec.reserve(currVVec.size());
			newLeftVVec.reserve(currVVec.size());
			onWhichC.reserve(currVVec.size());
			corners.reserve(currVVec.size());
			corners.emplace_back(0);
			onWhichC.emplace_back(currFc.V_types[currVVec[0]]);

			if (vOldToNew[currVVec[0]] != (uint32_t)-1)
				newLeftVVec.emplace_back(vOldToNew[currVVec[0]]);
			for (int j = 1; j < currVVec.size() - 1; ++j)
			{
				onWhichC.emplace_back(currFc.V_types[currVVec[j]]);
				if (currFc.V_types[currVVec[j]] == -1)
					corners.emplace_back(j);

				if (killVFlag[currVVec[j]])
					killedVVec.emplace_back(j);

				if (vOldToNew[currVVec[j]] != (uint32_t)-1)
					newLeftVVec.emplace_back(vOldToNew[currVVec[j]]);
			}
			corners.emplace_back(currVVec.size() - 1);
			onWhichC.emplace_back(currFc.V_types[currVVec[currVVec.size() - 1]]);
			if (vOldToNew[currVVec[currVVec.size() - 1]] != (uint32_t)-1)
				newLeftVVec.emplace_back(vOldToNew[currVVec[currVVec.size() - 1]]);

			if (killedVVec.empty())
				continue;

			bool haveMid = false;
			for (int j = 0; j < killedVVec.size(); ++j)
			{
				if (onWhichC[killedVVec[j]] == -1)
				{
					haveMid = true;
					break;
				}
			}
			if (!haveMid)
				continue;

			//std::memset(leftInNewVs, 0, qm.Vs_.size() * sizeof(bool));
			for (int j = 0; j < newLeftVVec.size(); ++j)
			{
				if (vNewToOld[newLeftVVec[j]] != (uint32_t)-1)
					leftInNewVs[vNewToOld[newLeftVVec[j]]] = true;
			}
			std::vector<uint32_t> realLeftVVec;
			realLeftVVec.reserve(currVVec.size());
			realLeftVVec.emplace_back(0);
			for (int j = 1; j < currVVec.size() - 1; ++j)
			{
				if (leftInNewVs[currVVec[j]])
					realLeftVVec.emplace_back(j);
			}
			realLeftVVec.emplace_back(currVVec.size() - 1);


			//更新各个段的信息
			std::vector<int> curves;
			std::vector<Eigen::Vector2d> ori_L, axa_L, cornerP;
			ori_L.reserve(corners.size() + 1);
			axa_L.reserve(corners.size() + 1);
			cornerP.reserve(corners.size() + 1);
			std::vector<double> realLength;
			realLength.reserve(corners.size() + 1);
			bool shouldCon0 = false;
			for (int j = 0; j < corners.size() - 1; ++j)
			{
				ori_L.emplace_back(qm.V_[currVVec[corners[j]]]);
				axa_L.emplace_back((qm.V_[currVVec[corners[j] + 1]] - qm.V_[currVVec[corners[j]]]).normalized());
				double len = 0;
				for (int k = corners[j]; k < corners[j + 1]; ++k)
				{
					len += (qm.V_[currVVec[k + 1]] - qm.V_[currVVec[k]]).norm();
				}
				realLength.emplace_back(len);
				if (j != 0)
				{
					cornerP.emplace_back(qm.V_[currVVec[corners[j]]]);
				}

				int currC = -1;
				if (corners[j + 1] - corners[j] != 1)
				{
					currC = onWhichC[corners[j] + 1];
				}
				else if (j == 0 && onWhichC[0] != -1 && onWhichC[1] == -1)
				{
					currC = onWhichC[0];
				}
				else if (j == corners.size() - 2 && onWhichC[currVVec.size() - 2] == -1 && onWhichC[currVVec.size() - 1] != -1)
				{
					currC = onWhichC[currVVec.size() - 1];
				}
				else
				{
					Eigen::Vector2d &cPos0 = qm.V_[currVVec[corners[j]]], &cPos1 = qm.V_[currVVec[corners[j + 1]]];
					uint32_t cor0 = (uint32_t)-1, cor1 = (uint32_t)-1;
					FindSourceCorner(cPos0, cor0);
					FindSourceCorner(cPos1, cor1);
					std::tuple<int, int> &corCur0 = mf.cornerCurves[cor0], &corCur1 = mf.cornerCurves[cor1];
					if (std::get<0>(corCur0) == std::get<0>(corCur1) || std::get<0>(corCur0) == std::get<1>(corCur1))
						currC = std::get<0>(corCur0);
					else if (std::get<1>(corCur0) == std::get<0>(corCur1) || std::get<1>(corCur0) == std::get<1>(corCur1))
						currC = std::get<1>(corCur0);
				}
				if (currC < 0)
				{
					//exit(-4396);
					shouldCon0 = true;
					break;
				}
				curves.emplace_back(currC);
			}
			if (shouldCon0)
			{
				delete[] ciVFlag;
				delete[] killVFlag;
				delete[] findVFlag;
				delete[] leftInNewVs;
				delete[] vFFlag;
				return false;
			}

			uint32_t allEdgeNum = realLeftVVec.size() - 1;
			if (cornerP.size() > realLeftVVec.size() - 2)
			{
				delete[] ciVFlag;
				delete[] killVFlag;
				delete[] findVFlag;
				delete[] leftInNewVs;
				delete[] vFFlag;
				return false;
			}

			double sumLength = std::accumulate(realLength.begin(), realLength.end(), 0.0);
			std::vector<uint32_t> targetNumVec;
			targetNumVec.reserve(realLength.size() + 1);
			uint32_t currSumENum = 0;
			for (int j = 0; j < realLength.size() - 1; ++j)
			{
				int currENum = std::round(allEdgeNum *realLength[j] / sumLength);
				if (currENum == 0)
					currENum = 1;
				targetNumVec.emplace_back(currENum);
				currSumENum += currENum;
			}
			if (allEdgeNum - currSumENum <= 0)
			{
				delete[] ciVFlag;
				delete[] killVFlag;
				delete[] findVFlag;
				delete[] leftInNewVs;
				return false;
			}
			targetNumVec.emplace_back(allEdgeNum - currSumENum);

			std::vector<uint32_t> newCorners;
			newCorners.reserve(cornerP.size());
			uint32_t newIdPos = 0;
			for (int j = 0; j < targetNumVec.size() - 1; ++j)
			{
				newIdPos += targetNumVec[j];
				newCorners.emplace_back(realLeftVVec[newIdPos]);
			}
			if (newIdPos >= realLeftVVec.size() - 1)
				exit(-3442);

			//处理度为2的点
			std::vector<bool> isCornerOcc(newCorners.size(), false);
			std::vector<uint32_t> va2V, corIdx;
			va2V.reserve(realLeftVVec.size());
			corIdx.reserve(realLeftVVec.size());
			for (int j = 0; j < realLeftVVec.size(); ++j)
			{
				int currV = currVVec[realLeftVVec[j]];
				if (qm.Vs_[currV].neighbor_fs.size() == 1)
				{
					int minDis = 1E10;
					int currIdx = -1;
					for (int k = 0; k < newCorners.size(); ++k)
					{
						int currDis = std::abs((int)newCorners[k] - (int)realLeftVVec[j]);
						if (currDis < minDis && !isCornerOcc[k])
						{
							minDis = currDis;
							currIdx = k;
						}
					}
					if (currIdx == -1)
					{
						delete[] ciVFlag;
						delete[] killVFlag;
						delete[] findVFlag;
						delete[] leftInNewVs;
						delete[] vFFlag;
						return false;
					}
					va2V.emplace_back(realLeftVVec[j]);
					corIdx.emplace_back(currIdx);
					isCornerOcc[currIdx] = true;
				}
			}
			if (!va2V.empty())
			{
				std::sort(corIdx.begin(), corIdx.end());
				for (int j = 0; j < corIdx.size(); ++j)
				{
					newCorners[corIdx[j]] = va2V[j];
				}
			}

			newLocalFc.V_types.resize(oi->localFc.V_types.size());
			uint32_t currInC = 0;
			uint32_t currCurveCount = 0, currCorner = 0;
			for (int j = 1; j < currVVec.size() - 1; ++j)
			{
				uint32_t localVId = oi->globalToLocalVVec[currVVec[j]];
				if (localVId == (uint32_t)-1 || vFFlag[localVId])
					continue;

				vFFlag[localVId] = true;
				if (currCorner < newCorners.size() && j == newCorners[currCorner])
				{
					currInC = 0;
					newLocalFc.ids_C.emplace_back(localVId);
					newLocalFc.C.emplace_back(cornerP[currCorner]);
					newLocalFc.V_types[localVId] = -1;
					++currCorner;
					++currCurveCount;
				}
				else
				{
					newLocalFc.ids_L.emplace_back(localVId);
					newLocalFc.on_which_L.emplace_back(curves[currCurveCount]);
					newLocalFc.axa_L.emplace_back(axa_L[currCurveCount]);
					newLocalFc.origin_L.emplace_back(ori_L[currCurveCount] + axa_L[currCurveCount] * (realLength[currCurveCount] / targetNumVec[currCurveCount])*(currInC + 1));
					newLocalFc.V_types[localVId] = curves[currCurveCount];
					++currInC;
				}
			}
		}

		for (int j = 0; j < oi->localFc.V_types.size(); ++j)
		{
			if (vFFlag[j])
				continue;

			newLocalFc.V_types[j] = oi->localFc.V_types[j];
		}
		for (int j = 0; j < oi->localFc.ids_C.size(); ++j)
		{
			if (vFFlag[oi->localFc.ids_C[j]])
				continue;

			newLocalFc.C.emplace_back(oi->localFc.C[j]);
			newLocalFc.ids_C.emplace_back(oi->localFc.ids_C[j]);
		}
		for (int j = 0; j < oi->localFc.ids_L.size(); ++j)
		{
			if (vFFlag[oi->localFc.ids_L[j]])
				continue;
			newLocalFc.ids_L.emplace_back(oi->localFc.ids_L[j]);
			newLocalFc.on_which_L.emplace_back(oi->localFc.on_which_L[j]);
			newLocalFc.origin_L.emplace_back(oi->localFc.origin_L[j]);
			newLocalFc.axa_L.emplace_back(oi->localFc.axa_L[j]);
		}

		oi->localFc = newLocalFc;

		delete[] ciVFlag;
		delete[] killVFlag;
		delete[] findVFlag;
		delete[] leftInNewVs;
		delete[] vFFlag;
		return true;
	}

	void Pipeline::UpdateOldNewAfterDelete(std::vector<uint32_t> &vOldToNew, std::vector<uint32_t> &vNewToOld, std::vector<uint32_t> &vNewToNewNew, std::vector<uint32_t> &vNewNewToNew)
	{
		std::vector<uint32_t> tempOldToNew(vOldToNew.size(), (uint32_t)-1), tempNewToOld;
		for (int i = 0; i < vOldToNew.size(); ++i)
		{
			uint32_t midV = vOldToNew[i];
			if (midV != (uint32_t)-1)
				tempOldToNew[i] = vNewToNewNew[midV];
		}
		for (int i = 0; i < vNewNewToNew.size(); ++i)
		{
			tempNewToOld.emplace_back(vNewToOld[vNewNewToNew[i]]);
		}
		vNewToOld = tempNewToOld;
		vOldToNew = tempOldToNew;
	}

	void Pipeline::UpdateFeatures(OptimizeInfo *oldOi, FeatureConstraints &oldFc, FeatureConstraints &newFc, std::vector<uint32_t> &vNewToOld, std::vector<uint32_t> &vNewNewToNew, std::vector<uint32_t> &vOldToNew, std::vector<uint32_t> &vNewToNewNew)
	{
		newFc.V_types.reserve(vNewNewToNew.size());
		uint32_t iniId;
		
		bool *vFlag;
		if (!vNewNewToNew.empty())
		{
			for (int i = 0; i < vNewNewToNew.size(); ++i)
			{
				iniId = vNewToOld[vNewNewToNew[i]];
				if (iniId == (uint32_t)-1)
					newFc.V_types.emplace_back(-2);
				else
					newFc.V_types.emplace_back(oldFc.V_types[iniId]);
			}
			vFlag = new bool[vNewNewToNew.size()];
			std::memset(vFlag, false, vNewNewToNew.size() * sizeof(bool));
		}
		else
		{
			for (int i = 0; i < vNewToOld.size(); ++i)
			{
				iniId = vNewToOld[i];
				if (iniId == (uint32_t)-1)
					newFc.V_types.emplace_back(-2);
				else
					newFc.V_types.emplace_back(oldFc.V_types[iniId]);
			}
			vFlag = new bool[vNewToOld.size()];
			std::memset(vFlag, false, vNewToOld.size() * sizeof(bool));
		}

		uint32_t midId, finalId;
		if (!vNewNewToNew.empty())
		{
			for (int i = 0; i < oldOi->localFc.ids_C.size(); ++i)	//先搞角点，角点固定了之后就不能动了
			{
				midId = vOldToNew[oldOi->localToGlobalVVec[oldOi->localFc.ids_C[i]]];
				if (midId == (uint32_t)-1)
					continue;
				finalId = vNewToNewNew[midId];
				if (finalId == (uint32_t)-1)
					continue;

				if (vFlag[finalId])
					continue;
				vFlag[finalId] = true;

				newFc.V_types[finalId] = -1;
				newFc.ids_C.emplace_back(finalId);
				newFc.C.emplace_back(oldOi->localFc.C[i]);
			}
			for (int i = 0; i < oldOi->localFc.ids_L.size(); ++i)
			{
				midId = vOldToNew[oldOi->localToGlobalVVec[oldOi->localFc.ids_L[i]]];
				if (midId == (uint32_t)-1)
					continue;
				finalId = vNewToNewNew[midId];
				if (finalId == (uint32_t)-1)
					continue;

				if (vFlag[finalId])
					continue;
				vFlag[finalId] = true;

				//newFc.V_types[finalId] = oldOi->localFc.on_which_L[i];
				newFc.ids_L.emplace_back(finalId);
				newFc.axa_L.emplace_back(oldOi->localFc.axa_L[i]);
				newFc.on_which_L.emplace_back(oldOi->localFc.on_which_L[i]);
				newFc.origin_L.emplace_back(oldOi->localFc.origin_L[i]);
			}

			for (int i = 0; i < oldFc.ids_C.size(); ++i)
			{
				midId = vOldToNew[oldFc.ids_C[i]];
				if (midId == (uint32_t)-1)
					continue;
				finalId = vNewToNewNew[midId];
				if (finalId == (uint32_t)-1)
					continue;

				if (vFlag[finalId])
					continue;
				vFlag[finalId] = true;

				newFc.V_types[finalId] = -1;
				newFc.ids_C.emplace_back(finalId);
				newFc.C.emplace_back(oldFc.C[i]);
			}
			for (int i = 0; i < oldFc.ids_L.size(); ++i)
			{
				midId = vOldToNew[oldFc.ids_L[i]];
				if (midId == (uint32_t)-1)
					continue;
				finalId = vNewToNewNew[midId];
				if (finalId == (uint32_t)-1)
					continue;

				if (vFlag[finalId])
					continue;
				vFlag[finalId] = true;

				//newFc.V_types[finalId] = oldFc.on_which_L[i];
				newFc.ids_L.emplace_back(finalId);
				newFc.axa_L.emplace_back(oldFc.axa_L[i]);
				newFc.on_which_L.emplace_back(oldFc.on_which_L[i]);
				newFc.origin_L.emplace_back(oldFc.origin_L[i]);
			}
		}
		else
		{
			for (int i = 0; i < oldOi->localFc.ids_C.size(); ++i)	//先搞角点，角点固定了之后就不能动了
			{
				finalId = vOldToNew[oldOi->localToGlobalVVec[oldOi->localFc.ids_C[i]]];
				
				if (finalId == (uint32_t)-1)
					continue;

				if (vFlag[finalId])
					continue;
				vFlag[finalId] = true;

				newFc.V_types[finalId] = -1;
				newFc.ids_C.emplace_back(finalId);
				newFc.C.emplace_back(oldOi->localFc.C[i]);
			}
			for (int i = 0; i < oldOi->localFc.ids_L.size(); ++i)
			{
				finalId = vOldToNew[oldOi->localToGlobalVVec[oldOi->localFc.ids_L[i]]];
				if (finalId == (uint32_t)-1)
					continue;

				if (vFlag[finalId])
					continue;
				vFlag[finalId] = true;

				//newFc.V_types[finalId] = oldOi->localFc.on_which_L[i];
				newFc.ids_L.emplace_back(finalId);
				newFc.axa_L.emplace_back(oldOi->localFc.axa_L[i]);
				newFc.on_which_L.emplace_back(oldOi->localFc.on_which_L[i]);
				newFc.origin_L.emplace_back(oldOi->localFc.origin_L[i]);
			}

			for (int i = 0; i < oldFc.ids_C.size(); ++i)
			{
				finalId = vOldToNew[oldFc.ids_C[i]];
				if (finalId == (uint32_t)-1)
					continue;

				if (vFlag[finalId])
					continue;
				vFlag[finalId] = true;

				newFc.V_types[finalId] = -1;
				newFc.ids_C.emplace_back(finalId);
				newFc.C.emplace_back(oldFc.C[i]);
			}
			for (int i = 0; i < oldFc.ids_L.size(); ++i)
			{
				finalId = vOldToNew[oldFc.ids_L[i]];
				if (finalId == (uint32_t)-1)
					continue;

				if (vFlag[finalId])
					continue;
				vFlag[finalId] = true;

				//newFc.V_types[finalId] = oldFc.on_which_L[i];
				newFc.ids_L.emplace_back(finalId);
				newFc.axa_L.emplace_back(oldFc.axa_L[i]);
				newFc.on_which_L.emplace_back(oldFc.on_which_L[i]);
				newFc.origin_L.emplace_back(oldFc.origin_L[i]);
			}
		}
		delete[] vFlag;
		vFlag = NULL;
	}

	void Pipeline::UpdateFeatures(OptimizeInfo *newOi, FeatureConstraints &oldFc, FeatureConstraints &newFc)
	{
		
		std::vector<bool> vFlag(vNewToOld_.size(), false);

		newFc.V_types.resize(vNewToOld_.size());
		FeatureConstraints &oiFc = newOi->localFc;
		for (int i = 0; i < oiFc.ids_C.size(); ++i)
		{
			uint32_t globalIdsL = newOi->localToGlobalVVec[oiFc.ids_C[i]];
			if (!vFlag[globalIdsL])
			{
				vFlag[globalIdsL] = true;
				newFc.V_types[globalIdsL] = oiFc.V_types[oiFc.ids_C[i]];
				newFc.ids_C.emplace_back(globalIdsL);
				newFc.C.emplace_back(oiFc.C[i]);
			}
		}
		for (int i = 0; i < oiFc.ids_L.size(); ++i)
		{
			uint32_t globalIdsL = newOi->localToGlobalVVec[oiFc.ids_L[i]];
			if (!vFlag[globalIdsL])
			{
				vFlag[globalIdsL] = true;
				newFc.V_types[globalIdsL] = oiFc.V_types[oiFc.ids_L[i]];
				newFc.ids_L.emplace_back(globalIdsL);
				newFc.on_which_L.emplace_back(oiFc.on_which_L[i]);
				newFc.axa_L.emplace_back(oiFc.axa_L[i]);
				newFc.origin_L.emplace_back(oiFc.origin_L[i]);
			}
		}

		for (int i = 0; i < vNewToOld_.size(); ++i)
		{
			if (!vFlag[i])
				newFc.V_types[i] = oldFc.V_types[vNewToOld_[i]];
		}

		for (int i = 0; i < oldFc.ids_C.size(); ++i)
		{
			uint32_t newIdsC = vOldToNew_[oldFc.ids_C[i]];
			if (newIdsC!=(uint32_t)-1 && !vFlag[newIdsC])
			{
				vFlag[newIdsC] = true;
				newFc.ids_C.emplace_back(newIdsC);
				newFc.C.emplace_back(oldFc.C[i]);
			}
		}
		for (int i = 0; i < oldFc.ids_L.size(); ++i)
		{
			uint32_t newIdsL = vOldToNew_[oldFc.ids_L[i]];
			if ( newIdsL!=(uint32_t)-1 && !vFlag[newIdsL])
			{
				vFlag[newIdsL] = true;
				newFc.ids_L.emplace_back(newIdsL);
				newFc.on_which_L.emplace_back(oldFc.on_which_L[i]);
				newFc.axa_L.emplace_back(oldFc.axa_L[i]);
				newFc.origin_L.emplace_back(oldFc.origin_L[i]);
			}
		}
	}

	void Pipeline::UpdateFeaturesAfterDelete(FeatureConstraints &fc, std::vector<uint32_t> &vNewToNewNew, std::vector<uint32_t> &vNewNewToNew)
	{
		FeatureConstraints newFc;
		newFc.V_types.resize(vNewNewToNew.size(), -2);
		for (int i = 0; i < vNewNewToNew.size(); ++i)
		{
			newFc.V_types[i] = fc.V_types[vNewNewToNew[i]];
		}

		for (int i = 0; i < fc.ids_C.size(); ++i)
		{
			uint32_t newIdsC = vNewToNewNew[fc.ids_C[i]];
			if (newIdsC != (uint32_t)-1)
			{
				newFc.ids_C.emplace_back(newIdsC);
				newFc.C.emplace_back(fc.C[i]);
			}
		}
		for (int i = 0; i < fc.ids_L.size(); ++i)
		{
			uint32_t newIdsL = vNewToNewNew[fc.ids_L[i]];
			if (newIdsL != (uint32_t)-1)
			{
				newFc.ids_L.emplace_back(newIdsL);
				newFc.axa_L.emplace_back(fc.axa_L[i]);
				newFc.on_which_L.emplace_back(fc.on_which_L[i]);
				newFc.origin_L.emplace_back(fc.origin_L[i]);
			}
		}
		fc = newFc;
	}

	void Pipeline::ProjectToSouceLine(OptimizeInfo *oi, MeshFeatures &mf)
	{
		std::vector<uint32_t> &currIdsC = oi->localFc.ids_C;
		for (int i = 0; i < currIdsC.size(); ++i)
		{
			oi->localVPos[currIdsC[i]] = oi->localFc.C[i];
		}
		std::vector<uint32_t> &currIdsL = oi->localFc.ids_L;
		for (int i = 0; i < currIdsL.size(); ++i)
		{
			const Eigen::Vector2d &currP = oi->localVPos[currIdsL[i]];
			Point segP(currP[0], currP[1], 0);
			uint32_t onWhichL = oi->localFc.on_which_L[i];
			pAndP = mf.aabbTrees[onWhichL]->closest_point_and_primitive(segP);

			Point closeP = pAndP.first;
			oi->localVPos[currIdsL[i]][0] = closeP[0]; oi->localVPos[currIdsL[i]][1] = closeP[1];
			oi->localFc.origin_L[i] = oi->localVPos[currIdsL[i]];

			int index = std::distance(mf.segList[onWhichL].begin(), pAndP.second);
			/*auto it = mf.segList[onWhichL].begin();
			for (int j = 0; j < index; ++j)
			{
				++it;
			}*/
			const Segment &currSeg = mf.segList[onWhichL][index];

			const Point &p0 = currSeg.vertex(0);
			const Point &p1 = currSeg.vertex(1);

			double currNorm = Eigen::Vector2d(p1[0] - p0[0], p1[1] - p0[1]).norm();
			oi->localFc.axa_L[i] = Eigen::Vector2d(p1[0] - p0[0], p1[1] - p0[1]).normalized();
		}

		std::function<void(const Point &, uint32_t, bool)> UpdateCurrCurveIdsL = [&](const Point &cornerP, uint32_t onWhichCur, bool isC0)
		{
			globalVec2d_[0] = cornerP[0];
			globalVec2d_[1] = cornerP[1];
			if (isC0)	//如果是segList的左边点
			{
				uint32_t closestRV = (uint32_t)-1;
				double closestRVNorm = 1E30;
				Eigen::Vector2d closestVec;
				std::vector<uint32_t> badVers, badVersIdsLId;
				badVers.reserve(oi->idsLDis.size());
				badVersIdsLId.reserve(oi->idsLDis.size());
				for (int i = 0; i < oi->idsLDis.size(); ++i)
				{
					uint32_t onWL = std::get<1>(oi->idsLDis[i]);
					uint32_t currV = std::get<2>(oi->idsLDis[i]);
					const Eigen::Vector2d &currVPos = oi->localVPos[currV];
					double currNorm = (currVPos - globalVec2d_).norm();
					if (onWL == onWhichCur && currNorm < CORNER_ERROR_THRE)
					{
						badVersIdsLId.emplace_back(std::get<3>(oi->idsLDis[i]));
						badVers.emplace_back(currV);
					}
					else if (onWL == onWhichCur && currNorm > CORNER_ERROR_THRE)
					{
						if (currNorm < closestRVNorm)
						{
							closestRV = currV;
							closestRVNorm = currNorm;
							closestVec = currVPos;
						}
					}
				}

				if (badVers.empty() || closestRV == (uint32_t)-1)
					return;

				double unitL = closestRVNorm / (badVers.size() + 1);
				const Eigen::Vector2d &unitVec = (closestVec - globalVec2d_).normalized();
				for (int j = 0; j < badVers.size(); ++j)
				{
					const Eigen::Vector2d &tempVec = globalVec2d_ + (j + 1) * unitL * unitVec;
					oi->localVPos[badVers[j]] = tempVec;
					oi->localFc.axa_L[badVersIdsLId[j]] = unitVec;
					oi->localFc.origin_L[badVersIdsLId[j]] = tempVec;
				}
			}
			else
			{
				uint32_t closestRV = (uint32_t)-1;
				double closestRVNorm = 1E30;
				Eigen::Vector2d closestVec;
				std::vector<uint32_t> badVers, badVersIdsLId;
				badVers.reserve(oi->idsLDis.size());
				badVersIdsLId.reserve(oi->idsLDis.size());
				for (int i = oi->idsLDis.size()-1; i >= 0; --i)
				{
					uint32_t onWL = std::get<1>(oi->idsLDis[i]);
					uint32_t currV = std::get<2>(oi->idsLDis[i]);
					const Eigen::Vector2d &currVPos = oi->localVPos[currV];
					double currNorm = (currVPos - globalVec2d_).norm();
					if (onWL == onWhichCur && currNorm < CORNER_ERROR_THRE)
					{
						badVersIdsLId.emplace_back(std::get<3>(oi->idsLDis[i]));
						badVers.emplace_back(currV);
					}
					else if (onWL == onWhichCur && currNorm > CORNER_ERROR_THRE)
					{
						if (currNorm < closestRVNorm)
						{
							closestRV = currV;
							closestRVNorm = currNorm;
							closestVec = currVPos;
						}
					}
				}

				if (badVers.empty() || closestRV == (uint32_t)-1)
					return;

				double unitL = closestRVNorm / (badVers.size() + 1);
				const Eigen::Vector2d &unitVec = (closestVec - globalVec2d_).normalized();
				for (int j = 0; j < badVers.size(); ++j)
				{
					const Eigen::Vector2d &tempVec = globalVec2d_ + (j + 1) * unitL * unitVec;
					oi->localVPos[badVers[j]] = tempVec;
					oi->localFc.axa_L[badVersIdsLId[j]] = unitVec;
					oi->localFc.origin_L[badVersIdsLId[j]] = tempVec;
				}
			}
		};
		//防止边上的点投影到corner上，从而造成退化
		bool *edgeFlag = new bool[mf.aabbTrees.size()];
		std::memset(edgeFlag, 0, mf.aabbTrees.size()*sizeof(bool));
		for (int i = 0; i < oi->localFc.ids_L.size(); ++i)
		{
			uint32_t onWhichL = oi->localFc.on_which_L[i];
			if (edgeFlag[onWhichL])
				continue;
			edgeFlag[onWhichL] = true;

			const Point &c0 = mf.segList[onWhichL][0].vertex(0), &c1 = mf.segList[onWhichL][mf.segList[onWhichL].size() - 1].vertex(1);
			UpdateCurrCurveIdsL(c0, onWhichL, true);
			UpdateCurrCurveIdsL(c1, onWhichL, false);
		}
		delete[] edgeFlag;
		edgeFlag = NULL;
	}

	void Pipeline::TestBuildAABBTrees(Point targetP)
	{
		std::vector<Segment> linesVec;
		std::vector<Point> pointVec(5);
		for (int i = 0; i < 5; ++i)
		{
			pointVec[i] = Point(0, i-2, 0);
		}
		for (int i = 0; i < 4; ++i)
		{
			linesVec.emplace_back(Segment(pointVec[i], pointVec[i + 1]));
		}


		SegmentTree aabbTree(linesVec.begin(), linesVec.end());
		aabbTree.accelerate_distance_queries();

		auto pAndP = aabbTree.closest_point_and_primitive(targetP);
		std::cout << pAndP.first << std::endl;
	}

	void Pipeline::BuildSegmentAABBTree(MeshFeatures &mf)
	{
		//std::vector<std::vector<Segment>>().swap(mf.segList);
		mf.segList.clear();
		for (int i = 0; i < mf.aabbTrees.size(); ++i)
		{
			if (mf.aabbTrees[i] != NULL)
			{
				delete mf.aabbTrees[i];
				mf.aabbTrees[i] == NULL;
			}
		}
		//std::vector<SegmentTree*>().swap(mf.aabbTrees);
		mf.aabbTrees.clear();

		for (int i = 0; i < mf.curveVs.size(); ++i)
		{
			std::vector<uint32_t> & currCurveVs = mf.curveVs[i];
			mf.segList.emplace_back(std::vector<Segment>());
			
			for (int j = 0; j < currCurveVs.size()-1; ++j)
			{
				Point p0(mf.bPoints[currCurveVs[j]][0], mf.bPoints[currCurveVs[j]][1], 0);
				Point p1(mf.bPoints[currCurveVs[j + 1]][0], mf.bPoints[currCurveVs[j + 1]][1], 0);
				mf.segList[mf.segList.size()-1].emplace_back(Segment(p0, p1));
			}
			SegmentTree *currST = new SegmentTree(mf.segList[i].begin(), mf.segList[i].end());
			currST->accelerate_distance_queries();
			mf.aabbTrees.emplace_back(currST);
		}
	}

	void Pipeline::BuildNewMeshOptimizeInfo(OptimizeInfo *oldOi, OptimizeInfo *newOi, QuadMesh *newQm)
	{
		if (oldOi == NULL || newOi == NULL)
		{
			std::cout << "Pipeline::BuildNewMeshOptimizeInfo" << std::endl;
			exit(-26);
		}

		std::vector<bool> localVFlag(oldOi->localToGlobalVVec.size(), false);
		std::vector<uint32_t> newLToOldL, oldLToNewL;
		newLToOldL.reserve(oldOi->localToGlobalVVec.size());
		oldLToNewL.resize(oldOi->localToGlobalVVec.size(), (uint32_t)-1);

		uint32_t currNewLId = 0;
		for (int i = 0; i < oldOi->localVsGroups.size(); ++i)
		{
			std::vector<uint32_t> &currVsGroup = oldOi->localVsGroups[i];
			newLToOldL.emplace_back(currVsGroup[0]);
			for (int j = 0; j < currVsGroup.size(); ++j)
			{
				if (oldOi->localFc.V_types[currVsGroup[j]] == -1)
				{
					newLToOldL[i] = currVsGroup[j];
					break;
				}
				else if (oldOi->localFc.V_types[currVsGroup[j]] >= 0)
					newLToOldL[i] = currVsGroup[j];
			}
			for (int j = 0; j < currVsGroup.size(); ++j)
			{
				localVFlag[currVsGroup[j]] = true;
				oldLToNewL[currVsGroup[j]] = i;
			}
		}
		for (int i = 0; i < oldOi->localToGlobalVVec.size(); ++i)
		{
			if (localVFlag[i] == false)
			{
				localVFlag[i] = true;
				oldLToNewL[i] = newLToOldL.size();
				newLToOldL.emplace_back(i);
			}
		}

		newOi->localToGlobalVVec.reserve(newLToOldL.size());
		for (int i = 0; i < newLToOldL.size(); ++i)
		{
			uint32_t newGlobalVId = vOldToNew_[oldOi->localToGlobalVVec[newLToOldL[i]]];
			if (newGlobalVId == (uint32_t)-1)
				exit(-936);
			newOi->localToGlobalVVec.emplace_back(newGlobalVId);
		}
		newOi->globalToLocalVVec.resize(newQm->Vs_.size(), (uint32_t)-1);
		for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
		{
			newOi->globalToLocalVVec[newOi->localToGlobalVVec[i]] = i;
		}
		for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
		{
			newOi->localVPos.emplace_back(oldOi->localVPos[newLToOldL[i]]);
		}
		for (int i = 0; i < oldOi->localTriIds.size(); ++i)
		{
			newOi->localTriIds.emplace_back(Eigen::Vector3i(oldLToNewL[oldOi->localTriIds[i][0]], oldLToNewL[oldOi->localTriIds[i][1]], oldLToNewL[oldOi->localTriIds[i][2]]));
		}
		newOi->localStanTriPos = oldOi->localStanTriPos;
		for (int i = 0; i < oldOi->localVsNotMove.size(); ++i)
		{
			newOi->localVsNotMove.emplace_back(oldLToNewL[oldOi->localVsNotMove[i]]);
		}
		newOi->localVsNotMovePos = oldOi->localVsNotMovePos;

		//fc.V_ids没有更新。
		FeatureConstraints &oldFc = oldOi->localFc, &newFc = newOi->localFc;
		for (int i = 0; i < newOi->localToGlobalVVec.size(); ++i)
		{
			newFc.V_types.emplace_back(oldFc.V_types[newLToOldL[i]]);
		}
		for (int i = 0; i < oldFc.ids_C.size(); ++i)
		{
			newFc.ids_C.emplace_back(oldLToNewL[oldFc.ids_C[i]]);
		}
		newFc.C = oldFc.C;
		for (int i = 0; i < oldFc.ids_L.size(); ++i)
		{
			uint32_t newLId = oldLToNewL[oldFc.ids_L[i]];
			if (newFc.V_types[newLId] >= 0)			//一个智障的错误：特征边上的点映过来可能变成corner点，但是没有加判断
			{
				newFc.ids_L.emplace_back(oldLToNewL[oldFc.ids_L[i]]);
				newFc.on_which_L.emplace_back(oldFc.on_which_L[i]);
				newFc.axa_L.emplace_back(oldFc.axa_L[i]);
				newFc.origin_L.emplace_back(oldFc.origin_L[i]);
			}
		}
	/*	newFc.on_which_L = oldFc.on_which_L;
		newFc.axa_L = oldFc.axa_L;
		newFc.origin_L = oldFc.origin_L;*/
	}

	void Pipeline::ComputeVNewToOld(std::vector<uint32_t> &vOldToNew, std::vector<uint32_t> &vNewToOld, FeatureConstraints &fcc)
	{
		std::vector<std::vector<uint32_t>> wholeNewToOld(vNewToOld.size());
		//std::vector<uint32_t>().swap(vNewToOld);
		vNewToOld.clear();
		for (int i = 0; i < vOldToNew.size(); ++i)
		{
			if (vOldToNew[i] != (uint32_t)- 1)
				wholeNewToOld[vOldToNew_[i]].emplace_back(i);
		}
		for (int i = 0; i < wholeNewToOld.size(); ++i)
		{
			if (wholeNewToOld[i].empty())
				exit(-27);
			vNewToOld.emplace_back(wholeNewToOld[i][0]);
			for (int j = 0; j < wholeNewToOld[i].size(); ++j)
			{
				if (fcc.V_types[wholeNewToOld[i][j]] == -1)
				{
					vNewToOld[i] = wholeNewToOld[i][j];
					break;
				}
				else if (fcc.V_types[wholeNewToOld[i][j]] >= 0)
					vNewToOld[i] = wholeNewToOld[i][j];
			}
		}
	}

	void Pipeline::SubdivideSheetsFixV2(QuadMesh &qm, BaseComplex &bc, FeatureConstraints &fcc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<RankingTuple3> &candidates)
	{
		if (qm.Fs_.size() > initialQuadNum_ * subdivideThre_)
			return;

		minusValue = initialQuadNum_ * (stopSubThre_ - subdivideThre_);
		std::vector<uint32_t> startEVec, fNumVec;
#if !USE_FIELD
		FindLongestOpenSheet(qm, startEVec, fNumVec);
#else
		//FindLongestOpenSheet(qm, startEVec, fNumVec);
		FindNearestOpenSheet(qm, startEVec, fNumVec);
#endif
		uint32_t currNum = 0;

		while (currNum < startEVec.size())
		{
			if (qm.Fs_.size() > initialQuadNum_ * stopSubThre_)
			{
				ExtractSheetsAndChords(qm, bc, sheets_, chords_);
#if !USE_FIELD
				RankingSheetsAndChords(sheets_, chords_, candidates);
#else
				RankingSheetsAndChordsByFields(qm, bc, sheets_, chords_, candidates);
#endif
				return;
			}


			if (SubdivideSheetFixV2(qm, bc, fcc, startEVec[currNum], fNumVec[currNum]))
			{
				currNum = 0;
#if !USE_FIELD
				FindLongestOpenSheet(qm, startEVec, fNumVec);
#else
				//FindLongestOpenSheet(qm, startEVec, fNumVec);
				FindNearestOpenSheet(qm, startEVec, fNumVec);
#endif
				continue;
			}
			else
			{
				++currNum;
			}
		}
		//for (int i = candidates.size()-1; i >=0; --i)
		//{
		//	if (std::get<1>(candidates[i]) != ElementType::SHEET || sheets[std::get<2>(candidates[i])].st!=SheetType::OPEN)
		//		continue;

		//	//SubdivideSheet(qm, bc, fcc, sheets[std::get<2>(candidates[i])], candidates);
		//	SubdivideSheet(qm, fcc);
		//	if (qm.Fs_.size() > initialQuadNum_ * subdivideThre_)
		//	{
		//		bc.ExtractBaseComplex(&qm);
		//		ExtractSheetsAndChords(qm, bc, sheets_, chords_);
		//		RankingSheetsAndChords(sheets_, chords_, candidates);
		//		return;
		//	}
		//}
		ExtractSheetsAndChords(qm, bc, sheets_, chords_);
#if !USE_FIELD
		RankingSheetsAndChords(sheets_, chords_, candidates);
#else
		RankingSheetsAndChordsByFields(qm, bc, sheets_, chords_, candidates);
#endif
	}

	void Pipeline::SubdivideSheets(QuadMesh &qm, BaseComplex &bc, FeatureConstraints &fcc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<RankingTuple3> &candidates)
	{
		if (qm.Fs_.size() > initialQuadNum_ * subdivideThre_)
			return;

		minusValue = initialQuadNum_ * (stopSubThre_ - subdivideThre_);
		std::vector<uint32_t> startEVec, fNumVec;
#if !USE_FIELD
		FindLongestOpenSheet(qm, startEVec, fNumVec);
#else
		//FindLongestOpenSheet(qm, startEVec, fNumVec);
		FindNearestOpenSheet(qm, startEVec, fNumVec);
#endif
		uint32_t currNum = 0;

		while (currNum < startEVec.size())
		{
			if (qm.Fs_.size() > initialQuadNum_ * stopSubThre_)
			{
				ExtractSheetsAndChords(qm, bc, sheets_, chords_);
#if !USE_FIELD
				RankingSheetsAndChords(sheets_, chords_, candidates);
#else
				RankingSheetsAndChordsByFields(qm, bc, sheets_, chords_, candidates);
#endif
				return;
			}


#if NEW_SUBDIVIDE
			if (SubdivideSheetFixV2(qm, bc, fcc, startEVec[currNum], fNumVec[currNum]))
#else
			if (SubdivideSheet(qm, bc, fcc, startEVec[currNum], fNumVec[currNum]))
#endif
			{
				currNum = 0;
#if !USE_FIELD
				FindLongestOpenSheet(qm, startEVec, fNumVec);
#else
				//FindLongestOpenSheet(qm, startEVec, fNumVec);
				FindNearestOpenSheet(qm, startEVec, fNumVec);
#endif
				continue;
			}
			else
			{
				++currNum;
			}
		}
		//for (int i = candidates.size()-1; i >=0; --i)
		//{
		//	if (std::get<1>(candidates[i]) != ElementType::SHEET || sheets[std::get<2>(candidates[i])].st!=SheetType::OPEN)
		//		continue;

		//	//SubdivideSheet(qm, bc, fcc, sheets[std::get<2>(candidates[i])], candidates);
		//	SubdivideSheet(qm, fcc);
		//	if (qm.Fs_.size() > initialQuadNum_ * subdivideThre_)
		//	{
		//		bc.ExtractBaseComplex(&qm);
		//		ExtractSheetsAndChords(qm, bc, sheets_, chords_);
		//		RankingSheetsAndChords(sheets_, chords_, candidates);
		//		return;
		//	}
		//}
		ExtractSheetsAndChords(qm, bc, sheets_, chords_);
#if !USE_FIELD
		RankingSheetsAndChords(sheets_, chords_, candidates);
#else
		RankingSheetsAndChordsByFields(qm, bc, sheets_, chords_, candidates);
#endif
	}

	bool Pipeline::SubdivideSheet(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, FeatureConstraints &fcc, BaseDataStructure::Sheet &sheet, std::vector<RankingTuple3> &candidates)
	{
		typedef std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> TupleUInt32_4;
		FeatureConstraints newFc = fcc;
		QuadMesh newQm;
		newQm.V_ = qm.V_;
		std::vector<TupleUInt32_4> finalFVs;
		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
			finalFVs.emplace_back(TupleUInt32_4(fvs[0], fvs[1], fvs[2], fvs[3]));
		}

		int diff = initialQuadNum_*subdivideThre_ - qm.Fs_.size();
		int splitOneAdd = sheet.qf_delete.size()/bc.Be_[sheet.b_middle_es[0]].es_link.size();
		
		int splitTimes = std::ceil(((double)diff) / ((double)splitOneAdd));

		//给所有子排序
		uint32_t startMidBE = sheet.b_middle_boundary_es[0];
		std::vector<uint32_t> startEsLink = bc.Be_[startMidBE].es_link;
		SortEdgeFromLongToShort(startEsLink, qm);
		if (splitTimes > startEsLink.size())
			splitTimes = startEsLink.size();

		std::function<void(uint32_t, uint32_t, uint32_t&, Eigen::Vector2d&, Eigen::Vector2d&)> vNewFcInfo = [&](uint32_t lVId, uint32_t rVId, uint32_t& curveId, Eigen::Vector2d& oriL, Eigen::Vector2d& axaL)
		{
			int lVType = fcc.V_types[lVId], rVType = fcc.V_types[rVId];
			if (lVType < 0 && rVType < 0)
			{
				std::cout << "Can't add point between two corners! " << std::endl;
				return false;
			}
			if (lVType > 0)
				curveId = lVType;
			else
				curveId = rVType;

			const Eigen::Vector2d &lVPos = qm.V_[lVId];
			const Eigen::Vector2d &rVPos = qm.V_[rVId];
			oriL = (lVPos + rVPos) * 0.5;
			axaL = (lVPos - rVPos).normalized();
		};
		uint32_t growNum = splitOneAdd, startVId = qm.Vs_.size();
		std::vector<Eigen::Vector2d> newVPoss;
		std::vector<uint32_t> newVIds, fIds, leftVIds, rightVIds;
		for (int i = 0; i < splitTimes; ++i)
		{
			uint32_t startE = startEsLink[i];
			GrowSplitSheet(startE, growNum, startVId, newVPoss, newVIds, fIds, leftVIds, rightVIds, qm);
			newQm.V_.insert(newQm.V_.end(), newVPoss.begin(), newVPoss.end());

			for (int j = 0; j < fIds.size(); ++j)
			{
				TupleUInt32_4 &currFvs = finalFVs[fIds[j]];
				std::get<0>(currFvs) = leftVIds[j];
				std::get<1>(currFvs) = newVIds[j];
				std::get<2>(currFvs) = newVIds[j + 1];
				std::get<3>(currFvs) = leftVIds[j + 1];

				finalFVs.emplace_back(TupleUInt32_4(newVIds[j], rightVIds[j], rightVIds[j + 1], newVIds[j + 1]));
			}

			//Build new fc
			uint32_t curveId = (uint32_t)-1;
			Eigen::Vector2d oriL, axaL;
			vNewFcInfo(leftVIds[0], rightVIds[0], curveId, oriL, axaL);
			newFc.V_types.emplace_back(curveId);
			newFc.ids_L.emplace_back(newVIds[0]);
			newFc.on_which_L.emplace_back(curveId);
			newFc.axa_L.emplace_back(axaL);
			newFc.origin_L.emplace_back(oriL);

			for (int j = 1; j < newVIds.size() - 1; ++j)
			{
				newFc.V_types.emplace_back(-2);
			}

			vNewFcInfo(leftVIds[leftVIds.size()-1], rightVIds[rightVIds.size()-1], curveId, oriL, axaL);
			newFc.V_types.emplace_back(curveId);
			newFc.ids_L.emplace_back(newVIds[newVIds.size()-1]);
			newFc.on_which_L.emplace_back(curveId);
			newFc.axa_L.emplace_back(axaL);
			newFc.origin_L.emplace_back(oriL);
		}

		std::function<void(QuadMesh &, std::vector<TupleUInt32_4> &)> buildQmFromTuple4 = [&](QuadMesh &nQm, std::vector<TupleUInt32_4> &fFvs)
		{
			QuadVertex qv;
			for (int i = 0; i < nQm.V_.size(); ++i)
			{
				qv.id = i;
				qv.bId = (uint32_t)-1;
				nQm.Vs_.emplace_back(qv);
			}

			for (int i = 0; i < fFvs.size(); ++i)
			{
				nQm.Fs_.emplace_back(QuadFace());
				QuadFace &currF = nQm.Fs_[i];
				currF.id = i;
				currF.bId = (uint32_t)-1;

				uint32_t vs[4];
				vs[0] = std::get<0>(fFvs[i]);
				vs[1] = std::get<1>(fFvs[i]);
				vs[2] = std::get<2>(fFvs[i]);
				vs[3] = std::get<3>(fFvs[i]);
				for (int j = 0; j < 4; ++j)
				{
					currF.vs.emplace_back(vs[j]);
					nQm.Vs_[vs[j]].neighbor_fs.emplace_back(i);
				}
			}
		};

		buildQmFromTuple4(newQm, finalFVs);
		if (!newQm.BuildConnectivity())
			return false;
		BaseComplex newBc;

		if (!isCornerBCV)
		{
#if !CORNER_BC
			newBc.ExtractBaseComplex(&newQm);
#else
			newBc.ExtractBaseComplex(&newQm, FEATURE_THRESHOLD);
#endif	
		}
		else
		{
			corners_.clear();
			for (int i = 0; i < newFc.ids_C.size(); ++i)
			{
				corners_.emplace_back(newFc.ids_C[i]);
			}
			newBc.ExtractBaseComplex(&newQm, corners_);
		}
		if (newBc.isError_)
			return false;
		if (!TopologyCheck(newQm, newBc))
		{
			std::cout << "Didn't pass topology check! " << std::endl;
			return false;
		}
		if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			std::cout << "Inverse Jacobi! " << std::endl;
			return false;
		}

#ifdef OUTPUT_MID_MESH
		QuadMeshIO qmi;
		qmi.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\subdivideMesh.obj");
#endif
		/*exit(0);*/

		qm = newQm;
		bc = newBc;
		fcc = newFc;
		ExtractSheetsAndChords(qm, bc, sheets_, chords_);
#if !USE_FIELD
		RankingSheetsAndChords(sheets_, chords_, candidates);
#else
		RankingSheetsAndChordsByFields(qm, bc, sheets_, chords_, candidates);
#endif

		return true;
	}

	void Pipeline::FindSubdivideSheet(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &startEVec, std::vector<uint32_t> &fNumVec, double aveLength)
	{
		startEVec.clear();
		fNumVec.clear();
		std::vector<bool> eFlag(qm.Es_.size(), false);	//标记一条边有没有找过
		bool *fFlag = new bool[qm.Fs_.size()];
		std::memset(fFlag, 0, qm.Fs_.size() * sizeof(bool));
		bool isOpenSheet = true;

		std::function<void(uint32_t, uint32_t, uint32_t &)> findOppsiteE = [&](uint32_t currF, uint32_t currE, uint32_t &oppE)
		{
			std::vector<uint32_t> &evs0 = qm.Es_[currE].vs;
			std::vector<uint32_t> &fes = qm.Fs_[currF].es;
			for (int i = 0; i < 4; ++i)
			{
				std::vector<uint32_t> &evs1 = qm.Es_[fes[i]].vs;
				if (evs0[0] != evs1[0] && evs0[0] != evs1[1] && evs0[1] != evs1[0] && evs0[1] != evs1[1])
				{
					oppE = fes[i];
					break;
				}
			}
		};
		std::function<bool(uint32_t preF, uint32_t preE, uint32_t &nextF, uint32_t &nextE)> findNextEF = [&](uint32_t preF, uint32_t preE, uint32_t &nextF, uint32_t &nextE) -> bool
		{
			if (qm.Es_[preE].boundary)
				return false;
			
			std::vector<uint32_t> &efs = qm.Es_[preE].neighbor_fs;
			if (efs[0] == preF)
				nextF = efs[1];
			else
				nextF = efs[0];

			findOppsiteE(nextF, preE, nextE);

			if (/*eFlag[nextE] ||*/ fFlag[nextF])
			{
				eFlag[nextE] = true;
				isOpenSheet = false;
				return false;
			}

			return true;
		};


		typedef std::tuple<double, uint32_t, uint32_t> TempTuple3;
		std::vector<TempTuple3> tempT3;
		for (int i = 0; i < qm.Es_.size(); ++i)
		{
			if (eFlag[i])
				continue;
			std::vector<uint32_t> &evs = qm.Es_[i].vs;
			double edgeLength = (qm.V_[evs[0]] - qm.V_[evs[1]]).norm();
			if (edgeLength < aveLength * 2)
			{
				eFlag[i] = true;
				continue;
			}

			isOpenSheet = true;
			std::memset(fFlag, 0, qm.Fs_.size() * sizeof(bool));
			uint32_t realStartE = (uint32_t)-1, finalFNum = 0;

			uint32_t startE = i, startF0 = qm.Es_[startE].neighbor_fs[0], startF1 = (uint32_t)-1;
			uint32_t preE = (uint32_t)-1, nextE = preE, preF = preE, nextF = preE;
			eFlag[startE] = true;
			fFlag[startF0] = true;
			if (!qm.Es_[startE].boundary)
			{
				startF1 = qm.Es_[startE].neighbor_fs[1];
				fFlag[startF1] = true;
			}

			findOppsiteE(startF0, startE, preE);
			preF = startF0;
			++finalFNum;

			if (eFlag[preE])
				continue;
			eFlag[preE] = true;

			while (findNextEF(preF, preE, nextF, nextE))
			{
				eFlag[nextE] = true; fFlag[nextF] = true;
				++finalFNum;
				preE = nextE; preF = nextF;
			}
			if (!isOpenSheet)
				continue;

			realStartE = preE;
			if (startF1 == (uint32_t)-1)
			{
				tempT3.emplace_back(TempTuple3(edgeLength, realStartE, finalFNum));
				continue;
			}

			findOppsiteE(startF1, startE, preE);
			preF = startF1;
			++finalFNum;

			if (eFlag[preE])
				continue;
			eFlag[preE] = true;

			while (findNextEF(preF, preE, nextF, nextE))
			{
				eFlag[nextE] = true; fFlag[nextF] = true;
				++finalFNum;
				preE = nextE; preF = nextF;
			}
			if (!isOpenSheet)
				continue;

			tempT3.emplace_back(TempTuple3(edgeLength, realStartE, finalFNum));
		}

		std::sort(tempT3.begin(), tempT3.end());
		std::reverse(tempT3.begin(), tempT3.end());
		for (int i = 0; i < tempT3.size(); ++i)
		{
			startEVec.emplace_back(std::get<1>(tempT3[i]));
			fNumVec.emplace_back(std::get<2>(tempT3[i]));
		}

		delete[] fFlag;
	}

	void Pipeline::SubdivideSheetByEdgeLength(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, FeatureConstraints &fcc, double aveLength)
	{
		/*if (qm.Fs_.size() > initialQuadNum_ * subdivideThre_)
			return;*/

		std::vector<uint32_t> startEVec, fNumVec;
//#if !USE_FIELD
//		FindLongestOpenSheet(qm, startEVec, fNumVec);
//#else
//		FindLongestOpenSheet(qm, startEVec, fNumVec);
//		//FindNearestOpenSheet(qm, startEVec, fNumVec);
//#endif
		FindSubdivideSheet(qm, startEVec, fNumVec, aveLength);
		uint32_t currNum = 0;

		while (currNum < startEVec.size())
		{
			/*if (qm.Fs_.size() > initialQuadNum_ * stopSubThre_)
			{
				ExtractSheetsAndChords(qm, bc, sheets_, chords_);
#if !USE_FIELD
				RankingSheetsAndChords(sheets_, chords_, candidates);
#else
				RankingSheetsAndChordsByFields(qm, sheets_, chords_, candidates_);
#endif
				return;
			}*/

#if NEW_SUBDIVIDE
			if (SubdivideSheetFixV2(qm, bc, fcc, startEVec[currNum], fNumVec[currNum]))
#else
			if (SubdivideSheet(qm, bc, fcc, startEVec[currNum], fNumVec[currNum]))
#endif
			{
				currNum = 0;
//#if !USE_FIELD
//				FindLongestOpenSheet(qm, startEVec, fNumVec);
//#else
//				FindLongestOpenSheet(qm, startEVec, fNumVec);
//				//FindNearestOpenSheet(qm, startEVec, fNumVec);
//#endif
				FindSubdivideSheet(qm, startEVec, fNumVec, aveLength);
				continue;
			}
			else
			{
				++currNum;
			}
		}

		/*ExtractSheetsAndChords(qm, bc, sheets_, chords_);
#if !USE_FIELD
		RankingSheetsAndChords(sheets_, chords_, candidates);
#else
		RankingSheetsAndChordsByFields(qm, sheets_, chords_, candidates_);
#endif*/
	}

	bool Pipeline::SubdivideSheetFixV2(QuadMesh &qm, BaseComplex &bc, FeatureConstraints &fcc, uint32_t startE, uint32_t growNum)
	{
		subdivideRings_ = minusValue / growNum + 1;
		typedef std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> TupleUInt32_4;
		FeatureConstraints newFc = fcc;
		QuadMesh newQm;
		newQm.V_ = qm.V_;
		std::vector<TupleUInt32_4> finalFVs;
		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
			finalFVs.emplace_back(TupleUInt32_4(fvs[0], fvs[1], fvs[2], fvs[3]));
		}

		std::function<bool(uint32_t, uint32_t, uint32_t&, Eigen::Vector2d&, Eigen::Vector2d&)> vNewFcInfo = [&](uint32_t lVId, uint32_t rVId, uint32_t& curveId, Eigen::Vector2d& oriL, Eigen::Vector2d& axaL)->bool
		{
			int lVType = fcc.V_types[lVId], rVType = fcc.V_types[rVId];
			if (lVType < 0 && rVType < 0)
			{
				if (lVType == -1 && rVType == -1)
				{
					Eigen::Vector2d &currP0 = qm.V_[lVId], &currP1 = qm.V_[rVId];
					uint32_t corner0 = (uint32_t)-1, corner1 = corner0;

					FindSourceCorner(currP0, corner0);
					FindSourceCorner(currP1, corner1);

					std::tuple<int, int> &cp0 = mf.cornerCurves[corner0], &cp1 = mf.cornerCurves[corner1];
					if (std::get<0>(cp0) == std::get<0>(cp1) || std::get<0>(cp0) == std::get<1>(cp1))
						curveId = std::get<0>(cp0);
					else if (std::get<1>(cp0) == std::get<0>(cp1) || std::get<1>(cp0) == std::get<1>(cp1))
						curveId = std::get<1>(cp0);
					else
						return false;
				}
				else
				{
					//std::cout << "Can't add point between two corners! " << std::endl;
					return false;
				}
			}
			else
			{
				if (lVType >= 0)
					curveId = lVType;
				else
					curveId = rVType;
			}

			const Eigen::Vector2d &lVPos = qm.V_[lVId];
			const Eigen::Vector2d &rVPos = qm.V_[rVId];
			oriL = (lVPos + rVPos) * 0.5;
			axaL = (lVPos - rVPos).normalized();
			return true;
		};
		uint32_t startVId = qm.Vs_.size();
		std::vector<Eigen::Vector2d> newVPoss;
		std::vector<uint32_t> newVIds, fIds, leftVIds, rightVIds;

		GrowSplitSheet(startE, growNum, startVId, newVPoss, newVIds, fIds, leftVIds, rightVIds, qm);
		newQm.V_.insert(newQm.V_.end(), newVPoss.begin(), newVPoss.end());

		std::vector<uint32_t> newFIds;
		newFIds.reserve(fIds.size());
		for (int j = 0; j < fIds.size(); ++j)
		{
			TupleUInt32_4 &currFvs = finalFVs[fIds[j]];
			std::get<0>(currFvs) = leftVIds[j];
			std::get<1>(currFvs) = newVIds[j];
			std::get<2>(currFvs) = newVIds[j + 1];
			std::get<3>(currFvs) = leftVIds[j + 1];

			newFIds.emplace_back(finalFVs.size());
			finalFVs.emplace_back(TupleUInt32_4(newVIds[j], rightVIds[j], rightVIds[j + 1], newVIds[j + 1]));
		}

		//Build new fc
		uint32_t curveId = (uint32_t)-1;
		Eigen::Vector2d oriL, axaL;
		if (!vNewFcInfo(leftVIds[0], rightVIds[0], curveId, oriL, axaL))
			return false;
		newFc.V_types.emplace_back(curveId);
		newFc.ids_L.emplace_back(newVIds[0]);
		newFc.on_which_L.emplace_back(curveId);
		newFc.axa_L.emplace_back(axaL);
		newFc.origin_L.emplace_back(oriL);

		for (int j = 1; j < newVIds.size() - 1; ++j)
		{
			newFc.V_types.emplace_back(-2);
		}


		if (!vNewFcInfo(leftVIds[leftVIds.size() - 1], rightVIds[rightVIds.size() - 1], curveId, oriL, axaL))
			return false;
		newFc.V_types.emplace_back(curveId);
		newFc.ids_L.emplace_back(newVIds[newVIds.size() - 1]);
		newFc.on_which_L.emplace_back(curveId);
		newFc.axa_L.emplace_back(axaL);
		newFc.origin_L.emplace_back(oriL);


		std::function<void(QuadMesh &, std::vector<TupleUInt32_4> &)> buildQmFromTuple4 = [&](QuadMesh &nQm, std::vector<TupleUInt32_4> &fFvs)
		{
			QuadVertex qv;
			for (int i = 0; i < nQm.V_.size(); ++i)
			{
				qv.id = i;
				qv.bId = (uint32_t)-1;
				nQm.Vs_.emplace_back(qv);
			}

			for (int i = 0; i < fFvs.size(); ++i)
			{
				nQm.Fs_.emplace_back(QuadFace());
				QuadFace &currF = nQm.Fs_[i];
				currF.id = i;
				currF.bId = (uint32_t)-1;

				uint32_t vs[4];
				vs[0] = std::get<0>(fFvs[i]);
				vs[1] = std::get<1>(fFvs[i]);
				vs[2] = std::get<2>(fFvs[i]);
				vs[3] = std::get<3>(fFvs[i]);
				for (int j = 0; j < 4; ++j)
				{
					currF.vs.emplace_back(vs[j]);
					nQm.Vs_[vs[j]].neighbor_fs.emplace_back(i);
				}
			}
		};

		buildQmFromTuple4(newQm, finalFVs);
		if (!newQm.BuildConnectivity())
			return false;
		if (!TopologyCheck(newQm))
			return false;
		/*if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			std::cout << "Inverse Jacobi! " << std::endl;
			return false;
		}*/
#ifdef OUTPUT_MID_MESH
		QuadMeshIO qmi;
		qmi.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\BeforeSubdivideMesh.obj");
#endif

#pragma region Local_Optimize
		OptimizeInfo newOi;
		std::vector<uint32_t> ltgIdsL;
		BuildSubdivideOptimizeInfoV2(newQm, newFc, fIds, newFIds, &newOi, ltgIdsL);

#ifndef  USE_SQUARE_SOLVER

		AESolver aes(&newOi, kdTree_, fields_);
		if (aes.isInverse)
			return false;
		ComputeIdsLDis(&newOi);
#if OPTIMIZE_FIELDS
		//ExtraFieldEnergyAlone efea(newQm, &newOi, kdTree_, fields_, 1, 1, newQm.ComputeAverageEdgeLength());
#endif
		for (int i = 0; i < AES_ITER; ++i)
		{
			//			aes.efe_.ExtractFields();
			aes.Optimize();
			aes.SetValueBack(&newOi);
			ProjectToSouceLine(&newOi, mf);
#if USE_3To0
			AdjustBoundaryPosInOptimization(newQm, &newOi);
#endif
			aes.SetValueFront(&newOi);
			if (aes.isInverse)
				return false;

#if OPTIMIZE_FIELDS
			//if (i != AES_ITER - 1)
			//{
				/*efea.Initialize();
				efea.SetValueFront(&newOi);
				efea.Optimize();
				efea.SetValueBack(&newOi);*/
				//}
#endif
		}
#else 

		AESolverSquare aes(&newOi);
		if (aes.isInverse)
			return false;
		for (int i = 0; i < AES_ITER; ++i)
		{
			aes.Optimize();
			aes.SetValueBack(&newOi);
			ProjectToSouceLine(&newOi, mf);
			aes.SetValueFront(&newOi);
			if (aes.isInverse)
				return false;
		}
#endif 

		for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
		{
			uint32_t newVId = newOi.localToGlobalVVec[i];
			newQm.V_[newVId] = newOi.localVPos[i];
		}
		if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			std::cout << "Subdivide Inverse Jacobi! " << std::endl;
			return false;
		}
		if (!ComputeHausdorff(backupQm_, newQm, hausdorffRadio, hausdorffThre))
			return false;
		for (int i = 0; i < ltgIdsL.size(); ++i)
		{
			uint32_t idsId = ltgIdsL[i];
			newFc.axa_L[idsId] = newOi.localFc.axa_L[i];
			newFc.origin_L[idsId] = newOi.localFc.origin_L[i];
		}
#pragma endregion

#ifdef OUTPUT_MID_MESH
		qmi.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\AfterSubdivideMesh.obj");
#endif

		BaseComplex newBc;
		if (!isCornerBCV)
		{
#if !CORNER_BC
			newBc.ExtractBaseComplex(&newQm);
#else
			newBc.ExtractBaseComplex(&newQm, FEATURE_THRESHOLD);
#endif	
		}
		else
		{
			corners_.clear();
			for (int i = 0; i < newFc.ids_C.size(); ++i)
			{
				corners_.emplace_back(newFc.ids_C[i]);
			}
			newBc.ExtractBaseComplex(&newQm, corners_);
		}
		if (newBc.isError_)
			return false;

		/*if (!TopologyCheck(newQm, newBc))
		{
			std::cout << "Didn't pass topology check! " << std::endl;
			return false;
		}*/

		/*exit(0);*/

		qm = newQm;
		bc = newBc;
		fcc = newFc;
		//ExtractSheetsAndChords(qm, bc, sheets_, chords_);
		//RankingSheetsAndChords(sheets_, chords_, candidates);

		return true;
	}

	bool Pipeline::SubdivideSheet(QuadMesh &qm, BaseComplex &bc, FeatureConstraints &fcc, uint32_t startE, uint32_t growNum)
	{
		subdivideRings_ = minusValue / growNum + 1;
		typedef std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> TupleUInt32_4;
		FeatureConstraints newFc = fcc;
		QuadMesh newQm;
		newQm.V_ = qm.V_;
		std::vector<TupleUInt32_4> finalFVs;
		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
			finalFVs.emplace_back(TupleUInt32_4(fvs[0], fvs[1], fvs[2], fvs[3]));
		}

		std::function<bool(uint32_t, uint32_t, uint32_t&, Eigen::Vector2d&, Eigen::Vector2d&)> vNewFcInfo = [&](uint32_t lVId, uint32_t rVId, uint32_t& curveId, Eigen::Vector2d& oriL, Eigen::Vector2d& axaL)->bool
		{
			int lVType = fcc.V_types[lVId], rVType = fcc.V_types[rVId];
			if (lVType < 0 && rVType < 0)
			{
				if (lVType == -1 && rVType == -1)
				{
					Eigen::Vector2d &currP0 = qm.V_[lVId], &currP1 = qm.V_[rVId];
					uint32_t corner0 = (uint32_t)-1, corner1 = corner0;

					FindSourceCorner(currP0, corner0);
					FindSourceCorner(currP1, corner1);

					std::tuple<int, int> &cp0 = mf.cornerCurves[corner0], &cp1= mf.cornerCurves[corner1];
					if (std::get<0>(cp0) == std::get<0>(cp1) || std::get<0>(cp0) == std::get<1>(cp1))
						curveId = std::get<0>(cp0);
					else if (std::get<1>(cp0) == std::get<0>(cp1) || std::get<1>(cp0) == std::get<1>(cp1))
						curveId = std::get<1>(cp0);
					else
						return false;
				}
				else
				{
					//std::cout << "Can't add point between two corners! " << std::endl;
					return false;
				}
			}
			else
			{
				if (lVType >= 0)
					curveId = lVType;
				else
					curveId = rVType;
			}

			const Eigen::Vector2d &lVPos = qm.V_[lVId];
			const Eigen::Vector2d &rVPos = qm.V_[rVId];
			oriL = (lVPos + rVPos) * 0.5;
			axaL = (lVPos - rVPos).normalized();
			return true;
		};
		uint32_t startVId = qm.Vs_.size();
		std::vector<Eigen::Vector2d> newVPoss;
		std::vector<uint32_t> newVIds, fIds, leftVIds, rightVIds;
		
		GrowSplitSheet(startE, growNum, startVId, newVPoss, newVIds, fIds, leftVIds, rightVIds, qm);
		newQm.V_.insert(newQm.V_.end(), newVPoss.begin(), newVPoss.end());

		std::vector<uint32_t> newFIds;
		newFIds.reserve(fIds.size());
		for (int j = 0; j < fIds.size(); ++j)
		{
			TupleUInt32_4 &currFvs = finalFVs[fIds[j]];
			std::get<0>(currFvs) = leftVIds[j];
			std::get<1>(currFvs) = newVIds[j];
			std::get<2>(currFvs) = newVIds[j + 1];
			std::get<3>(currFvs) = leftVIds[j + 1];

			newFIds.emplace_back(finalFVs.size());
			finalFVs.emplace_back(TupleUInt32_4(newVIds[j], rightVIds[j], rightVIds[j + 1], newVIds[j + 1]));
		}

		//Build new fc
		uint32_t curveId = (uint32_t)-1;
		Eigen::Vector2d oriL, axaL;
		if (!vNewFcInfo(leftVIds[0], rightVIds[0], curveId, oriL, axaL))
			return false;
		newFc.V_types.emplace_back(curveId);
		newFc.ids_L.emplace_back(newVIds[0]);
		newFc.on_which_L.emplace_back(curveId);
		newFc.axa_L.emplace_back(axaL);
		newFc.origin_L.emplace_back(oriL);

		for (int j = 1; j < newVIds.size() - 1; ++j)
		{
			newFc.V_types.emplace_back(-2);
		}

		
		if (!vNewFcInfo(leftVIds[leftVIds.size() - 1], rightVIds[rightVIds.size() - 1], curveId, oriL, axaL))
			return false;
		newFc.V_types.emplace_back(curveId);
		newFc.ids_L.emplace_back(newVIds[newVIds.size() - 1]);
		newFc.on_which_L.emplace_back(curveId);
		newFc.axa_L.emplace_back(axaL);
		newFc.origin_L.emplace_back(oriL);
		

		std::function<void(QuadMesh &, std::vector<TupleUInt32_4> &)> buildQmFromTuple4 = [&](QuadMesh &nQm, std::vector<TupleUInt32_4> &fFvs)
		{
			QuadVertex qv;
			for (int i = 0; i < nQm.V_.size(); ++i)
			{
				qv.id = i;
				qv.bId = (uint32_t)-1;
				nQm.Vs_.emplace_back(qv);
			}

			for (int i = 0; i < fFvs.size(); ++i)
			{
				nQm.Fs_.emplace_back(QuadFace());
				QuadFace &currF = nQm.Fs_[i];
				currF.id = i;
				currF.bId = (uint32_t)-1;

				uint32_t vs[4];
				vs[0] = std::get<0>(fFvs[i]);
				vs[1] = std::get<1>(fFvs[i]);
				vs[2] = std::get<2>(fFvs[i]);
				vs[3] = std::get<3>(fFvs[i]);
				for (int j = 0; j < 4; ++j)
				{
					currF.vs.emplace_back(vs[j]);
					nQm.Vs_[vs[j]].neighbor_fs.emplace_back(i);
				}
			}
		};

		buildQmFromTuple4(newQm, finalFVs);
		if (!newQm.BuildConnectivity())
			return false;
		if (!TopologyCheck(newQm))
			return false;
		/*if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			std::cout << "Inverse Jacobi! " << std::endl;
			return false;
		}*/
#ifdef OUTPUT_MID_MESH
		QuadMeshIO qmi;
		qmi.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\BeforeSubdivideMesh.obj");
#endif

#pragma region Local_Optimize
		OptimizeInfo newOi;
		std::vector<uint32_t> ltgIdsL;
		BuildSubdivideOptimizeInfo(newQm, newFc, fIds, newFIds, &newOi, ltgIdsL);

#ifndef  USE_SQUARE_SOLVER

		AESolver aes(&newOi, kdTree_, fields_);
		if (aes.isInverse)
			return false;
		ComputeIdsLDis(&newOi);
#if OPTIMIZE_FIELDS
		//ExtraFieldEnergyAlone efea(newQm, &newOi, kdTree_, fields_, 1, 1, newQm.ComputeAverageEdgeLength());
#endif
		for (int i = 0; i < AES_ITER; ++i)
		{
//			aes.efe_.ExtractFields();
			aes.Optimize();
			aes.SetValueBack(&newOi);
			ProjectToSouceLine(&newOi, mf);
#if USE_3To0
			AdjustBoundaryPosInOptimization(newQm, &newOi);
#endif
			aes.SetValueFront(&newOi);
			if (aes.isInverse)
				return false;

#if OPTIMIZE_FIELDS
			//if (i != AES_ITER - 1)
			//{
				/*efea.Initialize();
				efea.SetValueFront(&newOi);
				efea.Optimize();
				efea.SetValueBack(&newOi);*/
			//}
#endif
		}
#else 

		AESolverSquare aes(&newOi);
		if (aes.isInverse)
			return false;
		for (int i = 0; i < AES_ITER; ++i)
		{
			aes.Optimize();
			aes.SetValueBack(&newOi);
			ProjectToSouceLine(&newOi, mf);
			aes.SetValueFront(&newOi);
			if (aes.isInverse)
				return false;
		}
#endif 

		for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
		{
			uint32_t newVId = newOi.localToGlobalVVec[i];
			newQm.V_[newVId] = newOi.localVPos[i];
		}
		if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			std::cout << "Subdivide Inverse Jacobi! " << std::endl;
			return false;
		}
		if (!ComputeHausdorff(backupQm_, newQm, hausdorffRadio, hausdorffThre))
			return false;
		for (int i = 0; i < ltgIdsL.size(); ++i)
		{
			uint32_t idsId = ltgIdsL[i];
			newFc.axa_L[idsId] = newOi.localFc.axa_L[i];
			newFc.origin_L[idsId] = newOi.localFc.origin_L[i];
		}
#pragma endregion

#ifdef OUTPUT_MID_MESH
		qmi.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\AfterSubdivideMesh.obj");
#endif

		BaseComplex newBc;
		if (!isCornerBCV)
		{
#if !CORNER_BC
			newBc.ExtractBaseComplex(&newQm);
#else
			newBc.ExtractBaseComplex(&newQm, FEATURE_THRESHOLD);
#endif	
		}
		else
		{
			corners_.clear();
			for (int i = 0; i < newFc.ids_C.size(); ++i)
			{
				corners_.emplace_back(newFc.ids_C[i]);
			}
			newBc.ExtractBaseComplex(&newQm, corners_);
		}
		if (newBc.isError_)
			return false;

		/*if (!TopologyCheck(newQm, newBc))
		{
			std::cout << "Didn't pass topology check! " << std::endl;
			return false;
		}*/
		
		/*exit(0);*/

		qm = newQm;
		bc = newBc;
		fcc = newFc;
		//ExtractSheetsAndChords(qm, bc, sheets_, chords_);
		//RankingSheetsAndChords(sheets_, chords_, candidates);

		return true;
	}

	void Pipeline::GrowSplitSheet(uint32_t startE, uint32_t growNum, uint32_t &startVId, std::vector<Eigen::Vector2d> &newVPoss, std::vector<uint32_t> &newVIds, std::vector<uint32_t> &fIds, std::vector<uint32_t> &leftVIds, std::vector<uint32_t> &rightVIds, QuadMesh &qm)
	{
		std::vector<int> vLRFlag(qm.Vs_.size(), -1); //用来标记一个点是否属于左边或者右边， 如果左边则为0， 如果右边则为1，否则为-1
		//std::vector<Eigen::Vector2d>().swap(newVPoss);
		newVPoss.clear();
		newVPoss.resize(growNum + 1);		//这里是resize
		//std::vector<uint32_t>().swap(newVIds);
		newVIds.clear();
		newVIds.reserve(growNum + 1);
		//std::vector<uint32_t>().swap(fIds);
		fIds.clear();
		fIds.reserve(growNum);
		//std::vector<uint32_t>().swap(leftVIds);
		leftVIds.clear();
		leftVIds.resize(growNum + 1);
		//std::vector<uint32_t>().swap(rightVIds);
		rightVIds.clear();
		rightVIds.resize(growNum + 1);

		//startF
		std::vector<uint32_t> &startFs = qm.Es_[startE].neighbor_fs;
		if (startFs.size() != 1)
		{
			std::cout << "Error: Pipeline::GrowSplitSheet! " << std::endl;
			exit(-28);
		}
		uint32_t startF = startFs[0];

		//找第一个left 和 right。
		uint32_t leftV = qm.Es_[startE].vs[0], rightV = qm.Es_[startE].vs[1];
		std::function<void(uint32_t &, uint32_t &, uint32_t)> findLeftRightV = [&](uint32_t &lv, uint32_t &rv, uint32_t currF)
		{
			uint32_t lvId = (uint32_t)-1, rvId = (uint32_t)-1;
			for (int i = 0; i < 4; ++i)
			{
				if (qm.Fs_[currF].vs[i] == lv)
					lvId = i;
				else if (qm.Fs_[currF].vs[i] == rv)
					rvId = i;
			}
			if (lvId > rvId)
			{
				if (!(lvId == 3 && rvId == 0))
					std::swap(lv, rv);
			}
			else
			{
				if (lvId == 0 && rvId == 3)
					std::swap(lv, rv);
			}
		};
		std::function<void(uint32_t, Eigen::Vector2d&)> computeMidPos = [&](uint32_t currE, Eigen::Vector2d& midP)
		{
			uint32_t v0 = qm.Es_[currE].vs[0], v1 = qm.Es_[currE].vs[1];
			midP = (qm.V_[v0] + qm.V_[v1]) * 0.5;
		};
		findLeftRightV(leftV, rightV, startF);
		leftVIds[0] = leftV;
		rightVIds[0] = rightV;
		computeMidPos(startE, newVPoss[0]);
		newVIds.emplace_back(startVId++);

		std::vector<bool> fFindFlag(qm.Fs_.size(), false); //标记面有没有被找过。
		std::function<void(uint32_t, uint32_t&, uint32_t&)> findNextEF = [&](uint32_t preE, uint32_t& nextE, uint32_t& nextF)
		{
			std::vector<uint32_t> &eFs = qm.Es_[preE].neighbor_fs;
			if (eFs.size() == 1)
			{
				nextF = eFs[0];
			}
			else
			{
				if (fFindFlag[eFs[0]])
					nextF = eFs[1];
				else
					nextF = eFs[0];
			}
			fFindFlag[nextF] = true;

			std::vector<uint32_t> &preEVs = qm.Es_[preE].vs;
			std::vector<uint32_t> &fEs = qm.Fs_[nextF].es;
			for (int i = 0; i < fEs.size(); ++i)	//找下一个E.
			{
				if (fEs[i] == preE)
					continue;

				std::vector<uint32_t> &currEVs = qm.Es_[fEs[i]].vs;
				if (currEVs[0] != preEVs[0] && currEVs[0] != preEVs[1] && currEVs[1] != preEVs[0] && currEVs[1] != preEVs[1])
				{
					nextE = fEs[i];
					break;
				}
			}
		};

		std::vector<uint32_t> fvs;
		std::function<void(uint32_t, uint32_t, uint32_t, uint32_t&, uint32_t&)> findNewLRV = [&](uint32_t nextF, uint32_t preLv, uint32_t preRv, uint32_t &nextLv, uint32_t &nextRv)
		{
			std::vector<uint32_t> &plVVs = qm.Vs_[preLv].neighbor_vs;
			fvs = qm.Fs_[nextF].vs;
			std::vector<uint32_t> tempVVec;
			std::sort(plVVs.begin(), plVVs.end());
			std::sort(fvs.begin(), fvs.end());
			std::set_intersection(plVVs.begin(), plVVs.end(), fvs.begin(), fvs.end(), std::inserter(tempVVec, tempVVec.begin()));
			if (tempVVec[0] == preRv)
				nextLv = tempVVec[1];
			else
				nextLv = tempVVec[0];

			std::vector<uint32_t> &prVVs = qm.Vs_[preRv].neighbor_vs;
			tempVVec.clear();
			std::sort(prVVs.begin(), prVVs.end());
			std::set_intersection(prVVs.begin(), prVVs.end(), fvs.begin(), fvs.end(), std::inserter(tempVVec, tempVVec.begin()));
			if (tempVVec[0] == preLv)
				nextRv = tempVVec[1];
			else
				nextRv = tempVVec[0];
		};

		uint32_t preE = startE, nextE = (uint32_t)-1, nextF = (uint32_t)-1;
		for (int i = 0; i < growNum; ++i)
		{
			findNextEF(preE, nextE, nextF);
			fIds.emplace_back(nextF);
			computeMidPos(nextE, newVPoss[i+1]);
			newVIds.emplace_back(startVId++);

			uint32_t lastLv = leftVIds[i];
			uint32_t lastRv = rightVIds[i];
			findNewLRV(nextF, lastLv, lastRv, leftVIds[i + 1], rightVIds[i + 1]);

			preE = nextE;
		}
	}

	void Pipeline::SortEdgeFromLongToShort(std::vector<uint32_t> &es, BaseDataStructure::QuadMesh &qm)
	{
		std::vector<std::tuple<double, uint32_t>> lengthIdVec;
		for (int i = 0; i < es.size(); ++i)
		{
			uint32_t ev0 = qm.Es_[es[i]].vs[0];
			uint32_t ev1 = qm.Es_[es[i]].vs[1];
			lengthIdVec.emplace_back(std::tuple<double, uint32_t>((qm.V_[ev0] - qm.V_[ev1]).norm(), es[i]));
		}
		std::sort(lengthIdVec.begin(), lengthIdVec.end());

		es.clear();
		for (int i = 0; i < lengthIdVec.size(); ++i)
		{
			es.emplace_back(std::get<1>(lengthIdVec[i]));
		}
		std::reverse(es.begin(), es.end());
	}

	void Pipeline::FindNearestOpenSheet(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &startEVec, std::vector<uint32_t> &fNumVec)
	{
		//std::vector<uint32_t>().swap(startEVec);
		//std::vector<uint32_t>().swap(fNumVec);
		startEVec.clear();
		fNumVec.clear();
		std::vector<bool> eFlag(qm.Es_.size(), false);
		std::vector<bool> fFlag(qm.Fs_.size(), false);
		std::function<bool(uint32_t, uint32_t&, uint32_t&)> findNextEF = [&](uint32_t preE, uint32_t &nextE, uint32_t &nextF) -> bool
		{

			std::vector<uint32_t> &eFs = qm.Es_[preE].neighbor_fs;
			if (eFs.size() == 1)
				nextF = eFs[0];
			else
			{
				if (fFlag[eFs[0]])
					nextF = eFs[1];
				else
					nextF = eFs[0];
			}

			std::vector<uint32_t> &preEVs = qm.Es_[preE].vs;
			for (int i = 0; i < 4; ++i)
			{
				uint32_t currE = qm.Fs_[nextF].es[i];
				std::vector<uint32_t> &currEVs = qm.Es_[currE].vs;
				if (currEVs[0] != preEVs[0] && currEVs[0] != preEVs[1] && currEVs[1] != preEVs[0] && currEVs[1] != preEVs[1])
				{
					nextE = currE;
					break;
				}
			}

			if (qm.Es_[preE].boundary)
				return false;

			return true;
		};

		std::vector<std::tuple<double, uint32_t, uint32_t>> tuple3Vec;
		tuple3Vec.reserve(qm.Es_.size());
		std::vector<double> lengthVec, fieldVec;
		lengthVec.reserve(qm.Es_.size());
		fieldVec.reserve(qm.Es_.size());
		for (int i = 0; i < qm.Es_.size(); ++i)
		{
			if (eFlag[i] || !qm.Es_[i].boundary)
				continue;

			std::fill(fFlag.begin(), fFlag.end(), false);
			std::vector<uint32_t> currFacesVec;

			uint32_t preE = i, nextE = (uint32_t)-1, nextF = (uint32_t)-1;
			std::vector<uint32_t> esVec;
			eFlag[i] = true;
			esVec.emplace_back(i);
			bool isSelf = false;

			while (findNextEF(preE, nextE, nextF) || preE == i)
			{
				if (fFlag[nextF])
				{
					isSelf = true;
					break;	//说明碰到了之前碰到的四边形，也就是有自交，这种不要。
				}

				fFlag[nextF] = true;
				eFlag[nextE] = true;
				esVec.emplace_back(nextE);
				currFacesVec.emplace_back(nextF);
				preE = nextE;
			}
			if (isSelf)
				continue;

			double sumAngle = 0;
			for (int j = 0; j < currFacesVec.size(); ++j)
			{
				sumAngle += ComputeQuadFieldWeight(qm, currFacesVec[j]);
			}
			sumAngle /= currFacesVec.size();

			double aveLength = 0;
			for (int j = 0; j < esVec.size(); ++j)
			{
				uint32_t vId0 = qm.Es_[esVec[j]].vs[0];
				uint32_t vId1 = qm.Es_[esVec[j]].vs[1];
				aveLength += (qm.V_[vId0] - qm.V_[vId1]).norm();
			}
			aveLength /= esVec.size();

			lengthVec.emplace_back(aveLength); fieldVec.emplace_back(sumAngle);
			tuple3Vec.emplace_back(std::tuple<double, uint32_t, uint32_t>(sumAngle, i, esVec.size() - 1));
		}

		double maxL = *std::max_element(lengthVec.begin(), lengthVec.end());
		double minL = *std::min_element(lengthVec.begin(), lengthVec.end());
		double maxF = *std::max_element(fieldVec.begin(), fieldVec.end());
		double minF = *std::min_element(fieldVec.begin(), fieldVec.end());
		double lDis = maxL - minL;
		if (std::abs(lDis)<1E-5)
			lDis = 1;
		double fDis = maxF - minF;
		if (std::abs(fDis)<1E-5)
			fDis = 1;
		for (int j = 0; j < lengthVec.size(); ++j)
		{
			//lengthVec[j] = (1-(lengthVec[j] - minL) / lDis) * 1000;
			std::get<0>(tuple3Vec[j]) =  -lengthVec[j] + fieldVec[j] * 0.0003;
		}

		std::sort(tuple3Vec.begin(), tuple3Vec.end());
		for (int j = 0; j < tuple3Vec.size(); ++j)
		{
			startEVec.emplace_back(std::get<1>(tuple3Vec[j]));
			fNumVec.emplace_back(std::get<2>(tuple3Vec[j]));
		}
	}

	void Pipeline::FindLongestOpenSheet(QuadMesh &qm, std::vector<uint32_t> &startEVec, std::vector<uint32_t> &fNumVec)
	{
		//std::vector<uint32_t>().swap(startEVec);
		//std::vector<uint32_t>().swap(fNumVec);
		startEVec.clear();
		fNumVec.clear();
		std::vector<bool> eFlag(qm.Es_.size(), false);
		std::vector<bool> fFlag(qm.Fs_.size(), false);
		std::function<bool(uint32_t, uint32_t&, uint32_t&)> findNextEF = [&](uint32_t preE, uint32_t &nextE, uint32_t &nextF) -> bool
		{

			std::vector<uint32_t> &eFs = qm.Es_[preE].neighbor_fs;
			if (eFs.size() == 1)
				nextF = eFs[0];
			else
			{
				if (fFlag[eFs[0]])
					nextF = eFs[1];
				else
					nextF = eFs[0];
			}

			std::vector<uint32_t> &preEVs = qm.Es_[preE].vs;
			for (int i = 0; i < 4; ++i)
			{
				uint32_t currE = qm.Fs_[nextF].es[i];
				std::vector<uint32_t> &currEVs = qm.Es_[currE].vs;
				if (currEVs[0] != preEVs[0] && currEVs[0] != preEVs[1] && currEVs[1] != preEVs[0] && currEVs[1] != preEVs[1])
				{
					nextE = currE;
					break;
				}
			}

			if (qm.Es_[preE].boundary)
				return false;

			return true;
		};

		std::vector<std::tuple<double, uint32_t, uint32_t>> tuple3Vec;
		for (int i = 0; i < qm.Es_.size(); ++i)
		{
			if (eFlag[i] || !qm.Es_[i].boundary)
				continue;

			std::fill(fFlag.begin(), fFlag.end(), false);

			uint32_t preE = i, nextE = (uint32_t)-1, nextF = (uint32_t)-1;
			std::vector<uint32_t> esVec;
			eFlag[i] = true;
			esVec.emplace_back(i);
			bool isSelf = false;
			while (findNextEF(preE, nextE, nextF) || preE == i)
			{
				if (fFlag[nextF])
				{
					isSelf = true;
					break;	//说明碰到了之前碰到的四边形，也就是有自交，这种不要。
				}

				fFlag[nextF] = true;
				eFlag[nextE] = true;
				esVec.emplace_back(nextE);
				preE = nextE;
			}
			if (isSelf)
				continue;
			double aveLength = 0;
			for (int j = 0; j < esVec.size(); ++j)
			{
				uint32_t vId0 = qm.Es_[esVec[j]].vs[0];
				uint32_t vId1 = qm.Es_[esVec[j]].vs[1];
				aveLength += (qm.V_[vId0] - qm.V_[vId1]).norm();
			}
			aveLength /= esVec.size();
			tuple3Vec.emplace_back(std::tuple<double, uint32_t, uint32_t>(aveLength, i, esVec.size() - 1));
		}

		std::sort(tuple3Vec.begin(), tuple3Vec.end());
		std::reverse(tuple3Vec.begin(), tuple3Vec.end());
		for (int j = 0; j < tuple3Vec.size(); ++j)
		{
			startEVec.emplace_back(std::get<1>(tuple3Vec[j]));
			fNumVec.emplace_back(std::get<2>(tuple3Vec[j]));
		}
	}

	void Pipeline::BuildSubdivideOptimizeInfoV2(QuadMesh &qm, FeatureConstraints &fc, const std::vector<uint32_t> &oldFs, const std::vector<uint32_t> &newFs, OptimizeInfo *oi, std::vector<uint32_t> &localToGlobalIdsL)
	{
		bool *fFlag = new bool[qm.Fs_.size()];
		std::memset(fFlag, 0, qm.Fs_.size() * sizeof(bool));
		bool *vFlag = new bool[qm.Vs_.size()];
		std::memset(vFlag, 0, qm.Vs_.size() * sizeof(bool));
		bool *badFFlag = new bool[qm.Fs_.size()];
		std::memset(badFFlag, 0, qm.Fs_.size() * sizeof(bool));

		//找到所有的badF
		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			if (!qm.Fs_[i].boundary)
				continue;

			std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
			for (int j = 0; j < fvs.size(); ++j)
			{
				if (qm.Vs_[fvs[j]].boundary && qm.Vs_[fvs[j]].neighbor_fs.size() == 1)
				{
					badFFlag[i] = true;
					break;
				}
			}
		}

		double targetArea = initialArea_ / qm.Fs_.size();
		std::vector<uint32_t> localFs;
		//std::vector<uint32_t> localVs;
		//grow region
		for (int i = 0; i < oldFs.size(); ++i)
		{
			if (!fFlag[oldFs[i]] && !badFFlag[oldFs[i]])
			{
				fFlag[oldFs[i]] = true;
				localFs.emplace_back(oldFs[i]);
			}

			if (!fFlag[newFs[i]] && !badFFlag[newFs[i]])
			{
				fFlag[newFs[i]] = true;
				localFs.emplace_back(newFs[i]);
			}
		}

		std::vector<uint32_t> currFPool(localFs);
		std::vector<uint32_t> tempFPool;
		tempFPool.reserve(qm.Fs_.size());
		for (int i = 0; i < 4; ++i)
		{
			tempFPool.clear();
			for (int j = 0; j < currFPool.size(); ++j)
			{
				std::vector<uint32_t> &ffs = qm.Fs_[currFPool[j]].neighbor_fs;
				for (int k = 0; k < ffs.size(); ++k)
				{
					if (!fFlag[ffs[k]] && !badFFlag[ffs[k]])
					{
						fFlag[ffs[k]] = true;
						tempFPool.emplace_back(ffs[k]);
						localFs.emplace_back(ffs[k]);
					}
				}
			}
			tempFPool.swap(currFPool);
		}

		std::vector<uint32_t> currVPool;
		currVPool.reserve(qm.Vs_.size());
		for (int i = 0; i < localFs.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[localFs[i]].vs;
			for (int j = 0; j < fvs.size(); ++j)
			{
				if (!vFlag[fvs[j]])
				{
					vFlag[fvs[j]] = true;
					currVPool.emplace_back(fvs[j]);
				}
			}
		}

		//在找到所有局部点和局部面之后，开始build optimize info
		uint32_t newVId = 0;
		oi->globalToLocalVVec.resize(qm.Vs_.size(), (uint32_t)-1);
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if (vFlag[i])
			{
				oi->localVPos.emplace_back(qm.V_[i]);
				oi->globalToLocalVVec[i] = newVId++;
				oi->localToGlobalVVec.emplace_back(i);
			}
		}

		for (int i = 0; i < currVPool.size(); ++i)
		{
			std::vector<uint32_t> &vfs = qm.Vs_[currVPool[i]].neighbor_fs;
			bool isFix = false;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (!fFlag[vfs[j]])
				{
					isFix = true;
					break;
				}
			}
			if (isFix)
			{
				oi->localVsNotMove.emplace_back(oi->globalToLocalVVec[currVPool[i]]);
				oi->localVsNotMovePos.emplace_back(qm.V_[currVPool[i]]);
			}
		}

		std::vector<std::vector<Eigen::Vector2d>> tempVPosVec;
		Eigen::Vector3i tempVi;
		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			if (fFlag[i])
			{
				ComputeStandardTriangles(qm, i, targetArea, tempVPosVec);
				std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
				for (int j = 0; j < 4; ++j)
				{
					tempVi[0] = oi->globalToLocalVVec[fvs[quadTriTable[j][0]]];
					tempVi[1] = oi->globalToLocalVVec[fvs[quadTriTable[j][1]]];
					tempVi[2] = oi->globalToLocalVVec[fvs[quadTriTable[j][2]]];
					oi->localTriIds.emplace_back(tempVi);

					oi->localStanTriPos.emplace_back(tempVPosVec[j]);
				}
			}
		}

		FeatureConstraints &newFc = oi->localFc;
		newFc.V_types.reserve(oi->localToGlobalVVec.size());
		for (int i = 0; i < oi->localToGlobalVVec.size(); ++i)
		{
			uint32_t refeVId = oi->localToGlobalVVec[i];
			newFc.V_types.emplace_back(fc.V_types[refeVId]);
		}

		for (int i = 0; i < fc.ids_C.size(); ++i)
		{
			uint32_t newV = oi->globalToLocalVVec[fc.ids_C[i]];
			if (newV == (uint32_t)-1)
				continue;

			newFc.ids_C.emplace_back(newV);
			newFc.C.emplace_back(fc.C[i]);
		}

		localToGlobalIdsL.clear();
		for (int i = 0; i < fc.ids_L.size(); ++i)
		{
			uint32_t newV = oi->globalToLocalVVec[fc.ids_L[i]];
			if (newV == (uint32_t)-1)
				continue;

			localToGlobalIdsL.emplace_back(i);
			newFc.ids_L.emplace_back(newV);
			newFc.on_which_L.emplace_back(fc.on_which_L[i]);
			newFc.axa_L.emplace_back(fc.axa_L[i]);
			newFc.origin_L.emplace_back(fc.origin_L[i]);
		}

		delete[] fFlag;
		fFlag = NULL;
		delete[] vFlag;
		vFlag = NULL;
		delete[] badFFlag;
		badFFlag = NULL;
	}

	void Pipeline::BuildSubdivideOptimizeInfo(QuadMesh &qm, FeatureConstraints &fc, const std::vector<uint32_t> &oldFs, const std::vector<uint32_t> &newFs, OptimizeInfo *oi, std::vector<uint32_t> &localToGlobalIdsL)
	{
		std::vector<bool> fFlag(qm.Fs_.size(), false);
		std::vector<bool> vFlag(qm.Vs_.size(), false);

		std::vector<uint32_t> localFs;
		std::vector<double> fTargetArea(qm.Fs_.size(), -1);
		//std::vector<uint32_t> localVs;
		//grow region
		std::vector<uint32_t> tempFVec(2);
		for (int i = 0; i < oldFs.size(); ++i)
		{
			if (!fFlag[oldFs[i]])
			{
				fFlag[oldFs[i]] = true;
				localFs.emplace_back(oldFs[i]);
			}

			if (!fFlag[newFs[i]])
			{
				fFlag[newFs[i]] = true;
				localFs.emplace_back(newFs[i]);
			}
			tempFVec[0] = oldFs[i]; tempFVec[1] = newFs[i];
			double sumArea = ComputeFaceAreasSum(qm, tempFVec);
			fTargetArea[oldFs[i]] = sumArea; fTargetArea[newFs[i]] = sumArea;
		}

		std::vector<uint32_t> currVPool;
		for (int i = 0; i < localFs.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[localFs[i]].vs;
			for (int j = 0; j < 4; ++j)
			{
				if (!vFlag[fvs[j]])
				{
					vFlag[fvs[j]] = true;
					currVPool.emplace_back(fvs[j]);
				}
			}
		}

		//找它的4ring
		std::vector<uint32_t> tempVPool;
		for (int i = 0; i < 4; ++i)
		{
			tempVPool.clear();
			for (int j = 0; j < currVPool.size(); ++j)
			{
				/*std::vector<uint32_t> &vvs = qm.Vs_[currVPool[j]].neighbor_vs;
				for (int k = 0; k < vvs.size(); ++k)
				{
					if (!vFlag[vvs[k]])
					{
						vFlag[vvs[k]] = true;
						tempVPool.emplace_back(vvs[k]);
					}
				}*/
				std::vector<uint32_t> &vfs = qm.Vs_[currVPool[j]].neighbor_fs;
				for (int k = 0; k < vfs.size(); ++k)
				{
					std::vector<uint32_t> &vfvs = qm.Fs_[vfs[k]].vs;
					for (int m = 0; m < vfvs.size(); ++m)
					{
						if (!vFlag[vfvs[m]])
						{
							vFlag[vfvs[m]] = true;
							tempVPool.emplace_back(vfvs[m]);
						}
					}
					if (!fFlag[vfs[k]])
					{
						fFlag[vfs[k]] = true;
						localFs.emplace_back(vfs[k]);
						double sumArea = ComputeFaceArea(qm, vfs[k]);
						fTargetArea[vfs[k]] = sumArea;
					}
				}
			}
			currVPool.swap(tempVPool);
		}

		//在找到所有局部点和局部面之后，开始build optimize info
		uint32_t newVId = 0;
		oi->globalToLocalVVec.resize(qm.Vs_.size(), (uint32_t)-1);
		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			if (vFlag[i])
			{
				oi->localVPos.emplace_back(qm.V_[i]);
				oi->globalToLocalVVec[i] = newVId++;
				oi->localToGlobalVVec.emplace_back(i);
			}
		}

		for (int i = 0; i < currVPool.size(); ++i)
		{
			std::vector<uint32_t> &vfs = qm.Vs_[currVPool[i]].neighbor_fs;
			bool isFix = false;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (!fFlag[vfs[j]])
				{
					isFix = true;
					break;
				}
			}
			if (isFix)
			{
				oi->localVsNotMove.emplace_back(oi->globalToLocalVVec[currVPool[i]]);
				oi->localVsNotMovePos.emplace_back(qm.V_[currVPool[i]]);
			}
		}

		std::vector<std::vector<Eigen::Vector2d>> tempVPosVec;
		Eigen::Vector3i tempVi;
		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			if (fFlag[i])
			{
				ComputeStandardTriangles(qm, i, fTargetArea[i], tempVPosVec);
				std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
				for (int j = 0; j < 4; ++j)
				{
					tempVi[0] = oi->globalToLocalVVec[fvs[quadTriTable[j][0]]];
					tempVi[1] = oi->globalToLocalVVec[fvs[quadTriTable[j][1]]];
					tempVi[2] = oi->globalToLocalVVec[fvs[quadTriTable[j][2]]];
					oi->localTriIds.emplace_back(tempVi);

					oi->localStanTriPos.emplace_back(tempVPosVec[j]);
				}
			}
		}

		FeatureConstraints &newFc = oi->localFc;
		newFc.V_types.reserve(oi->localToGlobalVVec.size());
		for (int i = 0; i < oi->localToGlobalVVec.size(); ++i)
		{
			uint32_t refeVId = oi->localToGlobalVVec[i];
			newFc.V_types.emplace_back(fc.V_types[refeVId]);
		}

		for (int i = 0; i < fc.ids_C.size(); ++i)
		{
			uint32_t newV = oi->globalToLocalVVec[fc.ids_C[i]];
			if (newV == (uint32_t)-1)
				continue;

			newFc.ids_C.emplace_back(newV);
			newFc.C.emplace_back(fc.C[i]);
		}

		localToGlobalIdsL.clear();
		for (int i = 0; i < fc.ids_L.size(); ++i)
		{
			uint32_t newV = oi->globalToLocalVVec[fc.ids_L[i]];
			if (newV == (uint32_t)-1)
				continue;

			localToGlobalIdsL.emplace_back(i);
			newFc.ids_L.emplace_back(newV);
			newFc.on_which_L.emplace_back(fc.on_which_L[i]);
			newFc.axa_L.emplace_back(fc.axa_L[i]);
			newFc.origin_L.emplace_back(fc.origin_L[i]);
		}
	}

	double Pipeline::ComputeFaceArea(QuadMesh &qm, uint32_t faceId)
	{
		std::function<double(Eigen::Vector2d&, Eigen::Vector2d&)> crossLength = [&](Eigen::Vector2d& vec1, Eigen::Vector2d& vec2)->double
		{
			return std::abs(vec1[0] * vec2[1] - vec1[1] * vec2[0]);
		};

		std::vector<uint32_t> &fvs = qm.Fs_[faceId].vs;
		Eigen::Vector2d e0 = qm.V_[fvs[1]] - qm.V_[fvs[0]];
		Eigen::Vector2d e1 = qm.V_[fvs[1]] - qm.V_[fvs[2]];
		Eigen::Vector2d e2 = qm.V_[fvs[3]] - qm.V_[fvs[2]];
		Eigen::Vector2d e3 = qm.V_[fvs[3]] - qm.V_[fvs[0]];

		double sumArea = 0;
		sumArea += crossLength(e0, e3) / 2;
		sumArea += crossLength(e1, e2) / 2;
		/*sumArea += e0.cross(e3).norm() / 2;
		sumArea += e1.cross(e2).norm() / 2;*/
		return sumArea;
	}

	double Pipeline::ComputeFaceAreasSum(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &faceIds)
	{
		std::function<double(Eigen::Vector2d&, Eigen::Vector2d&)> crossLength = [&](Eigen::Vector2d& vec1, Eigen::Vector2d& vec2)->double
		{
			return std::abs(vec1[0] * vec2[1] - vec1[1] * vec2[0]);
		};

		double finalSum = 0;
		for (int i = 0; i < faceIds.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[faceIds[i]].vs;
			Eigen::Vector2d e0 = qm.V_[fvs[1]] - qm.V_[fvs[0]];
			Eigen::Vector2d e1 = qm.V_[fvs[1]] - qm.V_[fvs[2]];
			Eigen::Vector2d e2 = qm.V_[fvs[3]] - qm.V_[fvs[2]];
			Eigen::Vector2d e3 = qm.V_[fvs[3]] - qm.V_[fvs[0]];

			double sumArea = 0;
			sumArea += crossLength(e0, e3) / 2;
			sumArea += crossLength(e1, e2) / 2;
			/*sumArea += e0.cross(e3).norm() / 2;
			sumArea += e1.cross(e2).norm() / 2;*/
			finalSum += sumArea;
		}
		return finalSum;
	}

	void Pipeline::ComputeStandardTriangles(QuadMesh &qm, uint32_t faceId, double targetArea, std::vector<std::vector<Eigen::Vector2d>> &stanTris)
	{
		std::vector<uint32_t> &fvs = qm.Fs_[faceId].vs;
		double e0 = 0, e1 = 0;
#ifndef USE_SQUARE_AREA
		e0 += (qm.V_[fvs[1]] - qm.V_[fvs[0]]).norm();
		e0 += (qm.V_[fvs[3]] - qm.V_[fvs[2]]).norm();
		e0 /= 2;
		e1 += (qm.V_[fvs[2]] - qm.V_[fvs[1]]).norm();
		e1 += (qm.V_[fvs[3]] - qm.V_[fvs[0]]).norm();
		e1 /= 2;
#else
		e0 = 1;
		e1 = 1;
#endif

#ifndef USE_AVERAGE_AREA
		double radio = std::sqrt(targetArea / (e0*e1));
#else
		double radio = std::sqrt((initialArea_/qm.Fs_.size())/ (e0*e1));
#endif
		e0 *= radio; e1 *= radio;

		std::vector<Eigen::Vector2d> standVVec(4);
		standVVec[0][0] = 0; standVVec[0][1] = 0;
		standVVec[1][0] = e0; standVVec[1][1] = 0;
		standVVec[2][0] = e0; standVVec[2][1] = e1;
		standVVec[3][0] = 0; standVVec[3][1] = e1;

		stanTris.resize(4);
		for (int j = 0; j < 4; ++j)
		{
			stanTris[j].resize(3);
			stanTris[j][0] = standVVec[quadTriTable[j][0]];
			stanTris[j][1] = standVVec[quadTriTable[j][1]];
			stanTris[j][2] = standVVec[quadTriTable[j][2]];
		}
	}

	void Pipeline::FinalOptimizeWithCornerChangeFixV2(BaseDataStructure::QuadMesh &qm)
	{
		std::cout << "Begin Final Optimization!" << std::endl;
		QuadMesh newQm = qm;

		std::vector<uint32_t> faces;
		faces.reserve(newQm.Fs_.size());
		bool isGoodF = true;
		for (int i = 0; i < newQm.Fs_.size(); ++i)
		{
			isGoodF = true;
			if (newQm.Fs_[i].boundary)
			{
				std::vector<uint32_t> &fvs = newQm.Fs_[i].vs;
				for (int j = 0; j < fvs.size(); ++j)
				{
					if (newQm.Vs_[fvs[j]].boundary && newQm.Vs_[fvs[j]].neighbor_fs.size() == 1 && fc.V_types[fvs[j]] != -1)
					{
						isGoodF = false;
						break;
					}
				}
			}
			if (isGoodF)
				faces.emplace_back(i);
		}

		OptimizeInfo newOi;
		BuildNewOptimizeInfo(newQm, fc, faces, &newOi);

#if CHANGE_THE_CORNERS
		ChangeTheCornersInTwoCurves(newQm, &newOi);
#endif

		AESolverSquare aes(&newOi, mf, kdTree_, fields_, false);
		aes.ExtractBoundaryTriangles(&newQm);
		aes.SetQuadMesh(&newQm);
		if (aes.isInverse)
			return;
		ComputeIdsLDis(&newOi);
#if OPTIMIZE_FIELDS
		ExtraFieldEnergyAlone efea(newQm, &newOi, kdTree_, fields_, 1, 1, newQm.ComputeAverageEdgeLength());
#endif
		for (int i = 0; i < AES_ITER * 3; ++i)
		{
			aes.efe_.ExtractFields();
			aes.Optimize();
			aes.SetValueBack(&newOi);
			ProjectToSouceLine(&newOi, mf);
#if USE_3To0
			AdjustBoundaryPosInOptimization(newQm, &newOi);
#endif
			//for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
			//{
			//	//uint32_t newVId = newOi.localToGlobalVVec[i];
			//	newQm.V_[i] = newOi.localVPos[i];
			//}
			/*QuadMeshIO qmio;
			qmio.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\error8.obj");*/

			aes.SetValueFront(&newOi);
			if (aes.isInverse)
				return;

#if OPTIMIZE_FIELDS
			if (i < AES_ITER * 2)
			{
				efea.Initialize();
				efea.SetValueFront(&newOi);
				efea.Optimize();
				efea.SetValueBack(&newOi);
			}
#endif
		}

		for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
		{
			uint32_t newVId = newOi.localToGlobalVVec[i];
			newQm.V_[newVId] = newOi.localVPos[i];
		}
		/*QuadMeshIO qmio;
		qmio.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\error9.obj");*/
		if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			std::cout << "Inverse Jacobi! " << std::endl;
			return;
		}
		if (!ComputeHausdorff(backupQm_, newQm, hausdorffRadio, hausdorffThre))
			return;

		/*FeatureConstraints newFc;
		UpdateFeatures(&newOi, fc, newFc);*/
		UpdateFeaturesFromOi(fc, &newOi);

		std::cout << "Final Optimize OK" << std::endl;
		qm = newQm;
	}

	void Pipeline::FinalOptimizeWithCornerChange(BaseDataStructure::QuadMesh &qm)
	{
		std::cout << "Begin Final Optimization!" << std::endl;
		QuadMesh newQm = qm;

		OptimizeInfo newOi;
		newOi.globalToLocalVVec.reserve(newQm.Vs_.size());
		newOi.localToGlobalVVec.reserve(newQm.Vs_.size());
		for (int i = 0; i < newQm.Vs_.size(); ++i)
		{
			newOi.globalToLocalVVec.emplace_back(i);
			newOi.localToGlobalVVec.emplace_back(i);
			newOi.localVPos.emplace_back(newQm.V_[i]);
		}
		newOi.localFc = fc;

		std::vector<Eigen::Vector2d> standVVec(4);
		Eigen::Vector3i tempVi;
		standVVec[0][0] = 0; standVVec[0][1] = 0;
		for (int i = 0; i < newQm.Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fvs = newQm.Fs_[i].vs;
			double e0 = 0, e1 = 0;
#ifndef USE_SQUARE_AREA
			e0 += (newQm.V_[fvs[1]] - newQm.V_[fvs[0]]).norm();
			e0 += (newQm.V_[fvs[3]] - newQm.V_[fvs[2]]).norm();
			e0 /= 2;
			e1 += (newQm.V_[fvs[2]] - newQm.V_[fvs[1]]).norm();
			e1 += (newQm.V_[fvs[3]] - newQm.V_[fvs[0]]).norm();
			e1 /= 2;
#else
			e0 = 1;
			e1 = 1;
#endif

#ifndef USE_AVERAGE_AREA
			double referenceArea = ComputeFaceArea(newQm, i);
#else
			double referenceArea = initialArea_ / newQm.Fs_.size();
#endif

			double radio = std::sqrt(referenceArea / (e0*e1));
			e0 *= radio; e1 *= radio;

			standVVec[1][0] = e0; standVVec[1][1] = 0;
			standVVec[2][0] = e0; standVVec[2][1] = e1;
			standVVec[3][0] = 0; standVVec[3][1] = e1;

			std::vector<Eigen::Vector2d> tempVVec(3);
			for (int j = 0; j < 4; ++j)
			{
				tempVVec[0] = standVVec[quadTriTable[j][0]];
				tempVVec[1] = standVVec[quadTriTable[j][1]];
				tempVVec[2] = standVVec[quadTriTable[j][2]];
				newOi.localStanTriPos.emplace_back(tempVVec);

				tempVi[0] = fvs[quadTriTable[j][0]];
				tempVi[1] = fvs[quadTriTable[j][1]];
				tempVi[2] = fvs[quadTriTable[j][2]];
				newOi.localTriIds.emplace_back(tempVi);
			}
		}

#if CHANGE_THE_CORNERS
		ChangeTheCornersInTwoCurves(newQm, &newOi);
#endif

		AESolverSquare aes(&newOi, mf, kdTree_, fields_, false);
		aes.ExtractBoundaryTriangles(&newQm);
		if (aes.isInverse)
			return;
		ComputeIdsLDis(&newOi);
#if OPTIMIZE_FIELDS
		ExtraFieldEnergyAlone efea(newQm, &newOi, kdTree_, fields_, 1, 1, newQm.ComputeAverageEdgeLength());
#endif
		for (int i = 0; i < AES_ITER * 3; ++i)
		{
			aes.efe_.ExtractFields();
			aes.Optimize();
			aes.SetValueBack(&newOi);
			ProjectToSouceLine(&newOi, mf);
#if USE_3To0
			AdjustBoundaryPosInOptimization(newQm, &newOi);
#endif
			for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
			{
				//uint32_t newVId = newOi.localToGlobalVVec[i];
				newQm.V_[i] = newOi.localVPos[i];
			}
			/*QuadMeshIO qmio;
			qmio.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\error8.obj");*/

			aes.SetValueFront(&newOi);
			if (aes.isInverse)
				return;

#if OPTIMIZE_FIELDS
			if (i < AES_ITER*2)
			{
				efea.Initialize();
				efea.SetValueFront(&newOi);
				efea.Optimize();
				efea.SetValueBack(&newOi);
			}
#endif
		}

		for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
		{
			//uint32_t newVId = newOi.localToGlobalVVec[i];
			newQm.V_[i] = newOi.localVPos[i];
		}
		/*QuadMeshIO qmio;
		qmio.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\error9.obj");*/
		if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			std::cout << "Inverse Jacobi! " << std::endl;
			return;
		}
		if (!ComputeHausdorff(backupQm_, newQm, hausdorffRadio, hausdorffThre))
			return;

		/*FeatureConstraints newFc;
		UpdateFeatures(&newOi, fc, newFc);*/
		UpdateFeaturesFromOi(fc, &newOi);

		std::cout << "Final Optimize OK" << std::endl;
		qm = newQm;
	}

	void Pipeline::FinalOptimize(QuadMesh &qm)
	{
		std::cout << "Begin Final Optimization!" << std::endl;
		QuadMesh newQm = qm;

		OptimizeInfo newOi;
		newOi.globalToLocalVVec.reserve(newQm.Vs_.size());
		newOi.localToGlobalVVec.reserve(newQm.Vs_.size());
		for (int i = 0; i < newQm.Vs_.size(); ++i)
		{
			newOi.globalToLocalVVec.emplace_back(i);
			newOi.localToGlobalVVec.emplace_back(i);
			newOi.localVPos.emplace_back(newQm.V_[i]);
		}
		newOi.localFc = fc;

		std::vector<Eigen::Vector2d> standVVec(4);
		Eigen::Vector3i tempVi;
		standVVec[0][0] = 0; standVVec[0][1] = 0;
		for (int i = 0; i < newQm.Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fvs = newQm.Fs_[i].vs;
			double e0 = 0, e1 = 0;
#ifndef USE_SQUARE_AREA
			e0 += (newQm.V_[fvs[1]] - newQm.V_[fvs[0]]).norm();
			e0 += (newQm.V_[fvs[3]] - newQm.V_[fvs[2]]).norm();
			e0 /= 2;
			e1 += (newQm.V_[fvs[2]] - newQm.V_[fvs[1]]).norm();
			e1 += (newQm.V_[fvs[3]] - newQm.V_[fvs[0]]).norm();
			e1 /= 2;
#else
			e0 = 1;
			e1 = 1;
#endif

#ifndef USE_AVERAGE_AREA
			double referenceArea = ComputeFaceArea(newQm, i);
#else
			double referenceArea = initialArea_ / newQm.Fs_.size();
#endif

			double radio = std::sqrt(referenceArea / (e0*e1));
			e0 *= radio; e1 *= radio;

			standVVec[1][0] = e0; standVVec[1][1] = 0;
			standVVec[2][0] = e0; standVVec[2][1] = e1;
			standVVec[3][0] = 0; standVVec[3][1] = e1;

			std::vector<Eigen::Vector2d> tempVVec(3);
			for (int j = 0; j < 4; ++j)
			{
				tempVVec[0] = standVVec[quadTriTable[j][0]];
				tempVVec[1] = standVVec[quadTriTable[j][1]];
				tempVVec[2] = standVVec[quadTriTable[j][2]];
				newOi.localStanTriPos.emplace_back(tempVVec);

				tempVi[0] = fvs[quadTriTable[j][0]];
				tempVi[1] = fvs[quadTriTable[j][1]];
				tempVi[2] = fvs[quadTriTable[j][2]];
				newOi.localTriIds.emplace_back(tempVi);
			}
		}

//#if CHANGE_THE_CORNERS
//		ChangeTheCornersInTwoCurves(newQm, &newOi);
//#endif

		AESolverSquare aes(&newOi, mf, kdTree_, fields_, true);
		aes.ExtractBoundaryTriangles(&newQm);
		if (aes.isInverse)
			return;
		ComputeIdsLDis(&newOi);
#if OPTIMIZE_FIELDS
		ExtraFieldEnergyAlone efea(newQm, &newOi, kdTree_, fields_, 1, 1, newQm.ComputeAverageEdgeLength());
#endif
		for (int i = 0; i < AES_ITER * 3; ++i)
		{
			aes.efe_.ExtractFields();
			aes.Optimize();
			aes.SetValueBack(&newOi);
			ProjectToSouceLine(&newOi, mf);
#if USE_3To0
			AdjustBoundaryPosInOptimization(newQm, &newOi);
#endif
			//for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
			//{
			//	//uint32_t newVId = newOi.localToGlobalVVec[i];
			//	newQm.V_[i] = newOi.localVPos[i];
			//}
			//QuadMeshIO qmio;
			//qmio.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\error8.obj");

			aes.SetValueFront(&newOi);
			if (aes.isInverse)
				return;

#if OPTIMIZE_FIELDS
			if (i < AES_ITER*2)
			{
				efea.Initialize();
				efea.SetValueFront(&newOi);
				efea.Optimize();
				efea.SetValueBack(&newOi);
			}
#endif
		}

		for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
		{
			//uint32_t newVId = newOi.localToGlobalVVec[i];
			newQm.V_[i] = newOi.localVPos[i];
		}
		/*QuadMeshIO qmio;
		qmio.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\error9.obj");*/
		if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			std::cout << "Inverse Jacobi! " << std::endl;
			return;
		}
		if (!ComputeHausdorff(backupQm_, newQm, hausdorffRadio, hausdorffThre))
			return;

		/*FeatureConstraints newFc;
		UpdateFeatures(&newOi, fc, newFc);*/

		UpdateFeaturesFromOi(fc, &newOi);

		std::cout << "Final Optimize OK" << std::endl;
		qm = newQm;
	}

	void Pipeline::FinalOptimizeFixV2(QuadMesh &qm)
	{
		std::cout << "Begin Final Optimization Fix V2!" << std::endl;
		QuadMesh newQm = qm;

		std::vector<uint32_t> faces;
		faces.reserve(newQm.Fs_.size());
		bool isGoodF = true;
		for (int i = 0; i < newQm.Fs_.size(); ++i)
		{
			isGoodF = true;
			if (newQm.Fs_[i].boundary)
			{
				std::vector<uint32_t> &fvs = newQm.Fs_[i].vs;
				for (int j = 0; j < fvs.size(); ++j)
				{
					if (newQm.Vs_[fvs[j]].boundary && newQm.Vs_[fvs[j]].neighbor_fs.size() == 1 && fc.V_types[fvs[j]] != -1)
					{
						isGoodF = false;
						break;
					}
				}
			}
			if (isGoodF)
				faces.emplace_back(i);
		}
		/*for (int i = 0; i < newQm.Fs_.size(); ++i)
		{
			if (!newQm.Fs_[i].boundary)
				faces.emplace_back(i);
		}*/

		OptimizeInfo newOi;
		BuildNewOptimizeInfo(newQm, fc, faces, &newOi);
		/*newOi.globalToLocalVVec.reserve(newQm.Vs_.size());
		newOi.localToGlobalVVec.reserve(newQm.Vs_.size());
		for (int i = 0; i < newQm.Vs_.size(); ++i)
		{
			newOi.globalToLocalVVec.emplace_back(i);
			newOi.localToGlobalVVec.emplace_back(i);
			newOi.localVPos.emplace_back(newQm.V_[i]);
		}
		newOi.localFc = fc;

		std::vector<Eigen::Vector2d> standVVec(4);
		Eigen::Vector3i tempVi;
		standVVec[0][0] = 0; standVVec[0][1] = 0;
		for (int i = 0; i < newQm.Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fvs = newQm.Fs_[i].vs;
			double e0 = 0, e1 = 0;
#ifndef USE_SQUARE_AREA
			e0 += (newQm.V_[fvs[1]] - newQm.V_[fvs[0]]).norm();
			e0 += (newQm.V_[fvs[3]] - newQm.V_[fvs[2]]).norm();
			e0 /= 2;
			e1 += (newQm.V_[fvs[2]] - newQm.V_[fvs[1]]).norm();
			e1 += (newQm.V_[fvs[3]] - newQm.V_[fvs[0]]).norm();
			e1 /= 2;
#else
			e0 = 1;
			e1 = 1;
#endif

#ifndef USE_AVERAGE_AREA
			double referenceArea = ComputeFaceArea(newQm, i);
#else
			double referenceArea = initialArea_ / newQm.Fs_.size();
#endif

			double radio = std::sqrt(referenceArea / (e0*e1));
			e0 *= radio; e1 *= radio;

			standVVec[1][0] = e0; standVVec[1][1] = 0;
			standVVec[2][0] = e0; standVVec[2][1] = e1;
			standVVec[3][0] = 0; standVVec[3][1] = e1;

			std::vector<Eigen::Vector2d> tempVVec(3);
			for (int j = 0; j < 4; ++j)
			{
				tempVVec[0] = standVVec[quadTriTable[j][0]];
				tempVVec[1] = standVVec[quadTriTable[j][1]];
				tempVVec[2] = standVVec[quadTriTable[j][2]];
				newOi.localStanTriPos.emplace_back(tempVVec);

				tempVi[0] = fvs[quadTriTable[j][0]];
				tempVi[1] = fvs[quadTriTable[j][1]];
				tempVi[2] = fvs[quadTriTable[j][2]];
				newOi.localTriIds.emplace_back(tempVi);
			}
		}*/

		//#if CHANGE_THE_CORNERS
		//		ChangeTheCornersInTwoCurves(newQm, &newOi);
		//#endif

		AESolverSquare aes(&newOi, mf, kdTree_, fields_, true);
		aes.ExtractBoundaryTriangles(&newQm);
		if (aes.isInverse)
			return;
		ComputeIdsLDis(&newOi);
#if OPTIMIZE_FIELDS
		ExtraFieldEnergyAlone efea(newQm, &newOi, kdTree_, fields_, 1, 1, newQm.ComputeAverageEdgeLength());
#endif
		for (int i = 0; i < AES_ITER * 3; ++i)
		{
			aes.efe_.ExtractFields();
			aes.Optimize();
			aes.SetValueBack(&newOi);
			ProjectToSouceLine(&newOi, mf);
#if USE_3To0
			AdjustBoundaryPosInOptimization(newQm, &newOi);
#endif
			//for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
			//{
			//	//uint32_t newVId = newOi.localToGlobalVVec[i];
			//	newQm.V_[i] = newOi.localVPos[i];
			//}
			//QuadMeshIO qmio;
			//qmio.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\error8.obj");

			aes.SetValueFront(&newOi);
			if (aes.isInverse)
				return;

#if OPTIMIZE_FIELDS
			if (i < AES_ITER * 2)
			{
				efea.Initialize();
				efea.SetValueFront(&newOi);
				efea.Optimize();
				efea.SetValueBack(&newOi);
			}
#endif
		}

		for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
		{
			uint32_t newVId = newOi.localToGlobalVVec[i];
			newQm.V_[newVId] = newOi.localVPos[i];
		}
		/*QuadMeshIO qmio;
		qmio.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\error9.obj");*/
		if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			std::cout << "Inverse Jacobi! " << std::endl;
			return;
		}
		if (!ComputeHausdorff(backupQm_, newQm, hausdorffRadio, hausdorffThre))
			return;

		/*FeatureConstraints newFc;
		UpdateFeatures(&newOi, fc, newFc);*/

		UpdateFeaturesFromOi(fc, &newOi);

		std::cout << "Final Optimize OK" << std::endl;
		qm = newQm;
	}

	void Pipeline::FinalOptimizeInner(QuadMesh &qm)
	{
		std::cout << "Begin Final Optimization Fix V2!" << std::endl;
		QuadMesh newQm = qm;

		std::vector<uint32_t> faces;
		faces.reserve(newQm.Fs_.size());
		/*bool isGoodF = true;
		for (int i = 0; i < newQm.Fs_.size(); ++i)
		{
			isGoodF = true;
			if (newQm.Fs_[i].boundary)
			{
				std::vector<uint32_t> &fvs = newQm.Fs_[i].vs;
				for (int j = 0; j < fvs.size(); ++j)
				{
					if (newQm.Vs_[fvs[j]].boundary && newQm.Vs_[fvs[j]].neighbor_fs.size() == 1 && fc.V_types[fvs[j]] != -1)
					{
						isGoodF = false;
						break;
					}
				}
			}
			if (isGoodF)
				faces.emplace_back(i);
		}*/
		for (int i = 0; i < newQm.Fs_.size(); ++i)
		{
			if (newQm.Fs_[i].boundary)
				continue;

			std::vector<uint32_t> &fvs = newQm.Fs_[i].vs;
			bool isBF = false;
			for (int j = 0; j < fvs.size(); ++j)
			{
				if (newQm.Vs_[fvs[j]].boundary)
				{
					isBF = true;
					break;
				}
			}
			if (!isBF)
				faces.emplace_back(i);
		}

		OptimizeInfo newOi;
		BuildNewOptimizeInfo(newQm, fc, faces, &newOi);
		/*newOi.globalToLocalVVec.reserve(newQm.Vs_.size());
		newOi.localToGlobalVVec.reserve(newQm.Vs_.size());
		for (int i = 0; i < newQm.Vs_.size(); ++i)
		{
			newOi.globalToLocalVVec.emplace_back(i);
			newOi.localToGlobalVVec.emplace_back(i);
			newOi.localVPos.emplace_back(newQm.V_[i]);
		}
		newOi.localFc = fc;

		std::vector<Eigen::Vector2d> standVVec(4);
		Eigen::Vector3i tempVi;
		standVVec[0][0] = 0; standVVec[0][1] = 0;
		for (int i = 0; i < newQm.Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fvs = newQm.Fs_[i].vs;
			double e0 = 0, e1 = 0;
#ifndef USE_SQUARE_AREA
			e0 += (newQm.V_[fvs[1]] - newQm.V_[fvs[0]]).norm();
			e0 += (newQm.V_[fvs[3]] - newQm.V_[fvs[2]]).norm();
			e0 /= 2;
			e1 += (newQm.V_[fvs[2]] - newQm.V_[fvs[1]]).norm();
			e1 += (newQm.V_[fvs[3]] - newQm.V_[fvs[0]]).norm();
			e1 /= 2;
#else
			e0 = 1;
			e1 = 1;
#endif

#ifndef USE_AVERAGE_AREA
			double referenceArea = ComputeFaceArea(newQm, i);
#else
			double referenceArea = initialArea_ / newQm.Fs_.size();
#endif

			double radio = std::sqrt(referenceArea / (e0*e1));
			e0 *= radio; e1 *= radio;

			standVVec[1][0] = e0; standVVec[1][1] = 0;
			standVVec[2][0] = e0; standVVec[2][1] = e1;
			standVVec[3][0] = 0; standVVec[3][1] = e1;

			std::vector<Eigen::Vector2d> tempVVec(3);
			for (int j = 0; j < 4; ++j)
			{
				tempVVec[0] = standVVec[quadTriTable[j][0]];
				tempVVec[1] = standVVec[quadTriTable[j][1]];
				tempVVec[2] = standVVec[quadTriTable[j][2]];
				newOi.localStanTriPos.emplace_back(tempVVec);

				tempVi[0] = fvs[quadTriTable[j][0]];
				tempVi[1] = fvs[quadTriTable[j][1]];
				tempVi[2] = fvs[quadTriTable[j][2]];
				newOi.localTriIds.emplace_back(tempVi);
			}
		}*/

		//#if CHANGE_THE_CORNERS
		//		ChangeTheCornersInTwoCurves(newQm, &newOi);
		//#endif

		AESolverSquare aes(&newOi, mf, kdTree_, fields_, false);
		aes.SetQuadMesh(&newQm);
//		aes.ExtractBoundaryTriangles(&newQm);
		if (aes.isInverse)
			return;
		ComputeIdsLDis(&newOi);
#if OPTIMIZE_FIELDS
		//ExtraFieldEnergyAlone efea(newQm, &newOi, kdTree_, fields_, 1, 1, newQm.ComputeAverageEdgeLength());
#endif
		for (int i = 0; i < AES_ITER * 3; ++i)
		{
//			aes.efe_.ExtractFields();
			aes.Optimize();
			aes.SetValueBack(&newOi);
//			ProjectToSouceLine(&newOi, mf);
#if USE_3To0
			AdjustBoundaryPosInOptimization(newQm, &newOi);
#endif
			///*for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
			//{
			//	uint32_t newVId = newOi.localToGlobalVVec[i];
			//	newQm.V_[newVId] = newOi.localVPos[i];
			//}
			//QuadMeshIO qmio;
			//qmio.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\error8.obj");
			//if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
			//{
			//	std::cout << "Inverse Jacobi! " << std::endl;
			//	return;
			//}*/

			aes.SetValueFront(&newOi);
			if (aes.isInverse)
				return;

#if OPTIMIZE_FIELDS
			/*if (i < AES_ITER * 2)
			{
				efea.Initialize();
				efea.SetValueFront(&newOi);
				efea.Optimize();
				efea.SetValueBack(&newOi);
			}*/
#endif
		}

		for (int i = 0; i < newOi.localToGlobalVVec.size(); ++i)
		{
			uint32_t newVId = newOi.localToGlobalVVec[i];
			newQm.V_[newVId] = newOi.localVPos[i];
		}
		QuadMeshIO qmio;
		qmio.WriteQuadMesh(&newQm, "C:\\Users\\ChiZhang\\Desktop\\error9.obj");
		if (!newQm.JudgeQuadMeshJacobi() || newQm.minJacobi_ < MIN_ALLOWED_JACOBI)
		{
			//std::cout << newQm.minJacobi_ << std::endl;
			std::cout << "Inverse Jacobi! " << std::endl;
			return;
		}
		if (!ComputeHausdorff(backupQm_, newQm, hausdorffRadio, hausdorffThre))
			return;

		/*FeatureConstraints newFc;
		UpdateFeatures(&newOi, fc, newFc);*/

		//UpdateFeaturesFromOi(fc, &newOi);

		std::cout << "Final Optimize OK" << std::endl;
		qm = newQm;
	}

	void Pipeline::ComputeIdsLDis(OptimizeInfo *oi)
	{
		typedef std::tuple<double, uint32_t, uint32_t, uint32_t> TempTuple4;
		//std::vector<TempTuple4>().swap(oi->idsLDis);
		oi->idsLDis.clear();
		oi->idsLDis.reserve(oi->localFc.ids_L.size());
		for (int i = 0; i < oi->localFc.ids_L.size(); ++i)
		{
			uint32_t currL = oi->localFc.on_which_L[i];
			const Point &cornerP = mf.segList[currL][0].vertex(0);
			Point currP(oi->localVPos[oi->localFc.ids_L[i]][0], oi->localVPos[oi->localFc.ids_L[i]][1], 0);
			double tempDis = std::sqrt((currP - cornerP).squared_length());
			oi->idsLDis.emplace_back(TempTuple4(tempDis, currL, oi->localFc.ids_L[i], i));
		}
		std::sort(oi->idsLDis.begin(), oi->idsLDis.end());
	}

	void Pipeline::RefineQuadToTri(BaseDataStructure::QuadMesh &qm, MyMesh &triMesh)
	{
		triMesh.clear();
		std::vector<MyMesh::VertexHandle> vHandles;
		vHandles.reserve(qm.Vs_.size());
		std::vector<MyMesh::VertexHandle> fHandles;
		fHandles.reserve(3);

		for (int i = 0; i < qm.Vs_.size(); ++i)
		{
			vHandles.emplace_back(triMesh.add_vertex(MyMesh::Point(qm.V_[i][0], qm.V_[i][1], 0)));
		}

		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			Eigen::Vector2d avePos(0, 0);
			for (int j = 0; j < 4; ++j)
			{
				avePos += qm.V_[qm.Fs_[i].vs[j]];
			}
			avePos /= 4.0;

			MyMesh::VertexHandle tempV = triMesh.add_vertex(MyMesh::Point(avePos[0], avePos[1], 0));

			for (int j = 0; j < 4; ++j)
			{
				fHandles.clear();
				fHandles.emplace_back(tempV);
				fHandles.emplace_back(vHandles[qm.Fs_[i].vs[j]]);
				fHandles.emplace_back(vHandles[qm.Fs_[i].vs[(j+1)%4]]);
				triMesh.add_face(fHandles);
			}
		}

		triMesh.request_vertex_status();
		triMesh.request_edge_status();
		triMesh.request_face_status();

		triMesh.request_face_normals();
		triMesh.request_vertex_normals();
		triMesh.update_face_normals();
	}

	void Pipeline::BuildKDTree(std::vector<OpenMesh::Vec3d> &triMids)
	{
		if (kdTree_ != NULL)
		{
			delete kdTree_;
			kdTree_ = NULL;
		}

		ANNpointArray dataPts = annAllocPts(triMids.size(), 3);
		for (int i = 0; i < triMids.size(); ++i)
		{
			dataPts[i][0] = triMids[i][0];
			dataPts[i][1] = triMids[i][1];
			dataPts[i][2] = triMids[i][2];
		}
		kdTree_ = new ANNkd_tree(dataPts, triMids.size(), 3);
	}

	void Pipeline::RankingSheetsAndChordsByFields(BaseDataStructure::QuadMesh &qm, BaseDataStructure::BaseComplex &bc, std::vector<BaseDataStructure::Sheet> &sheets, std::vector<BaseDataStructure::Chord> &chords, std::vector<RankingTuple3> &candidates)
	{
		for (int i = 0; i < sheets.size(); ++i)
		{
			double sumAngle = 0;
			for (int j = 0; j < sheets[i].qf_delete.size(); ++j)
			{
				sumAngle += ComputeQuadFieldWeight(qm, sheets[i].qf_delete[j]);
			}
			
			sheets[i].weight = -sumAngle / sheets[i].qf_delete.size() + WEITH_WIDGET * sheets[i].weight;
		}

		for (int i = 0; i < chords.size(); ++i)
		{
			double sumAngle = 0;
			/*for (int j = 0; j < chords[i].qf_delete.size(); ++j)
			{
				sumAngle += ComputeQuadFieldWeight(qm, chords[i].qf_delete[j]);
			}*/
			for (int j = 0; j < bc.Bf_[chords[i].fid].fs_net.size(); ++j)
			{
				sumAngle += ComputeQuadFieldWeight(qm, bc.Bf_[chords[i].fid].fs_net[j]);
			}
			chords[i].weight = -sumAngle / bc.Bf_[chords[i].fid].fs_net.size() + WEITH_WIDGET * chords[i].weight;
		}

		//std::vector<RankingTuple3>().swap(candidates);
		candidates.clear();
		candidates.reserve(sheets.size() + chords.size());

		for (int i = 0; i < sheets.size(); ++i)
		{
			candidates.emplace_back(RankingTuple3(sheets[i].weight, ElementType::SHEET, sheets[i].id));
		}
		for (int i = 0; i < chords.size(); ++i)
		{
			candidates.emplace_back(RankingTuple3(chords[i].weight, ElementType::CHORD, chords[i].id));
		}
		std::sort(candidates.begin(), candidates.end());
	}

	double Pipeline::ComputeQuadFieldWeight(BaseDataStructure::QuadMesh &qm, uint32_t quadId)
	{
		const std::vector<uint32_t> &fvs = qm.Fs_[quadId].vs;
		const std::vector<uint32_t> &fes = qm.Fs_[quadId].es;

		Eigen::Vector2d midP(0, 0);
		for (int i = 0; i < 4; ++i)
		{
			midP += qm.V_[fvs[i]];
		}
		midP /= 4.0;

		ANNpoint ap = annAllocPt(3);
		ap[0] = midP[0]; ap[1] = midP[1]; ap[2] = 0;
		ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
		kdTree_->annkSearch(ap, 1, nnIdx, dists);
		const Eigen::Vector3d currStanField = fields_[nnIdx[0]];

		double currAngle = 0;
		for (int i = 0; i < 4; ++i)
		{
			const std::vector<uint32_t> &evs = qm.Es_[fes[i]].vs;
			currAngle += abs(ComputeEdgeFieldAngle(currStanField, qm.V_[evs[0]], qm.V_[evs[1]]));
		}
		currAngle /= 4;
		annDeallocPt(ap);
		delete[] nnIdx; delete[] dists;
		return currAngle;
	}

	double Pipeline::ComputeEdgeFieldAngle(const Eigen::Vector3d &field, const Eigen::Vector2d &ev0, const Eigen::Vector2d &ev1)
	{
		Eigen::Vector2d tVec0 = ev1 - ev0;
		Eigen::Vector2d tVec[4];
		tVec[0][0] = field[0]; tVec[0][1] = field[1];
		tVec0.normalize(); tVec[0].normalize();

		tVec[1][0] = -tVec[0][1]; tVec[1][1] = tVec[0][0];
		tVec[2][0] = -tVec[0][0]; tVec[2][1] = -tVec[0][1];
		tVec[3][0] = tVec[0][1]; tVec[3][1] = -tVec[0][0];

		double minAngle = 1000;
		for (int i = 0; i < 4; ++i)
		{
			if (tVec0[0] * tVec[i][0] + tVec0[1] * tVec[i][1] < 0)
				continue;

			double tempA = std::asin(tVec0[0] * tVec[i][1] - tVec0[1] * tVec[i][0]);
			if (std::abs(tempA) < std::abs(minAngle))
				minAngle = tempA;
		}

		return minAngle;
	}

	void Pipeline::DeleteValence2Vertices(QuadMesh &qm, FeatureConstraints &fcc, std::vector<uint32_t> &outFs)
	{
		const double offsetWeight = 0.15;
		std::function<bool(uint32_t, uint32_t, uint32_t, uint32_t&, uint32_t&)> findNextFE = [&](uint32_t startV, uint32_t preF, uint32_t preE, uint32_t &nextF, uint32_t &nextE)->bool
		{
			if (qm.Es_[preE].boundary)
				return false;

			if (qm.Es_[preE].neighbor_fs[0] == preF)
				nextF = qm.Es_[preE].neighbor_fs[1];
			else
				nextF = qm.Es_[preE].neighbor_fs[0];

			std::vector<uint32_t> &fes = qm.Fs_[nextF].es;
			for (int i = 0; i < 4; ++i)
			{
				if ((qm.Es_[fes[i]].vs[0] == startV || qm.Es_[fes[i]].vs[1] == startV) && fes[i] != preE)
				{
					nextE = fes[i];
					break;
				}
			}
			return true;
		};

		std::function<bool(uint32_t, uint32_t, uint32_t)> judgeRightEdge = [&](uint32_t currF, uint32_t currE, uint32_t targetV)->bool
		{
			uint32_t anotherV;
			if (qm.Es_[currE].vs[0] == targetV)
				anotherV = qm.Es_[currE].vs[1];
			else
				anotherV = qm.Es_[currE].vs[0];

			int tVId = qm.FindVInFId(targetV, currF);
			int aVId = qm.FindVInFId(anotherV, currF);
			if (tVId - aVId == 1 || (aVId == 3 && tVId == 0))
				return true;
			else
				return false;
		};

		int i = 0;
		while (i<qm.Vs_.size())
		{
			if (qm.Vs_[i].boundary && qm.Vs_[i].neighbor_fs.size() == 1 && fcc.V_types[i] != -1)
			{
				uint32_t badV = i, badF = qm.Vs_[i].neighbor_fs[0];
				int badVId = qm.FindVInFId(badV, badF);
				uint32_t startV = qm.Fs_[badF].vs[(badVId + 2) % 4];

				uint32_t oneSideE = (uint32_t)-1;
				std::vector<uint32_t> &fes = qm.Fs_[badF].es;
				for (int j = 0; j < 4; ++j)
				{
					if ((qm.Es_[fes[j]].vs[0] == startV || qm.Es_[fes[j]].vs[1] == startV) && judgeRightEdge(badF, fes[j], startV))
					{
						oneSideE = fes[j];
						break;
					}
				}
				uint32_t anSideE = (uint32_t)-1;
				for (int j = 0; j < 4; ++j)
				{
					if ((qm.Es_[fes[j]].vs[0] == startV || qm.Es_[fes[j]].vs[1] == startV) && fes[j]!=oneSideE)
					{
						anSideE = fes[j];
						break;
					}
				}

				std::vector<uint32_t> currFs;
				uint32_t preF, preE, nextF, nextE;
				//如果
				if (qm.Es_[oneSideE].boundary || qm.Es_[anSideE].boundary)
				{
					++i;
					continue;
				}
				else if (qm.Vs_[startV].boundary)
				{
					std::vector<bool> fFlag(qm.Fs_.size(), false);
					preF = badF; preE = oneSideE;
					while (findNextFE(startV, preF, preE, nextF, nextE))
					{
						fFlag[nextF] = true;
						currFs.emplace_back(nextF);
						preF = nextF; preE = nextE;
					}

					uint32_t anotherV;
					if (qm.Es_[preE].vs[0] == startV)
						anotherV = qm.Es_[preE].vs[1];
					else
						anotherV = qm.Es_[preE].vs[0];
					
					uint32_t rv, lv;
					if (qm.Es_[oneSideE].vs[0] == startV)
						rv = qm.Es_[oneSideE].vs[1];
					else
						rv = qm.Es_[oneSideE].vs[0];
					std::vector<uint32_t> &tempFvs = qm.Fs_[badF].vs;
					for (int j = 0; j < 4; ++j)
					{
						if (tempFvs[j] != startV && tempFvs[j] != badV && tempFvs[j] != rv)
						{
							lv = tempFvs[j];
							break;
						}
					}

					qm.V_[badV] = qm.V_[startV] + (qm.V_[rv] - qm.V_[lv])*offsetWeight*0.5;
					std::vector<BaseDataStructure::QuadFace> newFs;
					std::vector<BaseDataStructure::QuadVertex> newVs;
					newVs.reserve(qm.Vs_.size());
					newFs.reserve(qm.Fs_.size());
					for (int j = 0; j < qm.Vs_.size(); ++j)
					{
						BaseDataStructure::QuadVertex qv;
						qv.id = j; qv.bId = (uint32_t)-1;
						newVs.emplace_back(qv);
					}
					for (int j = 0; j < qm.Fs_.size(); ++j)
					{
						BaseDataStructure::QuadFace qf;
						qf.id = j; qf.bId = (uint32_t)-1;
						qf.vs.reserve(4);
						std::vector<uint32_t> &currFvs = qm.Fs_[j].vs;
						if (!fFlag[j] && j != badF)
						{
							qf.vs = currFvs;
							for (int k = 0; k < 4; ++k)
							{
								newVs[qf.vs[k]].neighbor_fs.emplace_back(j);
							}
						}
						else if (j == badF)
						{
							qf.vs.emplace_back(currFvs[(badVId + 1) % 4]);
							qf.vs.emplace_back(badV);
							qf.vs.emplace_back(currFvs[(badVId + 2) % 4]);
							qf.vs.emplace_back(currFvs[(badVId + 3) % 4]);
							for (int k = 0; k < 4; ++k)
							{
								newVs[qf.vs[k]].neighbor_fs.emplace_back(j);
							}
						}
						else
						{
							
							for (int k = 0; k < 4; ++k)
							{
								if (currFvs[k] == startV)
									qf.vs.emplace_back(badV);
								else
									qf.vs.emplace_back(currFvs[k]);
							}
							for (int k = 0; k < 4; ++k)
							{
								newVs[qf.vs[k]].neighbor_fs.emplace_back(j);
							}
						}
						newFs.emplace_back(qf);
					}
					qm.Es_.clear();
					qm.Vs_ = newVs;
					qm.Fs_ = newFs;
					qm.BuildConnectivity();
					/*QuadMeshIO qmi;
					qmi.WriteQuadMesh(&qm, "C:\\Users\\ChiZhang\\Desktop\\3210.obj");*/
				}
				else
				{
					std::vector<uint32_t> currEs;
					uint32_t startF;
					if (qm.Es_[oneSideE].neighbor_fs[0] == badF)
						startF = qm.Es_[oneSideE].neighbor_fs[1];
					else
						startF = qm.Es_[oneSideE].neighbor_fs[0];
					qm.FindOrderEdges(oneSideE, startF, startV, currEs, currFs);
					outFs.emplace_back(badF);
					if (!qm.SplitBadLine(badV, badF, startV, currEs[(currEs.size()-1) / 2], outFs))
					{
						++i;
						continue;
					}
					/*QuadMeshIO qmi;
					qmi.WriteQuadMesh(&qm, "C:\\Users\\ChiZhang\\Desktop\\3210.obj");*/
				}
				while (fcc.V_types.size() < qm.Vs_.size())
				{
					fcc.V_types.emplace_back(-2);
				}
				i = 0;

				//更新feature
				auto it = std::find(fcc.ids_L.begin(), fcc.ids_L.end(), badV);
				int vId = std::distance(fcc.ids_L.begin(), it);
				std::vector<uint32_t> &vvs = qm.Vs_[badV].neighbor_vs;
				std::vector<uint32_t> bvs; bvs.reserve(4);
				for (int j = 0; j < vvs.size(); ++j)
				{
					if (qm.Vs_[vvs[j]].boundary)
						bvs.emplace_back(vvs[j]);
				}
				if (bvs.size() < 2)
					exit(-12433);

				uint32_t currCurve;
				if (fcc.V_types[bvs[0]] != -1 || fcc.V_types[bvs[1]] != -1)
				{
					currCurve = fcc.V_types[bvs[0]] != -1 ? fcc.V_types[bvs[0]] : fcc.V_types[bvs[1]];
				}
				else
				{
					Eigen::Vector2d &cPos0 = qm.V_[bvs[0]], &cPos1 = qm.V_[bvs[1]];
					uint32_t cor0 = (uint32_t)-1, cor1 = (uint32_t)-1;
					FindSourceCorner(cPos0, cor0);
					FindSourceCorner(cPos1, cor1);
					std::tuple<int, int> &corCur0 = mf.cornerCurves[cor0], &corCur1 = mf.cornerCurves[cor1];
					if (std::get<0>(corCur0) == std::get<0>(corCur1) || std::get<0>(corCur0) == std::get<1>(corCur1))
						currCurve = std::get<0>(corCur0);
					else if (std::get<1>(corCur0) == std::get<0>(corCur1) || std::get<1>(corCur0) == std::get<1>(corCur1))
						currCurve = std::get<1>(corCur0);
				}
				fcc.V_types[badV] = currCurve;
				fcc.on_which_L[vId] = currCurve;
				fcc.origin_L[vId] = qm.V_[badV];
				fcc.axa_L[vId] = (qm.V_[bvs[0]] - qm.V_[bvs[1]]).normalized();
			}
			else
				++i;
		}
#ifdef OUTPUT_MID_MESH
		QuadMeshIO qmi;
		qmi.WriteQuadMesh(&qm, "C:\\Users\\ChiZhang\\Desktop\\321.obj");
#endif
	}


	void Pipeline::FindMovedFaces(BaseDataStructure::QuadMesh &qm, std::vector<uint32_t> &faces, FeatureConstraints &fcc)	//在优化过程中，把一些
	{
		std::vector<uint32_t> tempFaces, badFs;
		tempFaces.reserve(faces.size());
		badFs.reserve(qm.Fs_.size());

		bool *fFlag = new bool[qm.Fs_.size()];
		std::memset(fFlag, 0, qm.Fs_.size() * sizeof(bool));

		for (int i = 0; i < faces.size(); ++i)
		{
			fFlag[faces[i]] = true;
			bool isGoodF = true;
			if (qm.Fs_[faces[i]].boundary)
			{
				std::vector<uint32_t> &fvs = qm.Fs_[faces[i]].vs;
				for (int j = 0; j < 4; ++j)
				{
					if (qm.Vs_[fvs[j]].boundary && qm.Vs_[fvs[j]].neighbor_fs.size() == 1 && fcc.V_types[fvs[j]] != -1)
					{
						isGoodF = false;
						break;
					}
				}
			}
			if (isGoodF)
				tempFaces.emplace_back(faces[i]);
			else
				badFs.emplace_back(faces[i]);
		}
		if (faces.size() == tempFaces.size())
		{
			delete[] fFlag;
			fFlag = NULL;
			return;
		}

		faces = tempFaces;

		for (int i = 0; i < qm.Fs_.size(); ++i)
		{
			if (!fFlag[i] && qm.Fs_[i].boundary)
			{
				fFlag[i] = true;
				std::vector<uint32_t> &fvs = qm.Fs_[i].vs;
				for (int j = 0; j < 4; ++j)
				{
					if (qm.Vs_[fvs[j]].boundary && qm.Vs_[fvs[j]].neighbor_fs.size() == 1 && fcc.V_types[fvs[j]] != -1)
					{
						badFs.emplace_back(i);
						break;
					}
				}
			}
		}

		for (int i = 0; i < badFs.size(); ++i)
		{
			std::vector<uint32_t> &fvs = qm.Fs_[badFs[i]].vs;
			for (int j = 0; j < 4; ++j)
			{
				if (qm.Vs_[fvs[j]].boundary && qm.Vs_[fvs[j]].neighbor_fs.size() == 1 && fcc.V_types[fvs[j]] != -1)
				{
					ComputeV2BoundaryPos(qm, fvs[j], MIN_ALLOWED_JACOBI);
				}
			}
		}

		delete[] fFlag;
		fFlag = NULL;
	}
}

