#include "PolySplit.h"
#include "ExtraFunc.h"
#include "CheckPolyProperty.h"
#include <queue>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <functional>
#include <fstream>
#define OLD_LENGTH_THRE 0
#define SHORT_LENGTH_NOT_DO 1
#define FOR_QUADS_GEN 0
#define OUTPUT_SPLIT_LINES 1

const double SAMPLE_LENGTH = 0.1;
const double PI = 3.141592653;

double LARGER_ANGLE = 195.0*PI / 180.0;
const double FIELD_MATCHING_THRESHOLD = 10.0 * PI / 180.0;
const double FIELD_MATCHING_THRESHOLD_2 = 8.0 * PI / 180.0;
double DO_NOTHING_ANGLE_THRESHOLD = 15.0*PI / 180.0;
double AREA_THRESHOLD = 0.1;
double ANGLE_THRESHOLD = 45 * PI / 180.0;
double LENGTH_THRESHOLD = 0.007;
double LENGTH_THRESHOLD_0 = 0;

double AREA_THRESHOLD_2 = 0.003;
double BEFORE_ANGLE_MIN_0 = 240.0 * PI / 180.0;
double BEFORE_ANGLE_MAX_0 = 300.0 * PI / 180.0;
double BEFORE_ANGLE_MIN_1 = 330.0 * PI / 180.0;

double AFTER_ANGLE_MIN = 70.0 * PI / 180.0;
double AFTER_ANGLE_MAX = 110.0 * PI / 180.0;

double AFTER_ANGLE_MIN_THREE = 70.0 * PI / 180.0;
double AFTER_ANGLE_MAX_THREE = 110.0 * PI / 180.0;

double AREA_LENGTH_WEIGHT = 0.2;

double SHORT_NOT_DO_LENGTH = 0.007;

PolySplit::PolySplit()
{
}


PolySplit::~PolySplit()
{
	if (aabbTree_ != NULL)
		delete aabbTree_;
}

void PolySplit::BuildAABBTree(std::vector<Eigen::Vector2d> &inPoly)
{
	segList_.clear();
	segList_.reserve(inPoly.size());
	for (int i = 0; i < inPoly.size(); ++i)
	{
		Eigen::Vector2d &v0 = inPoly[i];
		Eigen::Vector2d &v1 = inPoly[(i + 1) % inPoly.size()];
		Point p0(v0[0], v0[1], 0), p1(v1[0], v1[1], 0);
		segList_.emplace_back(Segment(p0, p1));
	}
	aabbTree_ = new SegmentTree(segList_.begin(), segList_.end());
	aabbTree_->accelerate_distance_queries();
}

void PolySplit::InitializePoly(std::vector<Eigen::Vector2d> &inPoly)
{
	Vs_.clear();
	Vs_.reserve(inPoly.size());
	Es_.clear();
	Es_.reserve(inPoly.size() * 3);
	Fs_.clear();
	for (int i = 0; i < inPoly.size(); ++i)
	{
		PolyVertex pv;
		pv.id = i;
		pv.pos = inPoly[i];
		pv.neighbor_pfs.emplace_back(0);
		Vs_.emplace_back(pv);
	}
	PolyFace pf;
	pf.id = 0;
	//pf.pvs.reserve(inPoly.size());
	for (int i = 0; i < inPoly.size(); ++i)
	{
		pf.pvs.emplace_back(i);
		pf.pvs_new_add.emplace_back(false);
	}
	Fs_.emplace_back(pf);
	/*for (int i = 0; i < inPoly.size(); ++i)
	{
		PolyEdge pe;
		pe.id = i;
		pe.v0 = i;
		pe.v1 = (i + 1) % inPoly.size();
		pe.isBoundary = true;
		Es_.emplace_back(pe);

		Vs_[i].neighbor_pes.emplace_back(i);
		Vs_[(i + 1) % inPoly.size()].neighbor_pes.emplace_back(i);
	}*/

	sourcePolyArea_ = std::abs(ComputePolyArea(inPoly));
	diagLength_ = ComputePolyDiagLength(inPoly);
}

void PolySplit::InitializePoly(std::vector<Eigen::Vector2d> &inPoly, std::vector<bool> &inPolyAddNew)
{
	Vs_.clear();
	Vs_.reserve(inPoly.size());
	Es_.clear();
	Es_.reserve(inPoly.size() * 3);
	Fs_.clear();
	initialPoly_.reserve(inPoly.size());
	for (int i = 0; i < inPoly.size(); ++i)
	{
		PolyVertex pv;
		pv.id = i;
		pv.pos = inPoly[i];
		pv.neighbor_pfs.emplace_back(0);
		Vs_.emplace_back(pv);

		initialPoly_.emplace_back(i);
	}

	PolyFace pf;
	pf.id = 0;
	for (int i = 0; i < inPoly.size(); ++i)
	{
		pf.pvs.emplace_back(i);
		pf.pvs_new_add.emplace_back(inPolyAddNew[i]);
	}
	Fs_.emplace_back(pf);

	BuildAABBTree(inPoly);
	sourcePolyArea_ = std::abs(ComputePolyArea(inPoly));
	diagLength_ = ComputePolyDiagLength(inPoly);
}

double PolySplit::ComputePolyDiagLength(std::vector<Eigen::Vector2d> &poly)
{
	Eigen::Vector2d minP(1E20, 1E20), maxP(-1E20, -1E20);
	for (int i = 0; i < poly.size(); ++i)
	{
		if (poly[i][0] < minP[0])
			minP[0] = poly[i][0];
		if (poly[i][1] < minP[1])
			minP[1] = poly[i][1];

		if (poly[i][0] > maxP[0])
			maxP[0] = poly[i][0];
		if (poly[i][1] > maxP[1])
			maxP[1] = poly[i][1];
	}

	return (maxP - minP).norm();
}

double PolySplit::ComputePolyArea(std::vector<Eigen::Vector2d> &polyPos)
{
	double sumArea = 0;
	for (int i = 0; i < polyPos.size(); ++i)
	{
		Eigen::Vector2d vPos0 = polyPos[i], vPos1 = polyPos[(i + 1) % polyPos.size()];
		sumArea += (vPos0[0] + vPos1[0]) * (vPos1[1] - vPos0[1]);
	}
	return sumArea * 0.5;
}

double PolySplit::ComputePolyArea(std::vector<uint32_t> &polyVs)
{
	double sumArea = 0;
	for (int i = 0; i < polyVs.size(); ++i)
	{
		Eigen::Vector2d vPos0 = Vs_[polyVs[i]].pos, vPos1 = Vs_[polyVs[(i + 1) % polyVs.size()]].pos;
		sumArea += (vPos0[0] + vPos1[0]) * (vPos1[1] - vPos0[1]);
	}
	return sumArea * 0.5;
}

bool PolySplit::JudgePolyCCW(std::vector<uint32_t> &polyVs)
{
	double area = ComputePolyArea(polyVs);
	if (area > 0)
		return true;
	else
		return false;
}

bool PolySplit::JudgePolyCCW(std::vector<Eigen::Vector2d> &polyVPos)
{
	double area = ComputePolyArea(polyVPos);
	if (area > 0)
		return true;
	else
		return false;
}

void PolySplit::GetConcaveVertices(std::vector<uint32_t> &polyVs, std::vector<uint32_t> &ccVsIds)	//ccVsIds保存的是polyVs中的下标
{
	ccVsIds.clear();
	for (int i = 0; i < polyVs.size(); ++i)
	{
		int v0Id = polyVs[mmod(i-1, (int)polyVs.size())/*(i - 1) % polyVs.size()*/], v1Id = polyVs[i], v2Id = polyVs[(i + 1) % polyVs.size()];
		Eigen::Vector2d &pos0 = Vs_[v0Id].pos, &pos1 = Vs_[v1Id].pos, &pos2 = Vs_[v2Id].pos;
		double tempValue = (pos2[0] - pos1[0]) * (pos0[1] - pos1[1]) - (pos0[0] - pos1[0]) * (pos2[1] - pos1[1]);
		if (tempValue < 0)
			ccVsIds.emplace_back(i);
	}
}

void PolySplit::GetConcaveVerticesLarger(std::vector<uint32_t> &polyVs, std::vector<uint32_t> &ccVsIds, std::vector<bool> &polyAddNew)
{
	ccVsIds.clear();
	for (int i = 0; i < polyVs.size(); ++i)
	{
		if (polyAddNew[i])
			continue;

		int v0Id = polyVs[mmod(i - 1, (int)polyVs.size())/*(i - 1) % polyVs.size()*/], v1Id = polyVs[i], v2Id = polyVs[(i + 1) % polyVs.size()];
		Eigen::Vector2d &pos0 = Vs_[v0Id].pos, &pos1 = Vs_[v1Id].pos, &pos2 = Vs_[v2Id].pos;
		double tempValue = (pos2[0] - pos1[0]) * (pos0[1] - pos1[1]) - (pos0[0] - pos1[0]) * (pos2[1] - pos1[1]);
		double cosValue = (pos0 - pos1).normalized().dot((pos2 - pos1).normalized());
		double cosLarger = std::cos(LARGER_ANGLE);
		if (tempValue < 0 && cosValue > cosLarger)
			ccVsIds.emplace_back(i);
	}
}

void PolySplit::GetConcaveVerticesLarger(std::vector<uint32_t> &polyVs, std::vector<uint32_t> &ccVsIds)	//ccVsIds保存的是polyVs中的下标
{
	ccVsIds.clear();
	for (int i = 0; i < polyVs.size(); ++i)
	{
		int v0Id = polyVs[mmod(i - 1, (int)polyVs.size())/*(i - 1) % polyVs.size()*/], v1Id = polyVs[i], v2Id = polyVs[(i + 1) % polyVs.size()];
		Eigen::Vector2d &pos0 = Vs_[v0Id].pos, &pos1 = Vs_[v1Id].pos, &pos2 = Vs_[v2Id].pos;
		double tempValue = (pos2[0] - pos1[0]) * (pos0[1] - pos1[1]) - (pos0[0] - pos1[0]) * (pos2[1] - pos1[1]);
		double cosValue = (pos0 - pos1).normalized().dot((pos2 - pos1).normalized());
		double cosLarger = std::cos(LARGER_ANGLE);
		if (tempValue < 0 && cosValue > cosLarger)
			ccVsIds.emplace_back(i);
	}
}

bool PolySplit::CheckPointInsidePoly(std::vector<uint32_t> &poly, Eigen::Vector2d &midPos)
{
	std::vector<Eigen::Vector2d> polyVec;
	polyVec.reserve(poly.size());
	for (int i = 0; i < poly.size(); ++i)
	{
		polyVec.emplace_back(Vs_[poly[i]].pos);
	}

	if (CheckPolyProperty::CheckPointInsidePoly(polyVec, midPos))
		return true;
	else
		return false;
}

bool PolySplit::CheckPointOnBoundary(Eigen::Vector2d &pos)
{
	Point tPoint(pos[0], pos[1], 0);
	Point bPoint = aabbTree_->closest_point(tPoint);

	double dis = (Eigen::Vector2d(bPoint[0], bPoint[1]) - pos).norm();

	if (dis < 1E-5)
		return true;
	else
		return false;
}

int PolySplit::CheckPointPolyState(std::vector<uint32_t> &poly, Eigen::Vector2d &pos)
{
	std::vector<Eigen::Vector2d> polyVec;
	polyVec.reserve(poly.size());
	for (int i = 0; i < poly.size(); ++i)
	{
		polyVec.emplace_back(Vs_[poly[i]].pos);
	}

	int resu = CheckPolyProperty::CheckPointPolyState(polyVec, pos);
	if (resu == 0)
		return 0;
	else if (resu == 1)
		return 1;
	else
		return 2;
}

bool PolySplit::JudgeIfEdgeIntersectionOrOutside(std::vector<uint32_t> &polyVs, uint32_t v0, uint32_t v1)
{
	for (int i = 0; i < polyVs.size(); ++i)
	{
		if (JudgeIfEdgeIntersection(polyVs[i], polyVs[(i + 1) % polyVs.size()], v0, v1))
			return true;
	}
	
	Eigen::Vector2d midP = (Vs_[v0].pos + Vs_[v1].pos) / 2;
	std::vector<Eigen::Vector2d> polys;
	for (int i = 0; i < polyVs.size(); ++i)
	{
		polys.emplace_back(Vs_[polyVs[i]].pos);
	}

	if (!CheckPolyProperty::CheckPointInsidePoly(polys, midP))
		return true;
	else
		return false;
}

bool PolySplit::JudgeIfEdgeIntersection(Eigen::Vector2d &pos0, Eigen::Vector2d &pos1, Eigen::Vector2d &pos2, Eigen::Vector2d &pos3)
{
	double d0 = Cross2(pos1 - pos0, pos2 - pos0), d1 = Cross2(pos1 - pos0, pos3 - pos0), d2 = Cross2(pos3 - pos2, pos0 - pos2), d3 = Cross2(pos3 - pos2, pos1 - pos2);

	if (d0*d1 < 0 && d2*d3 < 0)
		return true;
	else
		return false;
}

bool PolySplit::JudgeIfOnSameLine(std::vector<bool> &polyAddNew, uint32_t v0Id, uint32_t v1Id)
{
	if (v0Id == v1Id)
		return false;

	bool hasOldV = false;
	int i = (v0Id+1)%polyAddNew.size();
	while (i != v1Id)
	{
		if (!polyAddNew[i])
		{
			hasOldV = true;
			break;
		}
		i = (i + 1) % polyAddNew.size();
	}
	if (!hasOldV)
		return true;

	hasOldV = false;
	i = (v1Id + 1) % polyAddNew.size();
	while (i != v0Id)
	{
		if (!polyAddNew[i])
		{
			hasOldV = true;
			break;
		}
		i = (i + 1) % polyAddNew.size();
	}
	if (!hasOldV)
		return true;

	return false;
}

bool PolySplit::JudgeIfOnSameLine(std::vector<uint32_t> &poly, uint32_t currId, Eigen::Vector2d &interPos)
{
	Eigen::Vector2d &pos0 = Vs_[poly[currId]].pos, &pos1 = Vs_[poly[mmod((int)currId - 1, (int)poly.size())]].pos, &pos2 = Vs_[poly[(currId + 1) % poly.size()]].pos;

	double cosValue0 = (interPos - pos0).normalized().dot((pos1 - pos0).normalized());
	double cosValue1 = (interPos - pos0).normalized().dot((pos2 - pos0).normalized());

	if (cosValue0 > std::cos(3.0 * PI / 180.0) || cosValue1 > std::cos(3.0*PI / 180.0))
		return true;

	return false;
}

bool PolySplit::JudgeIfEdgeIntersection(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3)
{
	if ((v0 == v2 && v1 == v3) || (v0 == v3 && v1 == v2))
		return true;
	
	if ((v0 == v2 && v1 != v3) || (v0 == v3 && v1 != v2) || (v1 == v2 && v0 != v3) || (v1 == v3 && v0 != v2))
	{
		return false;
	}

	Eigen::Vector2d &pos0 = Vs_[v0].pos, &pos1 = Vs_[v1].pos, &pos2 = Vs_[v2].pos, &pos3 = Vs_[v3].pos;
	double d0 = Cross2(pos1 - pos0, pos2 - pos0), d1 = Cross2(pos1 - pos0, pos3 - pos0), d2 = Cross2(pos3 - pos2, pos0 - pos2), d3 = Cross2(pos3 - pos2, pos1 - pos2);

	if (d0*d1 < 0 && d2*d3 < 0)
		return true;
	else
		return false;
}

bool PolySplit::JudgeIfEdgeIntersectionLarger(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3)
{
	if ((v0 == v2 && v1 == v3) || (v0 == v3 && v1 == v2))
		return true;

	if ((v0 == v2 && v1 != v3) || (v0 == v3 && v1 != v2) || (v1 == v2 && v0 != v3) || (v1 == v3 && v0 != v2))
	{
		return false;
	}

	Eigen::Vector2d &pos0 = Vs_[v0].pos, &pos1 = Vs_[v1].pos, &pos2 = Vs_[v2].pos, &pos3 = Vs_[v3].pos;
	double d0 = Cross2(pos1 - pos0, pos2 - pos0), d1 = Cross2(pos1 - pos0, pos3 - pos0), d2 = Cross2(pos3 - pos2, pos0 - pos2), d3 = Cross2(pos3 - pos2, pos1 - pos2);

	if (d0*d1 < 0 && d2*d3 < 0)
		return true;
	else
		return false;
}

void PolySplit::LineSample(uint32_t v0, uint32_t v1, std::vector<Eigen::Vector2d> &outSamples)
{
	Eigen::Vector2d &pos0 = Vs_[v0].pos, &pos1 = Vs_[v1].pos;
	double sl = SAMPLE_LENGTH;
	int num = (pos1 - pos0).norm() / sl;
	if (num < 10)
	{
		num = 10;
		sl = (pos1 - pos0).norm() / num;
	}
	Eigen::Vector2d unitVec = (pos1 - pos0).normalized();

	outSamples.clear();
	outSamples.reserve(num + 2);
	outSamples.emplace_back(pos0);
	for (int i = 0; i < num; ++i)
	{
		outSamples.emplace_back(pos0 + unitVec * (i + 1) * sl);
	}
	outSamples.emplace_back(pos1);
}

double PolySplit::LineFieldMatchAngle(uint32_t v0, uint32_t v1)
{
	std::vector<Eigen::Vector2d> outSamples;
	LineSample(v0, v1, outSamples);

	double aveAngle = 0;
	ANNpoint ap = annAllocPt(3);
	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
	for (int i = 0; i < outSamples.size(); ++i)
	{
		ap[0] = outSamples[i][0]; ap[1] = outSamples[i][1]; ap[2] = 0;
		fieldKdtree_->annkSearch(ap, 1, nnIdx, dists);
		const Eigen::Vector3d &currStanField = (*fields_)[nnIdx[0]];

		aveAngle += FieldMatchAngle(v0, v1, Eigen::Vector2d(currStanField[0], currStanField[1]));
	}
	aveAngle /= outSamples.size();
	annDeallocPt(ap);
	delete[] nnIdx; delete[] dists;

	return aveAngle;
}

double PolySplit::FieldMatchAngle(uint32_t v0, uint32_t v1, const Eigen::Vector2d &field)
{
	Eigen::Vector2d unitVec = (Vs_[v1].pos - Vs_[v0].pos).normalized();
	double angle = std::acos(field.normalized().dot(unitVec));
	angle = mmod(angle, PI/2);
	if (angle > PI / 4)
		return PI / 2 - angle;

	return angle;
}

void PolySplit::SplitPolyIntoTwo(std::vector<uint32_t> &polys, std::vector<bool> &polyAddNew, uint32_t v0Id, uint32_t v1Id, std::vector<uint32_t> &newPoly0, std::vector<uint32_t> &newPoly1, std::vector<bool> &newPAN0, std::vector<bool> &newPAN1)
{
	newPoly0.clear();
	newPoly1.clear();
	newPoly0.reserve(polys.size());
	newPoly1.reserve(polys.size());
	if (v1Id > v0Id)
	{
		for (int i = v0Id; i <= v1Id; ++i)
		{
			newPoly0.emplace_back(polys[i]);
			newPAN0.emplace_back(polyAddNew[i]);
		}
		for (int i = v1Id; i <= v0Id + polys.size(); ++i)
		{
			newPoly1.emplace_back(polys[i%polys.size()]);
			newPAN1.emplace_back(polyAddNew[i%polyAddNew.size()]);
		}
	}
	else
	{
		for (int i = v1Id; i <= v0Id; ++i)
		{
			newPoly1.emplace_back(polys[i]);
			newPAN1.emplace_back(polyAddNew[i]);
		}
		for (int i = v0Id; i <= v1Id + polys.size(); ++i)
		{
			newPoly0.emplace_back(polys[i%polys.size()]);
			newPAN0.emplace_back(polyAddNew[i%polyAddNew.size()]);
		}
	}
}

void PolySplit::SplitPolyIntoTwo(std::vector<uint32_t> &polys, uint32_t v0Id, uint32_t v1Id, std::vector<uint32_t> &newPoly0, std::vector<uint32_t> &newPoly1)
{
	/*if (v0Id > v1Id)
	{
		std::swap(v0Id, v1Id);
	}*/

	newPoly0.clear();
	newPoly1.clear();
	newPoly0.reserve(polys.size());
	newPoly1.reserve(polys.size());
	if (v1Id > v0Id)
	{
		for (int i = v0Id; i <= v1Id; ++i)
		{
			newPoly0.emplace_back(polys[i]);
		}
		for (int i = v1Id; i <= v0Id + polys.size(); ++i)
		{
			newPoly1.emplace_back(polys[i%polys.size()]);
		}
	}
	else
	{
		for (int i = v1Id; i <= v0Id; ++i)
		{
			newPoly1.emplace_back(polys[i]);
		}
		for (int i = v0Id; i <= v1Id + polys.size(); ++i)
		{
			newPoly0.emplace_back(polys[i%polys.size()]);
		}
	}
}

double PolySplit::ComputePolyAngle(std::vector<uint32_t> &polys, int vId)
{
	Eigen::Vector2d &pPre = Vs_[polys[mmod(vId-1, (int)polys.size())/*(vId - 1) % polys.size()*/]].pos, &pNow = Vs_[polys[vId]].pos, &pNext = Vs_[polys[(vId + 1) % polys.size()]].pos;
	Eigen::Vector2d vec0 = (pPre - pNow).normalized(), vec1 = (pNext - pNow).normalized();
	double angle = std::acos(vec0.dot(vec1));
	double sinV = Cross2(vec0, vec1);
	if (sinV > 0)
		angle = 2 * PI - angle;

	return angle;
}

double PolySplit::ComputePolyAngle(std::vector<Eigen::Vector2d> &polys, int vId)
{
	Eigen::Vector2d &pPre = polys[mmod(vId - 1, (int)polys.size())], &pNow = polys[vId], &pNext = polys[(vId + 1) % polys.size()];
	Eigen::Vector2d vec0 = (pPre - pNow).normalized(), vec1 = (pNext - pNow).normalized();
	double angle = std::acos(vec0.dot(vec1));
	double sinV = Cross2(vec0, vec1);
	if (sinV > 0)
		angle = 2 * PI - angle;

	return angle;
}

void PolySplit::SetScale(double radio, Eigen::Vector2d &minP)
{
	for (int i = 0; i < Vs_.size(); ++i)
	{
		Vs_[i].pos = minP + (Vs_[i].pos - minP)*radio;
	}
}

void PolySplit::BuildTopologyInfoFromVsFs()
{
	if (Vs_.empty() || Fs_.empty())
		return;
	Es_.clear();

	typedef std::tuple<uint32_t, uint32_t, uint32_t> TupleId3; //第一个顶点的id，第二个顶点的id（第一个小于第二个），f的编号
	std::vector<TupleId3> wholeInfoVec;
	uint32_t v0, v1;
	for (int i = 0; i < Fs_.size(); ++i)
	{
		std::vector<uint32_t> pvs;
		SetListToVec(pvs, Fs_[i].pvs);
		for (int j = 0; j < pvs.size(); ++j)
		{
			v0 = pvs[j]; v1 = pvs[(j + 1) % pvs.size()];
			if (v0 > v1)
				std::swap(v0, v1);
			wholeInfoVec.emplace_back(TupleId3(v0, v1, i));
		}
	}
	std::sort(wholeInfoVec.begin(), wholeInfoVec.end());

	int currNum = 0, currEId = 0;
	while (currNum < wholeInfoVec.size())
	{
		uint32_t v00 = std::get<0>(wholeInfoVec[currNum]), v01 = std::get<1>(wholeInfoVec[currNum]), f0 = std::get<2>(wholeInfoVec[currNum]), v10 = (uint32_t)-1, v11 = (uint32_t)-1, f1 = (uint32_t)-1;
		if (currNum != wholeInfoVec.size() - 1)
		{
			v10 = std::get<0>(wholeInfoVec[currNum + 1]);
			v11 = std::get<1>(wholeInfoVec[currNum + 1]);
			f1 = std::get<2>(wholeInfoVec[currNum + 1]);
		}

		if (currNum == wholeInfoVec.size() - 1 || !(v00 == v10 && v01 == v11))
		{
			PolyEdge pe;
			pe.id = currEId++;
			pe.isBoundary = true;
			pe.neighbor_f0 = f0;
			pe.neighbor_f1 = (uint32_t)-1;
			pe.v0 = v00;
			pe.v1 = v01;
			Es_.emplace_back(pe);

			Vs_[v00].neighbor_pes.emplace_back(currEId - 1);
			//Vs_[v00].neighbor_pfs.emplace_back(f0);
			Vs_[v01].neighbor_pes.emplace_back(currEId - 1);
			//Vs_[v01].neighbor_pfs.emplace_back(f0);

			//Fs_[f0].isBoundary = true;
			Fs_[f0].pes.emplace_back(currEId - 1);
			//Fs_[f0].pvs.emplace_back(v00);
			//Fs_[f0].pvs.emplace_back(v01);
			++currNum;
		}
		else
		{
			PolyEdge pe;
			pe.id = currEId++;
			pe.isBoundary = false;
			pe.neighbor_f0 = f0;
			pe.neighbor_f1 = f1;
			pe.v0 = v00;
			pe.v1 = v01;
			Es_.emplace_back(pe);

			Vs_[v00].neighbor_pes.emplace_back(currEId - 1);
			//Vs_[v00].neighbor_pfs.emplace_back(f0);
			//Vs_[v00].neighbor_pfs.emplace_back(f1);
			Vs_[v01].neighbor_pes.emplace_back(currEId - 1);
			//Vs_[v01].neighbor_pfs.emplace_back(f0);
			//Vs_[v00].neighbor_pfs.emplace_back(f1);

			Fs_[f0].pes.emplace_back(currEId - 1);
			//Fs_[f0].pvs.emplace_back(v00);
			//Fs_[f0].pvs.emplace_back(v01);
			Fs_[f1].pes.emplace_back(currEId - 1);
			//Fs_[f1].pvs.emplace_back(v00);
			//Fs_[f1].pvs.emplace_back(v01);

			Fs_[f0].neighbor_pfs.emplace_back(f1);
			Fs_[f1].neighbor_pfs.emplace_back(f0);
			currNum += 2;
		}
	}

	for (int i = 0; i < Fs_.size(); ++i)
	{
		std::vector<uint32_t> &fes = Fs_[i].pes;
		for (int j = 0; j < fes.size(); ++j)
		{
			if (Es_[fes[j]].isBoundary)
			{
				Fs_[i].isBoundary = true;
				break;
			}
			Fs_[i].isBoundary = false;
		}
	}

	for (int i = 0; i < Vs_.size(); ++i)
	{
		std::vector<uint32_t> &ves = Vs_[i].neighbor_pes;
		for (int j = 0; j < ves.size(); ++j)
		{
			if (Es_[ves[j]].isBoundary)
			{
				Vs_[i].isBoundary = true;
				break;
			}
			Vs_[i].isBoundary = false;
		}
	}
}

void PolySplit::BuildEFsInfo()
{
	//if (polysVec_.empty() || Vs_.empty())
	//	return;
	//Es_.clear();
	//Fs_.clear();
	//Fs_.reserve(polysVec_.size());

	//for (int i = 0; i < polysVec_.size(); ++i)
	//{
	//	PolyFace pf;
	//	pf.id = i;
	//	pf.isBoundary = false;
	//	Fs_.emplace_back(pf);
	//}

	//typedef std::tuple<uint32_t, uint32_t, uint32_t> TupleId3; //第一个顶点的id，第二个顶点的id（第一个小于第二个），f的编号。
	//std::vector<TupleId3> wholeInfoVec;
	//uint32_t v0, v1;
	//for (int i = 0; i < polysVec_.size(); ++i)
	//{
	//	for (int j = 0; j < polysVec_[i].size(); ++j)
	//	{
	//		v0 = polysVec_[i][j];
	//		v1 = polysVec_[i][(j + 1) % polysVec_[i].size()];
	//		if (v0 > v1)
	//			std::swap(v0, v1);

	//		wholeInfoVec.emplace_back(TupleId3(v0, v1, i));
	//	}
	//}

	//std::sort(wholeInfoVec.begin(), wholeInfoVec.end());
	//int currNum = 0;
	//int currEId = 0;
	//while (currNum < wholeInfoVec.size())
	//{
	//	uint32_t v00 = std::get<0>(wholeInfoVec[currNum]), v01 = std::get<1>(wholeInfoVec[currNum]), f0 = std::get<2>(wholeInfoVec[currNum]), v10 = (uint32_t)-1, v11 = (uint32_t)-1, f1 = (uint32_t)-1;
	//	if (currNum != wholeInfoVec.size() - 1)
	//	{
	//		v10 = std::get<0>(wholeInfoVec[currNum + 1]);
	//		v11 = std::get<1>(wholeInfoVec[currNum + 1]);
	//		f1 = std::get<2>(wholeInfoVec[currNum + 1]);
	//	}

	//	if (currNum == wholeInfoVec.size() - 1 || !(v00 == v10 && v01 == v11))
	//	{
	//		PolyEdge pe;
	//		pe.id = currEId++;
	//		pe.isBoundary = true;
	//		pe.neighbor_f0 = f0;
	//		pe.neighbor_f1 = (uint32_t)-1;
	//		pe.v0 = v00;
	//		pe.v1 = v01;
	//		Es_.emplace_back(pe);

	//		Vs_[v00].neighbor_pes.emplace_back(currEId - 1);
	//		Vs_[v00].neighbor_pfs.emplace_back(f0);
	//		Vs_[v01].neighbor_pes.emplace_back(currEId - 1);
	//		Vs_[v01].neighbor_pfs.emplace_back(f0);

	//		Fs_[f0].isBoundary = true;
	//		Fs_[f0].pes.emplace_back(currEId - 1);
	//		Fs_[f0].pvs.emplace_back(v00);
	//		Fs_[f0].pvs.emplace_back(v01);
	//		++currNum;
	//	}
	//	else
	//	{
	//		PolyEdge pe;
	//		pe.id = currEId++;
	//		pe.isBoundary = false;
	//		pe.neighbor_f0 = f0;
	//		pe.neighbor_f1 = f1;
	//		pe.v0 = v00;
	//		pe.v1 = v01;
	//		Es_.emplace_back(pe);

	//		Vs_[v00].neighbor_pes.emplace_back(currEId - 1);
	//		Vs_[v00].neighbor_pfs.emplace_back(f0);
	//		Vs_[v00].neighbor_pfs.emplace_back(f1);
	//		Vs_[v01].neighbor_pes.emplace_back(currEId - 1);
	//		Vs_[v01].neighbor_pfs.emplace_back(f0);
	//		Vs_[v00].neighbor_pfs.emplace_back(f1);

	//		Fs_[f0].pes.emplace_back(currEId - 1);
	//		Fs_[f0].pvs.emplace_back(v00);
	//		Fs_[f0].pvs.emplace_back(v01);
	//		Fs_[f1].pes.emplace_back(currEId - 1);
	//		Fs_[f1].pvs.emplace_back(v00);
	//		Fs_[f1].pvs.emplace_back(v01);

	//		Fs_[f0].neighbor_pfs.emplace_back(f1);
	//		Fs_[f1].neighbor_pfs.emplace_back(f0);
	//		currNum += 2;
	//	}
	//}
	
}

void PolySplit::TranslateToEigenPoly(std::vector<std::vector<Eigen::Vector2d>> &polys)
{
	/*polys.clear();
	std::vector<Eigen::Vector2d> currPoly;
	for (int i = 0; i < polysVec_.size(); ++i)
	{
		currPoly.clear();
		for (int j = 0; j < polysVec_[i].size(); ++j)
		{
			currPoly.emplace_back(Vs_[polysVec_[i][j]].pos);
		}
		polys.emplace_back(currPoly);
	}*/
}

void PolySplit::OutputPolysWithoutSelfIntersection(std::vector<Eigen::Vector2d> &inPoly, std::vector<bool> &inPolyAddNew)
{
	std::vector<uint32_t> pvs;
	std::vector<bool> pvsAddNew = inPolyAddNew;
	for (int i = 0; i < inPoly.size(); ++i)
	{
		pvs.emplace_back(i);
	}

	std::vector<uint32_t> oldVIds;
	oldVIds.reserve(pvs.size());
	for (int j = 0; j < pvsAddNew.size(); ++j)
	{
		if (!pvsAddNew[j])
			oldVIds.emplace_back(j);
	}
	if (oldVIds.empty())
		return;

	//std::vector<bool> backupPvsAddNew = pvsAddNew;
	for (int j = 0; j < oldVIds.size(); ++j)
	{
		for (int k = j + 2; k < oldVIds.size(); ++k)
		{
			uint32_t idJ = oldVIds[j], idJ1 = oldVIds[(j + 1) % oldVIds.size()], idK = oldVIds[k], idK1 = oldVIds[(k + 1) % oldVIds.size()];
			if (JudgeIfEdgeIntersection(pvs[idJ], pvs[idJ1], pvs[idK], pvs[idK1]))
			{
				int currId = idJ;
				while (currId != idJ1)
				{
					pvsAddNew[currId] = false;
					currId = (currId + 1) % pvs.size();
				}

				currId = idK;
				while (currId != idK1)
				{
					pvsAddNew[currId] = false;
					currId = (currId + 1) % pvs.size();
				}
			}
		}
	}
}

void PolySplit::OutputPolysWithoutSelfIntersection()
{
	for (int i = 0; i < Fs_.size(); ++i)
	{
		std::vector<uint32_t> pvs;
		std::vector<bool> pvsAddNew;
		SetListToVec(pvs, Fs_[i].pvs);
		SetListToVec(pvsAddNew, Fs_[i].pvs_new_add);

		std::vector<uint32_t> oldVIds;
		oldVIds.reserve(pvs.size());
		for (int j = 0; j < pvsAddNew.size(); ++j)
		{
			if (!pvsAddNew[j])
				oldVIds.emplace_back(j);
		}
		if (oldVIds.empty())
			return;

		//std::vector<bool> backupPvsAddNew = pvsAddNew;
		for (int j = 0; j < oldVIds.size(); ++j)
		{
			for (int k = j + 2; k < oldVIds.size(); ++k)
			{
				uint32_t idJ = oldVIds[j], idJ1 = oldVIds[(j+1)%oldVIds.size()], idK = oldVIds[k], idK1 = oldVIds[(k+1)%oldVIds.size()];
				if (JudgeIfEdgeIntersection(pvs[idJ], pvs[idJ1], pvs[idK], pvs[idK1]))
				{
					int currId = idJ;
					while (currId != idJ1)
					{
						pvsAddNew[currId] = false;
						currId = (currId + 1) % pvs.size();
					}

					currId = idK;
					while (currId != idK1)
					{
						pvsAddNew[currId] = false;
						currId = (currId + 1) % pvs.size();
					}
				}
			}
		}

		SetVecToList(pvsAddNew, Fs_[i].pvs_new_add);
	}
}

void PolySplit::TranslateToEigenPolyFromFs(std::vector<std::vector<Eigen::Vector2d>> &polys)
{
	polys.clear();
	std::vector<Eigen::Vector2d> currPoly;
	for (int i = 0; i < Fs_.size(); ++i)
	{
		currPoly.clear();
		std::list<uint32_t> &pvs = Fs_[i].pvs;
		std::list<bool> &pvsNewAdd = Fs_[i].pvs_new_add;
		auto itt = pvsNewAdd.begin();
		for (auto it = pvs.begin(); it != pvs.end(); ++it)
		{
			if (!(*itt))
				currPoly.emplace_back(Vs_[*it].pos);
			++itt;
		}
		polys.emplace_back(currPoly);
	}
}

bool PolySplit::CalcLinesIntersection2(Eigen::Vector2d &vPos0, Eigen::Vector2d &vPos1, Eigen::Vector2d &vPos2, Eigen::Vector2d &vec1, double &length0, double &length1, Eigen::Vector2d &interPos)
{
	Eigen::Vector2d vec0 = vPos1 - vPos0;

	if (std::abs(Cross2(vec0.normalized(), vec1.normalized())) < 1E-6)
		return false;

	double x0 = vPos0[0], y0 = vPos0[1], x1 = vPos2[0], y1 = vPos2[1], xvec0 = vec0[0], yvec0 = vec0[1], xvec1 = vec1[0], yvec1 = vec1[1];
	length1 = (x0 * yvec0 - y0 * xvec0 - x1 * yvec0 + y1 * xvec0) / (yvec0 * xvec1 - yvec1 * xvec0);
	length0 = (x1 * yvec1 - y1 * xvec1 - x0 * yvec1 + y0 * xvec1) / (xvec0 * yvec1 - yvec0 * xvec1);

	interPos = vPos0 + length0 * vec0;
	return true;
}

bool PolySplit::CalcLinesIntersection(Eigen::Vector2d &vPos0, Eigen::Vector2d &vPos1, Eigen::Vector2d &vPos2, Eigen::Vector2d &vPos3, double &length0, double &length1, Eigen::Vector2d &interPos)
{
	Eigen::Vector2d vec0 = vPos1 - vPos0;
	Eigen::Vector2d vec1 = vPos3 - vPos2;

	if (std::abs(Cross2(vec0.normalized(), vec1.normalized())) < 1E-6)
		return false;

	double x0 = vPos0[0], y0 = vPos0[1], x1 = vPos2[0], y1 = vPos2[1], xvec0 = vec0[0], yvec0 = vec0[1], xvec1 = vec1[0], yvec1 = vec1[1];
	length1 = (x0 * yvec0 - y0 * xvec0 - x1 * yvec0 + y1 * xvec0) / (yvec0 * xvec1 - yvec1 * xvec0);
	length0 = (x1 * yvec1 - y1 * xvec1 - x0 * yvec1 + y0 * xvec1) / (xvec0 * yvec1 - yvec0 * xvec1);

	interPos = vPos0 + length0 * vec0;
	return true;
}

bool PolySplit::GetPolyRayIntersection(std::vector<uint32_t> &poly, uint32_t v0Id, Eigen::Vector2d &vecDir, uint32_t &outV2Id, uint32_t &outV3Id, Eigen::Vector2d &interPos, std::vector<Eigen::Vector2d> &outPoly0, std::vector<Eigen::Vector2d> &outPoly1)
{
	double currlength = 1E20;
	bool isFoundGoodInter = false;
	uint32_t v2Id = (uint32_t)-1, v3Id = (uint32_t)-1;
	for (int i = 0; i < poly.size(); ++i)
	{
		uint32_t currV2 = poly[i], currV3 = poly[(i + 1) % poly.size()];
		if (v0Id == i || v0Id == (i + 1) % poly.size())
			continue;

		double length0 = 0, length1 = 0;
		Eigen::Vector2d currInterPos, currVecPos = Vs_[poly[v0Id]].pos + vecDir;
		bool isNotParallel = CalcLinesIntersection(Vs_[poly[v0Id]].pos, currVecPos, Vs_[currV2].pos, Vs_[currV3].pos, length0, length1, currInterPos);
		if (!isNotParallel || length0 >= 0 || length1 <= 0 || length1 >= 1)
			continue;

		if (std::abs(length0) < std::abs(currlength))
		{
			currlength = length0;
			v2Id = i;
			v3Id = (i + 1) % poly.size();
			isFoundGoodInter = true;
		}
	}

	if (!isFoundGoodInter)
		return false;
	else
	{
		outPoly0.clear();
		outPoly1.clear();
		outV2Id = v2Id;
		outV3Id = v3Id;
		interPos = Vs_[poly[v0Id]].pos + currlength * vecDir;

		/*auto it0 = std::find(poly.begin(), poly.end(), v0);
		int v0Id = std::distance(poly.begin(), it0);
		auto it2 = std::find(poly.begin(), poly.end(), v2);
		int v2Id = std::distance(poly.begin(), it2);
		auto it3 = std::find(poly.begin(), poly.end(), v3);
		int v3Id = std::distance(poly.begin(), it3);*/

		int i = v0Id;
		while (i != v2Id && i != v3Id)
		{
			outPoly0.emplace_back(Vs_[poly[i]].pos);
			i = (i + 1) % poly.size();
		}
		outPoly0.emplace_back(Vs_[poly[i]].pos);
		outPoly0.emplace_back(interPos);

		outPoly1.emplace_back(interPos);
		if (i == v2Id)
			i = v3Id;
		else
			i = v2Id;

		while (i != v0Id)
		{
			outPoly1.emplace_back(Vs_[poly[i]].pos);
			i = (i + 1) % poly.size();
		}
		outPoly1.emplace_back(Vs_[poly[i]].pos);
	}
	return true;
}

bool PolySplit::GetPolyRayIntersection(std::vector<uint32_t> &poly, uint32_t v0Id, uint32_t v1Id, uint32_t &outV2Id, uint32_t &outV3Id, Eigen::Vector2d &interPos, std::vector<Eigen::Vector2d> &outPoly0, std::vector<Eigen::Vector2d> &outPoly1)
{
	double currlength = 1E20;
	bool isFoundGoodInter = false;
	uint32_t v2Id = (uint32_t)-1, v3Id = (uint32_t)-1;
	for (int i = 0; i < poly.size(); ++i)
	{
		uint32_t currV2 = poly[i], currV3 = poly[(i + 1) % poly.size()];
		if (v0Id == i || v0Id == (i + 1) % poly.size() || v1Id == i || v1Id == (i + 1) % poly.size())
			continue;

		double length0 = 0, length1 = 0;
		Eigen::Vector2d currInterPos;
		bool isNotParallel = CalcLinesIntersection(Vs_[poly[v0Id]].pos, Vs_[poly[v1Id]].pos, Vs_[currV2].pos, Vs_[currV3].pos, length0, length1, currInterPos);
		if (!isNotParallel || length0 >= 0 || length1 <= 0 || length1 >= 1)
			continue;

		if (std::abs(length0) < std::abs(currlength))
		{
			currlength = length0;
			v2Id = i;
			v3Id = (i + 1) % poly.size();
			isFoundGoodInter = true;
		}
	}

	if (!isFoundGoodInter)
		return false;
	else
	{
		outPoly0.clear();
		outPoly1.clear();
		outV2Id = v2Id;
		outV3Id = v3Id;
		interPos = Vs_[poly[v0Id]].pos + currlength * (Vs_[poly[v1Id]].pos - Vs_[poly[v0Id]].pos);

		/*auto it0 = std::find(poly.begin(), poly.end(), v0);
		int v0Id = std::distance(poly.begin(), it0);
		auto it2 = std::find(poly.begin(), poly.end(), v2);
		int v2Id = std::distance(poly.begin(), it2);
		auto it3 = std::find(poly.begin(), poly.end(), v3);
		int v3Id = std::distance(poly.begin(), it3);*/

		int i = v0Id;
		while (i != v2Id && i != v3Id)
		{
			outPoly0.emplace_back(Vs_[poly[i]].pos);
			i = (i + 1) % poly.size();
		}
		outPoly0.emplace_back(Vs_[poly[i]].pos);
		outPoly0.emplace_back(interPos);

		outPoly1.emplace_back(interPos);
		if (i == v2Id)
			i = v3Id;
		else
			i = v2Id;

		while (i != v0Id)
		{
			outPoly1.emplace_back(Vs_[poly[i]].pos);
			i = (i + 1) % poly.size();
		}
		outPoly1.emplace_back(Vs_[poly[i]].pos);
	}
	return true;
}

void PolySplit::ComputePolyEdgesLength(std::vector<Eigen::Vector2d> &poly, std::vector<double> &edgesLength)
{
	edgesLength.clear();
	edgesLength.reserve(poly.size());
	for (int i = 0; i < poly.size(); ++i)
	{
		double currLength = (poly[(i + 1)%poly.size()] - poly[i]).norm();
		edgesLength.emplace_back(currLength);
	}
}

void PolySplit::ComputeShortestAndLongestEdge(std::vector<Eigen::Vector2d> &poly, std::vector<bool> &polyAddNew, double &minL, double &maxL)
{
	std::vector<Eigen::Vector2d> oldVertics;
	oldVertics.reserve(poly.size());

	for (int i = 0; i < polyAddNew.size(); ++i)
	{
		if (!polyAddNew[i])
			oldVertics.emplace_back(poly[i]);
	}

	if (oldVertics.empty())
		exit(-5347);
	std::vector<double> edgesLength;
	ComputePolyEdgesLength(oldVertics, edgesLength);
	if (edgesLength.empty())
		exit(-5348);

	minL = edgesLength[0];
	maxL = edgesLength[edgesLength.size() - 1];
}

void PolySplit::ComputeShortestAndLongestEdge(std::vector<Eigen::Vector2d> &poly, double &minL, double &maxL)
{
	std::vector<double> edgesLength;
	ComputePolyEdgesLength(poly, edgesLength);

	if (edgesLength.empty())
		exit(-3845);
	std::sort(edgesLength.begin(), edgesLength.end());
	minL = edgesLength[0];
	maxL = edgesLength[edgesLength.size() - 1];
}

bool PolySplit::ComputeVerticalIntersection(Eigen::Vector2d &v0, Eigen::Vector2d &v10, Eigen::Vector2d &v11, double &alpha, Eigen::Vector2d &interPos)
{
	Eigen::Vector2d vec1 = v11 - v10;
	Eigen::Vector2d vec0(vec1[1], -vec1[0]);
	vec0.normalize();

	double length1;
	CalcLinesIntersection2(v10, v11, v0, vec0, alpha, length1, interPos);
	if (alpha <= 0 || alpha >= 1)
		return false;

	return true;
}

bool PolySplit::JudgeIfVerticalIntersection(std::vector<uint32_t> &poly, uint32_t startVId, Eigen::Vector2d &interPos, uint32_t interV0Id, uint32_t interV1Id)
{
	for (int i = 0; i < poly.size(); ++i)
	{
		if (i == startVId || ((i + 1) % poly.size()) == startVId || (i == interV0Id && ((i + 1) % poly.size()) == interV1Id) || (i == interV1Id && ((i + 1) % poly.size()) == interV0Id))
			continue;

		Eigen::Vector2d startPos = Vs_[poly[startVId]].pos, interV0Pos = Vs_[poly[i]].pos, interV1Pos = Vs_[poly[(i + 1) % poly.size()]].pos;
		if (JudgeIfEdgeIntersection(startPos, interPos, interV0Pos, interV1Pos))
			return true;
	}
	return false;
}

void PolySplit::FindLinePolyNearIntersection(Eigen::Vector2d &vPos0, Eigen::Vector2d &vPos1, std::vector<Eigen::Vector2d> &poly, double &length0, double &length1, uint32_t &interId0, uint32_t &interId1, Eigen::Vector2d &interPos)
{
	length0 = 1E10; 
	length1 = 0;
	for (int i = 0; i < poly.size(); ++i)
	{
		double currLen0 = 0, currLen1 = 0;
		Eigen::Vector2d interPoss;
		if (CalcLinesIntersection(vPos0, vPos1, poly[i], poly[(i + 1) % poly.size()], currLen0, currLen1, interPoss))
		{
			if (currLen0 > 0 && currLen0 < length0 && currLen1<=1 && currLen1 >= 0)
			{
				length0 = currLen0;
				length1 = currLen1;
				interId0 = i; 
				interId1 = (i + 1) % poly.size();
				interPos = interPoss;
			}
		}
	}
}

void PolySplit::OutputInnerEdgesInfo(std::string &fileName, std::vector<uint32_t> &polyToTri, std::vector<Eigen::Vector2d> &sourcePoly)
{
	std::vector<uint32_t> vFlag(Vs_.size(), (uint32_t)-1);	//vFlag不为(uint32_t)-1表示这个点是某条内部线的端点
	int currId = 0;
	int innerENum = 0;
	for (int i = 0; i < Es_.size(); ++i)
	{
		if (Es_[i].isBoundary)
			continue;

		uint32_t v0 = Es_[i].v0, v1 = Es_[i].v1;
		
		if (vFlag[v0] == (uint32_t)-1)
			vFlag[v0] = currId++;
		if (vFlag[v1] == (uint32_t)-1)
			vFlag[v1] = currId++;

		++innerENum;
	}

	typedef std::tuple<uint32_t, uint32_t, double, Eigen::Vector2d, uint32_t, bool> Tuple5;	//最后一个true表示边界点,倒数第二个表示其原始id
	std::vector<Tuple5> outVsInfo;
	uint32_t sourcePointNum = sourcePoly.size();

	uint32_t boundVNum = 0, innerVNum = 0, currVId = 0;
	std::map<uint32_t, uint32_t> tempMap;
	for (int i = 0; i < Vs_.size(); ++i)
	{
		if (vFlag[i]==(uint32_t)-1)
			continue;
			
		//先考虑边界上的点，找它的局部位置
		//uint32_t currId = vFlag[i];
		if (Vs_[i].isBoundary)
		{	
			//如果这是一个和原始三角网格的顶点匹配的点
			if (i < sourcePointNum)
			{
				uint32_t v0 = polyToTri[i], v1 = polyToTri[(i + 1) % polyToTri.size()];
				double lambda = 0;
				outVsInfo.emplace_back( Tuple5(v0, v1, lambda, sourcePoly[i], i, true) );
				tempMap[i] = currVId++;
			}
			else  //否则找其和原始网格的交点
			{
				uint32_t anotherV = (uint32_t)-1;
				std::vector<uint32_t> &ves = Vs_[i].neighbor_pes;
				for (int j = 0; j < ves.size(); ++j)
				{
					if (!Es_[ves[j]].isBoundary)
					{
						if (Es_[ves[j]].v0 == i)
							anotherV = Es_[ves[j]].v1;
						else
							anotherV = Es_[ves[j]].v0;
						break;
					}
				}
				if (anotherV == (uint32_t)-1)
					exit(-35523);

				uint32_t currOutId0, currOutId1;
				double length0, length1;
				Eigen::Vector2d interPos;
				FindLinePolyNearIntersection(Vs_[anotherV].pos, Vs_[i].pos, sourcePoly, length0, length1, currOutId0, currOutId1, interPos);
				outVsInfo.emplace_back(Tuple5(polyToTri[currOutId0], polyToTri[currOutId1], length1, interPos, i, true));
				tempMap[i] = currVId++;
			}
			vFlag[i] = (uint32_t)-1;
		}
	}
	boundVNum = currVId;

	for (int i = 0; i < Vs_.size(); ++i)
	{
		if (vFlag[i] == (uint32_t)-1)
			continue;

		Eigen::Vector2d &currPos = Vs_[i].pos;

		for (int j = 0; j < outVsInfo.size(); ++j)
		{
			bool isBreak = false;
			for (int k = j + 1; k < outVsInfo.size(); ++k)
			{
				Eigen::Vector2d vPos0 = std::get<3>(outVsInfo[j]), vPos1 = std::get<3>(outVsInfo[k]);
				Eigen::Vector2d vVec(vPos0[1] - vPos1[1], vPos1[0] - vPos0[0]), interPos;
				vVec.normalize();
				double length0, length1;
				CalcLinesIntersection2(vPos0, vPos1, currPos, vVec, length0, length1, interPos);
				if (length0 < 1 && length0>0 && std::abs(length1) < 1E-8)
				{
					outVsInfo.emplace_back(Tuple5(j, k, length0, interPos, i, false));
					tempMap[i] = currVId++;
					isBreak = true;
					break;
				}
			}
			if (isBreak)
				break;
		}
	}
	innerVNum = currVId - boundVNum;

	std::ofstream ofs(fileName.c_str());
	ofs << boundVNum << std::endl;
	for (int i = 0; i < boundVNum; ++i)
	{
		ofs << std::get<0>(outVsInfo[i]) << " " << std::get<1>(outVsInfo[i]) << " " << std::get<2>(outVsInfo[i]) << std::endl;
	}
	ofs << innerVNum << std::endl;
	for (int i = 0; i < innerVNum; ++i)
	{
		ofs << std::get<0>(outVsInfo[i+boundVNum]) << " " << std::get<1>(outVsInfo[i+boundVNum]) << " " << std::get<2>(outVsInfo[i+boundVNum]) << std::endl;
	}
	ofs << innerENum << std::endl;
	for (int i = 0; i < Es_.size(); ++i)
	{
		if (Es_[i].isBoundary)
			continue;

		uint32_t v0 = Es_[i].v0, v1 = Es_[i].v1;
		ofs << tempMap[v0] << " " << tempMap[v1] << std::endl;
	}
}

void PolySplit::OutputInnerEdgesInfo(std::string &fileName)
{
	std::ofstream ofs(fileName.c_str());
	for (int i = 0; i < Es_.size(); ++i)
	{
		if (Es_[i].neighbor_f1 == (uint32_t)-1)
			continue;

		uint32_t v0 = Es_[i].v0, v1 = Es_[i].v1;
		ofs << Vs_[v0].pos[0] << " " << Vs_[v0].pos[1] << " " << Vs_[v1].pos[0] << " " << Vs_[v1].pos[1] << std::endl;
	}
}

void PolySplit::FindListInsertPos(std::list<uint32_t> &pvs, uint32_t v0, uint32_t v1, std::list<uint32_t>::iterator &it, int &diss)
{
	auto it0 = std::find(pvs.begin(), pvs.end(), v0), it1 = std::find(pvs.begin(), pvs.end(), v1);
	if (it0 == pvs.end() || it1 == pvs.end())
		exit(-908);

	int dis0 = std::distance(pvs.begin(), it0), dis1 = std::distance(pvs.begin(), it1);
	if (dis0 - dis1 == 1)
	{
		it = it0;
		diss = dis0;
	}
	else if (dis1 - dis0 == 1)
	{
		it = it1;
		diss = dis1;
	}
	else if ((dis0 == 0 && dis1 == pvs.size() - 1) || (dis1 == 0 && dis0 == pvs.size() - 1))
	{
		it = pvs.end();
		diss = pvs.size();
	}
	else
		exit(-909);
}

void PolySplit::SplitPolyIntoTwo(std::vector<uint32_t> &polys, uint32_t startVId, Eigen::Vector2d &interPos, uint32_t anV0Id, uint32_t anV1Id, std::vector<Eigen::Vector2d> &outPoly0, std::vector<Eigen::Vector2d> &outPoly1)
{
	outPoly0.clear();
	outPoly1.clear();

	/*auto it0 = std::find(poly.begin(), poly.end(), v0);
	int v0Id = std::distance(poly.begin(), it0);
	auto it2 = std::find(poly.begin(), poly.end(), v2);
	int v2Id = std::distance(poly.begin(), it2);
	auto it3 = std::find(poly.begin(), poly.end(), v3);
	int v3Id = std::distance(poly.begin(), it3);*/

	int i = startVId;
	while (i != anV0Id && i != anV1Id)
	{
		outPoly0.emplace_back(Vs_[polys[i]].pos);
		i = (i + 1) % polys.size();
	}
	outPoly0.emplace_back(Vs_[polys[i]].pos);
	outPoly0.emplace_back(interPos);

	outPoly1.emplace_back(interPos);
	if (i == anV0Id)
		i = anV1Id;
	else
		i = anV0Id;

	while (i != startVId)
	{
		outPoly1.emplace_back(Vs_[polys[i]].pos);
		i = (i + 1) % polys.size();
	}
	outPoly1.emplace_back(Vs_[polys[i]].pos);
}

int PolySplit::FindRealLargerOldVertexId(int currId, std::vector<bool> &polyAddNew)
{
	int i = currId;
	do
	{
		if (!polyAddNew[i] && i!=currId)
			return i;
		i = (i + 1) % polyAddNew.size();
	} while (i != currId);

	exit(-3454);
}

int PolySplit::FindLargerOldVertexId(int currId, std::vector<bool> &polyAddNew)
{
	int i = currId;
	do
	{
		if (!polyAddNew[i])
			return i;
		i = (i + 1) % polyAddNew.size();
	} while (i != currId);

	exit(-3452);
}

int PolySplit::FindRealSmallerOldVertexId(int currId, std::vector<bool> &polyAddNew)
{
	int i = currId;
	do
	{
		if (!polyAddNew[i] && i!=currId)
			return i;
		i = mmod(i - 1, (int)polyAddNew.size());
	} while (i != currId);

	exit(-3453);
}
int PolySplit::FindSmallerOldVertexId(int currId, std::vector<bool> &polyAddNew)
{
	int i = currId;
	do
	{
		if (!polyAddNew[i])
			return i;
		i = mmod(i - 1, (int)polyAddNew.size());
	} while (i != currId);

	exit(-3453);
}

void PolySplit::GetSubAddNew(int currId, int targetId, std::vector<bool> sourcePolyAddNew, std::vector<bool> &poly0AddNew, std::vector<bool> &poly1AddNew)
{

}

void PolySplit::SplitPoly(std::vector<Eigen::Vector2d> &inPoly, std::vector<bool> &inPolyAddNew, int polyNumThre, int resultNum)
{
	int numInEveryStep = int(std::pow((double)resultNum, 1.0 / 3) + 1);
	if (!JudgePolyCCW(inPoly))
	{
		std::reverse(inPoly.begin(), inPoly.end());
		std::reverse(inPolyAddNew.begin(), inPolyAddNew.end());
	}
	InitializePoly(inPoly, inPolyAddNew);

	currFaceId_.clear();

	std::queue<uint32_t> queueFaceId;
	queueFaceId.emplace(0);

	//第一步
	std::vector<uint32_t> ccVsIds;
	while (!queueFaceId.empty())
	{
		if (Fs_.size() >= polyNumThre)
		{
			BuildTopologyInfoFromVsFs();
			return;
		}
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35469);

		std::vector<TupleRankingLine> fieldRanking;
		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);

		std::vector<bool> concaveVFlag(currPolys.size(), false);
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			concaveVFlag[ccVsIds[i]] = true;
		}

		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			uint32_t v0 = currPolys[ccVsIds[i]];
#if SHORT_LENGTH_NOT_DO
			int iLargeId = FindRealLargerOldVertexId(ccVsIds[i], currPolyAddNew), iSmallId = FindRealSmallerOldVertexId(ccVsIds[i], currPolyAddNew);
			double length0 = (Vs_[currPolys[iLargeId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm(), length1 = (Vs_[currPolys[iSmallId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm();
			if (length0 < diagLength_ * SHORT_NOT_DO_LENGTH && length1 < diagLength_ * SHORT_NOT_DO_LENGTH)
				continue;
#endif
			for (int j = 0; j < currPolys.size(); ++j)
			{
				if (j == ccVsIds[i] || std::abs((int)(ccVsIds[i]) - (int)(j)) == 1 || currPolyAddNew[j])
					continue;

#if SHORT_LENGTH_NOT_DO
				int jLargeId = FindRealLargerOldVertexId(j, currPolyAddNew), jSmallId = FindRealSmallerOldVertexId(j, currPolyAddNew);
				double jlength0 = (Vs_[currPolys[jLargeId]].pos - Vs_[currPolys[j]].pos).norm(), jlength1 = (Vs_[currPolys[jSmallId]].pos - Vs_[currPolys[j]].pos).norm();
				if (jlength0 < diagLength_ * SHORT_NOT_DO_LENGTH && jlength1 < diagLength_ * SHORT_NOT_DO_LENGTH)
					continue;
#endif
				uint32_t v1 = currPolys[j];
				if (JudgeIfEdgeIntersectionOrOutside(currPolys, v0, v1))
					continue;

				if (JudgeIfOnSameLine(currPolyAddNew, ccVsIds[i], j))
					continue;

				double currFieldMatchAngle = LineFieldMatchAngle(v0, v1);
				bool isConcave = false;
				if (concaveVFlag[j])
				{
					currFieldMatchAngle -= 1E5;
					isConcave = true;
				}
				else
					isConcave = false;
				fieldRanking.emplace_back(TupleRankingLine(currFieldMatchAngle, v0, v1, ccVsIds[i], j, isConcave));
			}
		}

		if (fieldRanking.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		bool isFindGoodSplit = false;
		std::vector<uint32_t> newPoly0, newPoly1;
		std::sort(fieldRanking.begin(), fieldRanking.end());
		for (int i = 0; i < fieldRanking.size(); ++i)
		{
			uint32_t currSelect = i;
			int randomTillNum = std::min(int(fieldRanking.size()*2.0 / 3.0+1), numInEveryStep*2);
			if (randomTillNum == 0)
				randomTillNum = 1;

			if (currSelect < randomTillNum)
				currSelect = GetRandomNum(0, randomTillNum);
			double fieldValue = std::get<0>(fieldRanking[currSelect]);
			uint32_t v0 = std::get<1>(fieldRanking[currSelect]), v1 = std::get<2>(fieldRanking[currSelect]), v0Id = std::get<3>(fieldRanking[currSelect]), v1Id = std::get<4>(fieldRanking[currSelect]);
			bool isConcave = std::get<5>(fieldRanking[currSelect]);

			ANGLE_THRESHOLD = 75.0 * PI / 180.0;
			AREA_THRESHOLD = 0.003;
			/*if (isConcave)
			{
				ANGLE_THRESHOLD = 45.0 * PI / 180.0;
				AREA_THRESHOLD = 0.005;
			}
			else
			{
				ANGLE_THRESHOLD = 75.0 * PI / 180.0;
				AREA_THRESHOLD = 0.005;
			}*/

#if FOR_QUADS_GEN
			double fieldMatchV = fieldValue;
			if (fieldMatchV < 0)
				fieldMatchV += 1E5;
			if (fieldMatchV > FIELD_MATCHING_THRESHOLD)
				continue;
#else
			double fieldMatchV = fieldValue;
			if (fieldMatchV < 0)
				fieldMatchV += 1E5;
			if (fieldMatchV > FIELD_MATCHING_THRESHOLD_2)
				continue;
#endif
			double newLength = (Vs_[v0].pos - Vs_[v1].pos).norm();
#if OLD_LENGTH_THRE
			if (newLength < LENGTH_THRESHOLD * diagLength_)
				continue;
#else
			if (newLength < LENGTH_THRESHOLD_0 * diagLength_)
				continue;
#endif

			//SplitPolyIntoTwo(currPolys, v0Id, v1Id, newPoly0, newPoly1);
			std::vector<bool> newPAN0, newPAN1;
			SplitPolyIntoTwo(currPolys, currPolyAddNew, v0Id, v1Id, newPoly0, newPoly1, newPAN0, newPAN1);

			if (GetFalseNum(newPAN0) <= 3 || GetFalseNum(newPAN1) <= 3)
				continue;

			//是否不允许出现三角形？
			if (newPoly0.size() <= 3 || newPoly1.size() <= 3)
				continue;

#if !FOR_QUADS_GEN
			double poly0Area = std::abs(ComputePolyArea(newPoly0)), poly1Area = std::abs(ComputePolyArea(newPoly1));
			if (poly0Area < AREA_THRESHOLD * sourcePolyArea_ || poly1Area < AREA_THRESHOLD * sourcePolyArea_)
				continue;
#endif

			double angle0 = ComputePolyAngle(newPoly0, 0), angle1 = ComputePolyAngle(newPoly0, newPoly0.size() - 1), angle2 = ComputePolyAngle(newPoly1, 0), angle3 = ComputePolyAngle(newPoly1, newPoly1.size() - 1);
			if (angle0 < ANGLE_THRESHOLD || angle1 < ANGLE_THRESHOLD || angle2 < ANGLE_THRESHOLD || angle3 < ANGLE_THRESHOLD)
				continue;

			isFindGoodSplit = true;

#if OUTPUT_SPLIT_LINES
			splitEdges_.emplace_back(std::tuple<Eigen::Vector2d, Eigen::Vector2d>(Vs_[v0].pos, Vs_[v1].pos));
#endif

			//更新Vs_, Fs_拓扑信息
			Vs_[currPolys[v0Id]].neighbor_pfs.emplace_back(Fs_.size());
			Vs_[currPolys[v1Id]].neighbor_pfs.emplace_back(Fs_.size());
			Fs_[currFaceId].pvs.clear();
			Fs_[currFaceId].pvs_new_add.clear();

			int jj = v0Id;
			while (jj != v1Id)
			{
				Fs_[currFaceId].pvs.emplace_back(currPolys[jj]);
				if (jj == v0Id)
					Fs_[currFaceId].pvs_new_add.emplace_back(false);
				else
					Fs_[currFaceId].pvs_new_add.emplace_back(currPolyAddNew[jj]);
				jj = (jj + 1) % currPolys.size();
			}
			Fs_[currFaceId].pvs.emplace_back(currPolys[jj]);
			Fs_[currFaceId].pvs_new_add.emplace_back(false);

			for (int j = 0; j < newPoly1.size(); ++j)
			{
				if (newPoly1[j] != currPolys[v0Id] && newPoly1[j] != currPolys[v1Id])
				{
					std::vector<uint32_t> &vfs = Vs_[newPoly1[j]].neighbor_pfs;
					for (int k = 0; k < vfs.size(); ++k)
					{
						if (vfs[k] == currFaceId)
							vfs[k] = Fs_.size();
					}
				}
			}
			PolyFace pf;
			pf.id = Fs_.size();
			jj = v1Id;
			while (jj != v0Id)
			{
				pf.pvs.emplace_back(currPolys[jj]);
				if (jj == v1Id)
					pf.pvs_new_add.emplace_back(false);
				else
					pf.pvs_new_add.emplace_back(currPolyAddNew[jj]);
				jj = (jj + 1) % currPolys.size();
			}
			pf.pvs.emplace_back(currPolys[jj]);
			pf.pvs_new_add.emplace_back(false);
			Fs_.emplace_back(pf);
			break;
		}

		if (isFindGoodSplit)
		{
			queueFaceId.emplace(currFaceId);
			queueFaceId.emplace(Fs_.size() - 1);
		}
		else
		{
			currFaceId_.emplace_back(currFaceId);
		}
	}

	//第二步
	std::function<bool(double)> isAngleAfter = [&](double angle)->bool
	{
		if (angle > AFTER_ANGLE_MIN && angle < AFTER_ANGLE_MAX)
			return true;
		else
			return false;
	};

	typedef std::tuple<double, uint32_t, uint32_t, uint32_t, double, double> SplitTuple5;	//五个元素分别为分割后两个面积较小的那个；起始v, 分割点的左右两个v；分割点的pos
	std::vector<SplitTuple5> splitVec;
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();
	while (!queueFaceId.empty())
	{
		if (Fs_.size() >= polyNumThre)
		{
			BuildTopologyInfoFromVsFs();
			return;
		}
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35465);

		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		splitVec.reserve(ccVsIds.size());
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
#if SHORT_LENGTH_NOT_DO
			int iLargeId = FindRealLargerOldVertexId(ccVsIds[i], currPolyAddNew), iSmallId = FindRealSmallerOldVertexId(ccVsIds[i], currPolyAddNew);
			double length0 = (Vs_[currPolys[iLargeId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm(), length1 = (Vs_[currPolys[iSmallId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm();
			if (length0 < diagLength_ * SHORT_NOT_DO_LENGTH && length1 < diagLength_ * SHORT_NOT_DO_LENGTH)
				continue;
#endif
			for (int j = 0; j < 2; ++j)
			{
				int vPreId = mmod((int)ccVsIds[i] - 1, (int)currPolys.size()), vCurrId = ccVsIds[i];
				if (j == 1)
					vPreId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());

				double sourceAngle = ComputePolyAngle(currPolys, ccVsIds[i]);
				if ((sourceAngle<BEFORE_ANGLE_MIN_0 || sourceAngle>BEFORE_ANGLE_MAX_0) && sourceAngle < BEFORE_ANGLE_MIN_1)
					continue;

				/*if (sourceAngle > BEFORE_ANGLE_MIN_1)
				{
					AFTER_ANGLE_MIN = 60.0 * PI / 180.0;
					AFTER_ANGLE_MAX = 120.0 * PI / 180.0;
				}
				else
				{
					AFTER_ANGLE_MIN = 70.0 * PI / 180.0;
					AFTER_ANGLE_MAX = 110.0 * PI / 180.0;
				}*/
				AFTER_ANGLE_MIN = 70.0 * PI / 180.0;
				AFTER_ANGLE_MAX = 110.0 * PI / 180.0;

				if (!GetPolyRayIntersection(currPolys, vCurrId, vPreId, outV2Id_, outV3Id_, interPos_, outPoly0_, outPoly1_))
					continue;

				Eigen::Vector2d midPos = (Vs_[currPolys[vCurrId]].pos + interPos_) * 0.5;
				if (!CheckPointInsidePoly(currPolys, midPos))
					continue;
				//是否不允许出现三角形？
				if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
					continue;

				//直的点不算
				int realOutV2Id = FindSmallerOldVertexId(outV2Id_, currPolyAddNew), realOutV3Id = FindLargerOldVertexId(outV3Id_, currPolyAddNew);
				double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[realOutV2Id]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[realOutV3Id]].pos - interPos_).norm();
				std::vector<bool> poly0AddNew, poly1AddNew;
				GetSubAddNew(vCurrId, outV2Id_, outV3Id_, currPolyAddNew, poly0AddNew, poly1AddNew);
				if (GetFalseNum(poly0AddNew) <= 3 || GetFalseNum(poly1AddNew) <= 3)
					continue;

#if !FOR_QUADS_GEN
#if OLD_LENGTH_THRE
				if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
					continue;
#else
				double currLengthThre;
				//int tempResu = CheckPointPolyState(initialPoly_, interPos_);
				if (!CheckPointOnBoundary(interPos_))
					currLengthThre = LENGTH_THRESHOLD;
				else
					currLengthThre = LENGTH_THRESHOLD_0;
				if (newLength0 < LENGTH_THRESHOLD_0 * diagLength_ || newLength1 < currLengthThre * diagLength_ || newLength2 < currLengthThre * diagLength_)
					continue;
#endif
#endif

				double outArea0 = std::abs(ComputePolyArea(outPoly0_));
				double outArea1 = std::abs(ComputePolyArea(outPoly1_));
#if !FOR_QUADS_GEN
				if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
					continue;
#endif

				double angle0 = ComputePolyAngle(outPoly0_, outPoly0_.size() - 1), angle1 = ComputePolyAngle(outPoly1_, 0);
				if (!isAngleAfter(angle0) || !isAngleAfter(angle1))
					continue;

				double minEdge0 = 1.0, maxEdge0 = 1.0, minEdge1 = 1.0, maxEdge1 = 1.0;
				/*ComputeShortestAndLongestEdge(outPoly0_, minEdge0, maxEdge0);
				ComputeShortestAndLongestEdge(outPoly1_, minEdge1, maxEdge1);*/
				ComputeShortestAndLongestEdge(outPoly0_, poly0AddNew, minEdge0, maxEdge0);
				ComputeShortestAndLongestEdge(outPoly1_, poly1AddNew, minEdge1, maxEdge1);
				double minArea = std::min(outArea0, outArea1), maxArea = std::max(outArea0, outArea1);
				double ranking = minArea / maxArea + AREA_LENGTH_WEIGHT * 0.5 * (std::pow(minEdge0 / maxEdge0, 2) + std::pow(minEdge1 / maxEdge1, 2));

				splitVec.emplace_back(SplitTuple5(ranking, vCurrId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		int currSelect = 0;
		int randomTillNum = std::min(int(splitVec.size()*2.0 / 3.0 + 1), numInEveryStep*2);
		if (randomTillNum == 0)
			randomTillNum = 1;
		currSelect = GetRandomNum(0, randomTillNum);

		uint32_t startId = std::get<1>(splitVec[currSelect]);
		outV2Id_ = std::get<2>(splitVec[currSelect]);
		outV3Id_ = std::get<3>(splitVec[currSelect]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[currSelect]);
		pv.pos[1] = std::get<5>(splitVec[currSelect]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

#if OUTPUT_SPLIT_LINES
		splitEdges_.emplace_back(std::tuple<Eigen::Vector2d, Eigen::Vector2d>(pv.pos, Vs_[currPolys[startId]].pos));
#endif

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);


		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	//第三步
	std::function<bool(double)> isAngleAfter2 = [&](double angle)->bool
	{
		if (std::abs(angle) > 45.0*PI / 180.0)
			return true;
		else
			return false;
	};
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();

	while (!queueFaceId.empty())
	{
		if (Fs_.size() >= polyNumThre)
		{
			BuildTopologyInfoFromVsFs();
			return;
		}
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35464);

		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
#if SHORT_LENGTH_NOT_DO
			int iLargeId = FindRealLargerOldVertexId(ccVsIds[i], currPolyAddNew), iSmallId = FindRealSmallerOldVertexId(ccVsIds[i], currPolyAddNew);
			double length0 = (Vs_[currPolys[iLargeId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm(), length1 = (Vs_[currPolys[iSmallId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm();
			if (length0 < diagLength_ * SHORT_NOT_DO_LENGTH && length1 < diagLength_ * SHORT_NOT_DO_LENGTH)
				continue;
#endif
			uint32_t currId = ccVsIds[i];
			Eigen::Vector2d currPos = Vs_[currPolys[currId]].pos;

			for (int j = 0; j < currPolys.size(); ++j)
			{
				if (j == currId || ((j + 1) % currPolys.size()) == currId)
					continue;

				Eigen::Vector2d vPos0 = Vs_[currPolys[j]].pos, vPos1 = Vs_[currPolys[(j + 1) % currPolys.size()]].pos;
				double alpha = 0;
				Eigen::Vector2d interPos;
				if (!ComputeVerticalIntersection(currPos, vPos0, vPos1, alpha, interPos))
					continue;

				if (JudgeIfVerticalIntersection(currPolys, currId, interPos, j, (j + 1) % currPolys.size()))
					continue;

				Eigen::Vector2d midPos = (Vs_[currPolys[currId]].pos + interPos) * 0.5;
				if (!CheckPointInsidePoly(currPolys, midPos))
					continue;

				if (JudgeIfOnSameLine(currPolys, currId, interPos))
					continue;

				SplitPolyIntoTwo(currPolys, currId, interPos, j, (j + 1) % currPolys.size(), outPoly0_, outPoly1_);
				//是否不允许出现三角形？
				if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
					continue;
				outV2Id_ = j; outV3Id_ = (j + 1) % currPolys.size(); interPos_ = interPos;
				/*double newLength0 = (Vs_[currPolys[currId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[outV2Id_]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[outV3Id_]].pos - interPos_).norm();*/
				int realOutV2Id = FindSmallerOldVertexId(outV2Id_, currPolyAddNew), realOutV3Id = FindLargerOldVertexId(outV3Id_, currPolyAddNew);
				double newLength0 = (Vs_[currPolys[currId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[realOutV2Id]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[realOutV3Id]].pos - interPos_).norm();
				std::vector<bool> poly0AddNew, poly1AddNew;
				GetSubAddNew(currId, outV2Id_, outV3Id_, currPolyAddNew, poly0AddNew, poly1AddNew);
				if (GetFalseNum(poly0AddNew) <= 3 || GetFalseNum(poly1AddNew) <= 3)
					continue;

				//LENGTH_THRESHOLD = 0.0001;
#if !FOR_QUADS_GEN
#if OLD_LENGTH_THRE
				if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
					continue;
#else
				double currLengthThre;
				//int tempResu = CheckPointPolyState(initialPoly_, interPos_);
				if (!CheckPointOnBoundary(interPos_))
					currLengthThre = LENGTH_THRESHOLD;
				else
					currLengthThre = LENGTH_THRESHOLD_0;
				if (newLength0 < LENGTH_THRESHOLD_0 * diagLength_ || newLength1 < currLengthThre * diagLength_ || newLength2 < currLengthThre * diagLength_)
					continue;
#endif
#endif

#if !FOR_QUADS_GEN
				double outArea0 = std::abs(ComputePolyArea(outPoly0_));
				double outArea1 = std::abs(ComputePolyArea(outPoly1_));
				//AREA_THRESHOLD_2 = 0.0001;
				if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
					continue;
#endif

				double angle0 = ComputePolyAngle(outPoly0_, 0), angle1 = ComputePolyAngle(outPoly1_, outPoly1_.size() - 1);
				if (!isAngleAfter2(angle0) || !isAngleAfter2(angle1))
					continue;

				splitVec.emplace_back(SplitTuple5(std::min(angle0, angle1), currId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		int currSelect = 0;
		int randomTillNum = std::min(int(splitVec.size()*2.0 / 3.0 + 1), numInEveryStep*2);
		if (randomTillNum == 0)
			randomTillNum = 1;
		currSelect = GetRandomNum(0, randomTillNum);

		uint32_t startId = std::get<1>(splitVec[currSelect]);
		outV2Id_ = std::get<2>(splitVec[currSelect]);
		outV3Id_ = std::get<3>(splitVec[currSelect]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[currSelect]);
		pv.pos[1] = std::get<5>(splitVec[currSelect]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

#if OUTPUT_SPLIT_LINES
		splitEdges_.emplace_back(std::tuple<Eigen::Vector2d, Eigen::Vector2d>(pv.pos, Vs_[currPolys[startId]].pos));
#endif

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);

		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	//第四步
	std::function<bool(double)> isAngleAfter3 = [&](double angle)->bool
	{
		if (angle > AFTER_ANGLE_MIN_THREE && angle < AFTER_ANGLE_MAX_THREE)
			return true;
		else
			return false;
	};
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();

	while (!queueFaceId.empty())
	{
		if (Fs_.size() >= polyNumThre)
		{
			BuildTopologyInfoFromVsFs();
			return;
		}
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35464);

		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		splitVec.reserve(ccVsIds.size());
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
#if SHORT_LENGTH_NOT_DO
			int iLargeId = FindRealLargerOldVertexId(ccVsIds[i], currPolyAddNew), iSmallId = FindRealSmallerOldVertexId(ccVsIds[i], currPolyAddNew);
			double length0 = (Vs_[currPolys[iLargeId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm(), length1 = (Vs_[currPolys[iSmallId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm();
			if (length0 < diagLength_ * SHORT_NOT_DO_LENGTH && length1 < diagLength_ * SHORT_NOT_DO_LENGTH)
				continue;
#endif
			/*for (int j = 0; j < 2; ++j)
			{*/
			int vPreId = mmod((int)ccVsIds[i] - 1, (int)currPolys.size()), vCurrId = ccVsIds[i], vNextId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());
			/*if (j == 1)
				vPreId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());*/

				/*if (ccVsIds.size() == 1)
					DO_NOTHING_ANGLE_THRESHOLD = 20.0*PI / 180.0;
				else
					DO_NOTHING_ANGLE_THRESHOLD = 10.0*PI / 180.0;*/
					/*double concaveAngle0 = ComputePolyAngle(currPolys, ccVsIds[i]);
					if (std::abs(concaveAngle0 - PI) < DO_NOTHING_ANGLE_THRESHOLD)
						continue;*/

			Eigen::Vector2d aveAngleVec = (Vs_[currPolys[vPreId]].pos - Vs_[currPolys[vCurrId]].pos).normalized() + (Vs_[currPolys[vNextId]].pos - Vs_[currPolys[vCurrId]].pos).normalized();
			if (!GetPolyRayIntersection(currPolys, vCurrId, aveAngleVec, outV2Id_, outV3Id_, interPos_, outPoly0_, outPoly1_))
				continue;

			Eigen::Vector2d midPos = (Vs_[currPolys[vCurrId]].pos + interPos_) * 0.5;
			if (!CheckPointInsidePoly(currPolys, midPos))
				continue;
			//是否不允许出现三角形？
			if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
				continue;

			/*double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[outV2Id_]].pos - interPos_).norm(),
				newLength2 = (Vs_[currPolys[outV3Id_]].pos - interPos_).norm();*/
			int realOutV2Id = FindSmallerOldVertexId(outV2Id_, currPolyAddNew), realOutV3Id = FindLargerOldVertexId(outV3Id_, currPolyAddNew);
			double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[realOutV2Id]].pos - interPos_).norm(),
				newLength2 = (Vs_[currPolys[realOutV3Id]].pos - interPos_).norm();
			std::vector<bool> poly0AddNew, poly1AddNew;
			GetSubAddNew(vCurrId, outV2Id_, outV3Id_, currPolyAddNew, poly0AddNew, poly1AddNew);
			if (GetFalseNum(poly0AddNew) <= 3 || GetFalseNum(poly1AddNew) <= 3)
				continue;

#if !FOR_QUADS_GEN
#if OLD_LENGTH_THRE
			if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
				continue;
#else
			double currLengthThre;
			//int tempResu = CheckPointPolyState(initialPoly_, interPos_);
			if (!CheckPointOnBoundary(interPos_))
				currLengthThre = LENGTH_THRESHOLD;
			else
				currLengthThre = LENGTH_THRESHOLD_0;
			if (newLength0 < LENGTH_THRESHOLD_0 * diagLength_ || newLength1 < currLengthThre * diagLength_ || newLength2 < currLengthThre * diagLength_)
				continue;
#endif
#endif

			double outArea0 = std::abs(ComputePolyArea(outPoly0_));
			double outArea1 = std::abs(ComputePolyArea(outPoly1_));
#if !FOR_QUADS_GEN
			if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
				continue;
#endif

			double angle0 = ComputePolyAngle(outPoly0_, outPoly0_.size() - 1), angle1 = ComputePolyAngle(outPoly1_, 0);
			if (!isAngleAfter3(angle0) || !isAngleAfter3(angle1))
				continue;

			double minEdge0 = 1.0, maxEdge0 = 1.0, minEdge1 = 1.0, maxEdge1 = 1.0;
			/*ComputeShortestAndLongestEdge(outPoly0_, minEdge0, maxEdge0);
			ComputeShortestAndLongestEdge(outPoly1_, minEdge1, maxEdge1);*/
			ComputeShortestAndLongestEdge(outPoly0_, poly0AddNew, minEdge0, maxEdge0);
			ComputeShortestAndLongestEdge(outPoly1_, poly1AddNew, minEdge1, maxEdge1);
			double minArea = std::min(outArea0, outArea1), maxArea = std::max(outArea0, outArea1);
			double ranking = minArea / maxArea + AREA_LENGTH_WEIGHT * 0.5 * (std::pow(minEdge0 / maxEdge0, 2) + std::pow(minEdge1 / maxEdge1, 2));

			splitVec.emplace_back(SplitTuple5(ranking, vCurrId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			//}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		int currSelect = 0;
		int randomTillNum = std::min(int(splitVec.size()*2.0 / 3.0 + 1), numInEveryStep*2);
		if (randomTillNum == 0)
			randomTillNum = 1;
		currSelect = GetRandomNum(0, randomTillNum);

		uint32_t startId = std::get<1>(splitVec[currSelect]);
		outV2Id_ = std::get<2>(splitVec[currSelect]);
		outV3Id_ = std::get<3>(splitVec[currSelect]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[currSelect]);
		pv.pos[1] = std::get<5>(splitVec[currSelect]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

#if OUTPUT_SPLIT_LINES
		splitEdges_.emplace_back(std::tuple<Eigen::Vector2d, Eigen::Vector2d>(pv.pos, Vs_[currPolys[startId]].pos));
#endif

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);

		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	BuildTopologyInfoFromVsFs();
}

void PolySplit::SplitPoly(std::vector<Eigen::Vector2d> &inPoly, std::vector<bool> &inPolyAddNew, int polyNumThre)
{
	if (!JudgePolyCCW(inPoly))
	{
		std::reverse(inPoly.begin(), inPoly.end());
		std::reverse(inPolyAddNew.begin(), inPolyAddNew.end());
	}
	InitializePoly(inPoly, inPolyAddNew);

	currFaceId_.clear();

	std::queue<uint32_t> queueFaceId;
	queueFaceId.emplace(0);

	//第一步
	std::vector<uint32_t> ccVsIds;
	while (!queueFaceId.empty())
	{
		if (Fs_.size() >= polyNumThre)
		{
			BuildTopologyInfoFromVsFs();
			return;
		}
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35469);

		std::vector<TupleRankingLine> fieldRanking;
		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);

		std::vector<bool> concaveVFlag(currPolys.size(), false);
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			concaveVFlag[ccVsIds[i]] = true;
		}

		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			uint32_t v0 = currPolys[ccVsIds[i]];
#if SHORT_LENGTH_NOT_DO
			int iLargeId = FindRealLargerOldVertexId(ccVsIds[i], currPolyAddNew), iSmallId = FindRealSmallerOldVertexId(ccVsIds[i], currPolyAddNew);
			double length0 = (Vs_[currPolys[iLargeId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm(), length1 = (Vs_[currPolys[iSmallId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm();
			if (length0 < diagLength_ * SHORT_NOT_DO_LENGTH && length1 < diagLength_ * SHORT_NOT_DO_LENGTH)
				continue;
#endif
			for (int j = 0; j < currPolys.size(); ++j)
			{
				if (j == ccVsIds[i] || std::abs((int)(ccVsIds[i]) - (int)(j)) == 1 || currPolyAddNew[j])
					continue;

#if SHORT_LENGTH_NOT_DO
				int jLargeId = FindRealLargerOldVertexId(j, currPolyAddNew), jSmallId = FindRealSmallerOldVertexId(j, currPolyAddNew);
				double jlength0 = (Vs_[currPolys[jLargeId]].pos - Vs_[currPolys[j]].pos).norm(), jlength1 = (Vs_[currPolys[jSmallId]].pos - Vs_[currPolys[j]].pos).norm();
				if (jlength0 < diagLength_ * SHORT_NOT_DO_LENGTH && jlength1 < diagLength_ * SHORT_NOT_DO_LENGTH)
					continue;
#endif
				uint32_t v1 = currPolys[j];
				if (JudgeIfEdgeIntersectionOrOutside(currPolys, v0, v1))
					continue;

				if (JudgeIfOnSameLine(currPolyAddNew, ccVsIds[i], j))
					continue;

				double currFieldMatchAngle = LineFieldMatchAngle(v0, v1);
				bool isConcave = false;
				if (concaveVFlag[j])
				{
					currFieldMatchAngle -= 1E5;
					isConcave = true;
				}
				else
					isConcave = false;
				fieldRanking.emplace_back(TupleRankingLine(currFieldMatchAngle, v0, v1, ccVsIds[i], j, isConcave));
			}
		}

		if (fieldRanking.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		bool isFindGoodSplit = false;
		std::vector<uint32_t> newPoly0, newPoly1;
		std::sort(fieldRanking.begin(), fieldRanking.end());
		for (int i = 0; i < fieldRanking.size(); ++i)
		{
			double fieldValue = std::get<0>(fieldRanking[i]);
			uint32_t v0 = std::get<1>(fieldRanking[i]), v1 = std::get<2>(fieldRanking[i]), v0Id = std::get<3>(fieldRanking[i]), v1Id = std::get<4>(fieldRanking[i]);
			bool isConcave = std::get<5>(fieldRanking[i]);

			ANGLE_THRESHOLD = 75.0 * PI / 180.0;
			AREA_THRESHOLD = 0.003;
			/*if (isConcave)
			{
				ANGLE_THRESHOLD = 45.0 * PI / 180.0;
				AREA_THRESHOLD = 0.005;
			}
			else
			{
				ANGLE_THRESHOLD = 75.0 * PI / 180.0;
				AREA_THRESHOLD = 0.005;
			}*/

#if FOR_QUADS_GEN
			double fieldMatchV = fieldValue;
			if (fieldMatchV < 0)
				fieldMatchV += 1E5;
			if (fieldMatchV > FIELD_MATCHING_THRESHOLD)
				continue;
#else
			double fieldMatchV = fieldValue;
			if (fieldMatchV < 0)
				fieldMatchV += 1E5;
			if (fieldMatchV > FIELD_MATCHING_THRESHOLD_2)
				continue;
#endif
			double newLength = (Vs_[v0].pos - Vs_[v1].pos).norm();
#if OLD_LENGTH_THRE
			if (newLength < LENGTH_THRESHOLD * diagLength_)
				continue;
#else
			if (newLength < LENGTH_THRESHOLD_0 * diagLength_)
				continue;
#endif

			//SplitPolyIntoTwo(currPolys, v0Id, v1Id, newPoly0, newPoly1);
			std::vector<bool> newPAN0, newPAN1;
			SplitPolyIntoTwo(currPolys, currPolyAddNew, v0Id, v1Id, newPoly0, newPoly1, newPAN0, newPAN1);

			if (GetFalseNum(newPAN0) <= 3 || GetFalseNum(newPAN1) <= 3)
				continue;

			//是否不允许出现三角形？
			if (newPoly0.size() <= 3 || newPoly1.size() <= 3)
				continue;

#if !FOR_QUADS_GEN
			double poly0Area = std::abs(ComputePolyArea(newPoly0)), poly1Area = std::abs(ComputePolyArea(newPoly1));
			if (poly0Area < AREA_THRESHOLD * sourcePolyArea_ || poly1Area < AREA_THRESHOLD * sourcePolyArea_)
				continue;
#endif

			double angle0 = ComputePolyAngle(newPoly0, 0), angle1 = ComputePolyAngle(newPoly0, newPoly0.size() - 1), angle2 = ComputePolyAngle(newPoly1, 0), angle3 = ComputePolyAngle(newPoly1, newPoly1.size() - 1);
			if (angle0 < ANGLE_THRESHOLD || angle1 < ANGLE_THRESHOLD || angle2 < ANGLE_THRESHOLD || angle3 < ANGLE_THRESHOLD)
				continue;

			isFindGoodSplit = true;

#if OUTPUT_SPLIT_LINES
			splitEdges_.emplace_back(std::tuple<Eigen::Vector2d, Eigen::Vector2d>(Vs_[v0].pos, Vs_[v1].pos));
#endif

			//更新Vs_, Fs_拓扑信息
			Vs_[currPolys[v0Id]].neighbor_pfs.emplace_back(Fs_.size());
			Vs_[currPolys[v1Id]].neighbor_pfs.emplace_back(Fs_.size());
			Fs_[currFaceId].pvs.clear();
			Fs_[currFaceId].pvs_new_add.clear();

			int jj = v0Id;
			while (jj != v1Id)
			{
				Fs_[currFaceId].pvs.emplace_back(currPolys[jj]);
				if (jj == v0Id)
					Fs_[currFaceId].pvs_new_add.emplace_back(false);
				else
					Fs_[currFaceId].pvs_new_add.emplace_back(currPolyAddNew[jj]);
				jj = (jj + 1) % currPolys.size();
			}
			Fs_[currFaceId].pvs.emplace_back(currPolys[jj]);
			Fs_[currFaceId].pvs_new_add.emplace_back(false);

			for (int j = 0; j < newPoly1.size(); ++j)
			{
				if (newPoly1[j] != currPolys[v0Id] && newPoly1[j] != currPolys[v1Id])
				{
					std::vector<uint32_t> &vfs = Vs_[newPoly1[j]].neighbor_pfs;
					for (int k = 0; k < vfs.size(); ++k)
					{
						if (vfs[k] == currFaceId)
							vfs[k] = Fs_.size();
					}
				}
			}
			PolyFace pf;
			pf.id = Fs_.size();
			jj = v1Id;
			while (jj != v0Id)
			{
				pf.pvs.emplace_back(currPolys[jj]);
				if (jj == v1Id)
					pf.pvs_new_add.emplace_back(false);
				else
					pf.pvs_new_add.emplace_back(currPolyAddNew[jj]);
				jj = (jj + 1) % currPolys.size();
			}
			pf.pvs.emplace_back(currPolys[jj]);
			pf.pvs_new_add.emplace_back(false);
			Fs_.emplace_back(pf);
			break;
		}

		if (isFindGoodSplit)
		{
			queueFaceId.emplace(currFaceId);
			queueFaceId.emplace(Fs_.size() - 1);
		}
		else
		{
			currFaceId_.emplace_back(currFaceId);
		}
	}

	//第二步
	std::function<bool(double)> isAngleAfter = [&](double angle)->bool
	{
		if (angle > AFTER_ANGLE_MIN && angle < AFTER_ANGLE_MAX)
			return true;
		else
			return false;
	};

	typedef std::tuple<double, uint32_t, uint32_t, uint32_t, double, double> SplitTuple5;	//五个元素分别为分割后两个面积较小的那个；起始v, 分割点的左右两个v；分割点的pos
	std::vector<SplitTuple5> splitVec;
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();
	while (!queueFaceId.empty())
	{
		if (Fs_.size() >= polyNumThre)
		{
			BuildTopologyInfoFromVsFs();
			return;
		}
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35465);

		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		splitVec.reserve(ccVsIds.size());
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
#if SHORT_LENGTH_NOT_DO
			int iLargeId = FindRealLargerOldVertexId(ccVsIds[i], currPolyAddNew), iSmallId = FindRealSmallerOldVertexId(ccVsIds[i], currPolyAddNew);
			double length0 = (Vs_[currPolys[iLargeId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm(), length1 = (Vs_[currPolys[iSmallId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm();
			if (length0 < diagLength_ * SHORT_NOT_DO_LENGTH && length1 < diagLength_ * SHORT_NOT_DO_LENGTH)
				continue;
#endif
			for (int j = 0; j < 2; ++j)
			{
				int vPreId = mmod((int)ccVsIds[i] - 1, (int)currPolys.size()), vCurrId = ccVsIds[i];
				if (j == 1)
					vPreId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());

				double sourceAngle = ComputePolyAngle(currPolys, ccVsIds[i]);
				if ((sourceAngle<BEFORE_ANGLE_MIN_0 || sourceAngle>BEFORE_ANGLE_MAX_0) && sourceAngle < BEFORE_ANGLE_MIN_1)
					continue;

				/*if (sourceAngle > BEFORE_ANGLE_MIN_1)
				{
					AFTER_ANGLE_MIN = 60.0 * PI / 180.0;
					AFTER_ANGLE_MAX = 120.0 * PI / 180.0;
				}
				else
				{
					AFTER_ANGLE_MIN = 70.0 * PI / 180.0;
					AFTER_ANGLE_MAX = 110.0 * PI / 180.0;
				}*/
				AFTER_ANGLE_MIN = 70.0 * PI / 180.0;
				AFTER_ANGLE_MAX = 110.0 * PI / 180.0;

				if (!GetPolyRayIntersection(currPolys, vCurrId, vPreId, outV2Id_, outV3Id_, interPos_, outPoly0_, outPoly1_))
					continue;

				Eigen::Vector2d midPos = (Vs_[currPolys[vCurrId]].pos + interPos_) * 0.5;
				if (!CheckPointInsidePoly(currPolys, midPos))
					continue;
				//是否不允许出现三角形？
				if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
					continue;

				//直的点不算
				int realOutV2Id = FindSmallerOldVertexId(outV2Id_, currPolyAddNew), realOutV3Id = FindLargerOldVertexId(outV3Id_, currPolyAddNew);
				double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[realOutV2Id]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[realOutV3Id]].pos - interPos_).norm();
				std::vector<bool> poly0AddNew, poly1AddNew;
				GetSubAddNew(vCurrId, outV2Id_, outV3Id_, currPolyAddNew, poly0AddNew, poly1AddNew);
				if (GetFalseNum(poly0AddNew) <= 3 || GetFalseNum(poly1AddNew) <= 3)
					continue;

#if !FOR_QUADS_GEN
#if OLD_LENGTH_THRE
				if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
					continue;
#else
				double currLengthThre;
				//int tempResu = CheckPointPolyState(initialPoly_, interPos_);
				if (!CheckPointOnBoundary(interPos_))
					currLengthThre = LENGTH_THRESHOLD;
				else
					currLengthThre = LENGTH_THRESHOLD_0;
				if (newLength0 < LENGTH_THRESHOLD_0 * diagLength_ || newLength1 < currLengthThre * diagLength_ || newLength2 < currLengthThre * diagLength_)
					continue;
#endif
#endif

				double outArea0 = std::abs(ComputePolyArea(outPoly0_));
				double outArea1 = std::abs(ComputePolyArea(outPoly1_));
#if !FOR_QUADS_GEN
				if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
					continue;
#endif

				double angle0 = ComputePolyAngle(outPoly0_, outPoly0_.size() - 1), angle1 = ComputePolyAngle(outPoly1_, 0);
				if (!isAngleAfter(angle0) || !isAngleAfter(angle1))
					continue;

				double minEdge0 = 1.0, maxEdge0 = 1.0, minEdge1 = 1.0, maxEdge1 = 1.0;
				/*ComputeShortestAndLongestEdge(outPoly0_, minEdge0, maxEdge0);
				ComputeShortestAndLongestEdge(outPoly1_, minEdge1, maxEdge1);*/
				ComputeShortestAndLongestEdge(outPoly0_, poly0AddNew, minEdge0, maxEdge0);
				ComputeShortestAndLongestEdge(outPoly1_, poly1AddNew, minEdge1, maxEdge1);
				double minArea = std::min(outArea0, outArea1), maxArea = std::max(outArea0, outArea1);
				double ranking = minArea / maxArea + AREA_LENGTH_WEIGHT * 0.5 * (std::pow(minEdge0 / maxEdge0, 2) + std::pow(minEdge1 / maxEdge1, 2));

				splitVec.emplace_back(SplitTuple5(ranking, vCurrId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		uint32_t startId = std::get<1>(splitVec[0]);
		outV2Id_ = std::get<2>(splitVec[0]);
		outV3Id_ = std::get<3>(splitVec[0]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[0]);
		pv.pos[1] = std::get<5>(splitVec[0]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

#if OUTPUT_SPLIT_LINES
		splitEdges_.emplace_back(std::tuple<Eigen::Vector2d, Eigen::Vector2d>(pv.pos, Vs_[currPolys[startId]].pos));
#endif

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);


		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	//第三步
	std::function<bool(double)> isAngleAfter2 = [&](double angle)->bool
	{
		if (std::abs(angle) > 45.0*PI / 180.0)
			return true;
		else
			return false;
	};
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();

	while (!queueFaceId.empty())
	{
		if (Fs_.size() >= polyNumThre)
		{
			BuildTopologyInfoFromVsFs();
			return;
		}
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35464);

		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
#if SHORT_LENGTH_NOT_DO
			int iLargeId = FindRealLargerOldVertexId(ccVsIds[i], currPolyAddNew), iSmallId = FindRealSmallerOldVertexId(ccVsIds[i], currPolyAddNew);
			double length0 = (Vs_[currPolys[iLargeId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm(), length1 = (Vs_[currPolys[iSmallId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm();
			if (length0 < diagLength_ * SHORT_NOT_DO_LENGTH && length1 < diagLength_ * SHORT_NOT_DO_LENGTH)
				continue;
#endif
			uint32_t currId = ccVsIds[i];
			Eigen::Vector2d currPos = Vs_[currPolys[currId]].pos;

			for (int j = 0; j < currPolys.size(); ++j)
			{
				if (j == currId || ((j + 1) % currPolys.size()) == currId)
					continue;

				Eigen::Vector2d vPos0 = Vs_[currPolys[j]].pos, vPos1 = Vs_[currPolys[(j + 1) % currPolys.size()]].pos;
				double alpha = 0;
				Eigen::Vector2d interPos;
				if (!ComputeVerticalIntersection(currPos, vPos0, vPos1, alpha, interPos))
					continue;

				if (JudgeIfVerticalIntersection(currPolys, currId, interPos, j, (j + 1) % currPolys.size()))
					continue;

				Eigen::Vector2d midPos = (Vs_[currPolys[currId]].pos + interPos) * 0.5;
				if (!CheckPointInsidePoly(currPolys, midPos))
					continue;

				if (JudgeIfOnSameLine(currPolys, currId, interPos))
					continue;

				SplitPolyIntoTwo(currPolys, currId, interPos, j, (j + 1) % currPolys.size(), outPoly0_, outPoly1_);
				//是否不允许出现三角形？
				if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
					continue;
				outV2Id_ = j; outV3Id_ = (j + 1) % currPolys.size(); interPos_ = interPos;
				/*double newLength0 = (Vs_[currPolys[currId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[outV2Id_]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[outV3Id_]].pos - interPos_).norm();*/
				int realOutV2Id = FindSmallerOldVertexId(outV2Id_, currPolyAddNew), realOutV3Id = FindLargerOldVertexId(outV3Id_, currPolyAddNew);
				double newLength0 = (Vs_[currPolys[currId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[realOutV2Id]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[realOutV3Id]].pos - interPos_).norm();
				std::vector<bool> poly0AddNew, poly1AddNew;
				GetSubAddNew(currId, outV2Id_, outV3Id_, currPolyAddNew, poly0AddNew, poly1AddNew);
				if (GetFalseNum(poly0AddNew) <= 3 || GetFalseNum(poly1AddNew) <= 3)
					continue;

				//LENGTH_THRESHOLD = 0.0001;
#if !FOR_QUADS_GEN
#if OLD_LENGTH_THRE
				if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
					continue;
#else
				double currLengthThre;
				//int tempResu = CheckPointPolyState(initialPoly_, interPos_);
				if (!CheckPointOnBoundary(interPos_))
					currLengthThre = LENGTH_THRESHOLD;
				else
					currLengthThre = LENGTH_THRESHOLD_0;
				if (newLength0 < LENGTH_THRESHOLD_0 * diagLength_ || newLength1 < currLengthThre * diagLength_ || newLength2 < currLengthThre * diagLength_)
					continue;
#endif
#endif

#if !FOR_QUADS_GEN
				double outArea0 = std::abs(ComputePolyArea(outPoly0_));
				double outArea1 = std::abs(ComputePolyArea(outPoly1_));
				//AREA_THRESHOLD_2 = 0.0001;
				if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
					continue;
#endif

				double angle0 = ComputePolyAngle(outPoly0_, 0), angle1 = ComputePolyAngle(outPoly1_, outPoly1_.size() - 1);
				if (!isAngleAfter2(angle0) || !isAngleAfter2(angle1))
					continue;

				splitVec.emplace_back(SplitTuple5(std::min(angle0, angle1), currId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		uint32_t startId = std::get<1>(splitVec[0]);
		outV2Id_ = std::get<2>(splitVec[0]);
		outV3Id_ = std::get<3>(splitVec[0]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[0]);
		pv.pos[1] = std::get<5>(splitVec[0]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

#if OUTPUT_SPLIT_LINES
		splitEdges_.emplace_back(std::tuple<Eigen::Vector2d, Eigen::Vector2d>(pv.pos, Vs_[currPolys[startId]].pos));
#endif

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);

		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	//第四步
	std::function<bool(double)> isAngleAfter3 = [&](double angle)->bool
	{
		if (angle > AFTER_ANGLE_MIN_THREE && angle < AFTER_ANGLE_MAX_THREE)
			return true;
		else
			return false;
	};
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();

	while (!queueFaceId.empty())
	{
		if (Fs_.size() >= polyNumThre)
		{
			BuildTopologyInfoFromVsFs();
			return;
		}
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35464);

		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		splitVec.reserve(ccVsIds.size());
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
#if SHORT_LENGTH_NOT_DO
			int iLargeId = FindRealLargerOldVertexId(ccVsIds[i], currPolyAddNew), iSmallId = FindRealSmallerOldVertexId(ccVsIds[i], currPolyAddNew);
			double length0 = (Vs_[currPolys[iLargeId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm(), length1 = (Vs_[currPolys[iSmallId]].pos - Vs_[currPolys[ccVsIds[i]]].pos).norm();
			if (length0 < diagLength_ * SHORT_NOT_DO_LENGTH && length1 < diagLength_ * SHORT_NOT_DO_LENGTH)
				continue;
#endif
			/*for (int j = 0; j < 2; ++j)
			{*/
			int vPreId = mmod((int)ccVsIds[i] - 1, (int)currPolys.size()), vCurrId = ccVsIds[i], vNextId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());
			/*if (j == 1)
				vPreId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());*/

				/*if (ccVsIds.size() == 1)
					DO_NOTHING_ANGLE_THRESHOLD = 20.0*PI / 180.0;
				else
					DO_NOTHING_ANGLE_THRESHOLD = 10.0*PI / 180.0;*/
					/*double concaveAngle0 = ComputePolyAngle(currPolys, ccVsIds[i]);
					if (std::abs(concaveAngle0 - PI) < DO_NOTHING_ANGLE_THRESHOLD)
						continue;*/

			Eigen::Vector2d aveAngleVec = (Vs_[currPolys[vPreId]].pos - Vs_[currPolys[vCurrId]].pos).normalized() + (Vs_[currPolys[vNextId]].pos - Vs_[currPolys[vCurrId]].pos).normalized();
			if (!GetPolyRayIntersection(currPolys, vCurrId, aveAngleVec, outV2Id_, outV3Id_, interPos_, outPoly0_, outPoly1_))
				continue;

			Eigen::Vector2d midPos = (Vs_[currPolys[vCurrId]].pos + interPos_) * 0.5;
			if (!CheckPointInsidePoly(currPolys, midPos))
				continue;
			//是否不允许出现三角形？
			if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
				continue;

			/*double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[outV2Id_]].pos - interPos_).norm(),
				newLength2 = (Vs_[currPolys[outV3Id_]].pos - interPos_).norm();*/
			int realOutV2Id = FindSmallerOldVertexId(outV2Id_, currPolyAddNew), realOutV3Id = FindLargerOldVertexId(outV3Id_, currPolyAddNew);
			double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[realOutV2Id]].pos - interPos_).norm(),
				newLength2 = (Vs_[currPolys[realOutV3Id]].pos - interPos_).norm();
			std::vector<bool> poly0AddNew, poly1AddNew;
			GetSubAddNew(vCurrId, outV2Id_, outV3Id_, currPolyAddNew, poly0AddNew, poly1AddNew);
			if (GetFalseNum(poly0AddNew) <= 3 || GetFalseNum(poly1AddNew) <= 3)
				continue;

#if !FOR_QUADS_GEN
#if OLD_LENGTH_THRE
			if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
				continue;
#else
			double currLengthThre;
			//int tempResu = CheckPointPolyState(initialPoly_, interPos_);
			if (!CheckPointOnBoundary(interPos_))
				currLengthThre = LENGTH_THRESHOLD;
			else
				currLengthThre = LENGTH_THRESHOLD_0;
			if (newLength0 < LENGTH_THRESHOLD_0 * diagLength_ || newLength1 < currLengthThre * diagLength_ || newLength2 < currLengthThre * diagLength_)
				continue;
#endif
#endif

			double outArea0 = std::abs(ComputePolyArea(outPoly0_));
			double outArea1 = std::abs(ComputePolyArea(outPoly1_));
#if !FOR_QUADS_GEN
			if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
				continue;
#endif

			double angle0 = ComputePolyAngle(outPoly0_, outPoly0_.size() - 1), angle1 = ComputePolyAngle(outPoly1_, 0);
			if (!isAngleAfter3(angle0) || !isAngleAfter3(angle1))
				continue;

			double minEdge0 = 1.0, maxEdge0 = 1.0, minEdge1 = 1.0, maxEdge1 = 1.0;
			/*ComputeShortestAndLongestEdge(outPoly0_, minEdge0, maxEdge0);
			ComputeShortestAndLongestEdge(outPoly1_, minEdge1, maxEdge1);*/
			ComputeShortestAndLongestEdge(outPoly0_, poly0AddNew, minEdge0, maxEdge0);
			ComputeShortestAndLongestEdge(outPoly1_, poly1AddNew, minEdge1, maxEdge1);
			double minArea = std::min(outArea0, outArea1), maxArea = std::max(outArea0, outArea1);
			double ranking = minArea / maxArea + AREA_LENGTH_WEIGHT * 0.5 * (std::pow(minEdge0 / maxEdge0, 2) + std::pow(minEdge1 / maxEdge1, 2));

			splitVec.emplace_back(SplitTuple5(ranking, vCurrId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			//}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		uint32_t startId = std::get<1>(splitVec[0]);
		outV2Id_ = std::get<2>(splitVec[0]);
		outV3Id_ = std::get<3>(splitVec[0]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[0]);
		pv.pos[1] = std::get<5>(splitVec[0]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

#if OUTPUT_SPLIT_LINES
		splitEdges_.emplace_back(std::tuple<Eigen::Vector2d, Eigen::Vector2d>(pv.pos, Vs_[currPolys[startId]].pos));
#endif

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);

		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	BuildTopologyInfoFromVsFs();
}

void PolySplit::SplitPoly(std::vector<Eigen::Vector2d> &inPoly, std::vector<bool> &inPolyAddNew)
{
	if (!JudgePolyCCW(inPoly))
	{
		std::reverse(inPoly.begin(), inPoly.end());
		std::reverse(inPolyAddNew.begin(), inPolyAddNew.end());
	}
	InitializePoly(inPoly, inPolyAddNew);

	currFaceId_.clear();

	std::queue<uint32_t> queueFaceId;
	queueFaceId.emplace(0);

	//第一步
	std::vector<uint32_t> ccVsIds;
	while (!queueFaceId.empty())
	{
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35469);

		std::vector<TupleRankingLine> fieldRanking;
		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);

		std::vector<bool> concaveVFlag(currPolys.size(), false);
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			concaveVFlag[ccVsIds[i]] = true;
		}

		for (int i=0; i<ccVsIds.size(); ++i)
		{
			uint32_t v0 = currPolys[ccVsIds[i]];
			for (int j = 0; j < currPolys.size(); ++j)
			{
				if (j == ccVsIds[i] || std::abs((int)(ccVsIds[i]) - (int)(j)) == 1 || currPolyAddNew[j])
					continue;

				uint32_t v1 = currPolys[j];
				if (JudgeIfEdgeIntersectionOrOutside(currPolys, v0, v1))
					continue;

				if (JudgeIfOnSameLine(currPolyAddNew, ccVsIds[i], j))
					continue;

				double currFieldMatchAngle = LineFieldMatchAngle(v0, v1);
				bool isConcave = false;
				if (concaveVFlag[j])
				{
					currFieldMatchAngle -= 1E5;
					isConcave = true;
				}
				else
					isConcave = false;
				fieldRanking.emplace_back(TupleRankingLine(currFieldMatchAngle, v0, v1, ccVsIds[i], j, isConcave));
			}
		}

		if (fieldRanking.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		bool isFindGoodSplit = false;
		std::vector<uint32_t> newPoly0, newPoly1;
		std::sort(fieldRanking.begin(), fieldRanking.end());
		for (int i = 0; i < fieldRanking.size(); ++i)
		{
			double fieldValue = std::get<0>(fieldRanking[i]);
			uint32_t v0 = std::get<1>(fieldRanking[i]), v1 = std::get<2>(fieldRanking[i]), v0Id = std::get<3>(fieldRanking[i]), v1Id = std::get<4>(fieldRanking[i]);
			bool isConcave = std::get<5>(fieldRanking[i]);

			ANGLE_THRESHOLD = 75.0 * PI / 180.0;
			AREA_THRESHOLD = 0.003;
			/*if (isConcave)
			{
				ANGLE_THRESHOLD = 45.0 * PI / 180.0;
				AREA_THRESHOLD = 0.005;
			}
			else
			{
				ANGLE_THRESHOLD = 75.0 * PI / 180.0;
				AREA_THRESHOLD = 0.005;
			}*/

			/*if (fieldValue > FIELD_MATCHING_THRESHOLD)
				continue;*/
			double newLength = (Vs_[v0].pos - Vs_[v1].pos).norm();
			if (newLength < LENGTH_THRESHOLD * diagLength_)
				continue;

			//SplitPolyIntoTwo(currPolys, v0Id, v1Id, newPoly0, newPoly1);
			std::vector<bool> newPAN0, newPAN1;
			SplitPolyIntoTwo(currPolys, currPolyAddNew, v0Id, v1Id, newPoly0, newPoly1, newPAN0, newPAN1);

			if (GetFalseNum(newPAN0) <= 3 || GetFalseNum(newPAN1) <= 3)
				continue;

			//是否不允许出现三角形？
			if (newPoly0.size() <= 3 || newPoly1.size() <= 3)
				continue;

			double poly0Area = std::abs(ComputePolyArea(newPoly0)), poly1Area = std::abs(ComputePolyArea(newPoly1));
			if (poly0Area < AREA_THRESHOLD * sourcePolyArea_ || poly1Area < AREA_THRESHOLD * sourcePolyArea_)
				continue;

			double angle0 = ComputePolyAngle(newPoly0, 0), angle1 = ComputePolyAngle(newPoly0, newPoly0.size() - 1), angle2 = ComputePolyAngle(newPoly1, 0), angle3 = ComputePolyAngle(newPoly1, newPoly1.size() - 1);
			if (angle0 < ANGLE_THRESHOLD || angle1 < ANGLE_THRESHOLD || angle2 < ANGLE_THRESHOLD || angle3 < ANGLE_THRESHOLD)
				continue;

			isFindGoodSplit = true;

			//更新Vs_, Fs_拓扑信息
			Vs_[currPolys[v0Id]].neighbor_pfs.emplace_back(Fs_.size());
			Vs_[currPolys[v1Id]].neighbor_pfs.emplace_back(Fs_.size());
			Fs_[currFaceId].pvs.clear();
			Fs_[currFaceId].pvs_new_add.clear();

			int jj = v0Id;
			while (jj != v1Id)
			{
				Fs_[currFaceId].pvs.emplace_back(currPolys[jj]);
				if (jj == v0Id)
					Fs_[currFaceId].pvs_new_add.emplace_back(false);
				else
					Fs_[currFaceId].pvs_new_add.emplace_back(currPolyAddNew[jj]);
				jj = (jj + 1) % currPolys.size();
			}
			Fs_[currFaceId].pvs.emplace_back(currPolys[jj]);
			Fs_[currFaceId].pvs_new_add.emplace_back(false);

			for (int j = 0; j < newPoly1.size(); ++j)
			{
				if (newPoly1[j] != currPolys[v0Id] && newPoly1[j] != currPolys[v1Id])
				{
					std::vector<uint32_t> &vfs = Vs_[newPoly1[j]].neighbor_pfs;
					for (int k = 0; k < vfs.size(); ++k)
					{
						if (vfs[k] == currFaceId)
							vfs[k] = Fs_.size();
					}
				}
			}
			PolyFace pf;
			pf.id = Fs_.size();
			jj = v1Id;
			while (jj != v0Id)
			{
				pf.pvs.emplace_back(currPolys[jj]);
				if (jj == v1Id)
					pf.pvs_new_add.emplace_back(false);
				else
					pf.pvs_new_add.emplace_back(currPolyAddNew[jj]);
				jj = (jj + 1) % currPolys.size();
			}
			pf.pvs.emplace_back(currPolys[jj]);
			pf.pvs_new_add.emplace_back(false);
			Fs_.emplace_back(pf);
			break;
		}

		if (isFindGoodSplit)
		{
			queueFaceId.emplace(currFaceId);
			queueFaceId.emplace(Fs_.size() - 1);
		}
		else
		{
			currFaceId_.emplace_back(currFaceId);
		}
	}

	//第二步
	std::function<bool(double)> isAngleAfter = [&](double angle)->bool
	{
		if (angle > AFTER_ANGLE_MIN && angle < AFTER_ANGLE_MAX)
			return true;
		else
			return false;
	};

	typedef std::tuple<double, uint32_t, uint32_t, uint32_t, double, double> SplitTuple5;	//五个元素分别为分割后两个面积较小的那个；起始v, 分割点的左右两个v；分割点的pos
	std::vector<SplitTuple5> splitVec;
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();
	while (!queueFaceId.empty())
	{
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35465);

		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		splitVec.reserve(ccVsIds.size());
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				int vPreId = mmod((int)ccVsIds[i] - 1, (int)currPolys.size()), vCurrId = ccVsIds[i];
				if (j == 1)
					vPreId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());

				double sourceAngle = ComputePolyAngle(currPolys, ccVsIds[i]);
				if ((sourceAngle<BEFORE_ANGLE_MIN_0 || sourceAngle>BEFORE_ANGLE_MAX_0) && sourceAngle < BEFORE_ANGLE_MIN_1)
					continue;

				/*if (sourceAngle > BEFORE_ANGLE_MIN_1)
				{
					AFTER_ANGLE_MIN = 60.0 * PI / 180.0;
					AFTER_ANGLE_MAX = 120.0 * PI / 180.0;
				}
				else
				{
					AFTER_ANGLE_MIN = 70.0 * PI / 180.0;
					AFTER_ANGLE_MAX = 110.0 * PI / 180.0;
				}*/
				AFTER_ANGLE_MIN = 70.0 * PI / 180.0;
				AFTER_ANGLE_MAX = 110.0 * PI / 180.0;

				if (!GetPolyRayIntersection(currPolys, vCurrId, vPreId, outV2Id_, outV3Id_, interPos_, outPoly0_, outPoly1_))
					continue;

				Eigen::Vector2d midPos = (Vs_[currPolys[vCurrId]].pos + interPos_) * 0.5;
				if (!CheckPointInsidePoly(currPolys, midPos))
					continue;
				//是否不允许出现三角形？
				if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
					continue;

				//直的点不算
				int realOutV2Id = FindSmallerOldVertexId(outV2Id_, currPolyAddNew), realOutV3Id = FindLargerOldVertexId(outV3Id_, currPolyAddNew);
				double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[realOutV2Id]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[realOutV3Id]].pos - interPos_).norm();
				std::vector<bool> poly0AddNew, poly1AddNew;
				GetSubAddNew(vCurrId, outV2Id_, outV3Id_, currPolyAddNew, poly0AddNew, poly1AddNew);
				if (GetFalseNum(poly0AddNew) <= 3 || GetFalseNum(poly1AddNew) <= 3)
					continue;

				if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
					continue;

				double outArea0 = std::abs(ComputePolyArea(outPoly0_));
				double outArea1 = std::abs(ComputePolyArea(outPoly1_));
				if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
					continue;

				double angle0 = ComputePolyAngle(outPoly0_, outPoly0_.size() - 1), angle1 = ComputePolyAngle(outPoly1_, 0);
				if (!isAngleAfter(angle0) || !isAngleAfter(angle1))
					continue;

				double minEdge0 = 1.0, maxEdge0 = 1.0, minEdge1 = 1.0, maxEdge1 = 1.0;
				/*ComputeShortestAndLongestEdge(outPoly0_, minEdge0, maxEdge0);
				ComputeShortestAndLongestEdge(outPoly1_, minEdge1, maxEdge1);*/
				ComputeShortestAndLongestEdge(outPoly0_, poly0AddNew, minEdge0, maxEdge0);
				ComputeShortestAndLongestEdge(outPoly1_, poly1AddNew, minEdge1, maxEdge1);
				double minArea = std::min(outArea0, outArea1), maxArea = std::max(outArea0, outArea1);
				double ranking = minArea / maxArea + AREA_LENGTH_WEIGHT * 0.5 * (std::pow(minEdge0 / maxEdge0, 2) + std::pow(minEdge1 / maxEdge1, 2));

				splitVec.emplace_back(SplitTuple5(ranking, vCurrId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		uint32_t startId = std::get<1>(splitVec[0]);
		outV2Id_ = std::get<2>(splitVec[0]);
		outV3Id_ = std::get<3>(splitVec[0]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[0]);
		pv.pos[1] = std::get<5>(splitVec[0]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);


		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	//第三步
	std::function<bool(double)> isAngleAfter2 = [&](double angle)->bool
	{
		if (std::abs(angle) > 45.0*PI / 180.0)
			return true;
		else
			return false;
	};
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();

	while (!queueFaceId.empty())
	{
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35464);

		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			uint32_t currId = ccVsIds[i];
			Eigen::Vector2d currPos = Vs_[currPolys[currId]].pos;

			for (int j = 0; j < currPolys.size(); ++j)
			{
				if (j == currId || ((j + 1) % currPolys.size()) == currId)
					continue;

				Eigen::Vector2d vPos0 = Vs_[currPolys[j]].pos, vPos1 = Vs_[currPolys[(j + 1) % currPolys.size()]].pos;
				double alpha = 0;
				Eigen::Vector2d interPos;
				if (!ComputeVerticalIntersection(currPos, vPos0, vPos1, alpha, interPos))
					continue;

				if (JudgeIfVerticalIntersection(currPolys, currId, interPos, j, (j + 1) % currPolys.size()))
					continue;

				Eigen::Vector2d midPos = (Vs_[currPolys[currId]].pos + interPos) * 0.5;
				if (!CheckPointInsidePoly(currPolys, midPos))
					continue;

				if (JudgeIfOnSameLine(currPolys, currId, interPos))
					continue;

				SplitPolyIntoTwo(currPolys, currId, interPos, j, (j + 1) % currPolys.size(), outPoly0_, outPoly1_);
				//是否不允许出现三角形？
				if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
					continue;
				outV2Id_ = j; outV3Id_ = (j + 1) % currPolys.size(); interPos_ = interPos;
				/*double newLength0 = (Vs_[currPolys[currId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[outV2Id_]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[outV3Id_]].pos - interPos_).norm();*/
				int realOutV2Id = FindSmallerOldVertexId(outV2Id_, currPolyAddNew), realOutV3Id = FindLargerOldVertexId(outV3Id_, currPolyAddNew);
				double newLength0 = (Vs_[currPolys[currId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[realOutV2Id]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[realOutV3Id]].pos - interPos_).norm();
				std::vector<bool> poly0AddNew, poly1AddNew;
				GetSubAddNew(currId, outV2Id_, outV3Id_, currPolyAddNew, poly0AddNew, poly1AddNew);
				if (GetFalseNum(poly0AddNew) <= 3 || GetFalseNum(poly1AddNew) <= 3)
					continue;

				//LENGTH_THRESHOLD = 0.0001;
				if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
					continue;

				double outArea0 = std::abs(ComputePolyArea(outPoly0_));
				double outArea1 = std::abs(ComputePolyArea(outPoly1_));
				//AREA_THRESHOLD_2 = 0.0001;
				if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
					continue;

				double angle0 = ComputePolyAngle(outPoly0_, 0), angle1 = ComputePolyAngle(outPoly1_, outPoly1_.size() - 1);
				if (!isAngleAfter2(angle0) || !isAngleAfter2(angle1))
					continue;

				splitVec.emplace_back(SplitTuple5(std::min(angle0, angle1), currId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		uint32_t startId = std::get<1>(splitVec[0]);
		outV2Id_ = std::get<2>(splitVec[0]);
		outV3Id_ = std::get<3>(splitVec[0]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[0]);
		pv.pos[1] = std::get<5>(splitVec[0]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);

		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	//第四步
	std::function<bool(double)> isAngleAfter3 = [&](double angle)->bool
	{
		if (angle > AFTER_ANGLE_MIN_THREE && angle < AFTER_ANGLE_MAX_THREE)
			return true;
		else
			return false;
	};
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();

	while (!queueFaceId.empty())
	{
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);
		std::vector<bool> currPolyAddNew;
		SetListToVec(currPolyAddNew, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35464);

		GetConcaveVerticesLarger(currPolys, ccVsIds, currPolyAddNew);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		splitVec.reserve(ccVsIds.size());
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			/*for (int j = 0; j < 2; ++j)
			{*/
			int vPreId = mmod((int)ccVsIds[i] - 1, (int)currPolys.size()), vCurrId = ccVsIds[i], vNextId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());
			/*if (j == 1)
				vPreId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());*/

				/*if (ccVsIds.size() == 1)
					DO_NOTHING_ANGLE_THRESHOLD = 20.0*PI / 180.0;
				else
					DO_NOTHING_ANGLE_THRESHOLD = 10.0*PI / 180.0;*/
					/*double concaveAngle0 = ComputePolyAngle(currPolys, ccVsIds[i]);
					if (std::abs(concaveAngle0 - PI) < DO_NOTHING_ANGLE_THRESHOLD)
						continue;*/

			Eigen::Vector2d aveAngleVec = (Vs_[currPolys[vPreId]].pos - Vs_[currPolys[vCurrId]].pos).normalized() + (Vs_[currPolys[vNextId]].pos - Vs_[currPolys[vCurrId]].pos).normalized();
			if (!GetPolyRayIntersection(currPolys, vCurrId, aveAngleVec, outV2Id_, outV3Id_, interPos_, outPoly0_, outPoly1_))
				continue;

			Eigen::Vector2d midPos = (Vs_[currPolys[vCurrId]].pos + interPos_) * 0.5;
			if (!CheckPointInsidePoly(currPolys, midPos))
				continue;
			//是否不允许出现三角形？
			if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
				continue;

			/*double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[outV2Id_]].pos - interPos_).norm(),
				newLength2 = (Vs_[currPolys[outV3Id_]].pos - interPos_).norm();*/
			int realOutV2Id = FindSmallerOldVertexId(outV2Id_, currPolyAddNew), realOutV3Id = FindLargerOldVertexId(outV3Id_, currPolyAddNew);
			double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[realOutV2Id]].pos - interPos_).norm(),
				newLength2 = (Vs_[currPolys[realOutV3Id]].pos - interPos_).norm();
			std::vector<bool> poly0AddNew, poly1AddNew;
			GetSubAddNew(vCurrId, outV2Id_, outV3Id_, currPolyAddNew, poly0AddNew, poly1AddNew);
			if (GetFalseNum(poly0AddNew) <= 3 || GetFalseNum(poly1AddNew) <= 3)
				continue;

			if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
				continue;

			double outArea0 = std::abs(ComputePolyArea(outPoly0_));
			double outArea1 = std::abs(ComputePolyArea(outPoly1_));
			if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
				continue;

			double angle0 = ComputePolyAngle(outPoly0_, outPoly0_.size() - 1), angle1 = ComputePolyAngle(outPoly1_, 0);
			if (!isAngleAfter3(angle0) || !isAngleAfter3(angle1))
				continue;

			double minEdge0 = 1.0, maxEdge0 = 1.0, minEdge1 = 1.0, maxEdge1 = 1.0;
			/*ComputeShortestAndLongestEdge(outPoly0_, minEdge0, maxEdge0);
			ComputeShortestAndLongestEdge(outPoly1_, minEdge1, maxEdge1);*/
			ComputeShortestAndLongestEdge(outPoly0_, poly0AddNew, minEdge0, maxEdge0);
			ComputeShortestAndLongestEdge(outPoly1_, poly1AddNew, minEdge1, maxEdge1);
			double minArea = std::min(outArea0, outArea1), maxArea = std::max(outArea0, outArea1);
			double ranking = minArea / maxArea + AREA_LENGTH_WEIGHT * 0.5 * (std::pow(minEdge0 / maxEdge0, 2) + std::pow(minEdge1 / maxEdge1, 2));

			splitVec.emplace_back(SplitTuple5(ranking, vCurrId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			//}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		uint32_t startId = std::get<1>(splitVec[0]);
		outV2Id_ = std::get<2>(splitVec[0]);
		outV3Id_ = std::get<3>(splitVec[0]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[0]);
		pv.pos[1] = std::get<5>(splitVec[0]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);

		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	BuildTopologyInfoFromVsFs();
}

void PolySplit::ComputeTargetSubQuadsNum(std::vector<std::vector<Eigen::Vector2d>> &polys, uint32_t totalQuadsNum, std::vector<int> &subNum)
{
	std::vector<double> subPolyAreas;
	subPolyAreas.reserve(polys.size());
	subNum.reserve(polys.size());

	for (int i = 0; i < polys.size(); ++i)
	{
		double currArea = std::abs(ComputePolyArea(polys[i]));
		subPolyAreas.emplace_back(currArea);
	}

	for (int i = 0; i < polys.size(); ++i)
	{
		subNum.emplace_back((int)(totalQuadsNum * subPolyAreas[i] / sourcePolyArea_));
	}
}

void PolySplit::OutputPolysInfo(std::string &fileName)
{
	std::ofstream ofs(fileName.c_str());

	ofs << Vs_.size() << std::endl;

	for (int i = 0; i < Vs_.size(); ++i)
	{
		ofs << Vs_[i].pos[0] << " " << Vs_[i].pos[1] << std::endl;
	}

	ofs << Fs_.size() << std::endl;

	for (int i = 0; i < Fs_.size(); ++i)
	{
		ofs << Fs_[i].pvs.size() << " ";
		for (auto it = Fs_[i].pvs.begin(); it != Fs_[i].pvs.end(); ++it)
		{
			ofs << *it << " ";
		}
		ofs << std::endl;
	}
}

void PolySplit::OutputSplitEdgesInfo(std::string &fileName)
{
	if (splitEdges_.empty())
		return;

	std::ofstream ofs(fileName.c_str());

	ofs << splitEdges_.size() << std::endl;
	for (int i = 0; i < splitEdges_.size(); ++i)
	{
		Eigen::Vector2d &vPos0 = std::get<0>(splitEdges_[i]), &vPos1 = std::get<1>(splitEdges_[i]);
		ofs << vPos0[0] << " " << vPos0[1] << " " << vPos1[0] << " " << vPos1[1] << std::endl;
	}
}

void PolySplit::SplitPoly(std::vector<Eigen::Vector2d> &inPoly)
{
	if (!JudgePolyCCW(inPoly))
		std::reverse(inPoly.begin(), inPoly.end());
	InitializePoly(inPoly);

	currFaceId_.clear();
	
	std::queue<uint32_t> queueFaceId;
	queueFaceId.emplace(0);

	//第一步
	std::vector<uint32_t> ccVsIds;
	while (!queueFaceId.empty())
	{
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35463);

		std::vector<TupleRankingLine> fieldRanking;
		GetConcaveVerticesLarger(currPolys, ccVsIds);
		/*if (ccVsIds.size() >= 2)
			ANGLE_THRESHOLD = 45.0 * PI / 180.0;
		else
			ANGLE_THRESHOLD = 75.0 * PI / 180.0;*/
		
		/*std::vector<uint32_t> ccVsIdsLarger;
		GetConcaveVerticesLarger(currPolys, ccVsIdsLarger);*/
		std::vector<bool> concaveVFlag(currPolys.size(), false);
		for (int i = 0; i < ccVsIds/*Larger*/.size(); ++i)
		{
			concaveVFlag[ccVsIds/*Larger*/[i]] = true;
		}
		
		for (int i = 0; i < ccVsIds.size(); ++i)	//第一个点必须是concave点，但第二个点是所有点
		{
			/*double concaveAngle0 = ComputePolyAngle(currPolys, ccVsIds[i]);*/
			/*if (ccVsIds.size() == 1)
				DO_NOTHING_ANGLE_THRESHOLD = 20.0*PI / 180.0;
			else
				DO_NOTHING_ANGLE_THRESHOLD = 10.0*PI / 180.0;*/
			/*if (std::abs(concaveAngle0 - PI) < DO_NOTHING_ANGLE_THRESHOLD)
				continue;*/

			uint32_t v0 = currPolys[ccVsIds[i]];
			for (int j = 0; j < currPolys.size(); ++j)
			{
				if (j == ccVsIds[i] || std::abs((int)(ccVsIds[i]) - (int)(j)) == 1)
					continue;

				uint32_t v1 = currPolys[j];
				if (JudgeIfEdgeIntersectionOrOutside(currPolys, v0, v1))
					continue;

				double currFieldMatchAngle = LineFieldMatchAngle(v0, v1);
				bool isConcave = false;
				if (concaveVFlag[j])
				{
					currFieldMatchAngle -= 1E5;
					isConcave = true;
				}
				else
					isConcave = false;
				fieldRanking.emplace_back(TupleRankingLine(currFieldMatchAngle, v0, v1, ccVsIds[i], j, isConcave));
			}
		}
		if (fieldRanking.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		bool isFindGoodSplit = false;
		std::vector<uint32_t> newPoly0, newPoly1;
		std::sort(fieldRanking.begin(), fieldRanking.end());
		for (int i = 0; i < fieldRanking.size(); ++i)
		{
			double fieldValue = std::get<0>(fieldRanking[i]);
			uint32_t v0 = std::get<1>(fieldRanking[i]), v1 = std::get<2>(fieldRanking[i]), v0Id = std::get<3>(fieldRanking[i]), v1Id = std::get<4>(fieldRanking[i]);
			bool isConcave = std::get<5>(fieldRanking[i]);

			ANGLE_THRESHOLD = 75.0 * PI / 180.0;
			AREA_THRESHOLD = 0.003;
			/*if (isConcave)
			{
				ANGLE_THRESHOLD = 45.0 * PI / 180.0;
				AREA_THRESHOLD = 0.005;
			}
			else
			{
				ANGLE_THRESHOLD = 75.0 * PI / 180.0;
				AREA_THRESHOLD = 0.005;
			}*/

			/*if (fieldValue > FIELD_MATCHING_THRESHOLD)
				continue;*/
			double newLength = (Vs_[v0].pos - Vs_[v1].pos).norm();
			if (newLength < LENGTH_THRESHOLD * diagLength_)
				continue;

			SplitPolyIntoTwo(currPolys, v0Id, v1Id, newPoly0, newPoly1);

			//是否不允许出现三角形？
			if (newPoly0.size() <= 3 || newPoly1.size() <= 3)
				continue;

			double poly0Area = std::abs(ComputePolyArea(newPoly0)), poly1Area = std::abs(ComputePolyArea(newPoly1));
			if (poly0Area < AREA_THRESHOLD * sourcePolyArea_ || poly1Area < AREA_THRESHOLD * sourcePolyArea_)
				continue;

			double angle0 = ComputePolyAngle(newPoly0, 0), angle1 = ComputePolyAngle(newPoly0, newPoly0.size() - 1), angle2 = ComputePolyAngle(newPoly1, 0), angle3 = ComputePolyAngle(newPoly1, newPoly1.size() - 1);
			if (angle0 < ANGLE_THRESHOLD || angle1 < ANGLE_THRESHOLD || angle2 < ANGLE_THRESHOLD || angle3 < ANGLE_THRESHOLD)
				continue;

			isFindGoodSplit = true;

			//更新Vs_, Fs_拓扑信息
			Vs_[currPolys[v0Id]].neighbor_pfs.emplace_back(Fs_.size());
			Vs_[currPolys[v1Id]].neighbor_pfs.emplace_back(Fs_.size());
			Fs_[currFaceId].pvs.clear();
			Fs_[currFaceId].pvs_new_add.clear();
			for (int j = 0; j < newPoly0.size(); ++j)
			{
				Fs_[currFaceId].pvs.emplace_back(newPoly0[j]);
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			}

			for (int j = 0; j < newPoly1.size(); ++j)
			{
				if (newPoly1[j] != currPolys[v0Id] && newPoly1[j] != currPolys[v1Id])
				{
					std::vector<uint32_t> &vfs = Vs_[newPoly1[j]].neighbor_pfs;
					for (int k = 0; k < vfs.size(); ++k)
					{
						if (vfs[k] == currFaceId)
							vfs[k] = Fs_.size();
					}
				}
			}
			PolyFace pf;
			pf.id = Fs_.size();
			//pf.pvs.reserve(newPoly1.size());
			for (int j = 0; j < newPoly1.size(); ++j)
			{
				pf.pvs.emplace_back(newPoly1[j]);
				pf.pvs_new_add.emplace_back(false);
			}
			Fs_.emplace_back(pf);
			break;
		}

		if (isFindGoodSplit)
		{
			queueFaceId.emplace(currFaceId);
			queueFaceId.emplace(Fs_.size() - 1);
		}
		else
		{
			currFaceId_.emplace_back(currFaceId);
		}
	}

	//第二步
	std::function<bool(double)> isAngleAfter = [&](double angle)->bool
	{
		if (angle > AFTER_ANGLE_MIN && angle < AFTER_ANGLE_MAX)
			return true;
		else
			return false;
	};

	typedef std::tuple<double, uint32_t, uint32_t, uint32_t, double, double> SplitTuple5;	//五个元素分别为分割后两个面积较小的那个；起始v, 分割点的左右两个v；分割点的pos
	std::vector<SplitTuple5> splitVec;
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();
	while (!queueFaceId.empty())
	{
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35465);

		GetConcaveVertices(currPolys, ccVsIds);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		splitVec.reserve(ccVsIds.size());
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				int vPreId = mmod((int)ccVsIds[i] - 1, (int)currPolys.size()), vCurrId = ccVsIds[i];
				if (j == 1)
					vPreId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());

				double sourceAngle = ComputePolyAngle(currPolys, ccVsIds[i]);
				if ((sourceAngle<BEFORE_ANGLE_MIN_0 || sourceAngle>BEFORE_ANGLE_MAX_0) && sourceAngle < BEFORE_ANGLE_MIN_1)
					continue;

				/*if (sourceAngle > BEFORE_ANGLE_MIN_1)
				{
					AFTER_ANGLE_MIN = 60.0 * PI / 180.0;
					AFTER_ANGLE_MAX = 120.0 * PI / 180.0;
				}
				else
				{
					AFTER_ANGLE_MIN = 70.0 * PI / 180.0;
					AFTER_ANGLE_MAX = 110.0 * PI / 180.0;
				}*/
				AFTER_ANGLE_MIN = 70.0 * PI / 180.0;
				AFTER_ANGLE_MAX = 110.0 * PI / 180.0;

				if (!GetPolyRayIntersection(currPolys, vCurrId, vPreId, outV2Id_, outV3Id_, interPos_, outPoly0_, outPoly1_))
					continue;
				//是否不允许出现三角形？
				if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
					continue;

				double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[outV2Id_]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[outV3Id_]].pos - interPos_).norm();

				if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
					continue;

				double outArea0 = std::abs(ComputePolyArea(outPoly0_));
				double outArea1 = std::abs(ComputePolyArea(outPoly1_));
				if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
					continue;

				double angle0 = ComputePolyAngle(outPoly0_, outPoly0_.size() - 1), angle1 = ComputePolyAngle(outPoly1_, 0);
				if (!isAngleAfter(angle0) || !isAngleAfter(angle1))
					continue;

				double minEdge0 = 1.0, maxEdge0 = 1.0, minEdge1 = 1.0, maxEdge1 = 1.0;
				ComputeShortestAndLongestEdge(outPoly0_, minEdge0, maxEdge0);
				ComputeShortestAndLongestEdge(outPoly1_, minEdge1, maxEdge1);
				double minArea = std::min(outArea0, outArea1), maxArea = std::max(outArea0, outArea1);
				double ranking = minArea / maxArea + AREA_LENGTH_WEIGHT * 0.5 * (std::pow(minEdge0 / maxEdge0, 2) + std::pow(minEdge1 / maxEdge1, 2));

				splitVec.emplace_back(SplitTuple5(ranking, vCurrId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		uint32_t startId = std::get<1>(splitVec[0]);
		outV2Id_ = std::get<2>(splitVec[0]);
		outV3Id_ = std::get<3>(splitVec[0]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[0]);
		pv.pos[1] = std::get<5>(splitVec[0]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii==startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);

		
		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	//第三步
	std::function<bool(double)> isAngleAfter2 = [&](double angle)->bool
	{
		if (std::abs(angle) > 45.0*PI/180.0)
			return true;
		else
			return false;
	};
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();

	while (!queueFaceId.empty())
	{
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35464);

		GetConcaveVerticesLarger(currPolys, ccVsIds);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			uint32_t currId = ccVsIds[i];
			Eigen::Vector2d currPos = Vs_[currPolys[currId]].pos;

			for (int j = 0; j < currPolys.size(); ++j)
			{
				if (j == currId || ((j + 1) % currPolys.size()) == currId)
					continue;

				Eigen::Vector2d vPos0 = Vs_[currPolys[j]].pos, vPos1 = Vs_[currPolys[(j + 1) % currPolys.size()]].pos;
				double alpha = 0;
				Eigen::Vector2d interPos;
				if (!ComputeVerticalIntersection(currPos, vPos0, vPos1, alpha, interPos))
					continue;

				if (JudgeIfVerticalIntersection(currPolys, currId, interPos, j, (j + 1) % currPolys.size()))
					continue;

				SplitPolyIntoTwo(currPolys, currId, interPos, j, (j + 1) % currPolys.size(), outPoly0_, outPoly1_);
				//是否不允许出现三角形？
				if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
					continue;
				outV2Id_ = j; outV3Id_ = (j + 1) % currPolys.size(); interPos_ = interPos;
				double newLength0 = (Vs_[currPolys[currId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[outV2Id_]].pos - interPos_).norm(),
					newLength2 = (Vs_[currPolys[outV3Id_]].pos - interPos_).norm();

				if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
					continue;

				double outArea0 = std::abs(ComputePolyArea(outPoly0_));
				double outArea1 = std::abs(ComputePolyArea(outPoly1_));
				if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
					continue;

				double angle0 = ComputePolyAngle(outPoly0_, 0), angle1 = ComputePolyAngle(outPoly1_, outPoly1_.size() - 1);
				if (!isAngleAfter2(angle0) || !isAngleAfter2(angle1))
					continue;

				splitVec.emplace_back(SplitTuple5(std::min(angle0, angle1), currId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		uint32_t startId = std::get<1>(splitVec[0]);
		outV2Id_ = std::get<2>(splitVec[0]);
		outV3Id_ = std::get<3>(splitVec[0]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[0]);
		pv.pos[1] = std::get<5>(splitVec[0]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);

		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	//第四步
	std::function<bool(double)> isAngleAfter3 = [&](double angle)->bool
	{
		if (angle > AFTER_ANGLE_MIN_THREE && angle < AFTER_ANGLE_MAX_THREE)
			return true;
		else
			return false;
	};
	std::queue<uint32_t>().swap(queueFaceId);
	for (int i = 0; i < currFaceId_.size(); ++i)
	{
		queueFaceId.emplace(currFaceId_[i]);
	}
	currFaceId_.clear();

	while (!queueFaceId.empty())
	{
		uint32_t currFaceId = queueFaceId.front();
		queueFaceId.pop();
		std::vector<uint32_t> currPolys;
		SetListToVec(currPolys, Fs_[currFaceId].pvs);
		std::vector<bool> currNewAdd;
		SetListToVec(currNewAdd, Fs_[currFaceId].pvs_new_add);

		if (currPolys.size() == 3)
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}
		else if (currPolys.size() < 3)
			exit(-35464);

		GetConcaveVerticesLarger(currPolys, ccVsIds);
		if (ccVsIds.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		splitVec.clear();
		splitVec.reserve(ccVsIds.size());
		for (int i = 0; i < ccVsIds.size(); ++i)
		{
			/*for (int j = 0; j < 2; ++j)
			{*/
			int vPreId = mmod((int)ccVsIds[i] - 1, (int)currPolys.size()), vCurrId = ccVsIds[i], vNextId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());
			/*if (j == 1)
				vPreId = mmod((int)ccVsIds[i] + 1, (int)currPolys.size());*/

				/*if (ccVsIds.size() == 1)
					DO_NOTHING_ANGLE_THRESHOLD = 20.0*PI / 180.0;
				else
					DO_NOTHING_ANGLE_THRESHOLD = 10.0*PI / 180.0;*/
					/*double concaveAngle0 = ComputePolyAngle(currPolys, ccVsIds[i]);
					if (std::abs(concaveAngle0 - PI) < DO_NOTHING_ANGLE_THRESHOLD)
						continue;*/

			Eigen::Vector2d aveAngleVec = (Vs_[currPolys[vPreId]].pos - Vs_[currPolys[vCurrId]].pos).normalized() + (Vs_[currPolys[vNextId]].pos - Vs_[currPolys[vCurrId]].pos).normalized();
			if (!GetPolyRayIntersection(currPolys, vCurrId, aveAngleVec, outV2Id_, outV3Id_, interPos_, outPoly0_, outPoly1_))
				continue;
			//是否不允许出现三角形？
			if (outPoly0_.size() <= 3 || outPoly1_.size() <= 3)
				continue;

			double newLength0 = (Vs_[currPolys[vCurrId]].pos - interPos_).norm(), newLength1 = (Vs_[currPolys[outV2Id_]].pos - interPos_).norm(),
				newLength2 = (Vs_[currPolys[outV3Id_]].pos - interPos_).norm();

			if (newLength0 < LENGTH_THRESHOLD * diagLength_ || newLength1 < LENGTH_THRESHOLD * diagLength_ || newLength2 < LENGTH_THRESHOLD * diagLength_)
				continue;

			double outArea0 = std::abs(ComputePolyArea(outPoly0_));
			double outArea1 = std::abs(ComputePolyArea(outPoly1_));
			if (outArea0 < AREA_THRESHOLD_2 * sourcePolyArea_ || outArea1 < AREA_THRESHOLD_2 * sourcePolyArea_)
				continue;

			double angle0 = ComputePolyAngle(outPoly0_, outPoly0_.size() - 1), angle1 = ComputePolyAngle(outPoly1_, 0);
			if (!isAngleAfter3(angle0) || !isAngleAfter3(angle1))
				continue;

			double minEdge0 = 1.0, maxEdge0 = 1.0, minEdge1 = 1.0, maxEdge1 = 1.0;
			ComputeShortestAndLongestEdge(outPoly0_, minEdge0, maxEdge0);
			ComputeShortestAndLongestEdge(outPoly1_, minEdge1, maxEdge1);
			double minArea = std::min(outArea0, outArea1), maxArea = std::max(outArea0, outArea1);
			double ranking = minArea / maxArea + AREA_LENGTH_WEIGHT * 0.5 * (std::pow(minEdge0 / maxEdge0, 2) + std::pow(minEdge1 / maxEdge1, 2));

			splitVec.emplace_back(SplitTuple5(ranking, vCurrId, outV2Id_, outV3Id_, interPos_[0], interPos_[1]));
			//}
		}
		if (splitVec.empty())
		{
			currFaceId_.emplace_back(currFaceId);
			continue;
		}

		std::sort(splitVec.begin(), splitVec.end());
		std::reverse(splitVec.begin(), splitVec.end());

		uint32_t startId = std::get<1>(splitVec[0]);
		outV2Id_ = std::get<2>(splitVec[0]);
		outV3Id_ = std::get<3>(splitVec[0]);

		PolyVertex pv;
		pv.id = Vs_.size();
		pv.pos[0] = std::get<4>(splitVec[0]);
		pv.pos[1] = std::get<5>(splitVec[0]);
		pv.neighbor_pfs.emplace_back(currFaceId);
		pv.neighbor_pfs.emplace_back(Fs_.size());

		Vs_[currPolys[startId]].neighbor_pfs.emplace_back(Fs_.size());

		//获取这个T-junction的另一个面，并修改拓扑
		uint32_t outV2 = currPolys[outV2Id_], outV3 = currPolys[outV3Id_];
		std::vector<uint32_t> &neiPfsV2 = Vs_[outV2].neighbor_pfs, &neiPfsV3 = Vs_[outV3].neighbor_pfs;
		std::vector<uint32_t> interFs;
		std::sort(neiPfsV2.begin(), neiPfsV2.end()); std::sort(neiPfsV3.begin(), neiPfsV3.end());
		std::set_intersection(neiPfsV2.begin(), neiPfsV2.end(), neiPfsV3.begin(), neiPfsV3.end(), std::inserter(interFs, interFs.begin()));

		uint32_t anotherFace = (uint32_t)-1;
		for (int i = 0; i < interFs.size(); ++i)
		{
			if (interFs[i] != currFaceId)
			{
				anotherFace = interFs[i];
				break;
			}
		}
		if (anotherFace != (uint32_t)-1)
		{
			pv.neighbor_pfs.emplace_back(anotherFace);
			std::list<uint32_t>::iterator it;
			int dis = 0;
			FindListInsertPos(Fs_[anotherFace].pvs, outV2, outV3, it, dis);
			Fs_[anotherFace].pvs.insert(it, pv.id);
			auto itt = Fs_[anotherFace].pvs_new_add.begin();
			std::advance(itt, dis);
			Fs_[anotherFace].pvs_new_add.insert(itt, true);
		}
		Vs_.emplace_back(pv);

		//Build new Poly
		Fs_[currFaceId].pvs.clear();
		Fs_[currFaceId].pvs_new_add.clear();
		int ii = startId;
		while (ii != outV2Id_ && ii != outV3Id_)
		{
			Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
			if (ii == startId)
				Fs_[currFaceId].pvs_new_add.emplace_back(false);
			else
				Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
			ii = (ii + 1) % currPolys.size();
		}
		Fs_[currFaceId].pvs.emplace_back(currPolys[ii]);
		Fs_[currFaceId].pvs_new_add.emplace_back(currNewAdd[ii]);
		Fs_[currFaceId].pvs.emplace_back(pv.id);
		Fs_[currFaceId].pvs_new_add.emplace_back(false);

		PolyFace pf;
		pf.id = Fs_.size();
		pf.pvs.emplace_back(pv.id);
		pf.pvs_new_add.emplace_back(false);
		if (ii == outV2Id_)
			ii = outV3Id_;
		else
			ii = outV2Id_;

		while (ii != startId)
		{
			pf.pvs.emplace_back(currPolys[ii]);
			pf.pvs_new_add.emplace_back(currNewAdd[ii]);

			std::vector<uint32_t> &vfs = Vs_[currPolys[ii]].neighbor_pfs;
			for (int j = 0; j < vfs.size(); ++j)
			{
				if (vfs[j] == currFaceId)
					vfs[j] = Fs_.size();
			}
			ii = (ii + 1) % currPolys.size();
		}
		pf.pvs.emplace_back(currPolys[ii]);
		pf.pvs_new_add.emplace_back(false);
		Fs_.emplace_back(pf);

		queueFaceId.emplace(currFaceId);
		queueFaceId.emplace(Fs_.size() - 1);
	}

	BuildTopologyInfoFromVsFs();
}