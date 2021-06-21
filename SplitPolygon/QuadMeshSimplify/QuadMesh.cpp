#include "QuadMesh.h"
#include "AssistFunc.h"
#include <algorithm>
#include <functional>

namespace BaseDataStructure
{
	QuadMesh::QuadMesh()
	{
	}


	QuadMesh::~QuadMesh()
	{
	}

	bool QuadMesh::BuildConnectivity()
	{
		if (Vs_.empty() || Fs_.empty())
			return false;
		OritateFace();

		typedef std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t> ITuple5;

		//build Es
		if (!Es_.empty())
		{
			//std::vector<QuadEdge>().swap(Es_);
			Es_.clear();
			Es_.reserve(Fs_.size() * 4);
		}
		std::vector<ITuple5> repeatEdges;
		repeatEdges.reserve(Fs_.size() * 4);
		for (int i = 0; i < Fs_.size(); ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				uint32_t vIndex0 = Fs_[i].vs[j]; uint32_t vIndex1 = Fs_[i].vs[(j + 1) % 4];
				if (vIndex0 > vIndex1)		//make sure vIndex0 is smaller than vIndex1.
					std::swap(vIndex0, vIndex1);

				repeatEdges.emplace_back(ITuple5(vIndex0, vIndex1, i * 4 + j, i, j));
			}
		}

		std::sort(repeatEdges.begin(), repeatEdges.end());
		QuadEdge currQE;
		currQE.id = 0;
		currQE.bId = (uint32_t)-1;
		currQE.boundary = true;
		uint32_t v0Index = std::get<0>(repeatEdges[0]);
		uint32_t v1Index = std::get<1>(repeatEdges[0]);
		uint32_t fIndex = std::get<3>(repeatEdges[0]);
		currQE.vs.emplace_back(v0Index);		//e - v
		currQE.vs.emplace_back(v1Index);
		currQE.neighbor_fs.emplace_back(fIndex);		//e - f
		Vs_[v0Index].neighbor_es.emplace_back(0); Vs_[v1Index].neighbor_es.emplace_back(0);		//v - e
		Fs_[fIndex].es.emplace_back(0);	// f - e
		Es_.emplace_back(currQE);
		uint32_t eNum = 1;
		for (int i = 1; i < repeatEdges.size(); ++i)
		{
			if (std::get<0>(repeatEdges[i]) != std::get<0>(repeatEdges[i - 1]) || std::get<1>(repeatEdges[i]) != std::get<1>(repeatEdges[i - 1]))
			{
				currQE.id = eNum++;
				currQE.bId = (uint32_t)-1;
				currQE.boundary = true;
				uint32_t v0Index = std::get<0>(repeatEdges[i]);
				uint32_t v1Index = std::get<1>(repeatEdges[i]);
				uint32_t fIndex = std::get<3>(repeatEdges[i]);
				currQE.vs[0] = v0Index; currQE.vs[1] = v1Index;
				currQE.neighbor_fs.clear(); currQE.neighbor_fs.emplace_back(fIndex);
				Vs_[v0Index].neighbor_es.emplace_back(eNum-1); Vs_[v1Index].neighbor_es.emplace_back(eNum-1);		//v - e
				Fs_[fIndex].es.emplace_back(eNum-1);	// f - e
				Es_.emplace_back(currQE);
			}
			else
			{
				uint32_t currEId = eNum - 1;
				Es_[currEId].boundary = false;
				uint32_t fIndex = std::get<3>(repeatEdges[i]);
				Es_[currEId].neighbor_fs.emplace_back(fIndex);
				Fs_[fIndex].es.emplace_back(currEId);
			}
		}

		if (JudgeIfTopologyProblem())
		{
			std::cout << "Topology Problem!" << std::endl;
			return false;
			//std::cout << "Here!" << std::endl;
			//exit(-1);
		}

		for (int i = 0; i < Es_.size(); ++i)
		{
			uint32_t vIndex0 = Es_[i].vs[0];
			uint32_t vIndex1 = Es_[i].vs[1];
			uint32_t fIndex = Es_[i].neighbor_fs[0];
			//boundary vs & fs
			if (Es_[i].boundary)
			{
				Vs_[vIndex0].boundary = true;
				Vs_[vIndex1].boundary = true;
				Fs_[fIndex].boundary = true;
			}

			//v-neighbor_vs
			Vs_[vIndex0].neighbor_vs.emplace_back(vIndex1);
			Vs_[vIndex1].neighbor_vs.emplace_back(vIndex0);

			//f-neighbor_fs
			if (!Es_[i].boundary)
			{
				uint32_t fIndex1 = Es_[i].neighbor_fs[1];
				Fs_[fIndex].neighbor_fs.emplace_back(fIndex1);
				Fs_[fIndex1].neighbor_fs.emplace_back(fIndex);
			}
		}

		//e-neighbor_es
		for (int i = 0; i < Vs_.size(); ++i)
		{
			for (int j = 0; j < Vs_[i].neighbor_es.size(); ++j)
			{
				for (int k = 0; k < Vs_[i].neighbor_es.size(); ++k)
				{
					if (j != k)
					{
						Es_[Vs_[i].neighbor_es[j]].neighbor_es.emplace_back(Vs_[i].neighbor_es[k]);
					}
				}
			}
		}

		return true;
		//OrderElements();
	}

	void QuadMesh::ComputeMeshNormal()
	{
		for (int i = 0; i < Fs_.size(); ++i)
		{
			std::vector<uint32_t> &fvs = Fs_[i].vs;
			Eigen::Vector2d vPos0 = V_[fvs[0]];
			Eigen::Vector2d vPos1 = V_[fvs[1]];
			Eigen::Vector2d vPos2 = V_[fvs[2]];

			Eigen::Vector3d vec0(vPos2[0] - vPos1[0], vPos2[1] - vPos1[1], 0);
			Eigen::Vector3d vec1(vPos0[0] - vPos1[0], vPos0[1] - vPos1[1], 0);

			Eigen::Vector3d crossVec = vec0.cross(vec1);
			if (crossVec.norm() > 0.1)
			{
				meshNormal_ = crossVec;
				break;
			}
		}
	}

	double QuadMesh::ComputeAverageEdgeLength()
	{
		double sumLength = 0;
		for (int i = 0; i < Es_.size(); ++i)
		{
			uint32_t v0 = Es_[i].vs[0];
			uint32_t v1 = Es_[i].vs[1];
			sumLength += (V_[v0] - V_[v1]).norm();
		}
		return sumLength / Es_.size();
	}

	void QuadMesh::OritateFace()
	{
		if (Fs_.empty())
			return;

		uint32_t p0Index = Fs_[0].vs[0]; uint32_t p1Index = Fs_[0].vs[1]; uint32_t p3Index = Fs_[0].vs[3];
		Eigen::Vector3d p0(V_[p0Index][0], V_[p0Index][1], 0.0);
		Eigen::Vector3d p1(V_[p1Index][0], V_[p1Index][1], 0.0);
		Eigen::Vector3d p3(V_[p3Index][0], V_[p3Index][1], 0.0);

		Eigen::Vector3d normalVec = (p1 - p0).cross(p3 - p0).normalized();
		if (normalVec.dot(Eigen::Vector3d(0, 0, 1)) < 0)
		{
			for (int i = 0; i < Fs_.size(); ++i)
			{
				std::reverse(Fs_[i].vs.begin(), Fs_[i].vs.end());
			}
		}
	}

	bool QuadMesh::JudgeIfTopologyProblem() const
	{
		for (int i = 0; i < Vs_.size(); ++i)
		{
			if (Vs_[i].neighbor_fs.empty())
				return true;
		}
		return false;
	}

	void QuadMesh::ClearMesh()
	{
		/*std::vector<Eigen::Vector2d>().swap(V_);
		std::vector<QuadVertex>().swap(Vs_);
		std::vector<QuadEdge>().swap(Es_);
		std::vector<QuadFace>().swap(Fs_);*/
		V_.clear();
		Vs_.clear();
		Es_.clear();
		Fs_.clear();
	}

	//e - neighbor_es not do. If you want, please use v - es
	void QuadMesh::OrderElements()
	{
		std::function<uint32_t(uint32_t, uint32_t, const std::vector<uint32_t> &)> findCurrE = [&](uint32_t v0, uint32_t v1, const std::vector<uint32_t> &fes)->uint32_t
		{
			for (int j = 0; j < fes.size(); ++j)
			{
				uint32_t vIndex0 = Es_[fes[j]].vs[0];
				uint32_t vIndex1 = Es_[fes[j]].vs[1];

				if ((vIndex0 == v0 && vIndex1 == v1) || (vIndex0 == v1 && vIndex1 == v0))
					return fes[j];
			}
			exit(-1);
		};

		//Order f - es
		for (int i = 0; i < Fs_.size(); ++i)
		{
			const std::vector<uint32_t> &fvs = Fs_[i].vs;
			const std::vector<uint32_t> &fes = Fs_[i].es;

			uint32_t currE[4];
			for (int j = 0; j < fvs.size(); ++j)
			{
				uint32_t vIndex0 = fvs[j];
				uint32_t vIndex1 = fvs[(j + 1) % fvs.size()];
				currE[j] = findCurrE(vIndex0, vIndex1, fes);
			}
			for (int j = 0; j < Fs_[i].es.size(); ++j)
			{
				Fs_[i].es[j] = currE[j];
			}
		}

		//Order f - fs
		for (int i = 0; i < Fs_.size(); ++i)
		{
		std::vector<uint32_t> &fEs = Fs_[i].es;
		std::vector<uint32_t> tempFs;
		for (int j = 0; j < fEs.size(); ++j)
		{
			const std::vector<uint32_t> &eNeighborF = Es_[fEs[j]].neighbor_fs;
			if (eNeighborF.size() == 2)
			{
				if (eNeighborF[0] == i)
					tempFs.emplace_back(eNeighborF[1]);
				else
					tempFs.emplace_back(eNeighborF[0]);
			}
		}
		std::copy(tempFs.begin(), tempFs.end(), Fs_[i].neighbor_fs.begin());
		}
	}

	bool QuadMesh::JudgeQuadMeshJacobi()
	{
		if (Vs_.empty() || Fs_.empty())
			exit(-5);

		minJacobi_ = 10000;
		double startSign = 1;
		uint32_t v0, v1, v2;
		Eigen::Vector2d vec0, vec1;
		for (int i = 0; i < Fs_.size(); ++i)
		{
			v0 = Fs_[i].vs[quadTriTable[0][0]];
			v1 = Fs_[i].vs[quadTriTable[0][1]];
			v2 = Fs_[i].vs[quadTriTable[0][2]];
			vec0 = (V_[v2] - V_[v1]).normalized();
			vec1 = (V_[v0] - V_[v1]).normalized();
			startSign = vec0[0] * vec1[1] - vec0[1] * vec1[0];
			if (abs(startSign) < 0.1)
				continue;
			break;
		}

		for (int i = 0; i < Fs_.size(); ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				v0 = Fs_[i].vs[quadTriTable[j][0]];
				v1 = Fs_[i].vs[quadTriTable[j][1]];
				v2 = Fs_[i].vs[quadTriTable[j][2]];
				vec0 = (V_[v2] - V_[v1]).normalized();
				vec1 = (V_[v0] - V_[v1]).normalized();
				double currSign = vec0[0] * vec1[1] - vec0[1] * vec1[0];
				/*if (std::abs(currSign) < minJacobi_ && JudgeSign(startSign, currSign))
					minJacobi_ = std::abs(currSign);*/
				if (!JudgeSign(startSign, currSign))
				{
					//std::cout << i << std::endl;
					return false;
				}
				else if (std::abs(currSign) < minJacobi_)
					minJacobi_ = std::abs(currSign);

			}
		}
		return true;
	}

	bool QuadMesh::JudgeQuadMeshJacobiWithoutBoundary()
	{
		if (Vs_.empty() || Fs_.empty())
			exit(-5);

		minJacobi_ = 10000;
		double startSign = 1;
		uint32_t v0, v1, v2;
		Eigen::Vector2d vec0, vec1;
		for (int i = 0; i < Fs_.size(); ++i)
		{
			v0 = Fs_[i].vs[quadTriTable[0][0]];
			v1 = Fs_[i].vs[quadTriTable[0][1]];
			v2 = Fs_[i].vs[quadTriTable[0][2]];
			vec0 = (V_[v2] - V_[v1]).normalized();
			vec1 = (V_[v0] - V_[v1]).normalized();
			startSign = vec0[0] * vec1[1] - vec0[1] * vec1[0];
			if (abs(startSign) < 0.1)
				continue;
			break;
		}

		for (int i = 0; i < Fs_.size(); ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				v0 = Fs_[i].vs[quadTriTable[j][0]];
				v1 = Fs_[i].vs[quadTriTable[j][1]];
				v2 = Fs_[i].vs[quadTriTable[j][2]];
				if (Vs_[v1].boundary && Vs_[v1].neighbor_fs.size() == 1)
					continue;

				vec0 = (V_[v2] - V_[v1]).normalized();
				vec1 = (V_[v0] - V_[v1]).normalized();
				double currSign = vec0[0] * vec1[1] - vec0[1] * vec1[0];
				/*if (std::abs(currSign) < minJacobi_ && JudgeSign(startSign, currSign))
					minJacobi_ = std::abs(currSign);*/
				if (!JudgeSign(startSign, currSign))
					return false;
				else if (std::abs(currSign) < minJacobi_)
					minJacobi_ = std::abs(currSign);
			}
		}
		return true;
	}


	void QuadMesh::ComputeQuadMeshBBox(Eigen::Vector2d &minP, Eigen::Vector2d &maxP)
	{
		std::vector<Eigen::Vector2d> &qVs = V_;

		std::vector<double> vec0, vec1;
		vec0.reserve(qVs.size());
		vec1.reserve(qVs.size());
		for (int i = 0; i < qVs.size(); ++i)
		{
			vec0.emplace_back(qVs[i][0]);
			vec1.emplace_back(qVs[i][1]);
		}
		std::sort(vec0.begin(), vec0.end());
		std::sort(vec1.begin(), vec1.end());
		minP[0] = vec0[0]; minP[1] = vec1[0];
		maxP[0] = vec0[vec0.size() - 1]; maxP[1] = vec1[vec1.size() - 1];
	}

	void QuadMesh::ScaleQuadMesh(double targetBBox)
	{
		double averangeLength = ComputeAverageEdgeLength();
		Eigen::Vector2d minP, maxP;
		ComputeQuadMeshBBox(minP, maxP);

		double radio = targetBBox / averangeLength;

		for (int i = 0; i < V_.size(); ++i)
		{
			V_[i][0] = minP[0] + (V_[i][0] - minP[0])*radio;
			V_[i][1] = minP[1] + (V_[i][1] - minP[1])*radio;
		}
	}

	void QuadMesh::ScaleQuadMeshBack(double radio, Eigen::Vector2d &minP)
	{
		for (int i = 0; i < V_.size(); ++i)
		{
			V_[i][0] = (V_[i][0] - minP[0]) / radio + minP[0];
			V_[i][1] = (V_[i][1] - minP[1]) / radio + minP[1];
		}
	}
	
	void QuadMesh::ScaleQuadMesh(double targetLength, double &radio, Eigen::Vector2d &minP)
	{
		double averangeLength = ComputeAverageEdgeLength();
		Eigen::Vector2d maxP;
		ComputeQuadMeshBBox(minP, maxP);

		radio = targetLength / averangeLength;

		for (int i = 0; i < V_.size(); ++i)
		{
			V_[i][0] = minP[0] + (V_[i][0] - minP[0])*radio;
			V_[i][1] = minP[1] + (V_[i][1] - minP[1])*radio;
		}
	}

	bool QuadMesh::SplitBadLine(uint32_t badV, uint32_t badF, uint32_t startV, uint32_t startE, std::vector<uint32_t> &outFs)
	{
		const double offsetWeight = 0.15;
		tempEs.clear();
		tempVs.clear();
		tempOefs.clear();
		if (!FindSplitLine(startV, startE, tempEs, tempVs, tempOefs) || tempEs.empty())	//tempOefs的数量和除startV之外的内点的数量一致，换句话说，startV和边界点没有tempOefs.
			return false;

		/*int *vFlag = new int[Vs_.size()];
		std::memset(vFlag, 0, Vs_.size() * sizeof(int));*/
		int *fFlag = new int[Fs_.size()];
		std::memset(fFlag, 0, Fs_.size() * sizeof(int));	//1表示顺着es方向走的f（逆时针情形下），2表示逆着es方向走的f

		//先找起始边的左右两个面，后面的面的关系都可以从这两个面得到
		uint32_t f0 = Es_[startE].neighbor_fs[0], f1 = Es_[startE].neighbor_fs[1];
		std::function<int(uint32_t, uint32_t, uint32_t)> FindF1or2 = [&](uint32_t inF, uint32_t v0, uint32_t v1)->int
		{
			uint32_t v0Id, v1Id;
			std::vector<uint32_t> &fvs = Fs_[inF].vs;
			for (int i = 0; i < 4; ++i)
			{
				if (fvs[i] == v0)
					v0Id = i;
				else if (fvs[i] == v1)
					v1Id = i;
			}

			if (v0Id - v1Id == 1 || (v0Id == 0 && v1Id == 3))
				return 2;
			else if (v1Id - v0Id == 1 || (v0Id == 3 && v1Id == 0))
				return 1;
		};
		int f0State = FindF1or2(f0, tempVs[0], tempVs[1]);
		int f1State = FindF1or2(f1, tempVs[0], tempVs[1]);
		if (f0State == 2 && f1State == 1)
		{
			std::swap(f0State, f1State);
			std::swap(f0, f1);
		}
		fFlag[f0] = f0State;
		fFlag[f1] = f1State;

		//找初始点左右关联的面和点
		std::vector<uint32_t> currVs, currFs;
		std::function<bool(uint32_t, uint32_t, uint32_t&, uint32_t&)> findInitialNextVF = [&](uint32_t preE, uint32_t preF, uint32_t& nextE, uint32_t& nextF) -> bool
		{
			if (preF == badF)
				return false;

			std::vector<uint32_t> &fes = Fs_[preF].es;
			for (int i = 0; i < 4; ++i)
			{
				if (fes[i] != preE && (Es_[fes[i]].vs[0] == startV || Es_[fes[i]].vs[1] == startV))
				{
					nextE = fes[i];
					break;
				}
			}

			if (Es_[nextE].neighbor_fs[0] == preF)
				nextF = Es_[nextE].neighbor_fs[1];
			else
				nextF = Es_[nextE].neighbor_fs[0];

			return true;
		};

		uint32_t currVId = Vs_.size();
		typedef struct VInfoPair
		{
			uint32_t vId0, vId1;
			Eigen::Vector2d vPos0, vPos1;
		}VInfoPair;
		std::vector<VInfoPair> vVec;
		std::vector<uint32_t> changeVFs; //需要改变顶点的面
		uint32_t *vOldToNew = new uint32_t[Vs_.size()];
		std::memset(vOldToNew, (uint32_t)-1, sizeof(uint32_t)*Vs_.size());

		currVs.clear(); currFs.clear();
		currVs.emplace_back(tempVs[1]);
		uint32_t preE = startE, preF = f0, nextE = (uint32_t)-1, nextF = (uint32_t)-1;
		while (findInitialNextVF(preE, preF, nextE, nextF))
		{
			currFs.emplace_back(preF);
			fFlag[preF] = 1;
			if (Es_[nextE].vs[0] == startV)
				currVs.emplace_back(Es_[nextE].vs[1]);
			else
				currVs.emplace_back(Es_[nextE].vs[0]);
			preE = nextE;
			preF = nextF;
		}
		VInfoPair currVIP;
		currVIP.vId0 = startV;
		currVIP.vPos0 = V_[startV];
		for (int i = 0; i < currVs.size(); ++i)
		{
			currVIP.vPos0 += offsetWeight * (V_[currVs[i]] - V_[startV]);
		}
		currVs.clear(); currFs.clear();
		currVs.emplace_back(tempVs[1]);
		preE = startE, preF = f1, nextE = (uint32_t)-1, nextF = (uint32_t)-1;
		while (findInitialNextVF(preE, preF, nextE, nextF))
		{
			currFs.emplace_back(preF);
			changeVFs.emplace_back(preF);
			fFlag[preF] = 2;
			if (Es_[nextE].vs[0] == startV)
				currVs.emplace_back(Es_[nextE].vs[1]);
			else
				currVs.emplace_back(Es_[nextE].vs[0]);
			preE = nextE;
			preF = nextF;
		}
		vOldToNew[startV] = currVId;
		currVIP.vId1 = currVId++;
		currVIP.vPos1 = V_[startV];
		for (int i = 0; i < currVs.size(); ++i)
		{
			currVIP.vPos1 += offsetWeight * (V_[currVs[i]] - V_[startV]);
		}
		vVec.emplace_back(currVIP);

		for (int i = 1; i < tempVs.size() - 1; ++i)
		{
			uint32_t currEId = i - 1;
			uint32_t currNe = tempEs[i];
			std::vector<uint32_t> &orderEs = tempOefs[currEId].oEs;
			std::vector<uint32_t> &orderFs = tempOefs[currEId].oFs;
			auto it = std::find(orderEs.begin(), orderEs.end(), currNe);
			uint32_t oneSideFId = std::distance(orderEs.begin(), it);
			if (fFlag[orderFs[0]] == 2)
			{
				std::reverse(orderEs.begin(), orderEs.end());
				std::reverse(orderFs.begin(), orderFs.end());
				oneSideFId = orderFs.size() - 2 - oneSideFId;
			}

			VInfoPair vIP;
			vIP.vId0 = tempVs[i];
			vIP.vPos0 = V_[tempVs[i]] + offsetWeight * (V_[tempVs[i - 1]] - V_[tempVs[i]]);
			for (int j = 0; j <= oneSideFId; ++j)
			{
				uint32_t anotherV;
				if (Es_[orderEs[j]].vs[0] == tempVs[i])
					anotherV = Es_[orderEs[j]].vs[1];
				else
					anotherV = Es_[orderEs[j]].vs[0];
				vIP.vPos0 += offsetWeight * (V_[anotherV] - V_[tempVs[i]]);
			}
			vOldToNew[tempVs[i]] = currVId;
			vIP.vId1 = currVId++;
			vIP.vPos1 = V_[tempVs[i]] + offsetWeight * (V_[tempVs[i - 1]] - V_[tempVs[i]]);
			for (int j = oneSideFId; j < orderEs.size(); ++j)
			{
				uint32_t anotherV;
				if (Es_[orderEs[j]].vs[0] == tempVs[i])
					anotherV = Es_[orderEs[j]].vs[1];
				else
					anotherV = Es_[orderEs[j]].vs[0];
				vIP.vPos1 += offsetWeight * (V_[anotherV] - V_[tempVs[i]]);
			}
			vVec.emplace_back(vIP);
			for (int j = 0; j <= oneSideFId; ++j)
			{
				fFlag[orderFs[j]] = 1;
			}
			for (int j = oneSideFId + 1; j < orderFs.size(); ++j)
			{
				fFlag[orderFs[j]] = 2;
				changeVFs.emplace_back(orderFs[j]);
			}
		}

		//边界上的最后一个点。
		uint32_t finalV = tempVs[tempVs.size() - 1], finalE = tempEs[tempEs.size()-1];
		uint32_t finalE0 = (uint32_t)-1, finalE1 = (uint32_t)-1;
		std::function<bool(uint32_t, uint32_t, uint32_t&, uint32_t&)> findNVF = [&](uint32_t preE, uint32_t preF, uint32_t &nextE, uint32_t &nextF)->bool
		{
			if (Es_[preE].boundary)
				return false;

			std::vector<uint32_t> &fes = Fs_[preF].es;
			for (int i = 0; i < 4; ++i)
			{
				if (fes[i] != preE && (Es_[fes[i]].vs[0] == finalV || Es_[fes[i]].vs[1] == finalV))
				{
					nextE = fes[i];
					break;
				}
			}

			if (Es_[nextE].boundary)
				nextF = preF;
			else
			{
				if (Es_[nextE].neighbor_fs[0] == preF)
					nextF = Es_[nextE].neighbor_fs[1];
				else
					nextF = Es_[nextE].neighbor_fs[0];
			}
			
			return true;
		};

		vOldToNew[finalV] = badV;
		f0 = Es_[finalE].neighbor_fs[0];
		f1 = Es_[finalE].neighbor_fs[1];
		if (fFlag[f0] == 2)
			std::swap(f0, f1);
		changeVFs.emplace_back(f1);
		preE = finalE; preF = f0;
		while (findNVF(preE, preF, nextE, nextF))
		{
			fFlag[nextF] = 1;
			preE = nextE;
			preF = nextF;
		}
		preE = finalE; preF = f1;
		while (findNVF(preE, preF, nextE, nextF))
		{
			fFlag[nextF] = 2;
			changeVFs.emplace_back(nextF);
			preE = nextE;
			preF = nextF;
		}
		currVIP.vId0 = finalV;
		currVIP.vPos0 = V_[finalV];
		uint32_t anotherV;
		if (Es_[preE].vs[0] == finalV)
			anotherV = Es_[preE].vs[1];
		else
			anotherV = Es_[preE].vs[0];
		currVIP.vId1 = badV;
		currVIP.vPos1 = V_[finalV] + offsetWeight * (V_[anotherV] - V_[finalV]);
		vVec.emplace_back(currVIP);

		//构造新的qm
		std::vector<Eigen::Vector2d> newV = V_;
		std::vector<QuadVertex> newVs;
		std::vector<QuadFace> newFs;
		bool *fUsedFlag = new bool[Fs_.size()];
		std::memset(fUsedFlag, 0, Fs_.size() * sizeof(bool));
		
		for (int i = 0; i < Vs_.size(); ++i)
		{
			QuadVertex qv;
			qv.id = i;
			qv.bId = (uint32_t)-1;
			newVs.emplace_back(qv);
		}
		for (int i = 0; i < vVec.size()-1; ++i)
		{
			newV.emplace_back(Eigen::Vector2d());
			newV[vVec[i].vId0] = vVec[i].vPos0;
			newV[vVec[i].vId1] = vVec[i].vPos1;
			QuadVertex qv;
			qv.id = vVec[i].vId1;
			qv.bId = (uint32_t)-1;
			newVs.emplace_back(qv);
		}
		newV[vVec[vVec.size()-1].vId0] = vVec[vVec.size() - 1].vPos0;
		newV[vVec[vVec.size() - 1].vId1] = vVec[vVec.size() - 1].vPos1;

		std::sort(changeVFs.begin(), changeVFs.end());
		changeVFs.erase(std::unique(changeVFs.begin(), changeVFs.end()), changeVFs.end());
		for (int i = 0; i < changeVFs.size(); ++i)
		{
			fUsedFlag[changeVFs[i]] = true;
		}
		for (int i = 0; i < Fs_.size(); ++i)
		{
			QuadFace qf;
			qf.id = i;
			qf.bId = (uint32_t)-1;
			qf.vs.reserve(4);

			if (!fUsedFlag[i])
			{
				if (i != badF)
				{
					qf.vs = Fs_[i].vs;
				}
				else
				{
					int badVId = FindVInFId(badV, i);
					std::vector<uint32_t> &fvs = Fs_[i].vs;
					qf.vs.emplace_back(fvs[(badVId + 1) % 4]);
					qf.vs.emplace_back(vVec[0].vId1);
					qf.vs.emplace_back(fvs[(badVId + 2) % 4]);
					qf.vs.emplace_back(fvs[(badVId + 3) % 4]);
				}
			}
			else
			{
				std::vector<uint32_t> &fvs = Fs_[i].vs;
				for (int j = 0; j < 4; ++j)
				{
					if (vOldToNew[fvs[j]] != (uint32_t)-1)
						qf.vs.emplace_back(vOldToNew[fvs[j]]);
					else
						qf.vs.emplace_back(fvs[j]);
				}
			}

			for (int j = 0; j < 4; ++j)
			{
				newVs[qf.vs[j]].neighbor_fs.emplace_back(i);
			}
			newFs.emplace_back(qf);
		}
		//Add new faces to new quad mesh
		uint32_t currFId = newFs.size();
		for (int i = 0; i < vVec.size() - 1; ++i)
		{
			QuadFace qf;
			qf.id = currFId++;
			qf.bId = (uint32_t)-1;
			qf.vs.reserve(4);
			qf.vs.emplace_back(vVec[i].vId0);
			qf.vs.emplace_back(vVec[i].vId1);
			qf.vs.emplace_back(vVec[i + 1].vId1);
			qf.vs.emplace_back(vVec[i + 1].vId0);
			for (int j = 0; j < 4; ++j)
			{
				newVs[qf.vs[j]].neighbor_fs.emplace_back(currFId-1);
			}
			outFs.emplace_back(currFId - 1);
			newFs.emplace_back(qf);
		}
		V_ = newV;
		Vs_ = newVs;
		Fs_ = newFs;
		Es_.clear();
		BuildConnectivity();

		delete fUsedFlag;
		fUsedFlag = NULL;
		delete[] vOldToNew;
		vOldToNew = NULL;
		/*delete[] vFlag;
		vFlag = NULL;*/
		delete[] fFlag;
		fFlag = NULL;

		return true;
	}

	int QuadMesh::FindVInFId(uint32_t vId, uint32_t fId)
	{
		std::vector<uint32_t> &fvs = Fs_[fId].vs;
		for (int i = 0; i < 4; ++i)
		{
			if (fvs[i] == vId)
				return i;
		}
		return -1;
	}

	bool QuadMesh::FindSplitLine(uint32_t startV, uint32_t startE, std::vector<uint32_t> &es, std::vector<uint32_t> &vs, std::vector<OrderEFs> &oefs)
	{
		std::stack<SingularNode>().swap(snStack_);
		bool *vFlag = new bool[Vs_.size()];	//标记点有没有被找过
		std::memset(vFlag, 0, Vs_.size() * sizeof(bool));
		std::stack<SingularNode>().swap(snStack_);

		std::vector<std::vector<uint32_t>> subEs, subVs;
		std::vector<std::vector<OrderEFs>> subOefs;
		
		uint32_t subStartE = startE, subStartV = (uint32_t)-1;
		if (Es_[subStartE].vs[0] == startV)
			subStartV = Es_[subStartE].vs[1];
		else
			subStartV = Es_[subStartE].vs[0];

		vFlag[startV] = true;
		vs.emplace_back(startV);

		bool isBound = true;
		while (true)
		{
			std::vector<uint32_t> currEs, currVs;
			std::vector<OrderEFs> currOefs;
			bool isNotLoop = FindSubLine(subStartE, subStartV, currEs, currVs, vFlag, isBound, currOefs);
			if (isNotLoop)
			{
				if (isBound)
				{
					subEs.emplace_back(currEs);
					subVs.emplace_back(currVs);
					subOefs.emplace_back(currOefs);
					break;
				}
				else
				{
					snStack_.emplace(SingularNode());
					subEs.emplace_back(currEs);
					subVs.emplace_back(currVs);
					subOefs.emplace_back(currOefs);
					std::vector<uint32_t> tempFs;
					FindOrderEdges(currEs[currEs.size() - 1], currVs[currVs.size() - 1], snStack_.top().orderEs, tempFs);
					ReorderEs(snStack_.top().orderEs);
					snStack_.top().currE = currEs[currEs.size() - 1];
					snStack_.top().currEId = 0;
					snStack_.top().currV = currVs[currVs.size() - 1];
					snStack_.top().vValence = snStack_.top().orderEs.size() + 1;

					subStartE = snStack_.top().orderEs[0];
					if (Es_[subStartE].vs[0] == snStack_.top().currV)
						subStartV = Es_[subStartE].vs[1];
					else
						subStartV = Es_[subStartE].vs[0];
				}
			}
			else
			{
				if (subEs.empty() || !FindNewSubStart(snStack_, subStartE, subStartV))
				{
					delete[] vFlag;
					vFlag = NULL;
					return false;
				}

				while (subEs.size() > snStack_.size())
				{
					subEs.pop_back();
					subVs.pop_back();
					subOefs.pop_back();
				}
			}
		}

		for (int i = 0; i < subEs.size(); ++i)
		{
			for (int j = 0; j < subEs[i].size(); ++j)
			{
				es.emplace_back(subEs[i][j]);
				vs.emplace_back(subVs[i][j]);
			}
			for (int j = 0; j < subOefs[i].size(); ++j)
			{
				oefs.emplace_back(subOefs[i][j]);
			}
		}

		delete[] vFlag;
		vFlag = NULL;
		return true;
	}

	bool QuadMesh::FindNewSubStart(std::stack<SingularNode> &snStack, uint32_t &newSE, uint32_t &newSV)
	{
		while (true)
		{
			if (snStack.empty())
				return false;

			SingularNode &currSN = snStack.top();
			if (currSN.currEId < currSN.orderEs.size() - 1)
			{
				++currSN.currEId;
				newSE = currSN.orderEs[currSN.currEId];
				if (Es_[newSE].vs[0] == currSN.currV)
					newSV = Es_[newSE].vs[1];
				else
					newSV = Es_[newSE].vs[0];
				return true;
			}
			else
			{
				snStack.pop();
			}
		}
		return true;
	}

	void QuadMesh::ReorderEs(std::vector<uint32_t> &es)
	{
		if (es.empty())
			return;

		std::vector<uint32_t> newEs;
		newEs.reserve(es.size());
		int halfSize = es.size() / 2;
		if (es.size() % 2 == 0)
		{
			for (int i = 0; i < halfSize; ++i)
			{
				newEs.emplace_back(es[halfSize - 1 - i]);
				newEs.emplace_back(es[halfSize + i]);
			}
		}
		else
		{
			newEs.emplace_back(es[halfSize]);
			for (int i = 0; i < halfSize; ++i)
			{
				newEs.emplace_back(es[halfSize - 1 - i]);
				newEs.emplace_back(es[halfSize + 1 + i]);
			}
		}
		es.swap(newEs);
	}

	bool QuadMesh::FindSubLine(uint32_t startE, uint32_t startV, std::vector<uint32_t> &es, std::vector<uint32_t> &vs, bool *vFlag, bool &isBoundary, std::vector<OrderEFs> &oefs)
	{
		/*es.emplace_back(startE);
		vs.emplace_back(startV);*/

		//if (vFlag[startV] || qm.Vs_[startV].boundary)	//表示之前已经找过这个点了 或者是边界点
		//	return false;

		uint32_t preE = startE, preV = startV, nextE = (uint32_t)-1, nextV = (uint32_t)-1;
		std::vector<uint32_t> orderEs, orderFs;
		while (true)
		{
			if (vFlag[preV])
			{
				for (int i = 0; i < vs.size(); ++i)
				{
					vFlag[vs[i]] = false;
				}
				return false;
			}
			vFlag[preV] = true;

			if (Vs_[preV].boundary)
			{
				es.emplace_back(preE);
				vs.emplace_back(preV);
				isBoundary = true;
				return true;
			}
			FindOrderEdges(preE, preV, orderEs, orderFs);
			/*if (Vs_[preV].boundary)
			{
				es.emplace_back(preE);
				vs.emplace_back(preV);
				isBoundary = true;
				return true;
			}
			else*/ if (Vs_[preV].neighbor_es.size() != 4)
			{
				oefs.emplace_back(OrderEFs());
				oefs[oefs.size() - 1].oEs = orderEs;
				oefs[oefs.size() - 1].oFs = orderFs;
				es.emplace_back(preE);
				vs.emplace_back(preV);
				isBoundary = false;
				return true;
			}
			else
			{
				oefs.emplace_back(OrderEFs());
				oefs[oefs.size() - 1].oEs = orderEs;
				oefs[oefs.size() - 1].oFs = orderFs;
				es.emplace_back(preE);
				vs.emplace_back(preV);
			}

			nextE = orderEs[1];
			if (Es_[nextE].vs[0] == preV)
				nextV = Es_[nextE].vs[1];
			else
				nextV = Es_[nextE].vs[0];

			preV = nextV; preE = nextE;
		}
	}

	void QuadMesh::FindOrderEdges(uint32_t startE, uint32_t currV, std::vector<uint32_t> &orderEs, std::vector<uint32_t> &orderFs)
	{
		if (Vs_[currV].boundary)
		{
			std::cout << "QuadMesh::FindOrderEdges! " << std::endl;
			exit(-1323);
		}

		orderEs.clear();
		orderFs.clear();
		//orderEs.emplace_back(startE);
		tempEF = Es_[startE].neighbor_fs[0];
		orderFs.emplace_back(tempEF);

		std::function<bool(uint32_t, uint32_t, uint32_t&, uint32_t&)> findNextEF = [&](uint32_t preE, uint32_t preF, uint32_t &nextE, uint32_t &nextF)-> bool
		{
			std::vector<uint32_t> &fes = Fs_[preF].es;
			for (int i = 0; i < 4; ++i)
			{
				if (fes[i] != preE && (Es_[fes[i]].vs[0] == currV || Es_[fes[i]].vs[1] == currV))
				{
					nextE = fes[i];
					break;
				}
			}

			if (nextE == startE)
				return false;

			if (Es_[nextE].neighbor_fs[0] == preF)
				nextF = Es_[nextE].neighbor_fs[1];
			else
				nextF = Es_[nextE].neighbor_fs[0];

			return true;
		};
		
		uint32_t preE = startE, preF = tempEF, nextE = (uint32_t)-1, nextF = (uint32_t)-1;
		while (findNextEF(preE, preF, nextE, nextF))
		{
			orderEs.emplace_back(nextE);
			orderFs.emplace_back(nextF);
			preE = nextE;
			preF = nextF;
		}
	}

	void QuadMesh::FindOrderEdges(uint32_t startE, uint32_t startF, uint32_t currV, std::vector<uint32_t> &orderEs, std::vector<uint32_t> &orderFs)
	{
		if (Vs_[currV].boundary)
		{
			std::cout << "QuadMesh::FindOrderEdges! " << std::endl;
			exit(-1323);
		}

		orderEs.clear();
		orderFs.clear();
		//orderEs.emplace_back(startE);
		tempEF = startF;
		orderFs.emplace_back(tempEF);

		std::function<bool(uint32_t, uint32_t, uint32_t&, uint32_t&)> findNextEF = [&](uint32_t preE, uint32_t preF, uint32_t &nextE, uint32_t &nextF)-> bool
		{
			std::vector<uint32_t> &fes = Fs_[preF].es;
			for (int i = 0; i < 4; ++i)
			{
				if (fes[i] != preE && (Es_[fes[i]].vs[0] == currV || Es_[fes[i]].vs[1] == currV))
				{
					nextE = fes[i];
					break;
				}
			}

			if (nextE == startE)
				return false;

			if (Es_[nextE].neighbor_fs[0] == preF)
				nextF = Es_[nextE].neighbor_fs[1];
			else
				nextF = Es_[nextE].neighbor_fs[0];

			return true;
		};

		uint32_t preE = startE, preF = tempEF, nextE = (uint32_t)-1, nextF = (uint32_t)-1;
		while (findNextEF(preE, preF, nextE, nextF))
		{
			orderEs.emplace_back(nextE);
			orderFs.emplace_back(nextF);
			preE = nextE;
			preF = nextF;
		}
	}
}

