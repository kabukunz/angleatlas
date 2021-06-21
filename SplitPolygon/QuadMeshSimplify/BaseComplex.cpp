#include "BaseComplex.h"
#include "QuadMesh.h"
#include "QuadMeshIO.h"
#include "Pipeline.h"
#include <iostream>
#include <iterator>
#include <fstream>
#include <queue>
//using DataOperation::FeatureConstraints;

namespace BaseDataStructure
{

	BaseComplex::BaseComplex()
	{
	}


	BaseComplex::~BaseComplex()
	{
	}

	void BaseComplex::ExtractBCVerticesEdges(QuadMesh *quadMesh, double angleThre)
	{
		ClearBaseComplex();

		if (quadMesh == NULL)
		{
			std::cout << "BaseComplex::ExtractBaseComplex: quadMesh is empty. " << std::endl;
			return;
		}

		auto &Vs = quadMesh->Vs_;
		auto &Es = quadMesh->Es_;
		auto &Fs = quadMesh->Fs_;

		//Vertex flag.
		enum VertexState
		{
			NOT_FOUND,
			SINGULAR,
			NOT_SINGULAR_BUT_BV,
			ON_SINGULAR_EDGE
		};
		std::vector<VertexState> vertexFlag;
		vertexFlag.resize(Vs.size());
		std::fill(vertexFlag.begin(), vertexFlag.end(), VertexState::NOT_FOUND);

		//Find Base-Complex Vertices.
		Bv_.reserve(Vs.size());
		uint32_t fvNum = 0;
		for (int i = 0; i < Vs.size(); ++i)
		{
			if ((Vs[i].boundary && Vs[i].neighbor_fs.size() != 2)
				|| (!Vs[i].boundary&&Vs[i].neighbor_fs.size() != 4))
			{
				FrameVertex fv;
				fv.id = fvNum++;
				fv.qId = i;
				fv.boundary = Vs[i].boundary;
				fv.singular = true;
				Bv_.emplace_back(fv);

				Vs[i].bId = fvNum - 1;

				vertexFlag[i] = VertexState::SINGULAR;
			}

			//对于角点也给它搞成base complex的顶点试试
			if (Vs[i].boundary && vertexFlag[i]==VertexState::NOT_FOUND)
			{
				std::vector<uint32_t> &ves = Vs[i].neighbor_es;
				tempVs_[0] = (uint32_t)-1; tempVs_[1] = (uint32_t)-1;
				int currCount = 0;
				for (int j = 0; j < ves.size(); ++j)
				{
					if (Es[ves[j]].boundary)
					{
						if (Es[ves[j]].vs[0] == i)
							tempVs_[currCount++] = Es[ves[j]].vs[1];
						else
							tempVs_[currCount++] = Es[ves[j]].vs[0];
					}
					if (currCount == 2)
						break;
				}
				if (tempVs_[0] == (uint32_t)-1 || tempVs_[1] == (uint32_t)-1)
					continue;

				const Eigen::Vector2d &vec1 = quadMesh->V_[tempVs_[0]] - quadMesh->V_[i].normalized();
				const Eigen::Vector2d &vec2 = quadMesh->V_[tempVs_[1]] - quadMesh->V_[i].normalized();
				double cosValue = vec1.dot(vec2);
				if (cosValue > std::cos(angleThre))
				{
					FrameVertex fv;
					fv.id = fvNum++;
					fv.qId = i;
					fv.boundary = Vs[i].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);

					Vs[i].bId = fvNum - 1;

					vertexFlag[i] = VertexState::NOT_SINGULAR_BUT_BV;
				}
			}
		}

		std::vector<bool> vFlag(quadMesh->Vs_.size(), false); //标识一个vertex有没有被找过。
		std::function<bool(uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE)> findNextVE = [&](uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE) -> bool
		{
			std::vector<uint32_t> &currVs = quadMesh->Es_[preE].vs;
			if (currVs[0] == preV)
				nextV = currVs[1];
			else
				nextV = currVs[0];

			if (vFlag[nextV] == true)
				return false;

			std::vector<uint32_t> &nextEs = quadMesh->Vs_[nextV].neighbor_es;
			for (int i = 0; i < nextEs.size(); ++i)
			{
				if (quadMesh->Es_[nextEs[i]].boundary && nextEs[i] != preE)
				{
					nextE = nextEs[i];
					break;
				}
			}
			return true;
		};
		//找一条边界上边界点的函数
		std::function<bool(uint32_t startV, std::vector<uint32_t> &boundaryVs, std::vector<uint32_t> &bvsInThisCurve)> findSVInOneBoundary = [&](uint32_t startV, std::vector<uint32_t> &boundaryVs, std::vector<uint32_t> &bvsInThisCurve) -> bool
		{
			vFlag[startV] = true;
			std::vector<uint32_t> &ves = quadMesh->Vs_[startV].neighbor_es;
			uint32_t startE = (uint32_t)-1;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (quadMesh->Es_[ves[i]].boundary)
				{
					startE = ves[i];
					break;
				}
			}

			//std::vector<uint32_t> bvsInThisCurve;
			boundaryVs.emplace_back(startV);
			if (quadMesh->Vs_[startV].bId != (uint32_t)-1)
			{
				bvsInThisCurve.emplace_back(boundaryVs.size() - 1);
			}
			uint32_t preV = startV, preE = startE, nextV = (uint32_t)-1, nextE = (uint32_t)-1;
			while (findNextVE(preV, preE, nextV, nextE))
			{
				vFlag[nextV] = true;
				boundaryVs.emplace_back(nextV);
				if (quadMesh->Vs_[nextV].bId != (uint32_t)-1)
					bvsInThisCurve.emplace_back(boundaryVs.size() - 1);
				preV = nextV;
				preE = nextE;
			}

			if (bvsInThisCurve.empty() || bvsInThisCurve.size() == 1)
				return false;
			else
				return true;
		};

		//如果有些边界边上没有奇异点或者少，那么要人为加奇异点上去
		for (int i = 0; i < quadMesh->Vs_.size(); ++i)
		{
			if (!quadMesh->Vs_[i].boundary || vFlag[i])
				continue;

			std::vector<uint32_t> boundaryVs, bvsInThisCurve;
			if (!findSVInOneBoundary(i, boundaryVs, bvsInThisCurve))
			{
				if (bvsInThisCurve.empty())
				{
					uint32_t v0 = boundaryVs[0], v1 = boundaryVs[boundaryVs.size() / 2];
					FrameVertex fv;
					fv.id = fvNum++;
					fv.qId = v0;
					fv.boundary = Vs[v0].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);
					Vs[v0].bId = fvNum - 1;
					vertexFlag[v0] = VertexState::NOT_SINGULAR_BUT_BV;

					fv.id = fvNum++;
					fv.qId = v1;
					fv.boundary = Vs[v1].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);
					Vs[v1].bId = fvNum - 1;
					vertexFlag[v1] = VertexState::NOT_SINGULAR_BUT_BV;
				}
				else if (bvsInThisCurve.size() == 1)
				{
					uint32_t v0 = (bvsInThisCurve[0] + boundaryVs.size() / 2) % boundaryVs.size();
					FrameVertex fv;
					fv.id = fvNum++;
					fv.qId = v0;
					fv.boundary = Vs[v0].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);
					Vs[v0].bId = fvNum - 1;
					vertexFlag[v0] = VertexState::NOT_SINGULAR_BUT_BV;
				}
			}
		}

		std::function<bool(uint32_t, uint32_t, uint32_t &)> findNextEdge = [&](uint32_t preE, uint32_t preV, uint32_t &nextE)->bool
		{
			if (vertexFlag[preV] != VertexState::NOT_FOUND)
				return false;

			std::vector<uint32_t> &preEFs = Es[preE].neighbor_fs;
			std::vector<uint32_t> interFs;
			const std::vector<uint32_t> &currEs = Vs[preV].neighbor_es;
			for (int i = 0; i < currEs.size(); ++i)
			{
				if (currEs[i] == preE)
					continue;

				std::vector<uint32_t> &currEFs = Es[currEs[i]].neighbor_fs;
				interFs.clear();
				std::sort(preEFs.begin(), preEFs.end());
				std::sort(currEFs.begin(), currEFs.end());
				std::set_intersection(preEFs.begin(), preEFs.end(), currEFs.begin(), currEFs.end(), std::inserter(interFs, interFs.begin()));
				if (interFs.empty())
				{
					nextE = currEs[i];
					return true;
				}
			}

			return false;
		};

		//Find Base-Complex Edges.
		std::vector<bool> edgeFlag;	//edge flag to flag if the edge is found before.
		edgeFlag.resize(Es.size());
		std::fill(edgeFlag.begin(), edgeFlag.end(), false);
		uint32_t feNum = 0;
		for (int i = 0; i < Bv_.size(); ++i)
		{
			const QuadVertex &qv = Vs[Bv_[i].qId];
			const std::vector<uint32_t> &nes = qv.neighbor_es;
			for (int j = 0; j < nes.size(); ++j)
			{
				if (edgeFlag[nes[j]])
					continue;

				FrameEdge fe;
				Bv_[i].b_neighbor_es.emplace_back(feNum);
				fe.id = feNum++;
				Es[nes[j]].bId = feNum - 1;
				fe.boundary = Es[nes[j]].boundary;
				fe.b_vs.emplace_back(i);
				fe.vs_link.emplace_back(qv.id);
				//fe.es_link.emplace_back(nes[j]);
				edgeFlag[nes[j]] = true;

				uint32_t preE = nes[j], preV = 0, nextE = (uint32_t)-1;
				const std::vector<uint32_t> &eVs = Es[preE].vs;
				if (eVs[0] == qv.id)
					preV = eVs[1];
				else
					preV = eVs[0];
				while (findNextEdge(preE, preV, nextE))
				{
					edgeFlag[nextE] = true;
					vertexFlag[preV] = VertexState::ON_SINGULAR_EDGE;
					fe.vs_link.emplace_back(preV);
					fe.es_link.emplace_back(preE);
					Es[nextE].bId = feNum - 1;

					preE = nextE;
					auto &eVs = Es[nextE].vs;
					if (eVs[0] == preV)
						preV = eVs[1];
					else
						preV = eVs[0];
				}

				if (vertexFlag[preV] == VertexState::SINGULAR || vertexFlag[preV] == VertexState::NOT_SINGULAR_BUT_BV)
				{
					fe.vs_link.emplace_back(preV);
					fe.b_vs.emplace_back(Vs[preV].bId);
					fe.es_link.emplace_back(preE);
					Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(fe.id);
					Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(i);
					Bv_[i].b_neighbor_vs.emplace_back(Vs[preV].bId);

					Be_.emplace_back(fe);
				}
				else if (vertexFlag[preV] == VertexState::ON_SINGULAR_EDGE || (Vs[preV].boundary && !Es[preE].boundary))
				{
					/*if (vertexFlag[preV] != VertexState::ON_SINGULAR_EDGE)
					{
						double dde = 3432;
					}*/
					FrameVertex fv;
					fv.id = fvNum++;
					Vs[preV].bId = fvNum - 1;
					fv.qId = preV;
					fv.boundary = Vs[preV].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);

					fe.vs_link.emplace_back(preV);
					fe.b_vs.emplace_back(Vs[preV].bId);
					fe.es_link.emplace_back(preE);
					Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(fe.id);
					Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(i);
					Bv_[i].b_neighbor_vs.emplace_back(Vs[preV].bId);

					//split the source bc edge into two.
					if (vertexFlag[preV] == VertexState::ON_SINGULAR_EDGE)
					{
						auto &vEs = Vs[preV].neighbor_es;
						int seNum = 0;
						std::vector<uint32_t> ses;
						ses.resize(2);
						for (int k = 0; k < vEs.size(); ++k)
						{
							if (vEs[k] != preE && edgeFlag[vEs[k]] == true)
							{
								ses[seNum] = vEs[k];
								++seNum;
							}
						}
						if (seNum != 2)
						{
							std::cout << "BaseComplex::ExtractBaseComplex: Error in split! " << std::endl;
							isError_ = true;
							return;
							//exit(-1);
						}

						uint32_t bEIndex = Es[ses[0]].bId;
						auto esLink = Be_[bEIndex].es_link;
						auto vsLink = Be_[bEIndex].vs_link;
						Be_[bEIndex].es_link.clear();
						Be_[bEIndex].vs_link.clear();
						auto itE = std::find_first_of(esLink.begin(), esLink.end(), ses.begin(), ses.end());
						auto itV = std::find(vsLink.begin(), vsLink.end(), preV);
						Be_[bEIndex].es_link.resize(itE - esLink.begin() + 1);
						std::copy(esLink.begin(), itE + 1, Be_[bEIndex].es_link.begin());
						Be_[bEIndex].vs_link.resize(itV - vsLink.begin() + 1);
						std::copy(vsLink.begin(), itV + 1, Be_[bEIndex].vs_link.begin());

						uint32_t preFVIndex0 = Be_[bEIndex].b_vs[0], preFVIndex1 = Be_[bEIndex].b_vs[1];
						auto &ves = Bv_[preFVIndex1].b_neighbor_es;
						auto itEE = std::find(ves.begin(), ves.end(), bEIndex);
						*itEE = feNum;
						auto &vvs0 = Bv_[preFVIndex0].b_neighbor_vs;
						auto &vvs1 = Bv_[preFVIndex1].b_neighbor_vs;
						auto it = std::find(vvs0.begin(), vvs0.end(), preFVIndex1);
						*it = fvNum - 1;
						it = std::find(vvs1.begin(), vvs1.end(), preFVIndex0);
						*it = fvNum - 1;

						FrameEdge feNew;
						feNew.id = feNum++;
						feNew.b_vs.emplace_back(fvNum - 1);
						feNew.b_vs.emplace_back(Be_[bEIndex].b_vs[1]);
						feNew.boundary = Be_[bEIndex].boundary;
						feNew.es_link.resize(esLink.end() - itE - 1);
						std::copy(itE + 1, esLink.end(), feNew.es_link.begin());
						for (int k = 0; k < feNew.es_link.size(); ++k)
						{
							Es[feNew.es_link[k]].bId = feNew.id;
						}
						feNew.vs_link.resize(vsLink.end() - itV);
						std::copy(itV, vsLink.end(), feNew.vs_link.begin());

						Be_[bEIndex].b_vs[1] = fvNum - 1;
						Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(preFVIndex0);
						Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(preFVIndex1);
						Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(bEIndex);
						Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(feNew.id);
						Be_.emplace_back(fe);
						Be_.emplace_back(feNew);
					}
					else
					{
						Be_.emplace_back(fe);
					}

					vertexFlag[preV] = VertexState::NOT_SINGULAR_BUT_BV;
				}

			}
		}

		//build e--b-neighbor-es
		for (int i = 0; i < Bv_.size(); ++i)
		{
			auto &vEs = Bv_[i].b_neighbor_es;
			for (int j = 0; j < vEs.size(); ++j)
			{
				for (int k = 0; k < vEs.size(); ++k)
				{
					if (k == j)
						continue;
					Be_[vEs[j]].b_neighbor_es.emplace_back(vEs[k]);
				}
			}
		}
	}

	void BaseComplex::ExtractBCVerticesEdges(QuadMesh *quadMesh)
	{
		ClearBaseComplex();

		if (quadMesh == NULL)
		{
			std::cout << "BaseComplex::ExtractBaseComplex: quadMesh is empty. " << std::endl;
			return;
		}

		auto &Vs = quadMesh->Vs_;
		auto &Es = quadMesh->Es_;
		auto &Fs = quadMesh->Fs_;

		//Vertex flag.
		enum VertexState
		{
			NOT_FOUND,
			SINGULAR,
			NOT_SINGULAR_BUT_BV,
			ON_SINGULAR_EDGE
		};
		std::vector<VertexState> vertexFlag;
		vertexFlag.resize(Vs.size());
		std::fill(vertexFlag.begin(), vertexFlag.end(), VertexState::NOT_FOUND);

		//Find Base-Complex Vertices.
		Bv_.reserve(Vs.size());
		uint32_t fvNum = 0;
		for (int i = 0; i < Vs.size(); ++i)
		{
			if ((Vs[i].boundary && Vs[i].neighbor_fs.size() != 2)
				|| (!Vs[i].boundary&&Vs[i].neighbor_fs.size() != 4))
			{
				FrameVertex fv;
				fv.id = fvNum++;
				fv.qId = i;
				fv.boundary = Vs[i].boundary;
				fv.singular = true;
				Bv_.emplace_back(fv);

				Vs[i].bId = fvNum - 1;

				vertexFlag[i] = VertexState::SINGULAR;
			}
		}

		std::vector<bool> vFlag(quadMesh->Vs_.size(), false); //标识一个vertex有没有被找过。
		std::function<bool(uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE)> findNextVE = [&](uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE) -> bool
		{
			std::vector<uint32_t> &currVs = quadMesh->Es_[preE].vs;
			if (currVs[0] == preV)
				nextV = currVs[1];
			else
				nextV = currVs[0];

			if (vFlag[nextV] == true)
				return false;

			std::vector<uint32_t> &nextEs = quadMesh->Vs_[nextV].neighbor_es;
			for (int i = 0; i < nextEs.size(); ++i)
			{
				if (quadMesh->Es_[nextEs[i]].boundary && nextEs[i] != preE)
				{
					nextE = nextEs[i];
					break;
				}
			}
			return true;
		};
		//找一条边界上边界点的函数
		std::function<bool(uint32_t startV, std::vector<uint32_t> &boundaryVs, std::vector<uint32_t> &bvsInThisCurve)> findSVInOneBoundary = [&](uint32_t startV, std::vector<uint32_t> &boundaryVs, std::vector<uint32_t> &bvsInThisCurve) -> bool
		{
			vFlag[startV] = true;
			std::vector<uint32_t> &ves = quadMesh->Vs_[startV].neighbor_es;
			uint32_t startE = (uint32_t)-1;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (quadMesh->Es_[ves[i]].boundary)
				{
					startE = ves[i];
					break;
				}
			}

			//std::vector<uint32_t> bvsInThisCurve;
			boundaryVs.emplace_back(startV);
			if (quadMesh->Vs_[startV].bId != (uint32_t)-1)
			{
				bvsInThisCurve.emplace_back(boundaryVs.size()-1);
			}
			uint32_t preV = startV, preE = startE, nextV = (uint32_t)-1, nextE = (uint32_t)-1;
			while (findNextVE(preV, preE, nextV, nextE))
			{
				vFlag[nextV] = true;
				boundaryVs.emplace_back(nextV);
				if (quadMesh->Vs_[nextV].bId != (uint32_t)-1)
					bvsInThisCurve.emplace_back(boundaryVs.size()-1);
				preV = nextV;
				preE = nextE;
			}

			if (bvsInThisCurve.empty() || bvsInThisCurve.size()==1)
				return false;
			else
				return true;
		};

		//如果有些边界边上没有奇异点或者少，那么要人为加奇异点上去
		for (int i = 0; i < quadMesh->Vs_.size(); ++i)
		{
			if (!quadMesh->Vs_[i].boundary || vFlag[i])
				continue;

			std::vector<uint32_t> boundaryVs, bvsInThisCurve;
			if (!findSVInOneBoundary(i, boundaryVs, bvsInThisCurve))
			{
				if (bvsInThisCurve.empty())
				{
					uint32_t v0 = boundaryVs[0], v1 = boundaryVs[boundaryVs.size() / 2];
					FrameVertex fv;
					fv.id = fvNum++;
					fv.qId = v0;
					fv.boundary = Vs[v0].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);
					Vs[v0].bId = fvNum - 1;
					vertexFlag[v0] = VertexState::NOT_SINGULAR_BUT_BV;

					fv.id = fvNum++;
					fv.qId = v1;
					fv.boundary = Vs[v1].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);
					Vs[v1].bId = fvNum - 1;
					vertexFlag[v1] = VertexState::NOT_SINGULAR_BUT_BV;
				}
				else if (bvsInThisCurve.size() == 1)
				{
					uint32_t v0 = boundaryVs[(bvsInThisCurve[0] + boundaryVs.size()/2)%boundaryVs.size()];
					FrameVertex fv;
					fv.id = fvNum++;
					fv.qId = v0;
					fv.boundary = Vs[v0].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);
					Vs[v0].bId = fvNum - 1;
					vertexFlag[v0] = VertexState::NOT_SINGULAR_BUT_BV;
				}
			}  
		}

		std::function<bool(uint32_t, uint32_t, uint32_t &)> findNextEdge = [&](uint32_t preE, uint32_t preV, uint32_t &nextE)->bool
		{
			if (vertexFlag[preV] != VertexState::NOT_FOUND)
				return false;

			std::vector<uint32_t> &preEFs = Es[preE].neighbor_fs;
			std::vector<uint32_t> interFs;
			const std::vector<uint32_t> &currEs = Vs[preV].neighbor_es;
			for (int i = 0; i < currEs.size(); ++i)
			{
				if (currEs[i] == preE)
					continue;

				std::vector<uint32_t> &currEFs = Es[currEs[i]].neighbor_fs;
				interFs.clear();
				std::sort(preEFs.begin(), preEFs.end());
				std::sort(currEFs.begin(), currEFs.end());
				std::set_intersection(preEFs.begin(), preEFs.end(), currEFs.begin(), currEFs.end(), std::inserter(interFs, interFs.begin()));
				if (interFs.empty())
				{
					nextE = currEs[i];
					return true;
				}
			}
			
			return false;
		};

		//Find Base-Complex Edges.
		//std::vector<bool> edgeFlag ;	//edge flag to flag if the edge is found before.
		//edgeFlag.resize(Es.size());
		//std::fill(edgeFlag.begin(), edgeFlag.end(), false);
		bool *edgeFlag = new bool[Es.size()];

		uint32_t feNum = 0;
		std::vector<uint32_t> currVs;
		bool haveToRestart = false;

		std::function<void()> returnBegin = [&]()
		{
			std::memset(edgeFlag, 0, Es.size() * sizeof(bool));
			feNum = 0;
			haveToRestart = false;
			for (int i = 0; i < Bv_.size(); ++i)
			{
				Bv_[i].b_neighbor_es.clear();
				Bv_[i].b_neighbor_fs.clear();
				Bv_[i].b_neighbor_vs.clear();
			}
			for (int i = 0; i < Es.size(); ++i)
			{
				Es[i].bId = (uint32_t)-1;
			}
			for (int i = 0; i < Vs.size(); ++i)
			{
				if (vertexFlag[i] == VertexState::ON_SINGULAR_EDGE)
				{
					vertexFlag[i] = VertexState::NOT_FOUND;
				}
			}
			Be_.clear();
		};

		do
		{
			returnBegin();

			for (int i = 0; i < Bv_.size(); ++i)
			{
				const QuadVertex &qv = Vs[Bv_[i].qId];
				const std::vector<uint32_t> &nes = qv.neighbor_es;
				for (int j = 0; j < nes.size(); ++j)
				{
					if (edgeFlag[nes[j]])
						continue;

					FrameEdge fe;
					Bv_[i].b_neighbor_es.emplace_back(feNum);
					fe.id = feNum++;
					Es[nes[j]].bId = feNum - 1;
					fe.boundary = Es[nes[j]].boundary;
					fe.b_vs.emplace_back(i);
					fe.vs_link.emplace_back(qv.id);
					//fe.es_link.emplace_back(nes[j]);
					edgeFlag[nes[j]] = true;

					uint32_t preE = nes[j], preV = 0, nextE = (uint32_t)-1;
					const std::vector<uint32_t> &eVs = Es[preE].vs;
					if (eVs[0] == qv.id)
						preV = eVs[1];
					else
						preV = eVs[0];
					uint32_t startV = preV;
					currVs.clear();
					while (findNextEdge(preE, preV, nextE))
					{
						edgeFlag[nextE] = true;
						vertexFlag[preV] = VertexState::ON_SINGULAR_EDGE;
						currVs.emplace_back(preV);
						fe.vs_link.emplace_back(preV);
						fe.es_link.emplace_back(preE);
						Es[nextE].bId = feNum - 1;

						preE = nextE;
						auto &eVs = Es[nextE].vs;
						if (eVs[0] == preV)
							preV = eVs[1];
						else
							preV = eVs[0];
					}
					std::vector<uint32_t>::iterator it;
					if (!currVs.empty())
					{
						it = std::find(currVs.begin(), currVs.end(), preV);
					}
					if (!currVs.empty() && it!=currVs.end())
					{
						int dis = std::distance(it, currVs.end());
						if (currVs.size() < 3 || dis<2)
						{
							isError_ = true;
							delete[] edgeFlag;
							edgeFlag = NULL;
							return;
						}
						uint32_t newBVId = *(it+dis/2);

						haveToRestart = true;
						FrameVertex fv;
						fv.id = fvNum++;
						fv.qId = newBVId;
						fv.boundary = Vs[newBVId].boundary;
						fv.singular = false;
						Bv_.emplace_back(fv);
						Vs[i].bId = fvNum - 1;
						vertexFlag[i] = VertexState::NOT_SINGULAR_BUT_BV;
						break;
					}

					if (vertexFlag[preV] == VertexState::SINGULAR || vertexFlag[preV] == VertexState::NOT_SINGULAR_BUT_BV)
					{
						fe.vs_link.emplace_back(preV);
						fe.b_vs.emplace_back(Vs[preV].bId);
						fe.es_link.emplace_back(preE);
						Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(fe.id);
						Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(i);
						Bv_[i].b_neighbor_vs.emplace_back(Vs[preV].bId);

						Be_.emplace_back(fe);
					}
					else if (vertexFlag[preV] == VertexState::ON_SINGULAR_EDGE || (Vs[preV].boundary && !Es[preE].boundary))
					{
						/*if (vertexFlag[preV] != VertexState::ON_SINGULAR_EDGE)
						{
							double dde = 3432;
						}*/
						FrameVertex fv;
						fv.id = fvNum++;
						Vs[preV].bId = fvNum - 1;
						fv.qId = preV;
						fv.boundary = Vs[preV].boundary;
						fv.singular = false;
						Bv_.emplace_back(fv);

						fe.vs_link.emplace_back(preV);
						fe.b_vs.emplace_back(Vs[preV].bId);
						fe.es_link.emplace_back(preE);
						Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(fe.id);
						Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(i);
						Bv_[i].b_neighbor_vs.emplace_back(Vs[preV].bId);

						//split the source bc edge into two.
						if (vertexFlag[preV] == VertexState::ON_SINGULAR_EDGE)
						{
							auto &vEs = Vs[preV].neighbor_es;
							int seNum = 0;
							std::vector<uint32_t> ses;
							ses.resize(2);
							for (int k = 0; k < vEs.size(); ++k)
							{
								if (vEs[k] != preE && edgeFlag[vEs[k]] == true)
								{
									ses[seNum] = vEs[k];
									++seNum;
								}
							}
							if (seNum != 2)
							{
								std::cout << "BaseComplex::ExtractBaseComplex: Error in split! " << std::endl;
								isError_ = true;
								delete[] edgeFlag;
								edgeFlag = NULL;
								return;
								//exit(-1);
							}

							uint32_t bEIndex = Es[ses[0]].bId;
							auto esLink = Be_[bEIndex].es_link;
							auto vsLink = Be_[bEIndex].vs_link;
							Be_[bEIndex].es_link.clear();
							Be_[bEIndex].vs_link.clear();
							auto itE = std::find_first_of(esLink.begin(), esLink.end(), ses.begin(), ses.end());
							auto itV = std::find(vsLink.begin(), vsLink.end(), preV);
							Be_[bEIndex].es_link.resize(itE - esLink.begin() + 1);
							std::copy(esLink.begin(), itE + 1, Be_[bEIndex].es_link.begin());
							Be_[bEIndex].vs_link.resize(itV - vsLink.begin() + 1);
							std::copy(vsLink.begin(), itV + 1, Be_[bEIndex].vs_link.begin());

							uint32_t preFVIndex0 = Be_[bEIndex].b_vs[0], preFVIndex1 = Be_[bEIndex].b_vs[1];
							auto &ves = Bv_[preFVIndex1].b_neighbor_es;
							auto itEE = std::find(ves.begin(), ves.end(), bEIndex);
							*itEE = feNum;
							auto &vvs0 = Bv_[preFVIndex0].b_neighbor_vs;
							auto &vvs1 = Bv_[preFVIndex1].b_neighbor_vs;
							auto it = std::find(vvs0.begin(), vvs0.end(), preFVIndex1);
							*it = fvNum - 1;
							it = std::find(vvs1.begin(), vvs1.end(), preFVIndex0);
							*it = fvNum - 1;

							FrameEdge feNew;
							feNew.id = feNum++;
							feNew.b_vs.emplace_back(fvNum - 1);
							feNew.b_vs.emplace_back(Be_[bEIndex].b_vs[1]);
							feNew.boundary = Be_[bEIndex].boundary;
							feNew.es_link.resize(esLink.end() - itE - 1);
							std::copy(itE + 1, esLink.end(), feNew.es_link.begin());
							for (int k = 0; k < feNew.es_link.size(); ++k)
							{
								Es[feNew.es_link[k]].bId = feNew.id;
							}
							feNew.vs_link.resize(vsLink.end() - itV);
							std::copy(itV, vsLink.end(), feNew.vs_link.begin());

							Be_[bEIndex].b_vs[1] = fvNum - 1;
							Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(preFVIndex0);
							Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(preFVIndex1);
							Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(bEIndex);
							Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(feNew.id);
							Be_.emplace_back(fe);
							Be_.emplace_back(feNew);
						}
						else
						{
							Be_.emplace_back(fe);
						}

						vertexFlag[preV] = VertexState::NOT_SINGULAR_BUT_BV;
					}

				}
			}
		}
		while (haveToRestart);

		//build e--b-neighbor-es
		for (int i = 0; i < Bv_.size(); ++i)
		{
			auto &vEs = Bv_[i].b_neighbor_es;
			for (int j = 0; j < vEs.size(); ++j)
			{
				for (int k = 0; k < vEs.size(); ++k)
				{
					if (k == j)
						continue;
					Be_[vEs[j]].b_neighbor_es.emplace_back(vEs[k]);
				}
			}
		}

		delete[] edgeFlag;
		edgeFlag = NULL;

#pragma region Test_BC_1
		/*for (int i = 0; i < Be_.size(); ++i)
		{
			auto eEs = Be_[i].b_neighbor_es;
			std::sort(eEs.begin(), eEs.end());
			eEs.erase(std::unique(eEs.begin(), eEs.end()), eEs.end());
			if (eEs.size() != Be_[i].b_neighbor_es.size())
				std::cout << "Oh, no! " << std::endl;
		}

		for (int i = 0; i < Be_.size(); ++i)
		{
			auto &eEs = Be_[i].b_neighbor_es;
			auto vs0 = Be_[i].b_vs;
			for (int j = 0; j < eEs.size(); ++j)
			{
				auto vs1 = Be_[eEs[j]].b_vs;
				std::vector<uint32_t> result;
				std::sort(vs0.begin(), vs0.end());
				std::sort(vs1.begin(), vs1.end());
				std::set_intersection(vs0.begin(), vs0.end(), vs1.begin(), vs1.end(), std::inserter(result, result.begin()));
				if (result.empty())
					std::cout << "OOOh, no!" << std::endl;
			}
		}*/
#pragma endregion
	}

	void BaseComplex::ExtractBCVerticesEdges(QuadMesh *quadMesh, std::vector<uint32_t> &corners)
	{
		ClearBaseComplex();

		if (quadMesh == NULL)
		{
			std::cout << "BaseComplex::ExtractBaseComplex: quadMesh is empty. " << std::endl;
			return;
		}

		auto &Vs = quadMesh->Vs_;
		auto &Es = quadMesh->Es_;
		auto &Fs = quadMesh->Fs_;

		//Vertex flag.
		enum VertexState
		{
			NOT_FOUND,
			SINGULAR,
			NOT_SINGULAR_BUT_BV,
			ON_SINGULAR_EDGE
		};
		std::vector<VertexState> vertexFlag;
		vertexFlag.resize(Vs.size());
		std::fill(vertexFlag.begin(), vertexFlag.end(), VertexState::NOT_FOUND);

		//Find Base-Complex Vertices.
		Bv_.reserve(Vs.size());
		uint32_t fvNum = 0;
		for (int i = 0; i < Vs.size(); ++i)
		{
			if ((Vs[i].boundary && Vs[i].neighbor_fs.size() != 2)
				|| (!Vs[i].boundary&&Vs[i].neighbor_fs.size() != 4))
			{
				FrameVertex fv;
				fv.id = fvNum++;
				fv.qId = i;
				fv.boundary = Vs[i].boundary;
				fv.singular = true;
				Bv_.emplace_back(fv);

				Vs[i].bId = fvNum - 1;

				vertexFlag[i] = VertexState::SINGULAR;
			}
		}

		for (int i = 0; i < corners.size(); ++i)
		{
			uint32_t currId = corners[i];
			if (Vs[currId].bId == (uint32_t)-1)
			{
				FrameVertex fv;
				fv.id = fvNum++;
				fv.qId = currId;
				fv.boundary = Vs[currId].boundary;
				fv.singular = false;
				Bv_.emplace_back(fv);

				Vs[currId].bId = fvNum - 1;
				vertexFlag[currId] = VertexState::NOT_SINGULAR_BUT_BV;
			}
		}

		std::vector<bool> vFlag(quadMesh->Vs_.size(), false); //标识一个vertex有没有被找过。
		std::function<bool(uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE)> findNextVE = [&](uint32_t preV, uint32_t preE, uint32_t &nextV, uint32_t &nextE) -> bool
		{
			std::vector<uint32_t> &currVs = quadMesh->Es_[preE].vs;
			if (currVs[0] == preV)
				nextV = currVs[1];
			else
				nextV = currVs[0];

			if (vFlag[nextV] == true)
				return false;

			std::vector<uint32_t> &nextEs = quadMesh->Vs_[nextV].neighbor_es;
			for (int i = 0; i < nextEs.size(); ++i)
			{
				if (quadMesh->Es_[nextEs[i]].boundary && nextEs[i] != preE)
				{
					nextE = nextEs[i];
					break;
				}
			}
			return true;
		};
		//找一条边界上边界点的函数
		std::function<bool(uint32_t startV, std::vector<uint32_t> &boundaryVs, std::vector<uint32_t> &bvsInThisCurve)> findSVInOneBoundary = [&](uint32_t startV, std::vector<uint32_t> &boundaryVs, std::vector<uint32_t> &bvsInThisCurve) -> bool
		{
			vFlag[startV] = true;
			std::vector<uint32_t> &ves = quadMesh->Vs_[startV].neighbor_es;
			uint32_t startE = (uint32_t)-1;
			for (int i = 0; i < ves.size(); ++i)
			{
				if (quadMesh->Es_[ves[i]].boundary)
				{
					startE = ves[i];
					break;
				}
			}

			//std::vector<uint32_t> bvsInThisCurve;
			boundaryVs.emplace_back(startV);
			if (quadMesh->Vs_[startV].bId != (uint32_t)-1)
			{
				bvsInThisCurve.emplace_back(boundaryVs.size() - 1);
			}
			uint32_t preV = startV, preE = startE, nextV = (uint32_t)-1, nextE = (uint32_t)-1;
			while (findNextVE(preV, preE, nextV, nextE))
			{
				vFlag[nextV] = true;
				boundaryVs.emplace_back(nextV);
				if (quadMesh->Vs_[nextV].bId != (uint32_t)-1)
					bvsInThisCurve.emplace_back(boundaryVs.size() - 1);
				preV = nextV;
				preE = nextE;
			}

			if (bvsInThisCurve.empty() || bvsInThisCurve.size() == 1)
				return false;
			else
				return true;
		};

		//如果有些边界边上没有奇异点或者少，那么要人为加奇异点上去
		for (int i = 0; i < quadMesh->Vs_.size(); ++i)
		{
			if (!quadMesh->Vs_[i].boundary || vFlag[i])
				continue;

			std::vector<uint32_t> boundaryVs, bvsInThisCurve;
			if (!findSVInOneBoundary(i, boundaryVs, bvsInThisCurve))
			{
				if (bvsInThisCurve.empty())
				{
					uint32_t v0 = boundaryVs[0], v1 = boundaryVs[boundaryVs.size() / 2];
					FrameVertex fv;
					fv.id = fvNum++;
					fv.qId = v0;
					fv.boundary = Vs[v0].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);
					Vs[v0].bId = fvNum - 1;
					vertexFlag[v0] = VertexState::NOT_SINGULAR_BUT_BV;

					fv.id = fvNum++;
					fv.qId = v1;
					fv.boundary = Vs[v1].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);
					Vs[v1].bId = fvNum - 1;
					vertexFlag[v1] = VertexState::NOT_SINGULAR_BUT_BV;
				}
				else if (bvsInThisCurve.size() == 1)
				{
					uint32_t v0 = boundaryVs[(bvsInThisCurve[0] + boundaryVs.size() / 2) % boundaryVs.size()];
					FrameVertex fv;
					fv.id = fvNum++;
					fv.qId = v0;
					fv.boundary = Vs[v0].boundary;
					fv.singular = false;
					Bv_.emplace_back(fv);
					Vs[v0].bId = fvNum - 1;
					vertexFlag[v0] = VertexState::NOT_SINGULAR_BUT_BV;
				}
			}
		}

		std::function<bool(uint32_t, uint32_t, uint32_t &)> findNextEdge = [&](uint32_t preE, uint32_t preV, uint32_t &nextE)->bool
		{
			if (vertexFlag[preV] != VertexState::NOT_FOUND)
				return false;

			std::vector<uint32_t> &preEFs = Es[preE].neighbor_fs;
			std::vector<uint32_t> interFs;
			const std::vector<uint32_t> &currEs = Vs[preV].neighbor_es;
			for (int i = 0; i < currEs.size(); ++i)
			{
				if (currEs[i] == preE)
					continue;

				std::vector<uint32_t> &currEFs = Es[currEs[i]].neighbor_fs;
				interFs.clear();
				std::sort(preEFs.begin(), preEFs.end());
				std::sort(currEFs.begin(), currEFs.end());
				std::set_intersection(preEFs.begin(), preEFs.end(), currEFs.begin(), currEFs.end(), std::inserter(interFs, interFs.begin()));
				if (interFs.empty())
				{
					nextE = currEs[i];
					return true;
				}
			}

			return false;
		};

		//Find Base-Complex Edges.
		//std::vector<bool> edgeFlag ;	//edge flag to flag if the edge is found before.
		//edgeFlag.resize(Es.size());
		//std::fill(edgeFlag.begin(), edgeFlag.end(), false);
		bool *edgeFlag = new bool[Es.size()];

		uint32_t feNum = 0;
		std::vector<uint32_t> currVs;
		bool haveToRestart = false;

		std::function<void()> returnBegin = [&]()
		{
			std::memset(edgeFlag, 0, Es.size() * sizeof(bool));
			feNum = 0;
			haveToRestart = false;
			for (int i = 0; i < Bv_.size(); ++i)
			{
				Bv_[i].b_neighbor_es.clear();
				Bv_[i].b_neighbor_fs.clear();
				Bv_[i].b_neighbor_vs.clear();
			}
			for (int i = 0; i < Es.size(); ++i)
			{
				Es[i].bId = (uint32_t)-1;
			}
			for (int i = 0; i < Vs.size(); ++i)
			{
				if (vertexFlag[i] == VertexState::ON_SINGULAR_EDGE)
				{
					vertexFlag[i] = VertexState::NOT_FOUND;
				}
			}
			Be_.clear();
		};

		do
		{
			returnBegin();

			for (int i = 0; i < Bv_.size(); ++i)
			{
				const QuadVertex &qv = Vs[Bv_[i].qId];
				const std::vector<uint32_t> &nes = qv.neighbor_es;
				for (int j = 0; j < nes.size(); ++j)
				{
					if (edgeFlag[nes[j]])
						continue;

					FrameEdge fe;
					Bv_[i].b_neighbor_es.emplace_back(feNum);
					fe.id = feNum++;
					Es[nes[j]].bId = feNum - 1;
					fe.boundary = Es[nes[j]].boundary;
					fe.b_vs.emplace_back(i);
					fe.vs_link.emplace_back(qv.id);
					//fe.es_link.emplace_back(nes[j]);
					edgeFlag[nes[j]] = true;

					uint32_t preE = nes[j], preV = 0, nextE = (uint32_t)-1;
					const std::vector<uint32_t> &eVs = Es[preE].vs;
					if (eVs[0] == qv.id)
						preV = eVs[1];
					else
						preV = eVs[0];
					uint32_t startV = preV;
					currVs.clear();
					while (findNextEdge(preE, preV, nextE))
					{
						edgeFlag[nextE] = true;
						vertexFlag[preV] = VertexState::ON_SINGULAR_EDGE;
						currVs.emplace_back(preV);
						fe.vs_link.emplace_back(preV);
						fe.es_link.emplace_back(preE);
						Es[nextE].bId = feNum - 1;

						preE = nextE;
						auto &eVs = Es[nextE].vs;
						if (eVs[0] == preV)
							preV = eVs[1];
						else
							preV = eVs[0];
					}
					std::vector<uint32_t>::iterator it;
					if (!currVs.empty())
					{
						it = std::find(currVs.begin(), currVs.end(), preV);
					}
					if (!currVs.empty() && it != currVs.end())
					{
						int dis = std::distance(it, currVs.end());
						if (currVs.size() < 3 || dis < 2)
						{
							isError_ = true;
							delete[] edgeFlag;
							edgeFlag = NULL;
							return;
						}
						uint32_t newBVId = *(it + dis / 2);

						haveToRestart = true;
						FrameVertex fv;
						fv.id = fvNum++;
						fv.qId = newBVId;
						fv.boundary = Vs[newBVId].boundary;
						fv.singular = false;
						Bv_.emplace_back(fv);
						Vs[i].bId = fvNum - 1;
						vertexFlag[i] = VertexState::NOT_SINGULAR_BUT_BV;
						break;
					}

					if (vertexFlag[preV] == VertexState::SINGULAR || vertexFlag[preV] == VertexState::NOT_SINGULAR_BUT_BV)
					{
						fe.vs_link.emplace_back(preV);
						fe.b_vs.emplace_back(Vs[preV].bId);
						fe.es_link.emplace_back(preE);
						Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(fe.id);
						Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(i);
						Bv_[i].b_neighbor_vs.emplace_back(Vs[preV].bId);

						Be_.emplace_back(fe);
					}
					else if (vertexFlag[preV] == VertexState::ON_SINGULAR_EDGE || (Vs[preV].boundary && !Es[preE].boundary))
					{
						/*if (vertexFlag[preV] != VertexState::ON_SINGULAR_EDGE)
						{
							double dde = 3432;
						}*/
						FrameVertex fv;
						fv.id = fvNum++;
						Vs[preV].bId = fvNum - 1;
						fv.qId = preV;
						fv.boundary = Vs[preV].boundary;
						fv.singular = false;
						Bv_.emplace_back(fv);

						fe.vs_link.emplace_back(preV);
						fe.b_vs.emplace_back(Vs[preV].bId);
						fe.es_link.emplace_back(preE);
						Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(fe.id);
						Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(i);
						Bv_[i].b_neighbor_vs.emplace_back(Vs[preV].bId);

						//split the source bc edge into two.
						if (vertexFlag[preV] == VertexState::ON_SINGULAR_EDGE)
						{
							auto &vEs = Vs[preV].neighbor_es;
							int seNum = 0;
							std::vector<uint32_t> ses;
							ses.resize(2);
							for (int k = 0; k < vEs.size(); ++k)
							{
								if (vEs[k] != preE && edgeFlag[vEs[k]] == true)
								{
									ses[seNum] = vEs[k];
									++seNum;
								}
							}
							if (seNum != 2)
							{
								std::cout << "BaseComplex::ExtractBaseComplex: Error in split! " << std::endl;
								isError_ = true;
								delete[] edgeFlag;
								edgeFlag = NULL;
								return;
								//exit(-1);
							}

							uint32_t bEIndex = Es[ses[0]].bId;
							auto esLink = Be_[bEIndex].es_link;
							auto vsLink = Be_[bEIndex].vs_link;
							Be_[bEIndex].es_link.clear();
							Be_[bEIndex].vs_link.clear();
							auto itE = std::find_first_of(esLink.begin(), esLink.end(), ses.begin(), ses.end());
							auto itV = std::find(vsLink.begin(), vsLink.end(), preV);
							Be_[bEIndex].es_link.resize(itE - esLink.begin() + 1);
							std::copy(esLink.begin(), itE + 1, Be_[bEIndex].es_link.begin());
							Be_[bEIndex].vs_link.resize(itV - vsLink.begin() + 1);
							std::copy(vsLink.begin(), itV + 1, Be_[bEIndex].vs_link.begin());

							uint32_t preFVIndex0 = Be_[bEIndex].b_vs[0], preFVIndex1 = Be_[bEIndex].b_vs[1];
							auto &ves = Bv_[preFVIndex1].b_neighbor_es;
							auto itEE = std::find(ves.begin(), ves.end(), bEIndex);
							*itEE = feNum;
							auto &vvs0 = Bv_[preFVIndex0].b_neighbor_vs;
							auto &vvs1 = Bv_[preFVIndex1].b_neighbor_vs;
							auto it = std::find(vvs0.begin(), vvs0.end(), preFVIndex1);
							*it = fvNum - 1;
							it = std::find(vvs1.begin(), vvs1.end(), preFVIndex0);
							*it = fvNum - 1;

							FrameEdge feNew;
							feNew.id = feNum++;
							feNew.b_vs.emplace_back(fvNum - 1);
							feNew.b_vs.emplace_back(Be_[bEIndex].b_vs[1]);
							feNew.boundary = Be_[bEIndex].boundary;
							feNew.es_link.resize(esLink.end() - itE - 1);
							std::copy(itE + 1, esLink.end(), feNew.es_link.begin());
							for (int k = 0; k < feNew.es_link.size(); ++k)
							{
								Es[feNew.es_link[k]].bId = feNew.id;
							}
							feNew.vs_link.resize(vsLink.end() - itV);
							std::copy(itV, vsLink.end(), feNew.vs_link.begin());

							Be_[bEIndex].b_vs[1] = fvNum - 1;
							Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(preFVIndex0);
							Bv_[Vs[preV].bId].b_neighbor_vs.emplace_back(preFVIndex1);
							Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(bEIndex);
							Bv_[Vs[preV].bId].b_neighbor_es.emplace_back(feNew.id);
							Be_.emplace_back(fe);
							Be_.emplace_back(feNew);
						}
						else
						{
							Be_.emplace_back(fe);
						}

						vertexFlag[preV] = VertexState::NOT_SINGULAR_BUT_BV;
					}

				}
			}
		} while (haveToRestart);

		//build e--b-neighbor-es
		for (int i = 0; i < Bv_.size(); ++i)
		{
			auto &vEs = Bv_[i].b_neighbor_es;
			for (int j = 0; j < vEs.size(); ++j)
			{
				for (int k = 0; k < vEs.size(); ++k)
				{
					if (k == j)
						continue;
					Be_[vEs[j]].b_neighbor_es.emplace_back(vEs[k]);
				}
			}
		}

		delete[] edgeFlag;
		edgeFlag = NULL;
	}

	void BaseComplex::ExtractBaseComplex(QuadMesh *quadMesh, std::vector<uint32_t> &corners)
	{
		isError_ = false;
		ExtractBCVerticesEdges(quadMesh, corners);
		ExtractBCFaces(quadMesh);
	}

	void BaseComplex::ExtractBaseComplex(QuadMesh *quadMesh)
	{
		isError_ = false;
		ExtractBCVerticesEdges(quadMesh);
#pragma region Test_BC_2
		/*for (int i = 0; i < Bv_.size(); ++i)
		{
			std::cout << Bv_[i].qId << std::endl;
			for (int j = 0; j < Bv_[i].b_neighbor_vs.size(); ++j)
			{
				FrameVertex &currNV = Bv_[Bv_[i].b_neighbor_vs[j]];
				auto it = std::find(currNV.b_neighbor_vs.begin(), currNV.b_neighbor_vs.end(), i);
				if (it == currNV.b_neighbor_vs.end())
				{
					std::cout << "DDDSSS" << std::endl;
				}
			}
		}

		int count = 0;
		for (int i = 0; i < quadMesh->Vs_.size(); ++i)
		{
			if (quadMesh->Vs_[i].bId < 10000)
				++count;

		}
		if (count == Bv_.size())
			std::cout << std::endl << "Yes!" << std::endl;

		int eecount = 0;
		for (int i = 0; i < Be_.size(); ++i)
		{
			eecount += Be_[i].es_link.size();
		}
		int escount = 0;
		for (int i = 0; i < quadMesh->Es_.size(); ++i)
		{
			if (quadMesh->Es_[i].bId < 10000)
				++escount;
		}
		if (eecount == escount)
			std::cout << std::endl << "Yes!" << std::endl;


		for (int i = 0; i < Bv_.size(); ++i)
		{
			auto &vEs = Bv_[i].b_neighbor_es;
			for (int j = 0; j < vEs.size(); ++j)
			{
				auto &currE = Be_[vEs[j]];
				if (currE.b_vs[0] != i && currE.b_vs[1] != i)
				{
					std::cout << "NO! " << std::endl;
				}
			}
		}*/
#pragma endregion

#pragma region Output_BC_Edges
		//OutputBCEdges(quadMesh, "C:\\Users\\ChiZhang\\Desktop\\outEsLink.txt");
#pragma endregion
		ExtractBCFaces(quadMesh);
	}

	void BaseComplex::ExtractBaseComplex(QuadMesh *quadMesh, double angleThre)
	{
		isError_ = false;
		ExtractBCVerticesEdges(quadMesh);
		ExtractBCFaces(quadMesh);
	}

	void BaseComplex::ClearBaseComplex()
	{
		//std::vector<FrameVertex>().swap(Bv_);
		Bv_.clear();
		//std::vector<FrameEdge>().swap(Be_);
		Be_.clear();
		//std::vector<FrameFace>().swap(Bf_);
		Bf_.clear();
	}

	void BaseComplex::OutputBCEdges(QuadMesh *quadMesh, std::string ofileName)
	{
		std::ofstream ofs(ofileName.c_str());

		int sumCount = 0;
		for (int i = 0; i < Be_.size(); ++i)
		{
			sumCount += Be_[i].es_link.size();
		}
		ofs << sumCount << std::endl;
		for (int i = 0; i < Be_.size(); ++i)
		{
			auto esLink = Be_[i].es_link;
			for (int j = 0; j < esLink.size(); ++j)
			{
				uint32_t v0 = quadMesh->Es_[esLink[j]].vs[0];
				uint32_t v1 = quadMesh->Es_[esLink[j]].vs[1];
				const Eigen::Vector2d &vp0 = quadMesh->V_[v0];
				const Eigen::Vector2d &vp1 = quadMesh->V_[v1];
				ofs << vp0[0] << " " << vp0[1] << " " << vp1[0] << " " << vp1[1] << std::endl;
			}
		}
	}

	void BaseComplex::ExtractBCFaces(QuadMesh * quadMesh)
	{
		if (Be_.empty() || Bv_.empty())
		{
			return;
		}
		if (!Bf_.empty())
		{
			//std::vector<FrameFace>().swap(Bf_);
			Bf_.clear();
		}
		const std::vector<QuadVertex> &Vs = quadMesh->Vs_;
		const std::vector<QuadEdge> &Es = quadMesh->Es_;
		std::vector<QuadFace> &Fs = quadMesh->Fs_;
		std::vector<bool> fsFlag;
		fsFlag.resize(Fs.size());
		std::fill(fsFlag.begin(), fsFlag.end(), false);
		uint32_t ffCount = 0;
		std::vector<uint32_t> fAdjBE;
		std::vector<uint32_t> nextFaces;

		std::function<void(uint32_t, std::vector<uint32_t>&, std::vector<uint32_t>&)> findFaceAdjBEAndNewF = [&](uint32_t f, std::vector<uint32_t>& beId, std::vector<uint32_t>& nextFs)
		{
			beId.clear();
			nextFs.clear();
			const std::vector<uint32_t> &fes = Fs[f].es;
			for (int i = 0; i < fes.size(); ++i)
			{
				uint32_t currBE = Es[fes[i]].bId;
				if (currBE != (uint32_t)-1)
				{
					beId.emplace_back(currBE);
				}
				else
				{
					const std::vector<uint32_t> &eAdjFs = Es[fes[i]].neighbor_fs;
					uint32_t anotherFace = (eAdjFs[0] == f) ? eAdjFs[1] : eAdjFs[0];
					if (!fsFlag[anotherFace])
						nextFs.emplace_back(anotherFace);
				}
			}
		};
		
		std::vector<bool> facePoolFlag;
		facePoolFlag.resize(Fs.size());
		/*std::vector<bool> vsnetFlag;
		vsnetFlag.resize(Vs.size());
		std::vector<bool> besFlag;
		besFlag.resize(Be_.size());
		std::vector<bool> fvsFlag;
		fvsFlag.resize(Be_.size());*/
		bool *vsnetFlag = new bool[Vs.size()];
		bool *besFlag = new bool[Be_.size()];
		bool *fvsFlag = new bool[Be_.size()];
		for (int i = 0; i < Bv_.size(); ++i)
		{
			uint32_t currQVId = Bv_[i].qId;
			std::vector<uint32_t> currFs = Vs[currQVId].neighbor_fs;
			for (int j = 0; j < currFs.size(); ++j)
			{
				if (fsFlag[currFs[j]])
					continue;
				/*std::fill(vsnetFlag.begin(), vsnetFlag.end(), false);
				std::fill(besFlag.begin(), besFlag.end(), false);
				std::fill(fvsFlag.begin(), fvsFlag.end(), false);*/
				std::memset(vsnetFlag, 0, Vs.size() * sizeof(bool));
				std::memset(besFlag, 0, Be_.size() * sizeof(bool));
				std::memset(fvsFlag, 0, Be_.size() * sizeof(bool));
				std::vector<uint32_t> currBEs;
				FrameFace ff;
				ff.boundary = false;
				ff.id = ffCount++;
				std::queue<uint32_t> facesPool;
				facesPool.push(currFs[j]);
				while (!facesPool.empty())
				{
					uint32_t currFace = facesPool.front();
					facesPool.pop();
					if (Fs[currFace].boundary && ff.boundary==false)
						ff.boundary = true;
					fsFlag[currFace] = true;
					Fs[currFace].bId = ffCount - 1;
					ff.fs_net.emplace_back(currFace);
					for (int ii = 0; ii < 4; ++ii)
					{
						if (!vsnetFlag[Fs[currFace].vs[ii]])
						{
							vsnetFlag[Fs[currFace].vs[ii]] = true; 
							ff.vs_net.emplace_back(Fs[currFace].vs[ii]);
						}
					}

					findFaceAdjBEAndNewF(currFace, fAdjBE, nextFaces);
					if (fAdjBE.size() >= 2)
					{
						for (int ii = 0; ii < fAdjBE.size(); ++ii)
						{
							if (!besFlag[fAdjBE[ii]])
							{
								currBEs.emplace_back(fAdjBE[ii]);
								besFlag[fAdjBE[ii]] = true;
							}
						}
					}
					for (int k = 0; k < nextFaces.size(); ++k)
					{
						if (!facePoolFlag[nextFaces[k]])
						{
							facesPool.push(nextFaces[k]);
							facePoolFlag[nextFaces[k]] = true;
						}
					}
				}
				/*std::sort(ff.vs_net.begin(), ff.vs_net.end());
				ff.vs_net.erase(std::unique(ff.vs_net.begin(), ff.vs_net.end()), ff.vs_net.end());
				std::sort(currBEs.begin(), currBEs.end());
				currBEs.erase(std::unique(currBEs.begin(), currBEs.end()), currBEs.end());*/
				std::vector<uint32_t> fVs;
				fVs.reserve(currBEs.size() * 2);
				for (int k = 0; k < currBEs.size(); ++k)
				{
					for (int ii = 0; ii < 2; ++ii)
					{
						if (!fvsFlag[Be_[currBEs[k]].b_vs[ii]])
						{
							fVs.emplace_back(Be_[currBEs[k]].b_vs[ii]);
							fvsFlag[Be_[currBEs[k]].b_vs[ii]] = true;
						}
					}
				}
				/*std::sort(fVs.begin(), fVs.end());
				fVs.erase(std::unique(fVs.begin(), fVs.end()), fVs.end());*/
				if (currBEs.size() != 4)
				{
					/*QuadMeshIO qmi;
					qmi.WriteQuadMesh(quadMesh, "C:\\Users\\ChiZhang\\Desktop\\errorMesh.obj");
					OutputBCEdges(quadMesh, "C:\\Users\\ChiZhang\\Desktop\\errorMesh.txt");*/
					std::cout << "Error in finding Bf.b_neighbor. " << std::endl;
					//exit(-3);
					isError_ = true;
					return;
				}
				for (int k = 0; k < 4; ++k)
				{
					Be_[currBEs[k]].b_neighbor_fs.emplace_back(ffCount - 1);
					ff.b_es.emplace_back(currBEs[k]);
				}
				for (int k = 0; k < fVs.size(); ++k)
				{
					Bv_[fVs[k]].b_neighbor_fs.emplace_back(ffCount - 1);
					ff.b_vs.emplace_back(fVs[k]);
				}

				Bf_.emplace_back(ff);
			}
		}
		delete[] vsnetFlag;
		vsnetFlag = NULL;
		delete[] besFlag;
		besFlag = NULL;
		delete[] fvsFlag;
		fvsFlag = NULL;

		for (int i = 0; i < Be_.size(); ++i)
		{
			if (Be_[i].boundary)
				continue;
			if (Be_[i].b_neighbor_fs.size() != 2)
			{
				/*QuadMeshIO qmi;
				qmi.WriteQuadMesh(quadMesh, "C:\\Users\\ChiZhang\\Desktop\\errorMesh.obj");
				OutputBCEdges(quadMesh, "C:\\Users\\ChiZhang\\Desktop\\errorMesh.txt");*/
				std::cout << "Error in Here! " << std::endl;
				isError_ = true;
				//exit(0);
			}
			uint32_t efs0 = Be_[i].b_neighbor_fs[0];
			uint32_t efs1 = Be_[i].b_neighbor_fs[1];
			Bf_[efs0].b_neighbor_fs.emplace_back(efs1);
			Bf_[efs1].b_neighbor_fs.emplace_back(efs0);
		}
	}
}