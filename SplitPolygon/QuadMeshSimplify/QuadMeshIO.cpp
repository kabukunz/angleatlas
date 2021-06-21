#include "QuadMeshIO.h"
#include "QuadMesh.h"
#include <iostream>
#include <fstream>
#include <sstream>

namespace BaseDataStructure
{

	QuadMeshIO::QuadMeshIO()
	{
	}


	QuadMeshIO::~QuadMeshIO()
	{
	}

	void QuadMeshIO::ReadQuadMesh(QuadMesh *quadMesh, std::string fileName)
	{
		if (fileName.empty() || quadMesh==NULL)
		{
			std::cout << "QuadMeshIO::ReadQuadMesh: Empty file! " << std::endl;
			return;
		}
		quadMesh->ClearMesh();

		std::ifstream ifs(fileName.c_str());
		if (!ifs.is_open())
			return;

		std::istringstream iss;
		std::string line, word;
		double x, y, z;
		std::string fInfo[4];
		std::vector<std::string> strSplit;

		QuadVertex qv;
		uint32_t vIndex = 0, fIndex = 0;
		while (getline(ifs, line))
		{
			iss.clear();
			iss.str(line);
			iss >> word;
			if (word == "v")
			{
				iss >> x >> y >> z;
				quadMesh->V_.emplace_back(Eigen::Vector2d(x, y));
				qv.id = vIndex++;
				qv.bId = (uint32_t)-1;
				quadMesh->Vs_.emplace_back(qv);
			}
			else if (word == "f")
			{
				iss >> fInfo[0] >> fInfo[1] >> fInfo[2] >> fInfo[3];
				if (fInfo[0] == "" || fInfo[1] == "" || fInfo[2] == "" || fInfo[3] == "")
				{
					std::cout << "This is not a valid quad obj. " << std::endl;
					exit(0);
				}
				QuadFace qf;
				qf.id = fIndex++;
				qf.bId = (uint32_t)-1;
				for (int j = 0; j < 4; ++j)
				{
					strSplit.clear();
					SplitString(fInfo[j], strSplit, "/");
					int currIndex = std::stoi(strSplit[0]) - 1;
					qf.vs.emplace_back(currIndex);
					quadMesh->Vs_[currIndex].neighbor_fs.emplace_back(fIndex - 1);
				}
				quadMesh->Fs_.emplace_back(qf);
			}
		}
	}

	void QuadMeshIO::WriteQuadMesh(QuadMesh *quadMesh, std::string fileName)
	{
		std::ofstream ofs(fileName.c_str());
		for (int i = 0; i < quadMesh->V_.size(); ++i)
		{
			ofs << "v " << quadMesh->V_[i][0] << " " << quadMesh->V_[i][1] << " 0" << std::endl;
		}

		for (int i = 0; i < quadMesh->Fs_.size(); ++i)
		{
			ofs << "f " << quadMesh->Fs_[i].vs[0]+1 << " " << quadMesh->Fs_[i].vs[1]+1 << " " << quadMesh->Fs_[i].vs[2]+1 << " " << quadMesh->Fs_[i].vs[3]+1 << std::endl;
		}
	}

	void QuadMeshIO::SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
	{
		std::string::size_type pos1, pos2;
		pos2 = s.find(c);
		pos1 = 0;
		while (std::string::npos != pos2)
		{
			v.push_back(s.substr(pos1, pos2 - pos1));

			pos1 = pos2 + c.size();
			pos2 = s.find(c, pos1);
		}
		if (pos1 != s.length())
			v.push_back(s.substr(pos1));
	}
}