#include "GeoPolyReadWrite.h"
#include <fstream>
#include <sstream>
#include <map>

GeoPolyReadWrite::GeoPolyReadWrite()
{
}


GeoPolyReadWrite::~GeoPolyReadWrite()
{
}

void GeoPolyReadWrite::SplitString(const std::string &inStr, const std::string &key, std::vector<std::string> &outStrs)
{
	outStrs.clear();
	std::string currLeftStr = inStr;
	int currPos;
	currPos = currLeftStr.find(key);
	while (currPos != -1)
	{
		std::string currStr = currLeftStr.substr(0, currPos);
		outStrs.emplace_back(currStr);
		
		currLeftStr = currLeftStr.substr(currPos + key.size(), currLeftStr.size() - (currPos + key.size()));

		currPos = currLeftStr.find(key);
	}

	if (!currLeftStr.empty())
		outStrs.emplace_back(currLeftStr);
}

void GeoPolyReadWrite::ReadGeoToPoly(const std::string &fileName, std::vector<Eigen::Vector2d> &poly)
{
	poly.clear();

	std::ifstream ifs(fileName.c_str());
	std::string line;

	std::vector<Eigen::Vector2d> pointPos;
	std::vector<int> pointId;

	std::vector<Eigen::Vector2i> linePairs;
	std::vector<int> lineId;

	std::vector<int> lineLoop;

	while (std::getline(ifs, line))
	{
		if (!line.compare(0, 4, "Point", 0, 4))
		{
			int pos0 = line.find("(");
			int pos1 = line.find(")");
			int currId = std::stoi(line.substr(pos0 + 1, pos1 - pos0 - 1));
			pointId.emplace_back(currId);

			pos0 = line.find("{");
			pos1 = line.find("}");
			std::string dataStr = line.substr(pos0 + 1, pos1 - pos0 - 1);
			std::vector<std::string> outStr;
			SplitString(dataStr, ",", outStr);
			Eigen::Vector2d currPoss;
			currPoss[0] = std::stod(outStr[0]);
			currPoss[1] = std::stod(outStr[1]);
			pointPos.emplace_back(currPoss);
		}
		else if (!line.compare(0, 8, "Line Loop", 0, 8))
		{
			int pos0 = line.find("{");
			int pos1 = line.find("}");
			std::string dataStr = line.substr(pos0 + 1, pos1 - pos0 - 1);
			std::vector<std::string> outStr;
			SplitString(dataStr, ",", outStr);
			for (int j = 0; j < outStr.size(); ++j)
			{
				int llId = std::stoi(outStr[j]);
				lineLoop.emplace_back(llId);
			}
		}
		else if (!line.compare(0, 3, "Line", 0, 3))
		{
			int pos0 = line.find("(");
			int pos1 = line.find(")");
			int currId = std::stoi(line.substr(pos0 + 1, pos1 - pos0 - 1));
			lineId.emplace_back(currId);

			pos0 = line.find("{");
			pos1 = line.find("}");
			std::string dataStr = line.substr(pos0 + 1, pos1 - pos0 - 1);
			std::vector<std::string> outStr;
			SplitString(dataStr, ",", outStr);
			Eigen::Vector2i currLine;
			currLine[0] = std::stoi(outStr[0]);
			currLine[1] = std::stoi(outStr[1]);
			linePairs.emplace_back(currLine);
		}
	}

	int maxPId = *(std::max_element(pointId.begin(), pointId.end()));
	std::vector<Eigen::Vector2d> finalPPos;
	finalPPos.resize(maxPId + 1);
	for (int i = 0; i < pointId.size(); ++i)
	{
		finalPPos[pointId[i]] = pointPos[i];
	}

	int maxLId = *(std::max_element(lineId.begin(), lineId.end()));
	std::vector<Eigen::Vector2i> finalLPos;
	finalLPos.resize(maxLId + 1);
	for (int i = 0; i < lineId.size(); ++i)
	{
		finalLPos[lineId[i]] = linePairs[i];
	}

	for (int i = 0; i < lineLoop.size(); ++i)
	{
		poly.emplace_back(finalPPos[finalLPos[lineLoop[i]][0]]);
	}
}

void GeoPolyReadWrite::WritePolyToGeo(const std::string &fileName, std::vector<Eigen::Vector2d> &poly)
{
	if (poly.empty())
		return;

	std::ofstream ofs(fileName.c_str());
	uint32_t currNum = 1;
	ofs << "lc = 0.3;" << std::endl;
	for (int i = 0; i < poly.size(); ++i)
	{
		ofs << "Point(" << currNum++ << ") = {" << poly[i][0] << "," << poly[i][1] << ",0.0,lc};" << std::endl;
	}

	uint32_t pointNum = poly.size();
	currNum = 0;
	for (int i = 0; i < poly.size(); ++i)
	{
		ofs << "Line(" << currNum + 1 << ") = {" << (currNum%pointNum) + 1 << "," << ((currNum + 1) % pointNum) + 1 << "};" << std::endl;
		++currNum;
	}
	ofs << "Line Loop(1) = {";
	for (int i = 0; i < poly.size()-1; ++i)
	{
		ofs << i + 1 << ",";
	}
	ofs << poly.size() << "};" << std::endl;

	ofs << "Plane Surface(1) = {1};" << std::endl;
}

void GeoPolyReadWrite::ReadTxtToPoly(const std::string &fileName, std::vector<Eigen::Vector2d> &poly)
{
	poly.clear();

	std::ifstream ifs(fileName.c_str());
	
	int pointNum = 0;
	double value0 = 0, value1 = 0;
	ifs >> pointNum;
	for (int i = 0; i < pointNum; ++i)
	{
		ifs >> value0 >> value1;
		poly.emplace_back(Eigen::Vector2d(value0, value1));
	}
}

void GeoPolyReadWrite::ScalePoly(std::vector<Eigen::Vector2d> &poly, double targetDiagLength)
{
	Eigen::Vector2d minP(1E20, 1E20), maxP(-1E20, -1E20);
	for (int i = 0; i < poly.size(); ++i)
	{
		if (minP[0] > poly[i][0])
			minP[0] = poly[i][0];
		if (minP[1] > poly[i][1])
			minP[1] = poly[i][1];

		if (maxP[0] < poly[i][0])
			maxP[0] = poly[i][0];
		if (maxP[1] < poly[i][1])
			maxP[1] = poly[i][1];
	}

	double currDiagLength = (maxP - minP).norm();
	double radio = targetDiagLength / currDiagLength;

	for (int i = 0; i < poly.size(); ++i)
	{
		poly[i] = (poly[i] - minP) * radio + minP;
	}
}