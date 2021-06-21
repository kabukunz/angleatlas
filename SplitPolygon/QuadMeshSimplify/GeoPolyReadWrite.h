#pragma once
#include <string>
#include <Eigen/Sparse>
class GeoPolyReadWrite
{
public:
	GeoPolyReadWrite();
	~GeoPolyReadWrite();

	static void SplitString(const std::string &inStr, const std::string &key, std::vector<std::string> &outStrs);
	static void ReadGeoToPoly(const std::string &fileName, std::vector<Eigen::Vector2d> &poly);
	static void WritePolyToGeo(const std::string &fileName, std::vector<Eigen::Vector2d> &poly);

	static void ReadTxtToPoly(const std::string &fileName, std::vector<Eigen::Vector2d> &poly);

	static void ScalePoly(std::vector<Eigen::Vector2d> &poly, double targetDiagLength);
};

