#pragma once
#include <Eigen/Sparse>
#include <vector>

class CheckPolyProperty
{
public:
	CheckPolyProperty();
	~CheckPolyProperty();

	static bool CheckPolySimple(std::vector<Eigen::Vector2d> &poly);
	static bool CheckPointInsidePoly(std::vector<Eigen::Vector2d> &poly, Eigen::Vector2d &point);
	static int CheckPointPolyState(std::vector<Eigen::Vector2d> &poly, Eigen::Vector2d &point);
};

