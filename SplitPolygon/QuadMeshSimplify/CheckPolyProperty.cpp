#include "CheckPolyProperty.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <iostream>
#include "CheckPolyProperty.h"
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
using std::cout; using std::endl;


CheckPolyProperty::CheckPolyProperty()
{
}


CheckPolyProperty::~CheckPolyProperty()
{
}

bool CheckPolyProperty::CheckPolySimple(std::vector<Eigen::Vector2d> &poly)
{
	Point *points = new Point[poly.size()];
	for (int i = 0; i < poly.size(); ++i)
	{
		points[i] = Point(poly[i][0], poly[i][1]);
	}

	bool result = CGAL::is_simple_2(points, points + poly.size(), K());
	delete[] points;
	return result;
}

int CheckPolyProperty::CheckPointPolyState(std::vector<Eigen::Vector2d> &poly, Eigen::Vector2d &point)
{
	Point *points = new Point[poly.size()];
	for (int i = 0; i < poly.size(); ++i)
	{
		points[i] = Point(poly[i][0], poly[i][1]);
	}
	Point chP(point[0], point[1]);

	int result;
	if (CGAL::bounded_side_2(points, points + poly.size(), chP, K()) == CGAL::ON_UNBOUNDED_SIDE)
		result = 0;
	else if (CGAL::bounded_side_2(points, points + poly.size(), chP, K()) == CGAL::ON_BOUNDARY)
		result = 1;
	else
		result = 2;

	delete[] points;
	return result;
}

bool CheckPolyProperty::CheckPointInsidePoly(std::vector<Eigen::Vector2d> &poly, Eigen::Vector2d &point)
{
	Point *points = new Point[poly.size()];
	for (int i = 0; i < poly.size(); ++i)
	{
		points[i] = Point(poly[i][0], poly[i][1]);
	}
	Point chP(point[0], point[1]);

	bool result = true;
	if (CGAL::bounded_side_2(points, points + poly.size(), chP, K()) != CGAL::ON_UNBOUNDED_SIDE)
		result = true;
	else
		result = false;
	
	delete[] points;
	return result;
}