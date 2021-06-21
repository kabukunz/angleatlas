#pragma once
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>
typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef std::vector<Triangle>::iterator Triangle_Iterator;
typedef CGAL::AABB_triangle_primitive<K, Triangle_Iterator> primitive_;
typedef CGAL::AABB_traits<K, primitive_> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

typedef K::Segment_3 Segment;
typedef std::vector<Segment>::iterator Line_Iterator;
typedef CGAL::AABB_segment_primitive<K, Line_Iterator> Line_Primitive;
typedef CGAL::AABB_traits<K, Line_Primitive> Traits;
typedef CGAL::AABB_tree<Traits> SegmentTree;