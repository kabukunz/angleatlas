#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

namespace CommonFunctions
{
	inline int lfloor(double x) { return std::lround(std::floor(x)); }
	inline int lceil(double x) { return std::lround(std::ceil(x)); }

	template <int N>
	inline double vec_angle_acos(const OpenMesh::VectorT<double, N>& a, const OpenMesh::VectorT<double, N>& b)
	{
		return std::acos(std::max(-1.0, std::min(1.0, OpenMesh::dot(a.normalized(), b.normalized()))));
	}

	template <int N>
	inline double vec_angle_atan2(const OpenMesh::VectorT<double, N>& a, const OpenMesh::VectorT<double, N>& b)
	{
		return std::atan2(a[0] * b[1] - a[1] * b[0], a[0] * b[0] + a[1] * b[1]);
	}

	inline int period_id(int x, int size)
	{
		return (x < 0) ? (x + size) : (x >= size ? x - size : x);
	}

	template <int N>
	inline int get_tag(const OpenMesh::VectorT<double, N>& v)
	{
		uint b = (abs(v[0]) <= abs(v[1])) ? 1 : 0;
		uint a = (v[b] <= 0.0);

		return a * 2 + b;
	}
}

#endif