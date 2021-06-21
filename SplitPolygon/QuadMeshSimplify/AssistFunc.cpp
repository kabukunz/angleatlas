#include "AssistFunc.h"
#include <cmath>

namespace AssistFunctions
{

	AssistFunc::AssistFunc()
	{
	}


	AssistFunc::~AssistFunc()
	{
	}

	void AssistFunc::ComputeSVDDecomposition(double *inputMat, double *U, double &s1, double &s2, double *V)
	{
		double a = inputMat[0], b = inputMat[2], c = inputMat[1], d = inputMat[3];
		double theta = std::atan2(2 * a*c + 2 * b*d, a*a + b * b - c * c - d * d)*0.5;
		U[0] = std::cos(theta); U[1] = std::sin(theta); U[2] = -U[1]; U[3] = U[0];

		long double tempS1 = a * a + b * b + c * c + d * d;
		long double tempS2 = std::sqrt((a*a + b * b - c * c - d * d)*(a*a + b * b - c * c - d * d) + 4 * (a*c + b * d)*(a*c + b * d));
		s1 = std::sqrt((tempS1 + tempS2)*0.5);
		s2 = std::sqrt((tempS1 - tempS2)*0.5);
		long double tempValue = 4 * (long double(a)*long double(d) - long double(b) * long double(c)) * (long double(a)*long double(d) - long double(b) * long double(c));
		s2 = std::sqrt((tempValue/(tempS1 + tempS2))*0.5);
		//s2 = std::sqrt(a*a*0.5 - std::sqrt(4 * (a*c + b * d)* (a*c + b * d) + (a*a + b*b - c*c - d*d)*(a*a + b * b - c * c - d * d))*0.5 + b*b*0.5 + c*c*0.5 + d*d*0.5);
		/*if (s1 < 1E-6)
			s1 = 0;
		if (s2 < 1E-6)
			s2 = 0;*/

		double fye = atan2(2 * a*b + 2 * c*d, a*a - b * b + c * c - d * d)*0.5;
		double cosFye = std::cos(fye), sinFye = std::sin(fye);
		double s11 = (a * U[0] + c * U[1])*cosFye + (b*U[0] + d * U[1])*sinFye;
		double s22 = (a*U[1] - c * U[0])*sinFye + (-b * U[1] + d * U[0])*cosFye;


		double s11Sign, s22Sign;
		if (s11 > 0)
			s11Sign = 1;
		else if (s11 == 0)
			s11Sign = 0;
		else
			s11Sign = -1;

		if (s22 > 0)
			s22Sign = 1;
		else if (s22 == 0)
			s22Sign = 0;
		else
			s22Sign = -1;

		V[0] = s11Sign * cosFye; V[1] = s11Sign * sinFye; V[2] = -s22Sign * sinFye; V[3] = s22Sign * cosFye;
	}

	void AssistFunc::ComputeMinVPos(double x0, double y0, double x1, double y1, double theta, double &res0, double &res1)
	{
		double tempV = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
		if (tempV < 1E-15)
		{
			res0 = 1E30;
			res1 = 1E30;
			return;
		}
		double b = ((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0)) * std::cos(theta*0.5) * std::cos(theta*0.5) / (std::sin(theta) * std::sin(theta));

		res0 = -std::sqrt(-(tempV) *(tempV - 4 * b)) / (2 * tempV);
		res1 = std::sqrt(-(tempV) *(tempV - 4 * b)) / (2 * tempV);
	}

	double AssistFunc::ComputeAbsCross(const Eigen::Vector2d &vec0, const Eigen::Vector2d &vec1)
	{
		return std::abs(vec0[0] * vec1[1] - vec0[1] * vec1[0]);
	}
}