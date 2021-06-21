#ifndef FAST_MATH_H
#define FAST_MATH_H

namespace ig
{

	// The input must be in [0,pi/2].
	// max error sin0 = 1.7e-04, speed up = 4.0
	// max error sin1 = 1.9e-08, speed up = 2.8
	double FastSin0(double fAngle);
	double FastSin1(double fAngle);

	// The input must be in [0,pi/2]
	// max error cos0 = 1.2e-03, speed up = 4.5
	// max error cos1 = 6.5e-09, speed up = 2.8
	double FastCos0(double fAngle);
	double FastCos1(double fAngle);

	// The input must be in [0,pi/4].
	// max error tan0 = 8.1e-04, speed up = 5.6
	// max error tan1 = 1.9e-08, speed up = 3.4
	double FastTan0(double fAngle);
	double FastTan1(double fAngle);

	// The input must be in [0,1].
	// max error invsin0 = 6.8e-05, speed up = 7.5
	// max error invsin1 = 1.4e-07, speed up = 5.5
	double FastInvSin0(double fValue);
	double FastInvSin1(double fValue);

	// The input must be in [0,1].
	// max error invcos0 = 6.8e-05, speed up = 7.5
	// max error invcos1 = 1.4e-07, speed up = 5.7
	double FastInvCos0(double fValue);
	double FastInvCos1(double fValue);

	// The input must be in [-1,1]. 
	// max error invtan0 = 1.2e-05, speed up = 2.8
	// max error invtan1 = 2.3e-08, speed up = 1.8
	double FastInvTan0(double fValue);
	double FastInvTan1(double fValue);

	// A fast approximation to 1/sqrt.
	double FastInvSqrt(double fValue);

	// Fast approximations to exp(-x).  The input x must be in [0,infinity).
	// max error negexp0 = 0.00024, speed up = 25.4
	// max error negexp1 = 0.000024, speed up = 25.4
	// max error negexp2 = 0.0000024, speed up = 20.5
	// max error negexp3 = 0.00000025, speed up = 17.3
	double FastNegExp0(double fValue);
	double FastNegExp1(double fValue);
	double FastNegExp2(double fValue);
	double FastNegExp3(double fValue);
}
#endif


