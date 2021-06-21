#include "FastMath.h"
#include <cmath>

namespace ig
{
	const double HALF_PI = 1.5707963267948966192313216916398;

	//----------------------------------------------------------------------------
	double FastSin0(double fAngle)
	{
		double fASqr = fAngle*fAngle;
		double fResult = 7.61e-03;
		fResult *= fASqr;
		fResult -= 1.6605e-01;
		fResult *= fASqr;
		fResult += 1.0;
		fResult *= fAngle;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastSin1(double fAngle)
	{
		double fASqr = fAngle*fAngle;
		double fResult = -2.39e-08;
		fResult *= fASqr;
		fResult += 2.7526e-06;
		fResult *= fASqr;
		fResult -= 1.98409e-04;
		fResult *= fASqr;
		fResult += 8.3333315e-03;
		fResult *= fASqr;
		fResult -= 1.666666664e-01;
		fResult *= fASqr;
		fResult += 1.0;
		fResult *= fAngle;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastCos0(double fAngle)
	{
		double fASqr = fAngle*fAngle;
		double fResult = 3.705e-02;
		fResult *= fASqr;
		fResult -= 4.967e-01;
		fResult *= fASqr;
		fResult += 1.0;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastCos1(double fAngle)
	{
		double fASqr = fAngle*fAngle;
		double fResult = -2.605e-07;
		fResult *= fASqr;
		fResult += 2.47609e-05;
		fResult *= fASqr;
		fResult -= 1.3888397e-03;
		fResult *= fASqr;
		fResult += 4.16666418e-02;
		fResult *= fASqr;
		fResult -= 4.999999963e-01;
		fResult *= fASqr;
		fResult += 1.0;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastTan0(double fAngle)
	{
		double fASqr = fAngle*fAngle;
		double fResult = 2.033e-01;
		fResult *= fASqr;
		fResult += 3.1755e-01;
		fResult *= fASqr;
		fResult += 1.0;
		fResult *= fAngle;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastTan1(double fAngle)
	{
		double fASqr = fAngle*fAngle;
		double fResult = 9.5168091e-03;
		fResult *= fASqr;
		fResult += 2.900525e-03;
		fResult *= fASqr;
		fResult += 2.45650893e-02;
		fResult *= fASqr;
		fResult += 5.33740603e-02;
		fResult *= fASqr;
		fResult += 1.333923995e-01;
		fResult *= fASqr;
		fResult += 3.333314036e-01;
		fResult *= fASqr;
		fResult += 1.0;
		fResult *= fAngle;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastInvSin0(double fValue)
	{
		double fRoot = sqrt((1.0) - fValue);
		double fResult = -0.0187293;
		fResult *= fValue;
		fResult += 0.0742610;
		fResult *= fValue;
		fResult -= 0.2121144;
		fResult *= fValue;
		fResult += 1.5707288;
		fResult = HALF_PI - fRoot*fResult;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastInvSin1(double fValue)
	{
		double fRoot = sqrt(fabs((1.0) - fValue));
		double fResult = -0.0012624911;
		fResult *= fValue;
		fResult += 0.0066700901;
		fResult *= fValue;
		fResult -= 0.0170881256;
		fResult *= fValue;
		fResult += 0.0308918810;
		fResult *= fValue;
		fResult -= 0.0501743046;
		fResult *= fValue;
		fResult += 0.0889789874;
		fResult *= fValue;
		fResult -= 0.2145988016;
		fResult *= fValue;
		fResult += 1.5707963050;
		fResult = HALF_PI - fRoot*fResult;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastInvCos0(double fValue)
	{
		double fRoot = sqrt((1.0) - fValue);
		double fResult = -0.0187293;
		fResult *= fValue;
		fResult += 0.0742610;
		fResult *= fValue;
		fResult -= 0.2121144;
		fResult *= fValue;
		fResult += 1.5707288;
		fResult *= fRoot;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastInvCos1(double fValue)
	{
		double fRoot = sqrt(fabs((1.0) - fValue));
		double fResult = -0.0012624911;
		fResult *= fValue;
		fResult += 0.0066700901;
		fResult *= fValue;
		fResult -= 0.0170881256;
		fResult *= fValue;
		fResult += 0.0308918810;
		fResult *= fValue;
		fResult -= 0.0501743046;
		fResult *= fValue;
		fResult += 0.0889789874;
		fResult *= fValue;
		fResult -= 0.2145988016;
		fResult *= fValue;
		fResult += 1.5707963050;
		fResult *= fRoot;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastInvTan0(double fValue)
	{
		double fVSqr = fValue*fValue;
		double fResult = 0.0208351;
		fResult *= fVSqr;
		fResult -= 0.085133;
		fResult *= fVSqr;
		fResult += 0.180141;
		fResult *= fVSqr;
		fResult -= 0.3302995;
		fResult *= fVSqr;
		fResult += 0.999866;
		fResult *= fValue;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastInvTan1(double fValue)
	{
		double fVSqr = fValue*fValue;
		double fResult = 0.0028662257;
		fResult *= fVSqr;
		fResult -= 0.0161657367;
		fResult *= fVSqr;
		fResult += 0.0429096138;
		fResult *= fVSqr;
		fResult -= 0.0752896400;
		fResult *= fVSqr;
		fResult += 0.1065626393;
		fResult *= fVSqr;
		fResult -= 0.1420889944;
		fResult *= fVSqr;
		fResult += 0.1999355085;
		fResult *= fVSqr;
		fResult -= 0.3333314528;
		fResult *= fVSqr;
		fResult += 1.0;
		fResult *= fValue;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastNegExp0(double fValue)
	{
		double fResult = 0.0038278;
		fResult *= fValue;
		fResult += 0.0292732;
		fResult *= fValue;
		fResult += 0.2507213;
		fResult *= fValue;
		fResult += 1.0;
		fResult *= fResult;
		fResult *= fResult;
		fResult = (1.0) / fResult;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastNegExp1(double fValue)
	{
		double fResult = 0.00026695;
		fResult *= fValue;
		fResult += 0.00227723;
		fResult *= fValue;
		fResult += 0.03158565;
		fResult *= fValue;
		fResult += 0.24991035;
		fResult *= fValue;
		fResult += 1.0;
		fResult *= fResult;
		fResult *= fResult;
		fResult = (1.0) / fResult;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastNegExp2(double fValue)
	{
		double fResult = 0.000014876;
		fResult *= fValue;
		fResult += 0.000127992;
		fResult *= fValue;
		fResult += 0.002673255;
		fResult *= fValue;
		fResult += 0.031198056;
		fResult *= fValue;
		fResult += 0.250010936;
		fResult *= fValue;
		fResult += 1.0;
		fResult *= fResult;
		fResult *= fResult;
		fResult = (1.0) / fResult;
		return fResult;
	}
	//----------------------------------------------------------------------------
	double FastNegExp3(double fValue)
	{
		double fResult = 0.0000006906;
		fResult *= fValue;
		fResult += 0.0000054302;
		fResult *= fValue;
		fResult += 0.0001715620;
		fResult *= fValue;
		fResult += 0.0025913712;
		fResult *= fValue;
		fResult += 0.0312575832;
		fResult *= fValue;
		fResult += 0.2499986842;
		fResult *= fValue;
		fResult += 1.0;
		fResult *= fResult;
		fResult *= fResult;
		fResult = (1.0) / fResult;
		return fResult;
	}
	//----------------------------------------------------------------------------
}