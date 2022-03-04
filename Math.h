#ifndef _MATH_H_
#define _MATH_H_

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "mkl_cblas.h"
#include "mkl_lapacke.h"
#include "TypeDefine.h"
#include "Tensor.h"
#include "Point.h"
#include "Variable.h"
#include "EigenTensor.hpp"

namespace CfvMath
{

	// getMatrixGeneralInverse just operates with the right down (N-1,N-1) part of A
	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 4> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 4> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 7> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 7> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> &IA);

	// ScalarCfv::real getMoment(
	//	ScalarCfv::point p,
	//	ScalarCfv::point baryCenter,
	//	ScalarCfv::point scale,
	//	int index);

	// rO=1
	bool getMoment(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &A);
	bool getMomentRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &A);
	bool getMomentRBFB1(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &A,
		ScalarCfv::cellFieldData &cell); // adding 0
	bool getMomentRBFB1(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &A,
		ScalarCfv::cellFieldData &cell); // adding 1

	// rO=2
	bool getMoment(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> &A);

	bool getMomentRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &A);

	// rO=3
	bool getMoment(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> &A);

	const int mMapping[10] = {0, 1, 0, 2, 1, 0, 3, 2, 1, 0}; // 2 base derivative order index <-> 1 parameter
	const int nMapping[10] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3};

	int getCombination(int m, int n); // C^n_(m+n)
	int getFactorial(int m, int n);	  //(m+n)!

	ScalarCfv::real harmonic(ScalarCfv::real uL, ScalarCfv::real uR);

	// rO=1
	bool getDiffBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> &A);
	bool getDiffBaseValueRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 6> &A);
	bool getDiffBaseValueRBFB1(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getDiffBaseValueRBFB1(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 6> &A,
		ScalarCfv::cellFieldData &cell); // adding 1

	// rO=2
	bool getDiffBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> &A);
	bool getDiffBaseValueRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 6> &A);
	bool getDiffBaseValueRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 10> &A);
	// rO=3
	bool getDiffBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> &A);

	// rO=1
	bool getBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &A);
	bool getBaseValueRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &A);
	bool getBaseValueRBFB1(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &A,
		ScalarCfv::cellFieldData &cell); // adding 0
	bool getBaseValueRBFB1(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &A,
		ScalarCfv::cellFieldData &cell); // adding 1
	// rO=2
	bool getBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> &A);
	bool getBaseValueRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &A);
	// rO=3
	bool getBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> &A);

	template <typename T>
	T getSign(T a)
	{
		if (a <= 0)
		{
			return -1;
		}
		else
		{
			return 1;
		}
	}

	ScalarCfv::real W12(ScalarCfv::real u[], const int J);

	inline ScalarCfv::point GetFaceParam(ScalarCfv::cellType cellType, int ff, const ScalarCfv::point &faceP, bool inverse = false)
	{
		ScalarCfv::point pparaml;
		if (cellType == ScalarCfv::Quadrilateral)
			switch (ff)
			{
			case 1: // 12 left
				pparaml.y = inverse ? 1 - faceP.x : faceP.x;
				pparaml.x = 0;
				break;
			case 2: // 23 up
				pparaml.x = inverse ? 1 - faceP.x : faceP.x;
				pparaml.y = 1;
				break;
			case 3: // 34 right
				pparaml.y = inverse ? faceP.x : 1 - faceP.x;
				pparaml.x = 1;
				break;
			case 4: // 14 down
				pparaml.x = inverse ? faceP.x : 1 - faceP.x;
				pparaml.y = 0;
				break;

			default:
				std::cout << "ff is " << ff << std::endl;
				assert(false);
				break;
			}
		else
		{
			std::cout << "cell Type Wrong\n";
			assert(false);
		}
		return pparaml;
	}

	inline ScalarCfv::point getPoint(ScalarCfv::point p, ScalarCfv::cellFieldData &cell)
	{
		Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
		Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
										 {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
		Eigen::Vector2d pp = XiNj * Nj;
		return ScalarCfv::point(pp(0), pp(1));
	}
}
#endif