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


namespace CfvMath
{

	//getMatrixGeneralInverse just operates with the right down (N-1,N-1) part of A
	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> & A, 
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> & IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 4> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 4> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> & A,
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> & IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 7> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 7> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> & A,
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> & IA);

	//ScalarCfv::real getMoment(
	//	ScalarCfv::point p, 
	//	ScalarCfv::point baryCenter, 
	//	ScalarCfv::point scale, 
	//	int index);

	//rO=1
	bool getMoment(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> & A
		);
	bool getMomentRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &A);

	//rO=2
	bool getMoment(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> & A
		);

	bool getMomentRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &A);

	//rO=3
	bool getMoment(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> & A
		);

	const int mMapping[10] = { 0, 1, 0, 2, 1, 0, 3, 2, 1, 0 };//2 base derivative order index <-> 1 parameter
	const int nMapping[10] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };

	int getCombination(int  m, int  n);//C^n_(m+n)
	int getFactorial(int  m, int  n);//(m+n)!

	ScalarCfv::real harmonic(ScalarCfv::real uL, ScalarCfv::real uR);

	//rO=1
	bool getDiffBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> & moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> & A
		);
	bool getDiffBaseValueRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 6> &A);
	//rO=2
	bool getDiffBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> & moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> & A
		);
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
	//rO=3
	bool getDiffBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> & moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> & A
		);

		

	//rO=1
	bool getBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> & moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> & A
		);
	bool getBaseValueRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &A);
	//rO=2
	bool getBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> & moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> & A
		);
	bool getBaseValueRBFA1(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &A);
	//rO=3
	bool getBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> & moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> & A
		);

	template<typename T>
	T getSign(T  a) {
		if (a <= 0) {
			return -1;
		}
		else {
			return 1;
		}
	}

	ScalarCfv::real W12(ScalarCfv::real u[], const int J);

}
#endif