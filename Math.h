#ifndef _MATH_H_
#define _MATH_H_

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "TypeDefine.h"
#include "Tensor.h"
#include "Point.h"
#include "Variable.h"
#include "EigenTensor.hpp"

namespace CfvMath
{
	const ScalarCfv::real crbf = RBFB1_CRBF;
	// getMatrixGeneralInverse just operates with the right down (N-1,N-1) part of A
	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 4> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 4> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 5, 5> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 5, 5> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 7> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 7> &IA);

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> &A,
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> &IA);

	// template <unsigned mm, unsigned nn>
	// bool getMatrixGeneralInverse(
	// 	ScalarCfv::tensor2D<ScalarCfv::real, mm, nn> &A,
	// 	ScalarCfv::tensor2D<ScalarCfv::real, nn, mm> &IA)
	// {
	// 	int mrow = mm - 1;
	// 	int ncol = nn - 1;
	// 	MKL_INT m = mrow, n = ncol, lda = ncol, ldu = mrow, ldvt = ncol, info;
	// 	ScalarCfv::real temp_eps = 1e-16, alpha = 1.0, beta = 0.0;
	// 	ScalarCfv::real *a = new ScalarCfv::real[m * n];
	// 	ScalarCfv::real *inva = new ScalarCfv::real[n * m];
	// 	ScalarCfv::real *temp = new ScalarCfv::real[n * m];
	// 	MKL_INT min = (m < n) ? m : n;
	// 	ScalarCfv::real *superb = new ScalarCfv::real[min - 1];
	// 	ScalarCfv::real *s = new ScalarCfv::real[n];
	// 	ScalarCfv::real *splus = new ScalarCfv::real[n * m];
	// 	ScalarCfv::real *u = new ScalarCfv::real[ldu * m];
	// 	ScalarCfv::real *vt = new ScalarCfv::real[ldvt * n];

	// 	for (int ii = 0; ii < mrow; ++ii)
	// 	{
	// 		for (int jj = 0; jj < ncol; ++jj)
	// 		{
	// 			a[ii * n + jj] = A[ii + 1][jj + 1];
	// 		}
	// 	}

	// 	// SVD
	// 	info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m, n, a, lda, s, u, ldu, vt, ldvt, superb);
	// 	if (info > 0)
	// 	{
	// 		std::cout << "The algorithm computing SVD failed to converge." << std::endl;
	// 		exit(1);
	// 	}

	// 	for (int ii = 0; ii < ncol; ++ii)
	// 	{
	// 		for (int jj = 0; jj < mrow; ++jj)
	// 		{
	// 			splus[ii * mrow + jj] = 0.0;
	// 			if (jj == ii)
	// 			{
	// 				if (abs(s[ii]) > temp_eps)
	// 				{
	// 					splus[ii * mrow + jj] = 1.0 / s[ii];
	// 				}
	// 			}
	// 		}
	// 	}
	// 	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, n, alpha, vt, n, splus, m, beta, temp, m);
	// 	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, m, alpha, temp, m, u, m, beta, inva, m);

	// 	for (int ii = 0; ii < ncol; ++ii)
	// 	{
	// 		for (int jj = 0; jj < mrow; ++jj)
	// 		{
	// 			IA[ii + 1][jj + 1] = inva[ii * mrow + jj];
	// 		}
	// 	}

	// 	delete[] a;
	// 	delete[] inva;
	// 	delete[] temp;
	// 	delete[] superb;
	// 	delete[] s;
	// 	delete[] splus;
	// 	delete[] u;
	// 	delete[] vt;
	// 	a = NULL;
	// 	inva = NULL;
	// 	temp = NULL;
	// 	superb = NULL;
	// 	s = NULL;
	// 	splus = NULL;
	// 	u = NULL;
	// 	vt = NULL;

	// 	return true;
	// }

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

	bool getMomentRBFB1(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
		ScalarCfv::cellFieldData &cell); // adding 4

	bool getMomentRBFB1_7_3(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
		ScalarCfv::cellFieldData &cell); // adding 4

	bool getMomentRBFB1_5_3(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 5> &A,
		ScalarCfv::cellFieldData &cell); // adding 4

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

	bool getDiffBaseValueRBFB1(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 6> &A,
		ScalarCfv::cellFieldData &cell); // adding 4
	bool getDiffBaseValueRBFB1(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 3> &A,
		ScalarCfv::cellFieldData &cell); // adding 4

	bool getDiffBaseValueRBFB1_7_3(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 3> &A,
		ScalarCfv::cellFieldData &cell); // adding 4

	bool getDiffBaseValueRBFB1_5_3(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 5> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 5, 3> &A,
		ScalarCfv::cellFieldData &cell); // adding 4

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
	bool getBaseValueRBFB1(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getBaseValueRBFB1_7_3(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getBaseValueRBFB1_5_3(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 5> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 5> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

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
	/**
	 * EXTENDED:
	 *
	 */

	//////////////////////////
	bool getMomentRBFB1_4_3(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &A,
		ScalarCfv::cellFieldData &cell) // adding 1
		;
	bool getBaseValueRBFB1_4_3(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &A,
		ScalarCfv::cellFieldData &cell) // adding 1
		;

	bool getDiffBaseValueRBFB1_4_6(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 6> &A,
		ScalarCfv::cellFieldData &cell) // adding 1
		;
	bool getMomentRBFB1_4_6(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &A,
		ScalarCfv::cellFieldData &cell) // adding 1
		;
	bool getBaseValueRBFB1_4_6(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &A,
		ScalarCfv::cellFieldData &cell) // adding 1
		;

	bool getDiffBaseValueRBFB1_4_3(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 3> &A,
		ScalarCfv::cellFieldData &cell) // adding 1
		;
	bool getMomentRBFB1_7_6(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
		ScalarCfv::cellFieldData &cell) // adding 1
		;
	bool getBaseValueRBFB1_7_6(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
		ScalarCfv::cellFieldData &cell) // adding 1
		;

	bool getDiffBaseValueRBFB1_7_6(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 7, 6> &A,
		ScalarCfv::cellFieldData &cell) // adding 1
		;
	bool getMomentRBFB1_POLY(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getBaseValueRBFB1_POLY(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getDiffBaseValueRBFB1_POLY(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 3> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getMomentRBFB1_POLY(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 6> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getBaseValueRBFB1_POLY(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 6> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getDiffBaseValueRBFB1_POLY(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 6> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getMomentRBFB1_POLY(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 10> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getBaseValueRBFB1_POLY(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 10> &moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> &A,
		ScalarCfv::cellFieldData &cell); // adding 0

	bool getDiffBaseValueRBFB1_POLY(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 10> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> &A,
		ScalarCfv::cellFieldData &cell); // adding 0
	//////////////////////////

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

	inline Eigen::Matrix2d getIJacobi(ScalarCfv::point p, ScalarCfv::cellFieldData &cell)
	{
		assert(cell.cellType_ == ScalarCfv::Quadrilateral);
		Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
		Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
											 {-(1 - p.x), (1 - p.x), p.x, -p.x}};
		Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
										 {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
		Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose();
		return Jacobi.inverse();
	}

	inline Eigen::Matrix2d getFaceJacobi(ScalarCfv::point p, ScalarCfv::cellFieldData &cell, int iedge)
	{
		assert(cell.cellType_ == ScalarCfv::Quadrilateral);
		Eigen::Matrix2d dxijdmi; // dxij/dmi, m is parametric norm
		switch (iedge)
		{
		case 1:
			dxijdmi << -1, 0, 0, -1;
			break;
		case 2:
			dxijdmi << 0, 1, -1, 0;
			break;
		case 3:
			dxijdmi << 1, 0, 0, 1;
			break;
		case 4:
			dxijdmi << 0, -1, 1, 0;
			break;
		default:
			assert(false);
			break;
		}
		Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
		Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
											 {-(1 - p.x), (1 - p.x), p.x, -p.x}};
		Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
										 {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
		Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose();
		return dxijdmi * Jacobi;
	}

	// Cik += Aij Bjk
	template <class TA, class TB, class TC>
	void VVMatMat(TA &A, TB &B, TC &C,
				  int istart, int iend, int jstart, int jend, int kstart, int kend, bool clear = true)
	{
		if (clear)
			for (int i = istart; i < iend; i++)
				for (int k = kstart; k < kend; k++)
					C[i][k] = 0.0;
		for (int i = istart; i < iend; i++)
			for (int j = jstart; j < jend; j++)
				for (int k = kstart; k < kend; k++)
					C[i][k] += A[i][j] * B[j][k];
	}

	// ci += Aij * bj
	template <class TA, class TB, class TC>
	void VVMatVec(TA &A, TB &b, TC &c,
				  int istart, int iend, int jstart, int jend, bool clear = true)
	{
		if (clear)
			for (int i = istart; i < iend; i++)
				c[i] = 0.0;
		for (int i = istart; i < iend; i++)
			for (int j = jstart; j < jend; j++)
				c[i] += A[i][j] * b[j];
	}

	// ci += Aij * bj * alpha
	template <class TA, class TB, class TC>
	void VVMatVec(TA &A, TB &b, TC &c, ScalarCfv::real alpha,
				  int istart, int iend, int jstart, int jend, bool clear = true)
	{
		if (clear)
			for (int i = istart; i < iend; i++)
				c[i] = 0.0;
		for (int i = istart; i < iend; i++)
		{
			ScalarCfv::real inc = 0;
			for (int j = jstart; j < jend; j++)
				inc += A[i][j] * b[j];
			c[i] += inc * alpha;
		}
	}

	// single Eigen
	template <class TA, class TB, class TC>
	void VEMatVec(TA &A, const TB &b, TC &c,
				  int istart, int iend, int jstart, int jend, bool clear = true)
	{
		if (clear)
			for (int i = istart; i < iend; i++)
				c[i] = 0.0;
		for (int i = istart; i < iend; i++)
			for (int j = jstart; j < jend; j++)
				c[i] += A[i][j] * b(j);
	}

	// B = A
	template <class TA, class TB>
	void VVMatCopy(TA &A, TB &B, int istart, int iend, int jstart, int jend)
	{
		for (int i = istart; i < iend; i++)
			for (int j = jstart; j < jend; j++)
				B[i][j] = A[i][j];
	}

	// B = A
	template <class TA, class TB>
	void VVVecCopy(TA &A, TB &B, int istart, int iend)
	{
		for (int i = istart; i < iend; i++)
			B[i] = A[i];
	}

	// B = A // single Eigen
	template <class TA, class TB>
	void VEVecMatCopy(TA &A, TB &B, int brow, int istart, int iend)
	{
		for (int i = istart; i < iend; i++)
			B(brow,i) = A[i];
	}

	// B = alpha * A
	template <class TA, class TB>
	void VVVecCopy(TA &A, TB &B, ScalarCfv::real alpha, int istart, int iend)
	{
		for (int i = istart; i < iend; i++)
			B[i] = A[i] * alpha;
	}

	// B += alpha * A
	template <class TA, class TB>
	void VVVecAdd(TA &A, TB &B, ScalarCfv::real alpha, int istart, int iend)
	{
		for (int i = istart; i < iend; i++)
			B[i] += A[i] * alpha;
	}

	// Cik = sum_j -- Aij Wj Bkj
	template <class TA, class TB, class TC, class TW>
	void VVMatConjProd(TA &A, TB &B, TC &C, TW &W,
					   int istart, int iend, int jstart, int jend, int kstart, int kend, bool clear = true)
	{
		if (clear)
			for (int i = istart; i < iend; i++)
				for (int k = kstart; k < kend; k++)
					C[i][k] = 0.0;
		for (int i = istart; i < iend; i++)
			for (int j = jstart; j < jend; j++)
				for (int k = kstart; k < kend; k++)
					C[i][k] += A[i][j] * B[k][j] * W[j];
	}

	inline void RegularizeJacobian(Eigen::Matrix2d &dxjdxii)
	{
		Eigen::Vector2d x0 = dxjdxii(0, Eigen::all);
		Eigen::Vector2d x1 = dxjdxii(1, Eigen::all);
		Eigen::Vector2d x0n = x0 / x0.norm();
		Eigen::Vector2d x1n = x1 / x1.norm();
		Eigen::Vector2d xcn = x1n + x0n;
		xcn /= xcn.norm();
		Eigen::Vector2d xmn{xcn(1), -xcn(0)};
		Eigen::Vector2d x0n_new = xcn + xmn;
		x0n_new /= x0n_new.norm();
		Eigen::Vector2d x1n_new = xcn - xmn;
		x1n_new /= x1n_new.norm();
		if(x0.norm() < x1.norm())
			x0n_new = x0n, x1n_new << -x0n_new(1), x0n_new(0);
		else
			x1n_new = x1n, x0n_new << x1n_new(1), -x1n_new(0);

		// dxjdxii(0, Eigen::all) = x0n_new * x0.dot(x0n_new);
		// dxjdxii(1, Eigen::all) = x1n_new * x1.dot(x1n_new);
	}

	// inline void RegularizeJacobianSym(Eigen::Matrix2d &dxjdxii)
	// {
	// 	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(0.5*(dxjdxii+dxjdxii.transpose()));


	// }

	// // TEMPLATE BUT INSTANTIATED AT math.cpp
	void EigenLeastSquareInverse(const Eigen::MatrixXd &A, Eigen::MatrixXd &AI);


#ifndef RBFB1_GlobalPoly_ESC
	// Cik = sum_j -- Aij WjWj Bkj with tensored diffs
	template <class TA, class TB, class TC, class TW>
	void VVMatConjProd2DDiffCombine(TA &A, TB &B, TC &C, TW &W,
									int istart, int iend, int jstart, int jend, int kstart, int kend,
									ScalarCfv::faceFieldData &face,
									const ScalarCfv::point &n = ScalarCfv::point(0, 0),
									bool clear = true)
	{
		if (clear)
			for (int i = istart; i < iend; i++)
				for (int k = kstart; k < kend; k++)
					C[i][k] = 0.0;

		if (n.length() > 0.0)
		{
			// if (jstart == 0 && jend == 3)
			// {
			// 	for (int i = istart; i < iend; i++)
			// 		for (int k = kstart; k < kend; k++)
			// 		{
			// 			C[i][k] += A[i][0] * B[k][0] * W[0] * W[0];
			// 			ScalarCfv::real dn1A = A[i][1] * n.x + A[i][2] * n.y;
			// 			ScalarCfv::real dn1B = B[k][1] * n.x + B[k][2] * n.y;
			// 			C[i][k] += dn1A * dn1B * W[1] * W[1];
			// 		}
			// }
			// else if (jstart == 0 && jend == 6)
			// {
			// 	for (int i = istart; i < iend; i++)
			// 		for (int k = kstart; k < kend; k++)
			// 		{
			// 			C[i][k] += A[i][0] * B[k][0] * W[0] * W[0];
			// 			ScalarCfv::real dn1A = A[i][1] * n.x + A[i][2] * n.y;
			// 			ScalarCfv::real dn1B = B[k][1] * n.x + B[k][2] * n.y;
			// 			C[i][k] += dn1A * dn1B * W[1] * W[1];
			// 			ScalarCfv::real dn2A = A[i][4] * n.x * n.x + A[i][5] * n.y * n.x * 2 + A[i][5] * n.y * n.x * 2;
			// 			ScalarCfv::real dn2B = B[k][1] * n.x + B[k][2] * n.y;
			// 		}
			// }
			// else if (jstart == 0 && jend == 10)
			// {
			// }
			// else
			assert(false);
		}
		else
		{
			if (jstart == 0 && jend == 6)
			{
				const ScalarCfv::real epsR = 1e-20;
				ScalarCfv::real w2r;
				ScalarCfv::real n1;
				ScalarCfv::real n2;
				ScalarCfv::real t1;
				ScalarCfv::real t2;
				if (std::fabs(W[1]) >= epsR || std::fabs(W[2]) >= epsR)
				{
					w2r = std::fabs(W[1]) >= epsR ? W[3] / (W[1] * W[1]) : W[5] / (W[2] * W[2]);
					n1 = face.interFacialJacobi(0, 0);
					n2 = face.interFacialJacobi(0, 1);
					t1 = face.interFacialJacobi(1, 0);
					t2 = face.interFacialJacobi(1, 1);
				}
				else
				{
					w2r = 0.0;
					n1 = n2 = t1 = t2 = 0.0;
				}

				for (int i = istart; i < iend; i++)
					for (int k = kstart; k < kend; k++)
					{
						C[i][k] += A[i][0] * B[k][0] * W[0] * W[0];
						ScalarCfv::real csumA;
						ScalarCfv::real csumB;
						csumA = A[i][1] * n1 + A[i][2] * n2;
						csumB = B[k][1] * n1 + B[k][2] * n2;
						C[i][k] += csumA * csumB;
						csumA = A[i][1] * t1 + A[i][2] * t2;
						csumB = B[k][1] * t1 + B[k][2] * t2;
						C[i][k] += csumA * csumB;

						csumA = (A[i][3] * n1 * n1 +
								 A[i][4] * n1 * n2 * 2 +
								 A[i][5] * n2 * n2) *
								w2r;
						csumB = (B[k][3] * n1 * n1 +
								 B[k][4] * n1 * n2 * 2 +
								 B[k][5] * n2 * n2) *
								w2r;
						C[i][k] += csumA * csumB;
						csumA = (A[i][3] * n1 * t1 +
								 A[i][4] * n1 * t2 * 2 +
								 A[i][5] * n2 * t2) *
								w2r;
						csumB = (B[k][3] * n1 * t1 +
								 B[k][4] * n1 * t2 * 2 +
								 B[k][5] * n2 * t2) *
								w2r;
						C[i][k] += csumA * csumB * 2;
						csumA = (A[i][3] * t1 * t1 +
								 A[i][4] * t1 * t2 * 2 +
								 A[i][5] * t2 * t2) *
								w2r;
						csumB = (B[k][3] * t1 * t1 +
								 B[k][4] * t1 * t2 * 2 +
								 B[k][5] * t2 * t2) *
								w2r;
						C[i][k] += csumA * csumB;
					}
				// // MIXED:
				// for (int i = istart; i < iend; i++)
				// 	for (int k = kstart; k < kend; k++)
				// 	{
				// 		C[i][k] += A[i][0] * B[k][0] * W[0] * W[0];
				// 		ScalarCfv::real csumA = 0.0;
				// 		ScalarCfv::real csumB = 0.0;
				// 		ScalarCfv::real csumD = 0.0;
				// 		ScalarCfv::real csumW = 0.0;
				// 		for (int j = 1; j < 3; j++)
				// 		{
				// 			csumA += A[i][j] * (W[j]);
				// 			csumB += B[k][j] * (W[j]);
				// 			csumD += A[i][j] * B[k][j];
				// 			csumW += W[j] * W[j];
				// 		}
				// 		C[i][k] += 0.5 * (csumA * csumB + csumD * csumW);
				// 		// C[i][k] += csumA * csumB;
				// 		csumA = csumB = csumD = csumW = 0;
				// 		for (int j = 3; j < 6; j++)
				// 		{
				// 			csumA += A[i][j] * (W[j]);
				// 			csumB += B[k][j] * (W[j]);
				// 			csumD += A[i][j] * B[k][j];
				// 			csumW += W[j] * W[j];
				// 		}
				// 		csumW -= W[4] * W[4] * 0.5;
				// 		csumD += A[i][4] * B[k][4];
				// 		C[i][k] += csumA * csumB;
				// 		// C[i][k] += 0.5 * (csumA * csumB + csumD * csumW);
				// 		// C[i][k] = 0.67 * C[i][k] + 0.33 * (A[i][3] + A[i][5]) * (B[k][3] + B[k][5]) * (W[3] + W[5]) * (W[3] + W[5]);
				// 	}
			}
			else if (jstart == 0 && jend == 10)
			{
				// ScalarCfv::real w2r =
				// 	W[3] / (W[1] * W[1]) +
				// 	W[4] / (W[1] * W[2]) / 2. +
				// 	W[5] / (W[2] * W[2]);
				// w2r /= 3.;
				// ScalarCfv::real w3r =
				// 	W[6] / (W[1] * W[1] * W[1]) +
				// 	W[7] / (W[1] * W[1] * W[2]) / 3. +
				// 	W[8] / (W[1] * W[2] * W[2]) / 3. +
				// 	W[9] / (W[2] * W[2] * W[2]);
				// w3r /= 4.;
				const ScalarCfv::real epsR = 1e-20;
				// std::cout << W[1] << "\t" << W[2] << std::endl;
				// assert(std::fabs(W[1]) >= epsR || std::fabs(W[2]) >= epsR);
				ScalarCfv::real w2r, w3r;
				ScalarCfv::real n1;
				ScalarCfv::real n2;
				ScalarCfv::real t1;
				ScalarCfv::real t2;
				if (std::fabs(W[1]) >= epsR || std::fabs(W[2]) >= epsR)
				{
					w2r = std::fabs(W[1]) >= epsR ? W[3] / (W[1] * W[1]) : W[5] / (W[2] * W[2]);
					w3r = std::fabs(W[1]) >= epsR ? W[6] / (W[1] * W[1] * W[1]) : W[9] / (W[2] * W[2] * W[2]);
					n1 = face.interFacialJacobi(0, 0);
					n2 = face.interFacialJacobi(0, 1);
					t1 = face.interFacialJacobi(1, 0);
					t2 = face.interFacialJacobi(1, 1);
				}
				else
				{
					w2r = w3r = 0.0;
					n1 = n2 = t1 = t2 = 0.0;
				}
				// w2r *= 0.2;
				// w3r *= 0.1;

				// std::cout << w2r << "\t" << w3r << std::endl;

				// ScalarCfv::real n1 = W[1], n2 = W[2];
				// ScalarCfv::real t1 = -W[2], t2 = W[1];

				for (int i = istart; i < iend; i++)
					for (int k = kstart; k < kend; k++)
					{
						C[i][k] += A[i][0] * B[k][0] * W[0] * W[0];
						ScalarCfv::real csumA;
						ScalarCfv::real csumB;
						csumA = A[i][1] * n1 + A[i][2] * n2;
						csumB = B[k][1] * n1 + B[k][2] * n2;
						C[i][k] += csumA * csumB;
						csumA = A[i][1] * t1 + A[i][2] * t2;
						csumB = B[k][1] * t1 + B[k][2] * t2;
						C[i][k] += csumA * csumB;

						csumA = (A[i][3] * n1 * n1 +
								 A[i][4] * n1 * n2 * 2 +
								 A[i][5] * n2 * n2) *
								w2r;
						csumB = (B[k][3] * n1 * n1 +
								 B[k][4] * n1 * n2 * 2 +
								 B[k][5] * n2 * n2) *
								w2r;
						C[i][k] += csumA * csumB;
						csumA = (A[i][3] * n1 * t1 +
								 A[i][4] * n1 * t2 * 2 +
								 A[i][5] * n2 * t2) *
								w2r;
						csumB = (B[k][3] * n1 * t1 +
								 B[k][4] * n1 * t2 * 2 +
								 B[k][5] * n2 * t2) *
								w2r;
						C[i][k] += csumA * csumB * 2;
						csumA = (A[i][3] * t1 * t1 +
								 A[i][4] * t1 * t2 * 2 +
								 A[i][5] * t2 * t2) *
								w2r;
						csumB = (B[k][3] * t1 * t1 +
								 B[k][4] * t1 * t2 * 2 +
								 B[k][5] * t2 * t2) *
								w2r;
						C[i][k] += csumA * csumB;

						csumA = (A[i][6] * n1 * n1 * n1 +
								 A[i][7] * n1 * n1 * n2 * 3 +
								 A[i][8] * n1 * n2 * n2 * 3 +
								 A[i][9] * n2 * n2 * n2) *
								w3r;
						csumB = (B[k][6] * n1 * n1 * n1 +
								 B[k][7] * n1 * n1 * n2 * 3 +
								 B[k][8] * n1 * n2 * n2 * 3 +
								 B[k][9] * n2 * n2 * n2) *
								w3r;
						C[i][k] += csumA * csumB;
						csumA = (A[i][6] * n1 * n1 * t1 +
								 A[i][7] * n1 * n1 * t2 * 3 +
								 A[i][8] * n1 * n2 * t2 * 3 +
								 A[i][9] * n2 * n2 * t2) *
								w3r;
						csumB = (B[k][6] * n1 * n1 * t1 +
								 B[k][7] * n1 * n1 * t2 * 3 +
								 B[k][8] * n1 * n2 * t2 * 3 +
								 B[k][9] * n2 * n2 * t2) *
								w3r;
						C[i][k] += csumA * csumB * 3;
						csumA = (A[i][6] * n1 * t1 * t1 +
								 A[i][7] * n1 * t1 * t2 * 3 +
								 A[i][8] * n1 * t2 * t2 * 3 +
								 A[i][9] * n2 * t2 * t2) *
								w3r;
						csumB = (B[k][6] * n1 * t1 * t1 +
								 B[k][7] * n1 * t1 * t2 * 3 +
								 B[k][8] * n1 * t2 * t2 * 3 +
								 B[k][9] * n2 * t2 * t2) *
								w3r;
						C[i][k] += csumA * csumB * 3;
						csumA = (A[i][6] * t1 * t1 * t1 +
								 A[i][7] * t1 * t1 * t2 * 3 +
								 A[i][8] * t1 * t2 * t2 * 3 +
								 A[i][9] * t2 * t2 * t2) *
								w3r;
						csumB = (B[k][6] * t1 * t1 * t1 +
								 B[k][7] * t1 * t1 * t2 * 3 +
								 B[k][8] * t1 * t2 * t2 * 3 +
								 B[k][9] * t2 * t2 * t2) *
								w3r;
						C[i][k] += csumA * csumB;
					}

				// // Type Mix
				// for (int i = istart; i < iend; i++)
				// 	for (int k = kstart; k < kend; k++)
				// 	{
				// 		C[i][k] += A[i][0] * B[k][0] * W[0] * W[0];
				// 		ScalarCfv::real csumA = 0.0;
				// 		ScalarCfv::real csumB = 0.0;
				// 		ScalarCfv::real csumD = 0.0;
				// 		ScalarCfv::real csumW = 0.0;
				// 		for (int j = 1; j < 3; j++)
				// 		{
				// 			csumA += A[i][j] * (W[j]);
				// 			csumB += B[k][j] * (W[j]);
				// 			csumD += A[i][j] * B[k][j];
				// 			csumW += W[j] * W[j];
				// 		}
				// 		C[i][k] += 0.5 * (csumA * csumB + csumD * csumW);
				// 		// C[i][k] += csumA * csumB;
				// 		csumA = csumB = csumD = csumW = 0;
				// 		for (int j = 3; j < 6; j++)
				// 		{
				// 			csumA += A[i][j] * (W[j]);
				// 			csumB += B[k][j] * (W[j]);
				// 		}
				// 		C[i][k] += csumA * csumB;
				// 		csumA = csumB = csumD = csumW = 0;
				// 		for (int j = 6; j < 10; j++)
				// 		{
				// 			csumA += A[i][j] * (W[j]);
				// 			csumB += B[k][j] * (W[j]);
				// 		}
				// 		C[i][k] += csumA * csumB;
				// 		C[i][k] *= 0.95;
				// 	}

				// for (int i = istart; i < iend; i++)
				// 	for (int k = kstart; k < kend; k++)
				// 		for (int j = 0; j < 10; j++)
				// 			C[i][k] += 0.05 * A[i][j] * B[k][j] * W[j] * W[j];

				// //Type original
				// for (int i = istart; i < iend; i++)
				// 	for (int k = kstart; k < kend; k++)
				// 		for (int j = 0; j < 10; j++)
				// 			C[i][k] += A[i][j] * B[k][j] * W[j] * W[j];
			}
			else if (jstart == 0 && jend == 3)
				for (int i = istart; i < iend; i++)
					for (int k = kstart; k < kend; k++)
					{
						C[i][k] += A[i][0] * B[k][0] * W[0] * W[0];
						ScalarCfv::real csumA = 0.0;
						ScalarCfv::real csumB = 0.0;
						ScalarCfv::real csumD = 0.0;
						ScalarCfv::real csumW = 0.0;
						for (int j = 1; j < 3; j++)
						{
							csumA += A[i][j] * (W[j]);
							csumB += B[k][j] * (W[j]);
							csumD += A[i][j] * B[k][j];
							csumW += W[j] * W[j];
						}
						C[i][k] += 0.5 * (csumA * csumB + csumD * csumW); // OK for RBFB1+0 rO=1
																		  // C[i][k] += 0.5 * (csumA * csumB); // OK for RBFB1+0 rO=1
						C[i][k] = A[i][0] * B[k][0] * W[0] * W[0] +
								  (A[i][1] * W[1] + A[i][2] * W[2]) * (B[k][1] * W[1] + B[k][2] * W[2]) +
								  (-A[i][1] * W[2] + A[i][2] * W[1]) * (-B[k][1] * W[2] + B[k][2] * W[1]);
					}
			else
			{
				assert(false);
			}
		}
	}
#else
	template <class TA, class TB, class TC, class TW>
	void VVMatConjProd2DDiffCombine(TA &A, TB &B, TC &C, TW &W,
									int istart, int iend, int jstart, int jend, int kstart, int kend,
									ScalarCfv::faceFieldData &face,
									const ScalarCfv::point &n = ScalarCfv::point(0, 0),
									bool clear = true)
	{
		if (clear)
			for (int i = istart; i < iend; i++)
				for (int k = kstart; k < kend; k++)
					C[i][k] = 0.0;
		for (int i = istart; i < iend; i++)
			for (int j = jstart; j < jend; j++)
				for (int k = kstart; k < kend; k++)
					C[i][k] += A[i][j] * B[k][j] * W[j] * W[j];
	}
#endif

}
#endif