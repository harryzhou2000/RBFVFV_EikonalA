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
	const ScalarCfv::real crbf = 1;
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

	bool getDiffBaseValueRBFB1_4_3(
		ScalarCfv::point p,			 // parametric place
		ScalarCfv::point baryCenter, // dummy
		ScalarCfv::point scale,		 // dummy
		ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 4, 3> &A,
		ScalarCfv::cellFieldData &cell) // adding 1
		;
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

#ifndef RBFB1_GlobalPoly
	// Cik = sum_j -- Aij WjWj Bkj with tensored diffs
	template <class TA, class TB, class TC, class TW>
	void VVMatConjProd2DDiffCombine(TA &A, TB &B, TC &C, TW &W,
									int istart, int iend, int jstart, int jend, int kstart, int kend, bool clear = true)
	{
		if (clear)
			for (int i = istart; i < iend; i++)
				for (int k = kstart; k < kend; k++)
					C[i][k] = 0.0;
		if (jstart == 0 && jend == 6)
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
					C[i][k] += 0.5 * (csumA * csumB + csumD * csumW);
					// C[i][k] += csumA * csumB;
					csumA = csumB = csumD = csumW = 0;
					for (int j = 3; j < 6; j++)
					{
						csumA += A[i][j] * (W[j]);
						csumB += B[k][j] * (W[j]);
						csumD += A[i][j] * B[k][j];
						csumW += W[j] * W[j];
					}
					csumW -= W[4] * W[4] * 0.5;
					csumD += A[i][4] * B[k][4];
					C[i][k] += csumA * csumB;
					// C[i][k] += 0.5 * (csumA * csumB + csumD * csumW);
					// C[i][k] = 0.67 * C[i][k] + 0.33 * (A[i][3] + A[i][5]) * (B[k][3] + B[k][5]) * (W[3] + W[5]) * (W[3] + W[5]);
				}
		else if (jstart == 0 && jend == 10)
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
					C[i][k] += 0.5 * (csumA * csumB + csumD * csumW);
					// C[i][k] += csumA * csumB;
					csumA = csumB = csumD = csumW = 0;
					for (int j = 3; j < 6; j++)
					{
						csumA += A[i][j] * (W[j]);
						csumB += B[k][j] * (W[j]);
					}
					C[i][k] += csumA * csumB;
					csumA = csumB = csumD = csumW = 0;
					for (int j = 6; j < 10; j++)
					{
						csumA += A[i][j] * (W[j]);
						csumB += B[k][j] * (W[j]);
					}
					C[i][k] += csumA * csumB;
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
				}
		else
		{
			assert(false);
		}
	}
#else
	template <class TA, class TB, class TC, class TW>
	void VVMatConjProd2DDiffCombine(TA &A, TB &B, TC &C, TW &W,
									int istart, int iend, int jstart, int jend, int kstart, int kend, bool clear = true)
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