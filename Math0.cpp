#include "Math.h"
//20200220 ���������ļ����ڵ���-1����0�ſ�ʼ
namespace CfvMath
{

	const ScalarCfv::real crbf = 2.0;
	//ScalarCfv::real getMoment(
	//	ScalarCfv::point p, 
	//	ScalarCfv::point baryCenter, 
	//	ScalarCfv::point scale, 
	//	int index){
	//	ScalarCfv::point delta = p - baryCenter;
	//	delta.x = delta.x / scale.x;
	//	delta.y = delta.y / scale.y;
	//	switch (index){
	//	case 1:
	//		return delta.x;
	//		break;
	//	case 2:
	//		return delta.y;
	//		break;
	//	case 3:
	//		return std::pow(delta.x, 2);
	//		break;
	//	case 4:
	//		return delta.x*delta.y;
	//		break;
	//	case 5:
	//		return std::pow(delta.y, 2);
	//		break;
	//	case 6:
	//		return std::pow(delta.x, 3);
	//		break;
	//	case 7:
	//		return std::pow(delta.x, 2) * delta.y;
	//		break;
	//	case 8:
	//		return delta.x * std::pow(delta.y, 2);
	//		break;
	//	case 9:
	//		return std::pow(delta.y, 3);
	//		break;
	//	default:
	//		std::cout << "	error: moment index should from 1-9!" << std::endl;
	//		exit(1);
	//	}
	//}

	//rO=1
	bool getMoment(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> & A
		){
		ScalarCfv::point delta = p - baryCenter;
		delta.x = delta.x / scale.x;
		delta.y = delta.y / scale.y;

		A[1] = delta.x;
		A[2] = delta.y;
	
		return true;
	}

	//rO=2
	bool getMoment(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> & A
		){
	
		ScalarCfv::point delta = p - baryCenter;
		delta.x = delta.x / scale.x;
		delta.y = delta.y / scale.y;

		A[1] = delta.x;
		A[2] = delta.y;
	
		A[3] = std::pow(delta.x, 2);
		A[4] = delta.x*delta.y;
		A[5] = std::pow(delta.y, 2);
		return true;
	}

	//rO=3
	bool getMoment(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> & A
		){
	
		ScalarCfv::point delta = p - baryCenter;
		delta.x = delta.x / scale.x;
		delta.y = delta.y / scale.y;

		A[1] = delta.x;
		A[2] = delta.y;

		A[3] = std::pow(delta.x, 2);
		A[4] = delta.x*delta.y;
		A[5] = std::pow(delta.y, 2);

		A[6] = std::pow(delta.x, 3);
		A[7] = std::pow(delta.x, 2)*delta.y;
		A[8] = delta.x*std::pow(delta.y, 2);
		A[9] = std::pow(delta.y, 3);
		return true;
	}

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> & A,
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> & IA){

		int mrow = 3 - 1;
		int ncol = 3 - 1;
		MKL_INT m = mrow, n = ncol, lda = ncol, ldu = mrow, ldvt = ncol, info;
		ScalarCfv::real temp_eps = 1e-16, alpha = 1.0, beta = 0.0;
		ScalarCfv::real* a = new ScalarCfv::real[m*n];
		ScalarCfv::real* inva = new ScalarCfv::real[n*m];
		ScalarCfv::real* temp = new ScalarCfv::real[n*m];
		MKL_INT min = (m < n) ? m : n;
		ScalarCfv::real* superb = new ScalarCfv::real[min - 1];
		ScalarCfv::real* s = new ScalarCfv::real[n];
		ScalarCfv::real* splus = new ScalarCfv::real[n*m];
		ScalarCfv::real* u = new ScalarCfv::real[ldu*m];
		ScalarCfv::real* vt = new ScalarCfv::real[ldvt*n];

		for (int ii = 0; ii < mrow; ++ii){
			for (int jj = 0; jj < ncol; ++jj){
				a[ii*n + jj] = A[ii + 1][jj + 1];
				//a[ii*n + jj] = A[ii][jj];
			}
		}

		//SVD
		info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m, n, a, lda, s, u, ldu, vt, ldvt, superb);
		if (info > 0) {
			std::cout << "The algorithm computing SVD failed to converge." << std::endl;
			exit(1);
		}

		for (int ii = 0; ii < ncol; ++ii){
			for (int jj = 0; jj < mrow; ++jj){
				splus[ii*mrow + jj] = 0.0;
				if (jj == ii) {
					if (abs(s[ii])>temp_eps){
						splus[ii*mrow + jj] = 1.0 / s[ii];
					}
				}
			}
		}
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, n, alpha, vt, n, splus, m, beta, temp, m);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, m, alpha, temp, m, u, m, beta, inva, m);

		for (int ii = 0; ii < ncol; ++ii){
			for (int jj = 0; jj < mrow; ++jj){
				IA[ii + 1][jj + 1] = inva[ii*mrow + jj];
				//IA[ii][jj] = inva[ii*mrow + jj];
			}
		}

		delete[] a;
		delete[] inva;
		delete[] temp;
		delete[] superb;
		delete[] s;
		delete[] splus;
		delete[] u;
		delete[] vt;
		a = NULL;
		inva = NULL;
		temp = NULL;
		superb = NULL;
		s = NULL;
		splus = NULL;
		u = NULL;
		vt = NULL;

		return true;
	}

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> & A,
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> & IA){

		int mrow = 6 - 1;
		int ncol = 6 - 1;
		MKL_INT m = mrow, n = ncol, lda = ncol, ldu = mrow, ldvt = ncol, info;
		ScalarCfv::real temp_eps = 1e-16, alpha = 1.0, beta = 0.0;
		ScalarCfv::real* a = new ScalarCfv::real[m*n];
		ScalarCfv::real* inva = new ScalarCfv::real[n*m];
		ScalarCfv::real* temp = new ScalarCfv::real[n*m];
		MKL_INT min = (m < n) ? m : n;
		ScalarCfv::real* superb = new ScalarCfv::real[min - 1];
		ScalarCfv::real* s = new ScalarCfv::real[n];
		ScalarCfv::real* splus = new ScalarCfv::real[n*m];
		ScalarCfv::real* u = new ScalarCfv::real[ldu*m];
		ScalarCfv::real* vt = new ScalarCfv::real[ldvt*n];

		for (int ii = 0; ii < mrow; ++ii){
			for (int jj = 0; jj < ncol; ++jj){
				a[ii*n + jj] = A[ii + 1][jj + 1];
			}
		}

		//SVD
		info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m, n, a, lda, s, u, ldu, vt, ldvt, superb);
		if (info > 0) {
			std::cout << "The algorithm computing SVD failed to converge." << std::endl;
			exit(1);
		}

		for (int ii = 0; ii < ncol; ++ii){
			for (int jj = 0; jj < mrow; ++jj){
				splus[ii*mrow + jj] = 0.0;
				if (jj == ii) {
					if (abs(s[ii])>temp_eps){
						splus[ii*mrow + jj] = 1.0 / s[ii];
					}
				}
			}
		}
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, n, alpha, vt, n, splus, m, beta, temp, m);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, m, alpha, temp, m, u, m, beta, inva, m);

		for (int ii = 0; ii < ncol; ++ii){
			for (int jj = 0; jj < mrow; ++jj){
				IA[ii + 1][jj + 1] = inva[ii*mrow + jj];
			}
		}

		delete[] a;
		delete[] inva;
		delete[] temp;
		delete[] superb;
		delete[] s;
		delete[] splus;
		delete[] u;
		delete[] vt;
		a = NULL;
		inva = NULL;
		temp = NULL;
		superb = NULL;
		s = NULL;
		splus = NULL;
		u = NULL;
		vt = NULL;

		return true;
	}

	bool getMatrixGeneralInverse(
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> & A,
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> & IA){

		int mrow = 10 - 1;
		int ncol = 10 - 1;
		MKL_INT m = mrow, n = ncol, lda = ncol, ldu = mrow, ldvt = ncol, info;
		ScalarCfv::real temp_eps = 1e-12, alpha = 1.0, beta = 0.0;
		ScalarCfv::real* a = new ScalarCfv::real[m*n];
		ScalarCfv::real* inva = new ScalarCfv::real[n*m];
		ScalarCfv::real* temp = new ScalarCfv::real[n*m];
		MKL_INT min = (m < n) ? m : n;
		ScalarCfv::real* superb = new ScalarCfv::real[min - 1];
		ScalarCfv::real* s = new ScalarCfv::real[n];
		ScalarCfv::real* splus = new ScalarCfv::real[n*m];
		ScalarCfv::real* u = new ScalarCfv::real[ldu*m];
		ScalarCfv::real* vt = new ScalarCfv::real[ldvt*n];

		for (int ii = 0; ii < mrow; ++ii){
			for (int jj = 0; jj < ncol; ++jj){
				a[ii*n + jj] = A[ii + 1][jj + 1];
			}
		}

		//SVD
		info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m, n, a, lda, s, u, ldu, vt, ldvt, superb);
		if (info > 0) {
			std::cout << "The algorithm computing SVD failed to converge." << std::endl;
			exit(1);
		}

		for (int ii = 0; ii < ncol; ++ii){
			for (int jj = 0; jj < mrow; ++jj){
				splus[ii*mrow + jj] = 0.0;
				if (jj == ii) {
					if (abs(s[ii])>temp_eps){
						splus[ii*mrow + jj] = 1.0 / s[ii];
					}
				}
			}
		}
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, n, alpha, vt, n, splus, m, beta, temp, m);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, m, alpha, temp, m, u, m, beta, inva, m);

		for (int ii = 0; ii < ncol; ++ii){
			for (int jj = 0; jj < mrow; ++jj){
				IA[ii + 1][jj + 1] = inva[ii*mrow + jj];
			}
		}

		delete[] a;
		delete[] inva;
		delete[] temp;
		delete[] superb;
		delete[] s;
		delete[] splus;
		delete[] u;
		delete[] vt;
		a = NULL;
		inva = NULL;
		temp = NULL;
		superb = NULL;
		s = NULL;
		splus = NULL;
		u = NULL;
		vt = NULL;

		return true;
	}

	ScalarCfv::real harmonic(ScalarCfv::real uL, ScalarCfv::real uR){
	
		//return ((std::pow(uL, 2) + std::pow(uR, 2)) / (uL + uR));
		return ((std::pow(uL, 2) + std::pow(uR, 2)) / std::pow(uL + uR, 2));
	}

	int getCombination(int  m, int  n) { //C^n_(m+n)
		int numerator, denominator1, denominator2;
		if (m == 0 || n == 0) { return 1; }
		numerator = 1;
		denominator1 = 1;
		denominator2 = 1;
		for (int idx = 1; idx <= m + n; ++idx){
			numerator = numerator * idx;
		}
		for (int idx = 1; idx <= m; ++idx){
			denominator1 = denominator1 * idx;
		}
		for (int idx = 1; idx <= n; ++idx){
			denominator2 = denominator2 * idx;
		}
		return numerator / denominator1 / denominator2;
	}

	int getFactorial(int  m, int  n) {//(m+n)!
		int numerator;
		if (m == 0 && n == 0) { return 1; }
		numerator = 1;
		for (int idx = 1; idx <= m + n; ++idx){
			numerator = numerator * idx;
		}
		return numerator;
	}

	//rO=1
	bool getDiffBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> & moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> & A
		){
		//	|0		 & 0		 & 0		|
		//	|dphi_1	 & dphi_1/dx & dphi_1/dy|
		//	|dphi_2	 & dphi_2/dx & dphi_2/dy|
		ScalarCfv::point delta = p - baryCenter;
		delta.x = delta.x / scale.x;
		delta.y = delta.y / scale.y;

		A[1][0] = delta.x - moment[1];	A[1][1] = 1.0 / scale.x;	A[1][2] = 0.0;
		A[2][0] = delta.y - moment[2];	A[2][1] = 0.0;				A[2][2] = 1.0 / scale.y;
	
		return true;
	}

	//rO=2
	bool getDiffBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> & moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> & A
		){
		//	|0		& 0			& 0			& 0				& 0				& 0			   |
		//	|dphi_1 & dphi_1/dx & dphi_1/dy & dphi^2_1/dx^2 & dphi^2_1/dxdy & dphi^2_1/dy^2|
		//	|dphi_2 & dphi_2/dx & dphi_2/dy & dphi^2_2/dx^2 & dphi^2_2/dxdy & dphi^2_2/dy^2|
		//	|dphi_3 & dphi_3/dx & dphi_3/dy & dphi^2_3/dx^2 & dphi^2_3/dxdy & dphi^2_3/dy^2|
		//	|dphi_4 & dphi_4/dx & dphi_4/dy & dphi^2_4/dx^2 & dphi^2_4/dxdy & dphi^2_4/dy^2|
		//	|dphi_5 & dphi_5/dx & dphi_5/dy & dphi^2_5/dx^2 & dphi^2_5/dxdy & dphi^2_5/dy^2|
	
		ScalarCfv::point delta = p - baryCenter;
		delta.x = delta.x / scale.x;
		delta.y = delta.y / scale.y;

		A[1][0] = delta.x - moment[1]; A[1][1] = 1.0 / scale.x;		A[1][2] = 0.0;				A[1][3] = 0.0;	A[1][4] = 0.0;	A[1][5] = 0.0;
		A[2][0] = delta.y - moment[2]; A[2][1] = 0.0;				A[2][2] = 1.0 / scale.y;	A[2][3] = 0.0;	A[2][4] = 0.0;	A[2][5] = 0.0;

		A[3][0] = std::pow(delta.x, 2) - moment[3];		A[3][1] = 2.0 * delta.x / scale.x;	A[3][2] = 0.0;						A[3][3] = 2.0 / scale.x / scale.x;	A[3][4] = 0.0;						A[3][5] = 0.0;
		A[4][0] = delta.x*delta.y - moment[4];			A[4][1] = delta.y / scale.x;		A[4][2] = delta.x / scale.y;		A[4][3] = 0.0;						A[4][4] = 1.0 / scale.x / scale.y;	A[4][5] = 0.0;
		A[5][0] = std::pow(delta.y, 2) - moment[5];		A[5][1] = 0.0;						A[5][2] = 2.0 * delta.y / scale.y;	A[5][3] = 0.0;						A[5][4] = 0.0;						A[5][5] = 2.0 / scale.y / scale.y;

		return true;
	}
	//rO=3
	bool getDiffBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> & moment,
		ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> & A
		){
		//	|0		& 0			& 0			& 0				& 0				& 0				& 0				& 0				& 0					& 0			   |
		//	|dphi_1 & dphi_1/dx & dphi_1/dy & dphi^2_1/dx^2 & dphi^2_1/dxdy & dphi^2_1/dy^2 & dphi^3_1/dx^3 & dphi^3_1/dx^2dy & dphi^3_1/dxdy^2 & dphi^3_1/dy^3|
		//	|dphi_2 & dphi_2/dx & dphi_2/dy & dphi^2_2/dx^2 & dphi^2_2/dxdy & dphi^2_2/dy^2 & dphi^3_2/dx^3 & dphi^3_2/dx^2dy & dphi^3_2/dxdy^2 & dphi^3_2/dy^3|
		//	|dphi_3 & dphi_3/dx & dphi_3/dy & dphi^2_3/dx^2 & dphi^2_3/dxdy & dphi^2_3/dy^2 & dphi^3_3/dx^3 & dphi^3_3/dx^2dy & dphi^3_3/dxdy^2 & dphi^3_3/dy^3|
		//	|dphi_4 & dphi_4/dx & dphi_4/dy & dphi^2_4/dx^2 & dphi^2_4/dxdy & dphi^2_4/dy^2 & dphi^3_4/dx^3 & dphi^3_4/dx^2dy & dphi^3_4/dxdy^2 & dphi^3_4/dy^3|
		//	|dphi_5 & dphi_5/dx & dphi_5/dy & dphi^2_5/dx^2 & dphi^2_5/dxdy & dphi^2_5/dy^2 & dphi^3_5/dx^3 & dphi^3_5/dx^2dy & dphi^3_5/dxdy^2 & dphi^3_5/dy^3|
		//	|dphi_6 & dphi_6/dx & dphi_6/dy & dphi^2_6/dx^2 & dphi^2_6/dxdy & dphi^2_6/dy^2 & dphi^3_6/dx^3 & dphi^3_6/dx^2dy & dphi^3_6/dxdy^2 & dphi^3_6/dy^3|
		//	|dphi_7 & dphi_7/dx & dphi_7/dy & dphi^2_7/dx^2 & dphi^2_7/dxdy & dphi^2_7/dy^2 & dphi^3_7/dx^3 & dphi^3_7/dx^2dy & dphi^3_7/dxdy^2 & dphi^3_7/dy^3|
		//	|dphi_8 & dphi_8/dx & dphi_8/dy & dphi^2_8/dx^2 & dphi^2_8/dxdy & dphi^2_8/dy^2 & dphi^3_8/dx^3 & dphi^3_8/dx^2dy & dphi^3_8/dxdy^2 & dphi^3_8/dy^3|
		//	|dphi_9 & dphi_9/dx & dphi_9/dy & dphi^2_9/dx^2 & dphi^2_9/dxdy & dphi^2_9/dy^2 & dphi^3_9/dx^3 & dphi^3_9/dx^2dy & dphi^3_9/dxdy^2 & dphi^3_9/dy^3|
	
		ScalarCfv::point delta = p - baryCenter;
		delta.x = delta.x / scale.x;
		delta.y = delta.y / scale.y;

		A[1][0] = delta.x - moment[1]; A[1][1] = 1.0 / scale.x;		A[1][2] = 0.0;				A[1][3] = 0.0;	A[1][4] = 0.0;	A[1][5] = 0.0;	A[1][6] = 0.0;	A[1][7] = 0.0;	A[1][8] = 0.0;	A[1][9] = 0.0;
		A[2][0] = delta.y - moment[2]; A[2][1] = 0.0;				A[2][2] = 1.0 / scale.y;	A[2][3] = 0.0;	A[2][4] = 0.0;	A[2][5] = 0.0;	A[2][6] = 0.0;	A[2][7] = 0.0;	A[2][8] = 0.0;	A[2][9] = 0.0;

		A[3][0] = std::pow(delta.x, 2) - moment[3];	A[3][1] = 2.0 * delta.x / scale.x;	A[3][2] = 0.0;						A[3][3] = 2.0 / scale.x / scale.x;	A[3][4] = 0.0;						A[3][5] = 0.0;						A[3][6] = 0.0;	A[3][7] = 0.0;	A[3][8] = 0.0;	A[3][9] = 0.0;
		A[4][0] = delta.x*delta.y - moment[4];		A[4][1] = delta.y / scale.x;		A[4][2] = delta.x / scale.y;		A[4][3] = 0.0;						A[4][4] = 1.0 / scale.x / scale.y;	A[4][5] = 0.0;						A[4][6] = 0.0;	A[4][7] = 0.0;	A[4][8] = 0.0;	A[4][9] = 0.0;
		A[5][0] = std::pow(delta.y, 2) - moment[5];	A[5][1] = 0.0;						A[5][2] = 2.0 * delta.y / scale.y;	A[5][3] = 0.0;						A[5][4] = 0.0;						A[5][5] = 2.0 / scale.y / scale.y;	A[5][6] = 0.0;	A[5][7] = 0.0;	A[5][8] = 0.0;	A[5][9] = 0.0;

		A[6][0] = std::pow(delta.x, 3) - moment[6];	 A[6][1] = 3.0 * std::pow(delta.x, 2) / scale.x;	A[6][2] = 0.0;	A[6][3] = 6.0 * delta.x / scale.x / scale.x;
		A[6][4] = 0.0;	A[6][5] = 0.0;	A[6][6] = 6.0 / scale.x / scale.x / scale.x;	A[6][7] = 0.0;	A[6][8] = 0.0;	A[6][9] = 0.0;

		A[7][0] = std::pow(delta.x, 2)*delta.y - moment[7]; A[7][1] = 2.0 * delta.x * delta.y / scale.x;	A[7][2] = std::pow(delta.x, 2) / scale.y;		A[7][3] = 2.0 * delta.y / scale.x / scale.x;
		A[7][4] = 2.0 * delta.x / scale.x / scale.y;		A[7][5] = 0.0;									A[7][6] = 0.0;									A[7][7] = 2.0 / scale.x / scale.x / scale.y;	A[7][8] = 0.0;									A[7][9] = 0.0;

		A[8][0] = std::pow(delta.y, 2)*delta.x - moment[8]; A[8][1] = std::pow(delta.y, 2) / scale.x;		A[8][2] = 2.0 * delta.x * delta.y / scale.y;	A[8][3] = 0.0;
		A[8][4] = 2.0 * delta.y / scale.x / scale.y;		A[8][5] = 2.0 * delta.x / scale.y / scale.y;	A[8][6] = 0.0;									A[8][7] = 0.0;									A[8][8] = 2.0 / scale.x / scale.y / scale.y;	A[8][9] = 0.0;

		A[9][0] = std::pow(delta.y, 3) - moment[9];		A[9][1] = 0.0;	A[9][2] = 3.0 * std::pow(delta.y, 2) / scale.y;	A[9][3] = 0.0;
		A[9][4] = 0.0;	A[9][5] = 6.0 * delta.y / scale.y / scale.y;	A[9][6] = 0.0;	A[9][7] = 0.0;	A[9][8] = 0.0;	A[9][9] = 6.0 / scale.y / scale.y / scale.y;

		return true;
	}

	//rO=1
	bool getBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> & moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 3> & A
		){
		ScalarCfv::point delta = p - baryCenter;
		delta.x = delta.x / scale.x;
		delta.y = delta.y / scale.y;
		A[1] = delta.x - moment[1];
		A[2] = delta.y - moment[2];
		return true;
	}
	//rO=2
	bool getBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> & moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 6> & A
		){
		ScalarCfv::point delta = p - baryCenter;
		delta.x = delta.x / scale.x;
		delta.y = delta.y / scale.y;

		A[1] = delta.x - moment[1];
		A[2] = delta.y - moment[2];

		A[3] = std::pow(delta.x, 2) - moment[3];
		A[4] = delta.x * delta.y - moment[4];
		A[5] = std::pow(delta.y, 2) - moment[5];
		return true;
	}
	//rO=3
	bool getBaseValue(
		ScalarCfv::point p,
		ScalarCfv::point baryCenter,
		ScalarCfv::point scale,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> & moment,
		ScalarCfv::tensor1D<ScalarCfv::real, 10> & A
		){
		ScalarCfv::point delta = p - baryCenter;
		delta.x = delta.x / scale.x;
		delta.y = delta.y / scale.y;

		A[1] = delta.x - moment[1];
		A[2] = delta.y - moment[2];

		A[3] = std::pow(delta.x, 2) - moment[3];
		A[4] = delta.x * delta.y - moment[4];
		A[5] = std::pow(delta.y, 2) - moment[5];

		A[6] = std::pow(delta.x, 3) - moment[6];
		A[7] = std::pow(delta.x, 2) * delta.y - moment[7];
		A[8] = delta.x * std::pow(delta.y, 2) - moment[8];
		A[9] = std::pow(delta.y, 3) - moment[9];
		return true;
	}

	ScalarCfv::real W12(ScalarCfv::real u[], const int J){
		//ScalarCfv::real theta[n];

		ScalarCfv::real* theta = new ScalarCfv::real[J + 1];
		theta[0] = 1.0;
		for (int ii = 1; ii <= J; ++ii){
			theta[ii] = (u[ii] + 1e-12) / (u[0] + 1e-12);
		}
		ScalarCfv::real n = 10.0;// 10.0;
		ScalarCfv::real p = 4.0;//4.0
		ScalarCfv::real sumLocal1 = n;
		ScalarCfv::real sumLocal2 = n;
		for (int ii = 1; ii <= J; ++ii){
			sumLocal1 += std::pow(1.0 / theta[ii], (p - 1.0));
			sumLocal2 += std::pow(1.0 / theta[ii], p);
		}
	
		delete[] theta;
		theta = NULL;

		return u[0] * sumLocal1 / (sumLocal2 + 1e-12);
	}

}