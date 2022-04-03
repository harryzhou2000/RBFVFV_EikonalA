#ifndef _TYPE_DEFINE_H
#define _TYPE_DEFINE_H
#include <vector>
#include <assert.h>
#include <omp.h>

// #define USE_RBFA1
#define USE_RBFB1

#if defined USE_RBFB1
#define RBFB1_NO_GS

// #define RBFB1_CR_DIFFBOUND
// #define RBFB1_CR_DIFFDIM 0.1
// #define RBFB1_INCREMENT_CR
// #define USE_RBFB1_N
#define RBFB1_GlobalPoly
// #define RBFB1_GlobalPoly_ESC // for vfv functional research

// #define RBFB1_USE_UNITARY_MEANIJ
// #define RBFB1_USE_DELTA_Jacobi_INTERFACEJ
// #define RBFB1_USE_DELTA_NT_INTERFACEJ
#define RBFB1_USE_DELTA_INTERFACEJ	

#define RBFB1_RBF_USE_MEANIJ
#define RBFB1_POLY_USE_MEANIJ

#define RBFB1_AR0 4

// #define RBFB1_CR_INTERPOLATE

#define GRID_ROT 0

#define RBFB1GetMoment CfvMath::getMomentRBFB1_POLY
#define RBFB1GetBaseValue CfvMath::getBaseValueRBFB1_POLY
#define RBFB1GetDiffBaseValue CfvMath::getDiffBaseValueRBFB1_POLY

// #define RBFB1GetMoment CfvMath::getMomentRBFB1_7_6
// #define RBFB1GetBaseValue CfvMath::getBaseValueRBFB1_7_6
// #define RBFB1GetDiffBaseValue CfvMath::getDiffBaseValueRBFB1_7_6

// #define RBFB1GetMoment CfvMath::getMomentRBFB1_4_6
// #define RBFB1GetBaseValue CfvMath::getBaseValueRBFB1_4_6
// #define RBFB1GetDiffBaseValue CfvMath::getDiffBaseValueRBFB1_4_6

#define RBFB1GetMomentCR CfvMath::getMomentRBFB1_POLY
#define RBFB1GetBaseValueCR CfvMath::getBaseValueRBFB1_POLY
#define RBFB1GetDiffBaseValueCR CfvMath::getDiffBaseValueRBFB1_POLY

// #define RBFB1GetMomentCR CfvMath::getMomentRBFB1_Interp_15_3
// #define RBFB1GetBaseValueCR CfvMath::getBaseValueRBFB1_Interp_15_3
// #define RBFB1GetDiffBaseValueCR CfvMath::getDiffBaseValueRBFB1_Interp_15_3

// #define RBFB1GetMomentCR CfvMath::getMomentRBFB1_7_6
// #define RBFB1GetBaseValueCR CfvMath::getBaseValueRBFB1_7_6
// #define RBFB1GetDiffBaseValueCR CfvMath::getDiffBaseValueRBFB1_7_6

// #define RBFB1GetMomentCR CfvMath::getMomentRBFB1_4_6
// #define RBFB1GetBaseValueCR CfvMath::getBaseValueRBFB1_4_6
// #define RBFB1GetDiffBaseValueCR CfvMath::getDiffBaseValueRBFB1_4_6

/***
with VF=POLY2
unstable??
*/
// #define RBFB1GetMomentCR CfvMath::getMomentRBFB1_7_3
// #define RBFB1GetBaseValueCR CfvMath::getBaseValueRBFB1_7_3
// #define RBFB1GetDiffBaseValueCR CfvMath::getDiffBaseValueRBFB1_7_3

/***
with VF=POLY2
PHSpline c=1 works, but strange
*/
// #define RBFB1GetMomentCR CfvMath::getMomentRBFB1_5_3
// #define RBFB1GetBaseValueCR CfvMath::getBaseValueRBFB1_5_3
// #define RBFB1GetDiffBaseValueCR CfvMath::getDiffBaseValueRBFB1_5_3

/***
with VF=POLY2
MQ c=0.3 works
*/
// #define RBFB1GetMomentCR CfvMath::getMomentRBFB1_4_3
// #define RBFB1GetBaseValueCR CfvMath::getBaseValueRBFB1_4_3
// #define RBFB1GetDiffBaseValueCR CfvMath::getDiffBaseValueRBFB1_4_3


#define RBFB1_CRBF  0.3
#endif

constexpr int GLOBAL_NDOFS(int O) { return (O + 2) * (O + 1) / 2; }
constexpr int GLOBAL_NDIFFS(int O) { return (O + 2) * (O + 1) / 2; }

// constexpr int GLOBAL_NDOFS(int O) { return 6; }
// constexpr int GLOBAL_NDIFFS(int O) { return 6; }

// constexpr int GLOBAL_NDOFSCR(int O) { return 4; }
// constexpr int GLOBAL_NDIFFSCR(int O) { return 3; }

// constexpr int GLOBAL_NDOFSCR(int O) { return 15; }
// constexpr int GLOBAL_NDIFFSCR(int O) { return 3; }

constexpr int GLOBAL_NDOFSCR(int O) { return (O + 2) * (O + 1) / 2; }
constexpr int GLOBAL_NDIFFSCR(int O) { return (O + 2) * (O + 1) / 2; }
// #define TRIAL
// #define DEBUG_RBFB
// #define USE_RBF

#define rO 3 // reconstruction polynomal order {1,2,3}
#define fO 5 // face integral order {1,3,5}
#define vO 3 // volume integral order {1,2,3}
#define mO 1 // mesh order {1,2}

//-----------------------
#define O1mesh
//#define O2mesh
//-----------------------
#define orientFlagG -1
#define orientFlagF 1
//-----------------------
// #define IF_RESTRICT_RADIUS

namespace ScalarCfv
{
	enum boundaryType
	{
		Wall = -1,
		FarField = -2,
		Symmetric = -3,
		OutFlow = -4,
		InFlow = -5,
		Periodic = -10
	};
	enum boundaryCellType
	{
		InnerCell,
		WallCell,
		FarFieldCell,
		SymmetricCell,
		OutFlowCell,
		InFlowCell,
		PeriodicCell,
		MixedCell
	};
	enum timeMarchType
	{
		EX_RK3 = 0,
		EX_LUSGS = 1,
		IM_LUSGS = 2
	};

	enum cellType
	{
		Triangle = 3,
		Quadrilateral = 4
	};
	typedef double real;
	const real UNINITReal = 1e100;
	typedef std::vector<real> realVector;
	typedef std::vector<realVector> realVectorVector;

	template <typename T>
	T getMax(T a[], int idx)
	{
		T max = a[0];
		for (int ii = 1; ii < idx; ++ii)
		{
			if (max < a[ii])
			{
				max = a[ii];
			}
		}
		return max;
	}

	template <typename T>
	T getMin(T a[], int idx)
	{
		T min = a[0];
		for (int ii = 1; ii < idx; ++ii)
		{
			if (min > a[ii])
			{
				min = a[ii];
			}
		}
		return min;
	}

	template <typename T>
	T getMax(T a, T b)
	{
		return a > b ? a : b;
	}

	template <typename T>
	T getMin(T a, T b)
	{
		return a < b ? a : b;
	}

}  
  
 
// Linking issues
// #define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACKE
#ifdef USE_OPEN_BLAS
#include<cblas.h>
#include<lapacke.h>
#define MKL_INT lapack_int

#else
#include "mkl_cblas.h"
#include "mkl_lapacke.h"

#endif

#endif