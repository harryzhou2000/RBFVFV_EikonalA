#ifndef _TYPE_DEFINE_H
#define _TYPE_DEFINE_H
#include <vector>
#include <assert.h>

// #define USE_RBFA1
#define TRIAL
// #define DEBUG_RBFB
// #define USE_RBF

#define rO 1 //reconstruction polynomal order {1,2,3}
#define fO 5 //face integral order {1,3,5}
#define vO 3 //volume integral order {1,2,3}
#define mO 1 //mesh order {1,2}
//-----------------------
#define O1mesh
//#define O2mesh
//-----------------------
#define orientFlagG -1
#define orientFlagF 1
//-----------------------
#define IF_RESTRICT_RADIUS
namespace ScalarCfv
{
	enum boundaryType{
		Wall = -1,
		FarField = -2,
		Symmetric = -3,
		OutFlow = -4,
		InFlow = -5,
		Periodic = -10
	};
	enum boundaryCellType{
		InnerCell,
		WallCell,
		FarFieldCell,
		SymmetricCell,
		OutFlowCell,
		InFlowCell,
		PeriodicCell,
		MixedCell
	};
	enum timeMarchType{
		EX_RK3 = 0,
		EX_LUSGS = 1,
		IM_LUSGS = 2
	};

	enum cellType{
		Triangle = 3,
		Quadrilateral = 4
	};
	typedef double real;
	typedef std::vector<real> realVector;
	typedef std::vector<realVector> realVectorVector;

	template<typename T>
	T getMax(T a[], int idx) {
		T max = a[0];
		for (int ii = 1; ii<idx; ++ii){
			if (max < a[ii]){
				max = a[ii];
			}
		}
		return max;
	}    

	template<typename T>
	T getMin(T a[], int idx) {
		T min = a[0];
		for (int ii = 1; ii<idx; ++ii){
			if (min > a[ii]){
				min = a[ii];
			}
		}
		return min;
	}

	template<typename T>
	T getMax(T  a, T  b) {
		return a > b ? a : b;
	}

	template<typename T>
	T getMin(T a, T b) {
		return a < b ? a : b;
	}
}
#endif