#ifndef _VARIABLE_H
#define _VARIABLE_H
//#include "TypeDefine.h"
#include "GaussIntegral.h"
//#include "Point.h"
//#include "Tensor.h"
#include <cmath>

namespace ScalarCfv
{
	// template <unsigned int O>
	// class variable
	//{
	// public:
	//	variable() {};
	//	~variable() {};
	//	const static unsigned int NDOFS = (O + 2)*(O + 1) / 2 - 1;
	// private:
	//	real gama;
	//	real Cp;
	//	real Cv;
	//	real relaxFactor;
	//	real artificialViscousParameter;
	// public:
	//	std::vector<real> dT;
	//	std::vector<real> smoothIndicator;
	//	std::vector<real> artificialViscous;
	//	std::vector<tensor1D<NDOFS+1>> scalarVariableTn;
	//	std::vector<tensor1D<NDOFS+1>> scalarVariableTm;
	//	std::vector<tensor1D<NDOFS+1>> scalarVariableTnLimited;
	//	std::vector<tensor1D<NDOFS+1>> scalarVariableTnCR;
	//	//reconstruction Ax = b
	//	std::vector<tensor2D<real, NDOFS, NDOFS>> reconstructionMatrixPlus;
	//	std::vector<tensor2D<real, NDOFS, NDOFS>> reconstructionMatrix;
	//	std::vector<tensor2D<real, NDOFS, NDOFS>> reconstructionMatrixPlusCR;
	//	std::vector<tensor2D<real, NDOFS, NDOFS>> reconstructionMatrixCR;
	//	//boundary ??
	// };

	class cellFieldData : public cellGaussData
	{
	public:
		cellFieldData() : cellGaussData(){};
		~cellFieldData(){};
		// const static unsigned int NDOFS = (O + 2)*(O + 1) / 2; //������ֵ��

	public:
		unsigned int NDOFS;
		unsigned int NDOFSCR;
		unsigned int PG;

		real gama;
		real Cp;
		real Cv;
		real relaxFactor;
		real relaxFactorCR;
		real artificialViscousParameter;

	public:
		real dTPhysical;
		real dTPseudo;
		real smoothIndicator;
		real artificialViscous;

		// 20200324
		real lambdaCell;
		real timeMarchingRHSTn;

		// 20200404
		real phiFilter;

		// for multistep method,
		//[1] <-> the results from Tm
		//[2] <-> the results from Tm + (1)
		//[3] <-> the results from Tm + (2)
		std::vector<real> timeMarchingRHSRK;

		std::vector<real> baseMoment;
		std::vector<real> baseMomentCR;

		std::vector<real> scalarVariableTn; // current time step , inlcude the multistep time marching methods, Tn means the newest values
		std::vector<real> scalarVariableTm; // previous time step
		std::vector<real> scalarVariableTnLimited;
		std::vector<real> scalarVariableTnCR;
		std::vector<real> scalarVariableTnCRLimited;
		std::vector<real> deltaScalarVariableTn;
		std::vector<real> deltaScalarVariableTm;

		std::vector<std::vector<real>> vectorbij;
		std::vector<std::vector<real>> vectorAiiInversebij;
		std::vector<std::vector<real>> matrixAiiInverse;
		std::vector<std::vector<real>> matrixAii;
		std::vector<std::vector<std::vector<real>>> matrixBij;
		std::vector<std::vector<std::vector<real>>> matrixAiiInverseBij;
		std::vector<std::vector<real>> matrixAiiInverseCR;
		std::vector<std::vector<real>> matrixAiiCR;
		std::vector<std::vector<std::vector<real>>> matrixBijCR; // CR�ƺ�����Ҫ
		// std::vector<std::vector<real> > matrixAiInverseInverseCR;
		// std::vector<std::vector<real> > matrixAiCR;

		std::vector<std::pair<point, real>> parametricValue;

		std::vector<std::vector<std::vector<real>>> diffBaseValueData; // at gauss points
		std::vector<std::vector<std::vector<real>>> diffBaseValueDataCR; // at gauss points
		std::vector<std::vector<real>> diffBaseValueDataBary;
	};

	typedef std::vector<cellFieldData> cellFieldDataVector;

	class faceFieldData : public faceGaussData
	{
	public:
		faceFieldData() : faceGaussData(){};
		~faceFieldData(){};
		// unsigned int NDOFS = (O + 2)*(O + 1) / 2; //������ֵ��

	public:
		unsigned int NDOFS;
		unsigned int fPG;
		// tensor1D<NDOFS> flux;//??
		// tensor1D<NDOFS> scalarVariableTm; //previous time step
		real weightVF;
		std::vector<real> faceWeightVF;
		// std::vector<real> flux;
		real lambdaFace;

		std::vector<std::vector<std::vector<std::vector<real>>>> diffBaseValueData;// at gaussian points
		std::vector<std::vector<std::vector<std::vector<real>>>> diffBaseValueDataCR; // at gaussian points
		std::vector<std::vector<std::vector<real>>> diffBaseValueDataMid; 
	};
	typedef std::vector<faceFieldData> faceFieldDataVector;
}
#endif