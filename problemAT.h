#ifndef _PROBLEM_AT_H
#define _PROBLEM_AT_H
#include "FieldSolver.h"

namespace ScalarCfv
{
	class problemAT : public fieldsolver{

	public:
		problemAT(
			parameter* parameter,
			cellFieldDataVector* cellFieldData,
			faceFieldDataVector* faceFieldData,
			cellGaussDataVector* cellGaussData,
			faceGaussDataVector* faceGaussData,
			evolution<rO>* evolution,
			reconstruction<rO>* reconstruction
			) : fieldsolver(
			parameter,
			cellFieldData,
			faceFieldData,
			cellGaussData,
			faceGaussData,
			evolution,
			reconstruction
			){};
		~problemAT() {};

	public:


		virtual bool timeMarchingExplicitSSPRK(
			vertexVector* node,
			GaussIntegralCellO1Grid<vO>* gaussIntegralCell,
			GaussIntegralFaceO1Grid<fO>* gaussIntegralFace,
			const std::string & filenameInputBackup,
			const std::string & filenameOutputBackup,
			const std::string & filenameOutputSln,
			const std::string & filenameOutputResidual,
			const std::string & filenameOutputError);

		virtual bool getGRE(GaussIntegralCellO1Grid<vO>* gaussIntegralCell, real & L1, real & L2, real & L8);
		virtual bool exportGRE(
			int & flag,
			const int step,
			real L1, real L2, real L8,
			const std::string & filename);

	};
}

#endif