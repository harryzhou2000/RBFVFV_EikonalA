#ifndef _FIELDSOLVER_H
#define _FIELDSOLVER_H
//#include "TypeDefine.h"
//#include "Point.h"
//#include "Parameter.h"
//#include "Geometry.h"
//#include "GaussIntegral.h"
//#include "Variable.h"
#include "Reconstruction.h"
#include "Evolution.h"

namespace ScalarCfv
{
	class fieldsolver
	{
	public:
		fieldsolver(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			faceFieldDataVector *faceFieldData,
			cellGaussDataVector *cellGaussData,
			faceGaussDataVector *faceGaussData,
			evolution<rO> *evolution,
			reconstruction<rO> *reconstruction
		)
		{
			this->parameter_ = parameter;
			this->cellFieldData_ = cellFieldData;
			this->faceFieldData_ = faceFieldData;
			this->cellGaussData_ = cellGaussData;
			this->faceGaussData_ = faceGaussData;
			this->evolution_ = evolution;
			this->reconstruction_ = reconstruction;
		};
		~fieldsolver(){};

		parameter *parameter_;
		cellFieldDataVector *cellFieldData_;
		faceFieldDataVector *faceFieldData_;
		cellGaussDataVector *cellGaussData_;
		faceGaussDataVector *faceGaussData_;
		evolution<rO> *evolution_;
		reconstruction<rO> *reconstruction_;
		const static unsigned int NDOFS = (rO + 2) * (rO + 1) / 2;

	public:
		virtual bool allocateArray();
		virtual bool initializeArray();

		virtual bool timeMarchingExplicitSSPRK(
			vertexVector *node,
			GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
			GaussIntegralFaceO1Grid<fO> *gaussIntegralFace,
			const std::string &filenameInputBackup,
			const std::string &filenameOutputBackup,
			const std::string &filenameOutputSln,
			const std::string &filenameOutputResidual);

		virtual bool timeMarchingExplicit_DEBUG(
			vertexVector *node,
			GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
			GaussIntegralFaceO1Grid<fO> *gaussIntegralFace,
			const std::string &filenameInputBackup,
			const std::string &filenameOutputBackup,
			const std::string &filenameOutputSln,
			const std::string &filenameOutputResidual);

		virtual bool timeMarchingExplicitSSPRK(
			vertexVector *node,
			GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
			GaussIntegralFaceO2Grid<fO> *gaussIntegralFace,
			const std::string &filenameInputBackup,
			const std::string &filenameOutputBackup,
			const std::string &filenameOutputSln,
			const std::string &filenameOutputResidual);

		virtual bool timeMarchingExplicit_DEBUG(
			vertexVector *node,
			GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
			GaussIntegralFaceO2Grid<fO> *gaussIntegralFace,
			const std::string &filenameInputBackup,
			const std::string &filenameOutputBackup,
			const std::string &filenameOutputSln,
			const std::string &filenameOutputResidual);

		virtual bool readBackUp(
			const std::string &filename);

		virtual bool exportBackUp(
			const std::string &filename);

		virtual bool exportBackUp(
			const int ii,
			const std::string &filename);

		virtual bool exportCurrentSln(
			vertexVector *node,
			const int step,
			const std::string &filename);

		virtual bool getResidual(
			real &residual);

		virtual bool exportResidual(
			int &flag,
			const int step,
			real residual,
			const std::string &filename);
		//-------------------------------------------
		//get reconstruction
		//get evolution
		//post process
	};

	//eg.
	class circleproblem : public fieldsolver
	{

	public:
		circleproblem(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			faceFieldDataVector *faceFieldData,
			cellGaussDataVector *cellGaussData,
			faceGaussDataVector *faceGaussData,
			evolution<rO> *evolution,
			reconstruction<rO> *reconstruction) : fieldsolver(parameter,
															  cellFieldData,
															  faceFieldData,
															  cellGaussData,
															  faceGaussData,
															  evolution,
															  reconstruction){};
		~circleproblem(){};

	public:
		//bool getError();
	};

}
#endif