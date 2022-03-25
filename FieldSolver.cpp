#include "FieldSolver.h"
#include <cmath>
// 20200220 ���������ļ����ڵ���-1����0�ſ�ʼ //���޸� 20200326
namespace ScalarCfv
{
	bool fieldsolver::allocateArray()
	{

		std::cout << "allocating cell variables array..." << std::endl;
		cellFieldData_->resize(parameter_->cellNumber);
		// cellVector::iterator iterCell;
		cellGaussDataVector::iterator iterCellGaussData;
		cellFieldDataVector::iterator iterCellFieldData;

		for (iterCellGaussData = cellGaussData_->begin(); iterCellGaussData != cellGaussData_->end(); ++iterCellGaussData)
		{
			int ii = (*iterCellGaussData).index;
			iterCellFieldData = cellFieldData_->begin() + ii - 1;
			(*iterCellFieldData).index = (*iterCellGaussData).index;
			(*iterCellFieldData).index = (*iterCellGaussData).index;
			(*iterCellFieldData).isBoundaryCell = (*iterCellGaussData).isBoundaryCell;
			(*iterCellFieldData).isCountByboundaryCellType = (*iterCellGaussData).isCountByboundaryCellType;
			(*iterCellFieldData).cellType_ = (*iterCellGaussData).cellType_;
			(*iterCellFieldData).boundaryCellType_ = (*iterCellGaussData).boundaryCellType_;
			(*iterCellFieldData).volume = (*iterCellGaussData).volume;
			(*iterCellFieldData).baryCenter = (*iterCellGaussData).baryCenter;
			(*iterCellFieldData).cellSideOff = (*iterCellGaussData).cellSideOff;
			(*iterCellFieldData).lengthReference = (*iterCellGaussData).lengthReference;
			(*iterCellFieldData).cellNodeNumber = (*iterCellGaussData).cellNodeNumber;
			(*iterCellFieldData).cellFaceNumber = (*iterCellGaussData).cellFaceNumber;
			(*iterCellFieldData).cellCellNumber = (*iterCellGaussData).cellCellNumber;

			(*iterCellFieldData).parametricVolume = (*iterCellGaussData).parametricVolume;
			(*iterCellFieldData).lambdaCell = 0.0;

			(*iterCellFieldData).o2PointNumber = (*iterCellGaussData).o2PointNumber;
			// std::cout << "--------" << std::endl;
			// std::cout << (*iterCellFieldData).volume << "\t" << (*iterCellGaussData).volume << std::endl;

			if ((*iterCellFieldData).cellType_ == Triangle)
			{
				(*iterCellFieldData).PG = (*iterCellGaussData).PGTri[vO];
				(*iterCellFieldData).cellNodeNumber = (*iterCellGaussData).cellNodeNumber;
				(*iterCellFieldData).cellNode.resize((*iterCellGaussData).cellNodeNumber + 1);
				(*iterCellFieldData).gaussPairVector_.resize((*iterCellGaussData).PGTri[vO]);
				(*iterCellFieldData).parametricValueTri.resize((*iterCellGaussData).PGTri[vO]);
				(*iterCellFieldData).parametricValue.resize((*iterCellGaussData).PGTri[vO]);
				for (int jj = 0; jj < (*iterCellGaussData).PGTri[vO]; ++jj)
				{
					(*iterCellFieldData).gaussPairVector_[jj].JacobiCof = (*iterCellGaussData).gaussPairVector_[jj].JacobiCof;
					(*iterCellFieldData).gaussPairVector_[jj].normalVector = (*iterCellGaussData).gaussPairVector_[jj].normalVector;
					(*iterCellFieldData).gaussPairVector_[jj].p = (*iterCellGaussData).gaussPairVector_[jj].p;

					(*iterCellFieldData).parametricValueTri[jj].first = (*iterCellGaussData).parametricValueTri[jj].first;
					(*iterCellFieldData).parametricValueTri[jj].second = (*iterCellGaussData).parametricValueTri[jj].second;

					(*iterCellFieldData).parametricValue[jj].first = (*iterCellGaussData).parametricValueTri[jj].first;
					(*iterCellFieldData).parametricValue[jj].second = (*iterCellGaussData).parametricValueTri[jj].second;
				}
			}
			else if ((*iterCellFieldData).cellType_ == Quadrilateral)
			{
				(*iterCellFieldData).PG = (*iterCellGaussData).PGQuad[vO];
				(*iterCellFieldData).cellNodeNumber = (*iterCellGaussData).cellNodeNumber;
				(*iterCellFieldData).cellNode.resize((*iterCellGaussData).cellNodeNumber + 1);
				(*iterCellFieldData).gaussPairVector_.resize((*iterCellGaussData).PGQuad[vO]);
				(*iterCellFieldData).parametricValueQuad.resize((*iterCellGaussData).PGQuad[vO]);
				(*iterCellFieldData).parametricValue.resize((*iterCellGaussData).PGQuad[vO]);
				for (int jj = 0; jj < (*iterCellGaussData).PGQuad[vO]; ++jj)
				{
					(*iterCellFieldData).gaussPairVector_[jj].JacobiCof = (*iterCellGaussData).gaussPairVector_[jj].JacobiCof;
					(*iterCellFieldData).gaussPairVector_[jj].normalVector = (*iterCellGaussData).gaussPairVector_[jj].normalVector;
					(*iterCellFieldData).gaussPairVector_[jj].p = (*iterCellGaussData).gaussPairVector_[jj].p;

					(*iterCellFieldData).parametricValueQuad[jj].first = (*iterCellGaussData).parametricValueQuad[jj].first;
					(*iterCellFieldData).parametricValueQuad[jj].second = (*iterCellGaussData).parametricValueQuad[jj].second;

					(*iterCellFieldData).parametricValue[jj].first = (*iterCellGaussData).parametricValueQuad[jj].first;
					(*iterCellFieldData).parametricValue[jj].second = (*iterCellGaussData).parametricValueQuad[jj].second;
				}
			}

			//	std::cout << "------------------------" << std::endl;
			//	std::cout << "PG: " << ii << std::endl;
			//	std::cout << (*iterCellFieldData).gaussPairVector_[1].p.x << "\t" << (*iterCellFieldData).gaussPairVector_[1].p.y
			//		<< "\t" << (*iterCellFieldData).gaussPairVector_[9].p.x << "\t" << (*iterCellFieldData).gaussPairVector_[9].p.y << std::endl;
			//	std::cout << (*iterCellFieldData).gaussPairVector_[1].JacobiCof << "\t" << (*iterCellFieldData).gaussPairVector_[9].JacobiCof << std::endl;
			//	std::cout << phi[ii][4] << "\t" << phi[ii][5] << "\t" << phi[ii][6] << std::endl;
			//	std::cout << phi[ii][7] << "\t" << phi[ii][8] << "\t" << phi[ii][9] << std::endl;

			//	system("pause");

			for (int jj = 1; jj < (*iterCellFieldData).cellNodeNumber + 1; ++jj)
			{
				(*iterCellFieldData).cellNode[jj].first = (*iterCellGaussData).cellNode[jj].first;
				(*iterCellFieldData).cellNode[jj].second = (*iterCellGaussData).cellNode[jj].second;
			}

			(*iterCellFieldData).cellFaceIndex.resize((*iterCellFieldData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterCellFieldData).cellFaceNumber + 1; ++jj)
			{
				(*iterCellFieldData).cellFaceIndex[jj] = (*iterCellGaussData).cellFaceIndex[jj];
				(*iterCellFieldData).cellFaceIndex[jj] = (*iterCellGaussData).cellFaceIndex[jj];
			}

			(*iterCellFieldData).cellCellIndex.resize((*iterCellFieldData).cellCellNumber + 1);
			for (int jj = 0; jj < (*iterCellFieldData).cellCellNumber + 1; ++jj)
			{
				(*iterCellFieldData).cellCellIndex[jj] = (*iterCellGaussData).cellCellIndex[jj];
				(*iterCellFieldData).cellCellIndex[jj] = (*iterCellGaussData).cellCellIndex[jj];
			}

			(*iterCellFieldData).cellFaceSideOff.resize((*iterCellFieldData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterCellFieldData).cellFaceNumber + 1; ++jj)
			{
				(*iterCellFieldData).cellFaceSideOff[jj] = (*iterCellGaussData).cellFaceSideOff[jj];
				(*iterCellFieldData).cellFaceSideOff[jj] = (*iterCellGaussData).cellFaceSideOff[jj];
			}

			(*iterCellFieldData).cellFaceWeight.resize((*iterCellFieldData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterCellFieldData).cellFaceNumber + 1; ++jj)
			{
				(*iterCellFieldData).cellFaceWeight[jj] = (*iterCellGaussData).cellFaceWeight[jj];
				(*iterCellFieldData).cellFaceWeight[jj] = (*iterCellGaussData).cellFaceWeight[jj];
			}

			// iterCellFieldData NDOFS should conform with that of recstruciton
			// (*iterCellFieldData).NDOFS = (rO + 2)*(rO + 1) / 2;
			(*iterCellFieldData).NDOFS = reconstruction_->ReturnNDOFS();
			(*iterCellFieldData).NDOFSCR = reconstruction_->ReturnNDOFSCR();
			// std::cout << (*iterCellFieldData).NDOFS << std::endl;
			// exit(1);

			(*iterCellFieldData).baseMoment.resize((*iterCellFieldData).NDOFS);
			(*iterCellFieldData).baseMomentCR.resize((*iterCellFieldData).NDOFSCR);

			(*iterCellFieldData).scalarVariableTn.resize((*iterCellFieldData).NDOFS);
			(*iterCellFieldData).scalarVariableTm.resize((*iterCellFieldData).NDOFS);
			(*iterCellFieldData).scalarVariableTnLimited.resize((*iterCellFieldData).NDOFS);
			(*iterCellFieldData).scalarVariableTnCR.resize((*iterCellFieldData).NDOFSCR);
			(*iterCellFieldData).deltaScalarVariableTn.resize((*iterCellFieldData).NDOFS);
			(*iterCellFieldData).deltaScalarVariableTm.resize((*iterCellFieldData).NDOFS);
			(*iterCellFieldData).scalarVariableTnCRLimited.resize((*iterCellFieldData).NDOFSCR);

			(*iterCellFieldData).matrixAiiInverse.resize((*iterCellFieldData).NDOFS);
			(*iterCellFieldData).matrixAii.resize((*iterCellFieldData).NDOFS);
			(*iterCellFieldData).matrixAiiInverseCR.resize((*iterCellFieldData).NDOFSCR);
			(*iterCellFieldData).matrixAiiCR.resize((*iterCellFieldData).NDOFSCR);

			for (int jj = 0; jj < static_cast<int>((*iterCellFieldData).NDOFS); ++jj)
			{
				(*iterCellFieldData).matrixAiiInverse[jj].resize((*iterCellFieldData).NDOFS);
				(*iterCellFieldData).matrixAii[jj].resize((*iterCellFieldData).NDOFS);
			}
			for (int jj = 0; jj < static_cast<int>((*iterCellFieldData).NDOFSCR); ++jj)
			{
				(*iterCellFieldData).matrixAiiInverseCR[jj].resize((*iterCellFieldData).NDOFSCR);
				(*iterCellFieldData).matrixAiiCR[jj].resize((*iterCellFieldData).NDOFSCR);
			}

			// cellFaceNumber->cellCellNumber
			(*iterCellFieldData).vectorbij.resize((*iterCellFieldData).cellCellNumber + 1);
			(*iterCellFieldData).vectorAiiInversebij.resize((*iterCellFieldData).cellCellNumber + 1);
			(*iterCellFieldData).matrixBij.resize((*iterCellFieldData).cellCellNumber + 1);
			(*iterCellFieldData).matrixAiiInverseBij.resize((*iterCellFieldData).cellCellNumber + 1);
			(*iterCellFieldData).matrixBijCR.resize((*iterCellFieldData).cellCellNumber + 1);

			for (int jj = 1; jj < (*iterCellFieldData).cellCellNumber + 1; ++jj)
			{
				(*iterCellFieldData).vectorbij[jj].resize((*iterCellFieldData).NDOFS);
				(*iterCellFieldData).vectorAiiInversebij[jj].resize((*iterCellFieldData).NDOFS);
				(*iterCellFieldData).matrixBij[jj].resize((*iterCellFieldData).NDOFS);
				(*iterCellFieldData).matrixAiiInverseBij[jj].resize((*iterCellFieldData).NDOFS);
				(*iterCellFieldData).matrixBijCR[jj].resize((*iterCellFieldData).NDOFSCR);
				for (int kk = 0; kk < static_cast<int>((*iterCellFieldData).NDOFS); ++kk)
				{
					(*iterCellFieldData).matrixBij[jj][kk].resize((*iterCellFieldData).NDOFS);
					(*iterCellFieldData).matrixAiiInverseBij[jj][kk].resize((*iterCellFieldData).NDOFS);
				}
				for (int kk = 0; kk < static_cast<int>((*iterCellFieldData).NDOFSCR); ++kk)
				{
					(*iterCellFieldData).matrixBijCR[jj][kk].resize((*iterCellFieldData).NDOFSCR);
				}
			}
			(*iterCellFieldData).timeMarchingRHSRK.resize(parameter_->nStepTimeMarching + 1); //ָ���빫ʽ��Ӧ��

			iterCellFieldData->diffBaseValueData.resize(iterCellFieldData->PG);
			for (auto &i : iterCellFieldData->diffBaseValueData)
			{
				i.resize(reconstruction_->ReturnNDOFS());
				for (auto &j : i)
					j.resize(reconstruction_->ReturnNDIFFS(), UNINITReal);
			}

			iterCellFieldData->diffBaseValueDataCR.resize(iterCellFieldData->PG);
			for (auto &i : iterCellFieldData->diffBaseValueDataCR)
			{
				i.resize(reconstruction_->ReturnNDOFSCR());
				for (auto &j : i)
					j.resize(reconstruction_->ReturnNDIFFSCR(), UNINITReal);
			}

			iterCellFieldData->diffBaseValueDataBary.resize(reconstruction_->ReturnNDOFS());
			for (auto &i : iterCellFieldData->diffBaseValueDataBary)
				i.resize(reconstruction_->ReturnNDIFFS(), UNINITReal);
		}

		std::cout << " ..cell variables array has been allocated." << std::endl;
		std::cout << "allocating face variables array..." << std::endl;
		faceFieldData_->resize(parameter_->faceNumber);
		// cellVector::iterator iterCell;
		faceGaussDataVector::iterator iterFaceGaussData;
		faceFieldDataVector::iterator iterFaceFieldData;

		for (iterFaceGaussData = faceGaussData_->begin(); iterFaceGaussData != faceGaussData_->end(); ++iterFaceGaussData)
		{
			int ii = (*iterFaceGaussData).index;
			iterFaceFieldData = faceFieldData_->begin() + ii - 1;
			(*iterFaceFieldData).index = (*iterFaceGaussData).index;
			(*iterFaceFieldData).faceProperty = (*iterFaceGaussData).faceProperty;
			(*iterFaceFieldData).faceNodeNumber = (*iterFaceGaussData).faceNodeNumber;
			(*iterFaceFieldData).faceCellNumber = (*iterFaceGaussData).faceCellNumber;
			(*iterFaceFieldData).spectralRadius = (*iterFaceGaussData).spectralRadius;
			(*iterFaceFieldData).area = (*iterFaceGaussData).area;
			(*iterFaceFieldData).areaReference = (*iterFaceGaussData).areaReference;
			(*iterFaceFieldData).normalVector = (*iterFaceGaussData).normalVector;

			(*iterFaceFieldData).sideOff = (*iterFaceGaussData).sideOff;
			(*iterFaceFieldData).faceNodeNumber = (*iterFaceGaussData).faceNodeNumber;
			(*iterFaceFieldData).faceNode.resize((*iterFaceFieldData).faceNodeNumber + 1);

			(*iterFaceFieldData).lambdaFace = 0.0;

			(*iterFaceFieldData).o2PointNumber = (*iterFaceGaussData).o2PointNumber;

			// std::cout << "--------" << std::endl;
			// std::cout << (*iterFaceFieldData).area << "\t" << (*iterFaceGaussData).area << std::endl;
			// std::cout << (*iterFaceFieldData).normalVector.x << "\t" << (*iterFaceFieldData).normalVector.y
			//	<< "\t" << (*iterFaceGaussData).normalVector.x << "\t" << (*iterFaceGaussData).normalVector.y << std::endl;
			(*iterFaceFieldData).faceCellIndex.resize(2 + 1);
			(*iterFaceFieldData).faceCellIndex[1] = (*iterFaceGaussData).faceCellIndex[1];
			(*iterFaceFieldData).faceCellIndex[2] = (*iterFaceGaussData).faceCellIndex[2];

			for (int jj = 1; jj < (*iterFaceFieldData).faceNodeNumber + 1; ++jj)
			{
				(*iterFaceFieldData).faceNode[jj].first = (*iterFaceGaussData).faceNode[jj].first;
				(*iterFaceFieldData).faceNode[jj].second = (*iterFaceGaussData).faceNode[jj].second;
			}

			(*iterFaceFieldData).parametricValue.resize(static_cast<int>((*iterFaceGaussData).PG[fO]));
			for (int jj = 0; jj < static_cast<int>((*iterFaceGaussData).PG[fO]); ++jj)
			{
				(*iterFaceFieldData).parametricValue[jj].first = (*iterFaceGaussData).parametricValue[jj].first;
				(*iterFaceFieldData).parametricValue[jj].second = (*iterFaceGaussData).parametricValue[jj].second;
			}

			(*iterFaceFieldData).parametricArea = (*iterFaceGaussData).parametricArea;

			(*iterFaceFieldData).gaussPairVector_.resize(static_cast<int>((*iterFaceGaussData).PG[fO]));
			for (int jj = 0; jj < static_cast<int>((*iterFaceGaussData).PG[fO]); ++jj)
			{
				(*iterFaceFieldData).gaussPairVector_[jj].JacobiCof = (*iterFaceGaussData).gaussPairVector_[jj].JacobiCof;
				(*iterFaceFieldData).gaussPairVector_[jj].normalVector = (*iterFaceGaussData).gaussPairVector_[jj].normalVector;
				(*iterFaceFieldData).gaussPairVector_[jj].p = (*iterFaceGaussData).gaussPairVector_[jj].p;
			}

			(*iterFaceFieldData).fPG = (*iterFaceGaussData).PG[fO];
			(*iterFaceFieldData).NDOFS = reconstruction_->ReturnNDIFFS();
			(*iterFaceFieldData).faceWeightVF.resize((*iterFaceFieldData).NDOFS);

			(*iterFaceFieldData).diffBaseValueData.resize(2);
			for (auto &i : (*iterFaceFieldData).diffBaseValueData)
			{
				i.resize((*iterFaceGaussData).PG[fO]);
				for (auto &j : i)
				{
					j.resize(reconstruction_->ReturnNDOFS());
					for (auto &k : j)
						k.resize(reconstruction_->ReturnNDIFFS(), UNINITReal);
				}
			}
			(*iterFaceFieldData).diffBaseValueDataCR.resize(2);
			for (auto &i : (*iterFaceFieldData).diffBaseValueDataCR)
			{
				i.resize((*iterFaceGaussData).PG[fO]);
				for (auto &j : i)
				{
					j.resize(reconstruction_->ReturnNDOFSCR());
					for (auto &k : j)
						k.resize(reconstruction_->ReturnNDIFFSCR(), UNINITReal);
				}
			}
			(*iterFaceFieldData).diffBaseValueDataMid.resize(2);
			for (auto &j : (*iterFaceFieldData).diffBaseValueDataMid)
			{
				j.resize(reconstruction_->ReturnNDOFS());
				for (auto &k : j)
					k.resize(reconstruction_->ReturnNDIFFS(), UNINITReal);
			}

			// flux?
			// debug
			//	std::cout << "----------------------------------------------" << std::endl;
			//	std::cout << (*iterFaceFieldData).index << std::endl;
			//	std::cout << (*iterFaceFieldData).faceNode[1].first << "\t" << (*iterFaceFieldData).faceNode[2].first << std::endl;
			//	std::cout << (*iterFaceFieldData).normalVector.x << "\t" << (*iterFaceFieldData).normalVector.y << std::endl;
		}
		std::cout << " ..face variables array has been allocated." << std::endl;
		return true;
	}

	bool fieldsolver::initializeArray()
	{
		//	real Lb = -10.0;
		//	real Rb = 20.0;
		//	real Db = -10.0;
		//	real Ub = 10.0;
		//	real pi = 3.1415926535898;
		//	cellFieldDataVector::iterator iterCellFieldData;
		//	for (iterCellFieldData = cellFieldData_->begin(); iterCellFieldData != cellFieldData_->end(); ++iterCellFieldData){
		//		(*iterCellFieldData).scalarVariableTn[0] = 0.0;
		//		(*iterCellFieldData).artificialViscousParameter = 1.0 / std::pow((rO + 1), 2) * parameter_->cofAV;
		//	}

		real Lb = -5.0;
		real Rb = 5.0;
		real Db = -5.0;
		real Ub = 5.0;
		real pi = 3.1415926535898;
		cellFieldDataVector::iterator iterCellFieldData;
		for (iterCellFieldData = cellFieldData_->begin(); iterCellFieldData != cellFieldData_->end(); ++iterCellFieldData)
		{
			//	if ((*iterCellFieldData).baryCenter.x <= 0.0
			//		&& (*iterCellFieldData).baryCenter.x <= (*iterCellFieldData).baryCenter.y
			//		&& -(*iterCellFieldData).baryCenter.x >= (*iterCellFieldData).baryCenter.y){
			//		(*iterCellFieldData).scalarVariableTn[0] = (*iterCellFieldData).baryCenter.x - Lb;
			//	}
			//	else if ((*iterCellFieldData).baryCenter.y <= 0.0
			//		&& (*iterCellFieldData).baryCenter.y <= (*iterCellFieldData).baryCenter.x
			//		&& -(*iterCellFieldData).baryCenter.y >= (*iterCellFieldData).baryCenter.x){
			//		(*iterCellFieldData).scalarVariableTn[0] = (*iterCellFieldData).baryCenter.y - Db;
			//	}
			//	else if ((*iterCellFieldData).baryCenter.x > 0.0
			//		&& -(*iterCellFieldData).baryCenter.x <= (*iterCellFieldData).baryCenter.y
			//		&& (*iterCellFieldData).baryCenter.x >= (*iterCellFieldData).baryCenter.y){
			//		(*iterCellFieldData).scalarVariableTn[0] = Rb - (*iterCellFieldData).baryCenter.x;
			//	}
			//	else if ((*iterCellFieldData).baryCenter.y > 0.0
			//		&& -(*iterCellFieldData).baryCenter.y <= (*iterCellFieldData).baryCenter.x
			//		&& (*iterCellFieldData).baryCenter.y >= (*iterCellFieldData).baryCenter.x){
			//		(*iterCellFieldData).scalarVariableTn[0] = Ub - (*iterCellFieldData).baryCenter.y;
			//	}
			(*iterCellFieldData).scalarVariableTn[0] = 0.0;
			(*iterCellFieldData).artificialViscousParameter = std::pow((rO + 1), 2) * parameter_->cofAV;
		}
		return true;
	}

	bool fieldsolver::timeMarchingExplicitSSPRK(
		vertexVector *node,
		GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO1Grid<fO> *gaussIntegralFace,
		const std::string &filenameInputBackup,
		const std::string &filenameOutputBackup,
		const std::string &filenameOutputSln,
		const std::string &filenameOutputResidual)
	{

		real residual;
		int IOflag = 0;
		std::ofstream outBoundCheck(ScalarCfv::fineOut_BndCheck);

		allocateArray();
		if (parameter_->isContinue == true)
		{
			readBackUp(filenameInputBackup);
		}
		else
		{
			initializeArray();
		}

		reconstruction_->initBaseMomentAndRelaxFactor(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			gaussIntegralCell);
		reconstruction_->initFaceWeight(
			parameter_,
			cellFieldData_,
			faceFieldData_);
		reconstruction_->initReconstructionMatrixAndVector(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace);
		reconstruction_->initReconstructionMatrixAndVectorCR(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace);
		exportCurrentSln(node, 0, filenameOutputSln);
		// main loop
		parameter_->refTn = 0.0; // global time marching
		for (int iStep = parameter_->nStart; iStep <= parameter_->nEnd; ++iStep)
		{
			evolution_->getTimeStep(
				parameter_,
				cellFieldData_,
				faceFieldData_);
			if (parameter_->isLocalTimeMarching == false)
			{
				parameter_->refTn += parameter_->refDeltaTn;
			}

			for (int irk = 1; irk <= parameter_->nStepTimeMarching; ++irk)
			{
				for (int isubrk = 1; isubrk <= parameter_->nInnerStepTimeMarching; ++isubrk)
				{
					for (int irec = 1; irec <= 1; irec++)
					reconstruction_->excuteReconstruction(
						parameter_,
						cellFieldData_,
						cellGaussData_,
						faceFieldData_,
						faceGaussData_,
						gaussIntegralCell,
						gaussIntegralFace);
					reconstruction_->excuteReconstructionCR(
						parameter_,
						cellFieldData_,
						cellGaussData_,
						faceFieldData_,
						faceGaussData_,
						gaussIntegralCell,
						gaussIntegralFace);
					// deal with AV?
					evolution_->getSourceTerm(
						parameter_,
						cellFieldData_,
						faceFieldData_,
						gaussIntegralCell);

					//------------------------------

					evolution_->getArtificialViscosityTerm(
						parameter_,
						cellFieldData_,
						faceFieldData_,
						gaussIntegralFace);

					//------------------------------

					if (parameter_->nStepTimeMarching == 3)
					{
						evolution_->excuteInnerLoopExplicitSSPRKStep3(
							irk,
							cellFieldData_);
					}
					else if (parameter_->nStepTimeMarching == 4)
					{
						evolution_->excuteInnerLoopExplicitSSPRKStep4(
							irk,
							cellFieldData_);
					}
					// evolution_->applyBoundaryLimiter(
					//	parameter_,
					//	cellFieldData_);
				}

				if (iStep % parameter_->nScreenOutput == 0)
				{
					std::cout << "	RK step:"
							  << "\t"
							  << "\t" << irk << "has been done." << std::endl;
				}
			}
			// evolution_->updateArtificalViscosity(
			//	iStep,
			//	1.0 / std::pow((rO + 1), 2) * parameter_->cofAV,
			//	cellFieldData_);
			// post -> residual
			getResidual(residual);
			exportResidual(
				IOflag,
				iStep,
				residual,
				filenameOutputResidual);

			// evolution_->excuteFilter(parameter_,
			//	cellFieldData_);

			if (iStep % parameter_->nScreenOutput == 0)
			{
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "Current step:"
						  << "\t"
						  << "\t" << iStep << std::endl;
			}
			if (iStep % parameter_->nFileOutput == 0)
			{
				exportCurrentSln(node, iStep, filenameOutputSln);
				outBoundCheck << "Step" << iStep << std::endl;
				reconstruction_->getBoundaryValueCheck(
					parameter_,
					cellFieldData_,
					cellGaussData_,
					faceFieldData_,
					faceGaussData_,
					gaussIntegralCell,
					gaussIntegralFace, outBoundCheck);
			}

			if (iStep % parameter_->nBackUp == 0)
			{
				// exportBackUp(filenameOutputBackup);
				exportBackUp(iStep, filenameOutputBackup);
			}
		}
		outBoundCheck.close();
		return true;
	}

	bool fieldsolver::timeMarchingExplicitSSPRK(
		vertexVector *node,
		GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO2Grid<fO> *gaussIntegralFace,
		const std::string &filenameInputBackup,
		const std::string &filenameOutputBackup,
		const std::string &filenameOutputSln,
		const std::string &filenameOutputResidual)
	{

		real residual;
		int IOflag = 0;

		allocateArray();
		if (parameter_->isContinue == true)
		{
			readBackUp(filenameInputBackup);
		}
		else
		{
			initializeArray();
		}

		reconstruction_->initBaseMomentAndRelaxFactor(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			gaussIntegralCell);
		reconstruction_->initFaceWeight(
			parameter_,
			cellFieldData_,
			faceFieldData_);
		reconstruction_->initReconstructionMatrixAndVector(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace);
		reconstruction_->initReconstructionMatrixAndVectorCR(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace);
		exportCurrentSln(node, 0, filenameOutputSln);
		// main loop
		parameter_->refTn = 0.0; // global time marching
		for (int iStep = parameter_->nStart; iStep <= parameter_->nEnd; ++iStep)
		{
			evolution_->getTimeStep(
				parameter_,
				cellFieldData_,
				faceFieldData_);
			if (parameter_->isLocalTimeMarching == false)
			{
				parameter_->refTn += parameter_->refDeltaTn;
			}

			for (int irk = 1; irk <= parameter_->nStepTimeMarching; ++irk)
			{
				for (int isubrk = 1; isubrk <= parameter_->nInnerStepTimeMarching; ++isubrk)
				{
					// for (int irec = 1; irec <= 3; irec++)
					reconstruction_->excuteReconstruction(
						parameter_,
						cellFieldData_,
						cellGaussData_,
						faceFieldData_,
						faceGaussData_,
						gaussIntegralCell,
						gaussIntegralFace);
					reconstruction_->excuteReconstructionCR(
						parameter_,
						cellFieldData_,
						cellGaussData_,
						faceFieldData_,
						faceGaussData_,
						gaussIntegralCell,
						gaussIntegralFace);
					// deal with AV?
					evolution_->getSourceTerm(
						parameter_,
						cellFieldData_,
						faceFieldData_,
						gaussIntegralCell);
					if (parameter_->nStepTimeMarching == 3)
					{
						evolution_->excuteInnerLoopExplicitSSPRKStep3(
							irk,
							cellFieldData_);
					}
					else if (parameter_->nStepTimeMarching == 4)
					{
						evolution_->excuteInnerLoopExplicitSSPRKStep4(
							irk,
							cellFieldData_);
					}
				}

				if (iStep % parameter_->nScreenOutput == 0)
				{
					std::cout << "	RK step:"
							  << "\t"
							  << "\t" << irk << "has been done." << std::endl;
				}
			}
			evolution_->updateArtificalViscosity(
				iStep,
				1.0 / std::pow((rO + 1), 2) * parameter_->cofAV,
				cellFieldData_);
			// post -> residual
			getResidual(residual);
			exportResidual(
				IOflag,
				iStep,
				residual,
				filenameOutputResidual);

			if (iStep % parameter_->nScreenOutput == 0)
			{
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "Current step:"
						  << "\t"
						  << "\t" << iStep << std::endl;
			}
			if (iStep % parameter_->nFileOutput == 0)
			{
				exportCurrentSln(node, iStep, filenameOutputSln);
			}

			if (iStep % parameter_->nBackUp == 0)
			{
				exportBackUp(filenameOutputBackup);
			}
		}
		return true;
	}

	bool fieldsolver::timeMarchingExplicit_DEBUG(
		vertexVector *node,
		GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO1Grid<fO> *gaussIntegralFace,
		const std::string &filenameInputBackup,
		const std::string &filenameOutputBackup,
		const std::string &filenameOutputSln,
		const std::string &filenameOutputResidual)
	{

		real residual;
		int IOflag = 0;
		cellFieldDataVector::iterator iterCell;

		allocateArray();
		if (parameter_->isContinue == true)
		{
			readBackUp(filenameInputBackup);
		}
		else
		{
			initializeArray();
		}

		reconstruction_->initBaseMomentAndRelaxFactor(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			gaussIntegralCell);
		reconstruction_->initFaceWeight(
			parameter_,
			cellFieldData_,
			faceFieldData_);
		reconstruction_->initReconstructionMatrixAndVector(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace);
		reconstruction_->initReconstructionMatrixAndVectorCR(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace);
		exportCurrentSln(node, 0, filenameOutputSln);
		//#define RECON_STEADY

#ifdef RECON_STEADY

		for (int ii = 1; ii < 2; ++ii)
		{
			reconstruction_->excuteReconstruction(
				parameter_,
				cellFieldData_,
				cellGaussData_,
				faceFieldData_,
				faceGaussData_,
				gaussIntegralCell,
				gaussIntegralFace);
		}

		reconstruction_->excuteReconstructionCR(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace);

		exportCurrentSln(node, 1, filenameOutputSln);
		system("pause");
#endif

		// main loop
		parameter_->refTn = 0.0; // global time marching
		for (int iStep = parameter_->nStart; iStep <= parameter_->nEnd; ++iStep)
		{
			evolution_->getTimeStep(
				parameter_,
				cellFieldData_,
				faceFieldData_);
			if (parameter_->isLocalTimeMarching == false)
			{
				parameter_->refTn += parameter_->refDeltaTn;
			}

			reconstruction_->excuteReconstruction(
				parameter_,
				cellFieldData_,
				cellGaussData_,
				faceFieldData_,
				faceGaussData_,
				gaussIntegralCell,
				gaussIntegralFace);

			// evolution_->excuteWBAPLimiter4th(
			//	parameter_,
			//	cellFieldData_);

			reconstruction_->excuteReconstructionCR(
				parameter_,
				cellFieldData_,
				cellGaussData_,
				faceFieldData_,
				faceGaussData_,
				gaussIntegralCell,
				gaussIntegralFace);

			// evolution_->excuteWBAPLimiter4thCR(
			//	parameter_,
			//	cellFieldData_);

			// deal with AV?
			evolution_->getSourceTerm(
				parameter_,
				cellFieldData_,
				faceFieldData_,
				gaussIntegralCell);

			//...
			for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
			{
				(*iterCell).scalarVariableTm[0] = (*iterCell).scalarVariableTn[0];
			}
			for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
			{
				(*iterCell).scalarVariableTn[0] = (*iterCell).scalarVariableTm[0] +
												  (*iterCell).dTPhysical * (*iterCell).timeMarchingRHSTn;
			}

			evolution_->updateArtificalViscosity(
				iStep,
				1.0 / std::pow((rO + 1), 2) * parameter_->cofAV,
				cellFieldData_);
			// post -> residual
			getResidual(residual);
			exportResidual(IOflag, iStep, residual, filenameOutputResidual);

			if (iStep % parameter_->nScreenOutput == 0)
			{
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "Current step:"
						  << "\t"
						  << "\t" << iStep << std::endl;
			}
			if (iStep % parameter_->nFileOutput == 0)
			{
				exportCurrentSln(node, iStep, filenameOutputSln);
			}
			if (iStep % parameter_->nBackUp == 0)
			{
				exportBackUp(filenameOutputBackup);
			}
		}
		return true;
	}

	bool fieldsolver::timeMarchingExplicit_DEBUG(
		vertexVector *node,
		GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO2Grid<fO> *gaussIntegralFace,
		const std::string &filenameInputBackup,
		const std::string &filenameOutputBackup,
		const std::string &filenameOutputSln,
		const std::string &filenameOutputResidual)
	{

		real residual;
		int IOflag = 0;
		cellFieldDataVector::iterator iterCell;

		allocateArray();
		if (parameter_->isContinue == true)
		{
			readBackUp(filenameInputBackup);
		}
		else
		{
			initializeArray();
		}

		reconstruction_->initBaseMomentAndRelaxFactor(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			gaussIntegralCell);
		reconstruction_->initFaceWeight(
			parameter_,
			cellFieldData_,
			faceFieldData_);
		reconstruction_->initReconstructionMatrixAndVector(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace);
		reconstruction_->initReconstructionMatrixAndVectorCR(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace);
		exportCurrentSln(node, 0, filenameOutputSln);
		//#define RECON_STEADY

#ifdef RECON_STEADY

		for (int ii = 1; ii < 2; ++ii)
		{
			reconstruction_->excuteReconstruction(
				parameter_,
				cellFieldData_,
				cellGaussData_,
				faceFieldData_,
				faceGaussData_,
				gaussIntegralCell,
				gaussIntegralFace);
		}

		reconstruction_->excuteReconstructionCR(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace);

		exportCurrentSln(node, 1, filenameOutputSln);
		system("pause");
#endif

		// main loop
		parameter_->refTn = 0.0; // global time marching
		for (int iStep = parameter_->nStart; iStep <= parameter_->nEnd; ++iStep)
		{
			evolution_->getTimeStep(
				parameter_,
				cellFieldData_,
				faceFieldData_);
			if (parameter_->isLocalTimeMarching == false)
			{
				parameter_->refTn += parameter_->refDeltaTn;
			}

			reconstruction_->excuteReconstruction(
				parameter_,
				cellFieldData_,
				cellGaussData_,
				faceFieldData_,
				faceGaussData_,
				gaussIntegralCell,
				gaussIntegralFace);

			// evolution_->excuteWBAPLimiter4th(
			//	parameter_,
			//	cellFieldData_);

			reconstruction_->excuteReconstructionCR(
				parameter_,
				cellFieldData_,
				cellGaussData_,
				faceFieldData_,
				faceGaussData_,
				gaussIntegralCell,
				gaussIntegralFace);

			// evolution_->excuteWBAPLimiter4thCR(
			//	parameter_,
			//	cellFieldData_);

			// deal with AV?
			evolution_->getSourceTerm(
				parameter_,
				cellFieldData_,
				faceFieldData_,
				gaussIntegralCell);

			//...
			for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
			{
				(*iterCell).scalarVariableTm[0] = (*iterCell).scalarVariableTn[0];
			}
			for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
			{
				(*iterCell).scalarVariableTn[0] = (*iterCell).scalarVariableTm[0] +
												  (*iterCell).dTPhysical * (*iterCell).timeMarchingRHSTn;
			}

			evolution_->updateArtificalViscosity(
				iStep,
				1.0 / std::pow((rO + 1), 2) * parameter_->cofAV,
				cellFieldData_);
			// post -> residual
			getResidual(residual);
			exportResidual(IOflag, iStep, residual, filenameOutputResidual);

			if (iStep % parameter_->nScreenOutput == 0)
			{
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "Current step:"
						  << "\t"
						  << "\t" << iStep << std::endl;
			}
			if (iStep % parameter_->nFileOutput == 0)
			{
				exportCurrentSln(node, iStep, filenameOutputSln);
			}
			if (iStep % parameter_->nBackUp == 0)
			{
				exportBackUp(filenameOutputBackup);
			}
		}
		return true;
	}

	bool fieldsolver::readBackUp(
		const std::string &filename)
	{

		std::ifstream fileIn;
		std::string fileName = filename;

		fileIn.open(fileName.c_str());
		if (!fileIn)
		{
			std::cerr << "	error: unable to open \"" << fileName << "\"\n";
			exit(1);
		}

		cellFieldDataVector::iterator iterCell;
		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			fileIn >> (*iterCell).scalarVariableTn[0];
		}
		return true;
	};

	bool fieldsolver::exportBackUp(
		const std::string &filename)
	{

		std::ofstream fileOut;
		std::string fileName = filename;
		fileOut.open(fileName.c_str());
		if (!fileOut)
		{
			std::cerr << "	error: unable to open \"" << fileName << "\"\n";
			exit(1);
		}

		cellFieldDataVector::iterator iterCell;
		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			fileOut << (*iterCell).scalarVariableTn[0] << "\t";
		}

		return true;
	}

	bool fieldsolver::exportBackUp(
		const int ii,
		const std::string &filename)
	{

		// std::ofstream fileOut;
		// std::string fileName = filename;
		// fileOut.open(fileName.c_str());
		// if (!fileOut){
		//	std::cerr << "	error: unable to open \"" << fileName << "\"\n";
		//	exit(1);
		// }

		std::ofstream fileOut;
		std::string fileName = filename;
		std::stringstream ss;
		ss << ii;
		fileName += "_";
		fileName += ss.str();
		fileName += ".sav";
		fileOut.open(fileName.c_str());
		if (!fileOut)
		{
			std::cerr << "		error: unable to open \"" << fileName << "\"\n";
			exit(1);
		}

		cellFieldDataVector::iterator iterCell;
		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			fileOut << (*iterCell).scalarVariableTn[0] << "\t";
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				fileOut << (*iterCell).scalarVariableTn[kk] << "\t";
			}
		}

		return true;
	}

	bool fieldsolver::exportCurrentSln(
		vertexVector *node,
		const int step,
		const std::string &filename)
	{
		//#define DEBUG

#ifdef DEBUG
		int lineControl = 5;

		std::ofstream fileOut;
		std::string fileName = filename;
		std::stringstream ss;
		ss << step;
		fileName += "_";
		fileName += ss.str();
		fileName += ".plt";
		fileOut.open(fileName.c_str());
		if (!fileOut)
		{
			std::cerr << "		error: unable to open \"" << fileName << "\"\n";
			exit(1);
		}
		std::cout << "	in Tecplot format, to file : " << filename << "..." << std::endl;
		fileOut << "VARIABLES = \"x\", \"y\", \"f\", \"m1\", \"m2\", \"m3\", \"m4\", \"m5\", \"m6\", \"m7\", \"m8\", \"m9\"\n";
		fileOut << "ZONE N =" << parameter_->nodeNumber << ","
				<< "E=" << parameter_->cellNumber << ","
				<< "VARLOCATION=([1-2]=NODAL,[3-12]=CELLCENTERED)" << std::endl;
		fileOut << ", "
				<< "DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL" << std::endl;
		fileOut << std::setprecision(15);

		vertexVector::iterator iterNode;
		for (iterNode = node->begin(); iterNode != node->end(); ++iterNode)
		{
			fileOut << (*iterNode).nodePhysical.x << "\t";
			if (((*iterNode).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';
		for (iterNode = node->begin(); iterNode != node->end(); ++iterNode)
		{
			fileOut << (*iterNode).nodePhysical.y << "\t";
			if (((*iterNode).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';
		cellFieldDataVector::iterator iterCell;
		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			real f = (*iterCell).timeMarchingRHSTn;
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			real f = (*iterCell).scalarVariableTnCR[1];
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			real f = (*iterCell).scalarVariableTnCR[2];
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			real f = (*iterCell).scalarVariableTnCR[3];
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			real f = (*iterCell).scalarVariableTnCR[4];
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			real f = (*iterCell).scalarVariableTnCR[5];
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			real f = (*iterCell).scalarVariableTnCR[6];
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			real f = (*iterCell).scalarVariableTnCR[7];
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			real f = (*iterCell).scalarVariableTnCR[8];
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			real f = (*iterCell).scalarVariableTnCR[9];
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			if ((*iterCell).cellType_ == Triangle)
			{
				for (int jj = 1; jj < (*iterCell).cellNodeNumber + 1; ++jj)
				{
					fileOut << (*iterCell).cellNode[jj].first << " ";
				}
				fileOut << (*iterCell).cellNode[(*iterCell).cellNodeNumber].first;
				fileOut << '\n';
			}
			else if ((*iterCell).cellType_ == Quadrilateral)
			{
				for (int jj = 1; jj < (*iterCell).cellNodeNumber + 1; ++jj)
				{
					fileOut << (*iterCell).cellNode[jj].first << " ";
				}
				fileOut << '\n';
			}
		}
		fileOut << std::;
		fileOut.close();

#else
		int lineControl = 5;

		std::ofstream fileOut;
		std::string fileName = filename;
		std::stringstream ss;
		ss << step;
		fileName += "_";
		fileName += ss.str();
		fileName += ".plt";
		fileOut.open(fileName.c_str());
		if (!fileOut)
		{
			std::cerr << "		error: unable to open \"" << fileName << "\"\n";
			exit(1);
		}
		std::cout << "	exporting sln in Tecplot format, to file : " << fileName << "..." << '\n';

#ifndef USE_RBFB1_N
		fileOut << "VARIABLES = \"x\", \"y\", \"sln\", \"smooth\", \"dx\", \"dy\"\n";
		fileOut << "ZONE N =" << parameter_->nodeNumber << ","
				<< "E=" << parameter_->cellNumber << ","
				<< "VARLOCATION=([1-2]=NODAL,[3-6]=CELLCENTERED)" << '\n';
#else
		fileOut << "VARIABLES = \"x\", \"y\", \"sln\", \"smooth\", \"dx\", \"dy\"\n";
		fileOut << "ZONE N =" << parameter_->nodeNumber << ","
				<< "E=" << parameter_->cellNumber << ","
				<< "VARLOCATION=([1-3]=NODAL,[4-6]=CELLCENTERED)" << '\n';
#endif
		fileOut << ", "
				<< "DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL" << '\n';
		fileOut << std::setprecision(15);

		vertexVector::iterator iterNode;
		for (iterNode = node->begin(); iterNode != node->end(); ++iterNode)
		{
			fileOut << (*iterNode).nodePhysical.x << "\t";
			if (((*iterNode).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';
		for (iterNode = node->begin(); iterNode != node->end(); ++iterNode)
		{
			fileOut << (*iterNode).nodePhysical.y << "\t";
			if (((*iterNode).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';
		cellFieldDataVector::iterator iterCell;

		std::vector<real> nodeSln(node->size(), 0.0);
		std::vector<real> nodeCount(node->size(), 0.0);
		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			// base moment i
			tensor1D<real, NDOFS> momentI;
			// base i
			tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
			point p = (*iterCell).baryCenter;
			point baryCenterI = (*iterCell).baryCenter;
			point scaleI = (*iterCell).lengthReference;
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				momentI[kk] = (*iterCell).baseMoment[kk];
			}

#ifndef USE_RBFB1_N

#ifdef USE_RBFB1
			assert(iterCell->cellType_ == Quadrilateral);
			RBFB1GetDiffBaseValue(
				point(0.5, 0.5),
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI,
				*iterCell);
#else
			CfvMath::getDiffBaseValue(
				p,
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI);
#endif

#else

			RBFB1GetDiffBaseValue(
				point(0, 0),
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI,
				*iterCell);
			nodeCount[iterCell->cellNode[1].first - 1] += 1;
			nodeSln[iterCell->cellNode[1].first - 1] += (*iterCell).scalarVariableTn[0];
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				nodeSln[iterCell->cellNode[1].first - 1] += matrixDiffBaseI[kk][0] * (*iterCell).scalarVariableTn[kk];
			}

			RBFB1GetDiffBaseValue(
				point(0, 1),
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI,
				*iterCell);
			nodeCount[iterCell->cellNode[2].first - 1] += 1;
			nodeSln[iterCell->cellNode[2].first - 1] += (*iterCell).scalarVariableTn[0];
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				nodeSln[iterCell->cellNode[2].first - 1] += matrixDiffBaseI[kk][0] * (*iterCell).scalarVariableTn[kk];
			}

			RBFB1GetDiffBaseValue(
				point(1, 0),
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI,
				*iterCell);
			nodeCount[iterCell->cellNode[3].first - 1] += 1;
			nodeSln[iterCell->cellNode[3].first - 1] += (*iterCell).scalarVariableTn[0];
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				nodeSln[iterCell->cellNode[3].first - 1] += matrixDiffBaseI[kk][0] * (*iterCell).scalarVariableTn[kk];
			}

			RBFB1GetDiffBaseValue(
				point(1, 1),
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI,
				*iterCell);
			nodeCount[iterCell->cellNode[4].first - 1] += 1;
			nodeSln[iterCell->cellNode[4].first - 1] += (*iterCell).scalarVariableTn[0];
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				nodeSln[iterCell->cellNode[4].first - 1] += matrixDiffBaseI[kk][0] * (*iterCell).scalarVariableTn[kk];
			}
#endif

			real f = (*iterCell).scalarVariableTn[0];
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				f += matrixDiffBaseI[kk][0] * (*iterCell).scalarVariableTn[kk];
			}
#ifndef USE_RBFB1_N
			fileOut << f << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
#endif
		}
#ifdef USE_RBFB1_N
		for (int i = 0; i < nodeSln.size(); i++)
		{
			fileOut << nodeSln[i] / nodeCount[i] << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
#endif
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			// base moment i
			tensor1D<real, NDOFS> momentI;
			// base i
			tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
			point p = (*iterCell).baryCenter;
			point baryCenterI = (*iterCell).baryCenter;
			point scaleI = (*iterCell).lengthReference;
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				momentI[kk] = (*iterCell).baseMoment[kk];
			}
#ifndef USE_RBFB1
			CfvMath::getDiffBaseValue(
				p,
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI);
#else
			assert(iterCell->cellType_ == Quadrilateral);
			RBFB1GetDiffBaseValue(
				point(0.5, 0.5),
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI,
				*iterCell);
#endif

			point dPhi;
			dPhi.setZero();
			real tempX, tempY;
			tempX = 0.0;
			tempY = 0.0;
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				tempX += matrixDiffBaseI[kk][1] * (*iterCell).scalarVariableTn[kk];
				tempY += matrixDiffBaseI[kk][2] * (*iterCell).scalarVariableTn[kk];
			}
			dPhi.x = tempX;
			dPhi.y = tempY;

			fileOut << std::fabs(1.0 - std::pow(dPhi.length(), 2)) << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		//------------------
		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			// base moment i
			tensor1D<real, NDOFS> momentI;
			// base i
			tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
			point p = (*iterCell).baryCenter;
			point baryCenterI = (*iterCell).baryCenter;
			point scaleI = (*iterCell).lengthReference;
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				momentI[kk] = (*iterCell).baseMoment[kk];
			}
#ifndef USE_RBFB1
			CfvMath::getDiffBaseValue(
				p,
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI);
#else
			assert(iterCell->cellType_ == Quadrilateral);
			RBFB1GetDiffBaseValue(
				point(0.5, 0.5),
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI,
				*iterCell);
#endif

			real dfdx = 0.0;
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				dfdx += matrixDiffBaseI[kk][1] * (*iterCell).scalarVariableTn[kk];
			}

			fileOut << dfdx << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			// base moment i
			tensor1D<real, NDOFS> momentI;
			// base i
			tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
			point p = (*iterCell).baryCenter;
			point baryCenterI = (*iterCell).baryCenter;
			point scaleI = (*iterCell).lengthReference;
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				momentI[kk] = (*iterCell).baseMoment[kk];
			}
#ifndef USE_RBFB1
			CfvMath::getDiffBaseValue(
				p,
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI);
#else
			assert(iterCell->cellType_ == Quadrilateral);
			RBFB1GetDiffBaseValue(
				point(0.5, 0.5),
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI,
				*iterCell);

#endif

			real dfdy = 0.0;
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				dfdy += matrixDiffBaseI[kk][2] * (*iterCell).scalarVariableTn[kk];
			}

			fileOut << dfdy << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << '\n';
		//------------------

		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			if ((*iterCell).cellType_ == Triangle)
			{
				for (int jj = 1; jj < (*iterCell).cellNodeNumber + 1 - (*iterCell).o2PointNumber; ++jj)
				{
					fileOut << (*iterCell).cellNode[jj].first << " ";
				}
				fileOut << (*iterCell).cellNode[(*iterCell).cellNodeNumber - (*iterCell).o2PointNumber].first;
				fileOut << '\n';
			}
			else if ((*iterCell).cellType_ == Quadrilateral)
			{
				for (int jj = 1; jj < (*iterCell).cellNodeNumber + 1 - (*iterCell).o2PointNumber; ++jj)
				{
					fileOut << (*iterCell).cellNode[jj].first << " ";
				}
				fileOut << '\n';
			}
		}
		fileOut << '\n';
		fileOut.close();
#endif
		return true;
	};

	bool fieldsolver::getResidual(
		real &residual)
	{
#ifdef IF_RESTRICT_RADIUS
		residual = 0.0;
		real volTemp = 0.0;
		real R = 12.0;
		cellFieldDataVector::iterator iterCell;
		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			if ((*iterCell).baryCenter.length() < R)
			{
				continue;
			}
			residual += std::fabs((*iterCell).scalarVariableTn[0] - (*iterCell).scalarVariableTm[0]) * (*iterCell).volume;
			volTemp += (*iterCell).volume;
		}
		residual /= volTemp;
#else
		residual = 0.0;
		real volTemp = 0.0;
		cellFieldDataVector::iterator iterCell;
		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell)
		{
			residual += std::fabs((*iterCell).scalarVariableTn[0] - (*iterCell).scalarVariableTm[0]) * (*iterCell).volume;
			volTemp += (*iterCell).volume;
		}
		residual /= volTemp;
#endif

		return true;
	};

	bool fieldsolver::exportResidual(
		int &flag,
		const int step,
		real residual,
		const std::string &filename)
	{

		std::ofstream fileOut;
		fileOut.open(filename.c_str(), std::ios_base::app);
		if (!fileOut)
		{
			std::cerr << "		error: unable to open \"" << filename << "\"\n";
			exit(1);
		}
		if (flag == 0)
		{
			fileOut << "VARIABLES = \"Steps\", \"L1(Res)\"\n";
		}
		fileOut << step << "\t" << residual << std::endl;
		if (isnan(residual))
		{
			std::cout << "\
 | \\ | |     /\\     | \\ | |           |  _ \\  | |       / __ \\  \\ \\        / / | |  | | |  __ \\    | | | | | |\n\
 |  \\| |    /  \\    |  \\| |           | |_) | | |      | |  | |  \\ \\  /\\  / /  | |  | | | |__) |   | | | | | |\n\
 | . ` |   / /\\ \\   | . ` |           |  _ <  | |      | |  | |   \\ \\/  \\/ /   | |  | | |  ___/    | | | | | |\n\
 | |\\  |  / ____ \\  | |\\  |           | |_) | | |____  | |__| |    \\  /\\  /    | |__| | | |        |_| |_| |_|\n\
 |_| \\_| /_/    \\_\\ |_| \\_|           |____/  |______|  \\____/      \\/  \\/      \\____/  |_|        (_) (_) (_)\n\
";
			exit(-5);
		}
		flag += 1;
		return true;
	};

}