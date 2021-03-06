#ifndef _EVOLUTION_H
#define _EVOLUTION_H
//#include "TypeDefine.h"
//#include "Point.h"
//#include "Parameter.h"
//#include "Geometry.h"
//#include "GaussIntegral.h"
#include "Variable.h"
#include "Parameter.h"
#include <math.h>

namespace ScalarCfv
{
	template <unsigned int O>
	class evolution
	{
	public:
		evolution(){};
		~evolution(){};

		const static unsigned int NDOFS = (O + 2) * (O + 1) / 2;
		const static unsigned int NDIFFS = (O + 2) * (O + 1) / 2;

		// const static unsigned int NDOFS = (O + 2) * (O + 1) / 2 + 1;
		// const static unsigned int NDIFFS = (O + 2 + 1) * (O + 1 + 1) / 2;

	public:
		virtual bool getFaceFlux(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			faceFieldDataVector *faceFieldData)
		{
			return true;
		};
		virtual bool getSourceTerm(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			faceFieldDataVector *faceFieldData,
			GaussIntegralCellO1Grid<vO> *gaussIntegralCell); //
		virtual bool getSourceTerm(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			faceFieldDataVector *faceFieldData,
			GaussIntegralCellO2Grid<vO> *gaussIntegralCell); //

		virtual bool getArtificialViscosityTerm(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			faceFieldDataVector *faceFieldData,
			GaussIntegralFaceO1Grid<fO> *gaussIntegralFace); //

		virtual bool getTimeStep(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			faceFieldDataVector *faceFieldData);

		virtual bool excuteInnerLoopExplicitSSPRKStep3( // time marching order 3
			const int iStep,
			cellFieldDataVector *cellFieldData); //��չʱoverride ��дtimemarch
		virtual bool updateArtificalViscosity(
			const int iStep,
			const real refValue,
			cellFieldDataVector *cellFieldData);

		virtual bool excuteFilter(
			parameter *parameter,
			cellFieldDataVector *cellFieldData); //

		virtual bool excuteWBAPLimiter4th(
			parameter *parameter,
			cellFieldDataVector *cellFieldData); //

		virtual bool excuteWBAPLimiter4thCR(
			parameter *parameter,
			cellFieldDataVector *cellFieldData); //

		virtual bool excuteInnerLoopExplicitSSPRKStep4( // time marching order 4
			const int iStep,
			cellFieldDataVector *cellFieldData); //��չʱoverride ��дtimemarch
		virtual bool excuteInnerLoopImplicitSDIRK()
		{
			return true;
		};

		virtual bool applyBoundaryLimiter( // time marching order 4
			parameter *parameter,
			cellFieldDataVector *cellFieldData); //��չʱoverride ��дtimemarch
	};

	//------------------------------------------------------------------------------------------
	template <unsigned int O>
	bool evolution<O>::getTimeStep(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		faceFieldDataVector *faceFieldData)
	{

		//		std::cout << "computing time step..." << std::endl;
		cellFieldDataVector::iterator iterCellFieldData_;
		cellFieldDataVector::iterator iterCellFieldData;
		faceFieldDataVector::iterator iterFaceFieldData;

		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			(*iterCellFieldData).lambdaCell = 0.0;
		}
		for (iterFaceFieldData = faceFieldData->begin(); iterFaceFieldData != faceFieldData->end(); ++iterFaceFieldData)
		{
			// use low order
			int cl = (*iterFaceFieldData).faceCellIndex[1];
			int cr = (*iterFaceFieldData).faceCellIndex[2];
			point p1 = (*iterFaceFieldData).faceNode[1].second;
			point p2 = (*iterFaceFieldData).faceNode[2].second;
			point pMid = 0.5 * (p1 + p2);
			point normalVector = (*iterFaceFieldData).normalVector;
			point uNV = (1.0 / normalVector.length()) * normalVector;
			// diff base matrix i
			tensor2D<real, NDOFS, NDOFS> matrixDiffBaseI;
			// base moment i
			tensor1D<real, NDOFS> momentI;
			iterCellFieldData_ = cellFieldData->begin() + cl - 1;
			point baryCenterI = (*iterCellFieldData_).baryCenter;
			point scaleI = (*iterCellFieldData_).lengthReference;
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
			}
			CfvMath::getDiffBaseValue(
				pMid,
				baryCenterI,
				scaleI,
				momentI,
				matrixDiffBaseI);
			real dxL = 0.0;
			real dyL = 0.0;
			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
			{
				dxL += matrixDiffBaseI[kk][1] * (*iterCellFieldData_).scalarVariableTn[kk]; // based on VR
				dyL += matrixDiffBaseI[kk][2] * (*iterCellFieldData_).scalarVariableTn[kk];
			}
			point uL;
			uL.x = dxL;
			uL.y = dyL;
			if (cr > 0)
			{
				// diff base matrix j
				tensor2D<real, NDOFS, NDOFS> matrixDiffBaseJ;
				// base moment j
				tensor1D<real, NDOFS> momentJ;
				iterCellFieldData_ = cellFieldData->begin() + cr - 1;
				point baryCenterJ = (*iterCellFieldData_).baryCenter;
				point scaleJ = (*iterCellFieldData_).lengthReference;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					momentJ[kk] = (*iterCellFieldData_).baseMoment[kk];
				}
				CfvMath::getDiffBaseValue(
					pMid + (*iterFaceFieldData).sideOff, //���ڱ߽�����
					baryCenterJ,
					scaleJ,
					momentJ,
					matrixDiffBaseJ);
				real dxR = 0.0;
				real dyR = 0.0;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					dxR += matrixDiffBaseJ[kk][1] * (*iterCellFieldData_).scalarVariableTn[kk]; // based on VR
					dyR += matrixDiffBaseJ[kk][2] * (*iterCellFieldData_).scalarVariableTn[kk];
				}
				point uR;
				uR.x = dxR;
				uR.y = dyR;
				point uMid = 0.5 * (uL + uR);
				(*iterFaceFieldData).lambdaFace = (*iterFaceFieldData).area * (std::fabs(getInnerProduct(uMid, uNV)) + 1.0);
			}
			else
			{
				point uR;
				uR.setZero();
				point uMid = 0.5 * (uL + uR);
				(*iterFaceFieldData).lambdaFace = (*iterFaceFieldData).area * (std::fabs(getInnerProduct(uMid, uNV)) + 1.0);
			}

			iterCellFieldData_ = cellFieldData->begin() + cl - 1;
			(*iterCellFieldData_).lambdaCell += (*iterFaceFieldData).lambdaFace;
			if (cr > 0)
			{
				iterCellFieldData_ = cellFieldData->begin() + cr - 1;
				(*iterCellFieldData_).lambdaCell += (*iterFaceFieldData).lambdaFace;
			}
		}
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			(*iterCellFieldData).dTPhysical = parameter->CFL * (*iterCellFieldData).volume / ((*iterCellFieldData).lambdaCell + parameter->EPS);
			(*iterCellFieldData).dTPseudo = (*iterCellFieldData).dTPhysical;

			// debug
			//(*iterCellFieldData).dTPhysical = 1e-2;
			//(*iterCellFieldData).dTPseudo = 1e-2;
		}
		if (parameter->isLocalTimeMarching == false)
		{
			// real dTMin = 1e6;
			parameter->refDeltaTn = 1e6;
			for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
			{
				parameter->refDeltaTn = ((*iterCellFieldData).dTPseudo < parameter->refDeltaTn) ? (*iterCellFieldData).dTPseudo : parameter->refDeltaTn;
			}
			for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
			{
				(*iterCellFieldData).dTPhysical = parameter->refDeltaTn;
				(*iterCellFieldData).dTPseudo = parameter->refDeltaTn;
			}
		}
		//		std::cout << " ..time step has been computed." << std::endl;
		return true;
	}

	//------------------------------------------------------------------------------------------

	template <unsigned int O>
	bool evolution<O>::excuteFilter(
		parameter *parameter,
		cellFieldDataVector *cellFieldData)
	{

		real ratio = 1.0;
		real filterWeightI;
		real filterWeightJ;
		real sumLocal;
		cellFieldDataVector::iterator iterCellFieldData;
		cellFieldDataVector::iterator iterCellFieldData_;
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			if ((*iterCellFieldData).boundaryCellType_ == FarFieldCell)
			{
				filterWeightI = ratio / ((*iterCellFieldData).cellCellNumber + ratio);
				filterWeightJ = 1.0 / ((*iterCellFieldData).cellCellNumber + ratio);
				sumLocal = filterWeightI * (*iterCellFieldData).scalarVariableTn[0];
				//(*iterCellFieldData).phiFilter
				for (int ii = 1; ii <= (*iterCellFieldData).cellCellNumber; ++ii)
				{
					int cr = (*iterCellFieldData).cellCellIndex[ii];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					sumLocal += filterWeightJ * (*iterCellFieldData_).scalarVariableTn[0];
				}
				(*iterCellFieldData).phiFilter = sumLocal;
			}
		}

		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			//||
			//	(*iterCellFieldData).boundaryCellType_ == SymmetricCell

			if ((*iterCellFieldData).boundaryCellType_ == FarFieldCell)
			{
				(*iterCellFieldData).scalarVariableTn[0] = (*iterCellFieldData).phiFilter;
			}
		}
		return true;
	}

	//------------------------------------------------------------------------------------------

	template <unsigned int O>
	bool evolution<O>::excuteWBAPLimiter4th(
		parameter *parameter,
		cellFieldDataVector *cellFieldData)
	{

		real dPhi1[5]; // 4 + 1   0<->i
		real dPhi2[5];

		real dPhi3[5];
		real dPhi4[5];
		real dPhi5[5];

		real dPhi6[5];
		real dPhi7[5];
		real dPhi8[5];
		real dPhi9[5];

		cellFieldDataVector::iterator iterCellFieldData;
		cellFieldDataVector::iterator iterCellFieldData_;
		// 4th
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			dPhi6[0] = (*iterCellFieldData).scalarVariableTn[6];
			dPhi7[0] = (*iterCellFieldData).scalarVariableTn[7];
			dPhi8[0] = (*iterCellFieldData).scalarVariableTn[8];
			dPhi9[0] = (*iterCellFieldData).scalarVariableTn[9];
			for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber; ++ii)
			{
				int cr = (*iterCellFieldData).cellCellIndex[ii];
				iterCellFieldData_ = cellFieldData->begin() + cr - 1;
				real rx = (*iterCellFieldData).lengthReference.x / (*iterCellFieldData_).lengthReference.x;
				real ry = (*iterCellFieldData).lengthReference.y / (*iterCellFieldData_).lengthReference.y;
				dPhi6[ii] = (*iterCellFieldData_).scalarVariableTn[6] * std::pow(rx, 3);
				dPhi7[ii] = (*iterCellFieldData_).scalarVariableTn[7] * std::pow(rx, 2) * ry;
				dPhi8[ii] = (*iterCellFieldData_).scalarVariableTn[8] * rx * std::pow(ry, 2);
				dPhi9[ii] = (*iterCellFieldData_).scalarVariableTn[9] * std::pow(ry, 3);
			}

			(*iterCellFieldData).scalarVariableTnLimited[6] = CfvMath::W12(dPhi6, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnLimited[7] = CfvMath::W12(dPhi7, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnLimited[8] = CfvMath::W12(dPhi8, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnLimited[9] = CfvMath::W12(dPhi9, (*iterCellFieldData).cellCellNumber - 1);
		}

		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			(*iterCellFieldData).scalarVariableTn[6] = (*iterCellFieldData).scalarVariableTnLimited[6];
			(*iterCellFieldData).scalarVariableTn[7] = (*iterCellFieldData).scalarVariableTnLimited[7];
			(*iterCellFieldData).scalarVariableTn[8] = (*iterCellFieldData).scalarVariableTnLimited[8];
			(*iterCellFieldData).scalarVariableTn[9] = (*iterCellFieldData).scalarVariableTnLimited[9];
		}
		// 3rd
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			dPhi3[0] = (*iterCellFieldData).scalarVariableTn[3];
			dPhi4[0] = (*iterCellFieldData).scalarVariableTn[4];
			dPhi5[0] = (*iterCellFieldData).scalarVariableTn[5];
			for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber; ++ii)
			{
				int cr = (*iterCellFieldData).cellCellIndex[ii];
				iterCellFieldData_ = cellFieldData->begin() + cr - 1;
				real rx = (*iterCellFieldData).lengthReference.x / (*iterCellFieldData_).lengthReference.x;
				real ry = (*iterCellFieldData).lengthReference.y / (*iterCellFieldData_).lengthReference.y;
				real deltaX = ((*iterCellFieldData).baryCenter.x - ((*iterCellFieldData_).baryCenter.x + (*iterCellFieldData).cellFaceSideOff[ii].x)) / (*iterCellFieldData_).lengthReference.x;
				real deltaY = ((*iterCellFieldData).baryCenter.y - ((*iterCellFieldData_).baryCenter.y + (*iterCellFieldData).cellFaceSideOff[ii].y)) / (*iterCellFieldData_).lengthReference.y;
				dPhi3[ii] = ((*iterCellFieldData_).scalarVariableTn[3] + 3.0 * (*iterCellFieldData_).scalarVariableTn[6] * deltaX + (*iterCellFieldData_).scalarVariableTn[7] * deltaY) * std::pow(rx, 2);
				dPhi4[ii] = ((*iterCellFieldData_).scalarVariableTn[4] + 2.0 * (*iterCellFieldData_).scalarVariableTn[7] * deltaX + 2.0 * (*iterCellFieldData_).scalarVariableTn[8] * deltaY) * rx * ry;
				dPhi5[ii] = ((*iterCellFieldData_).scalarVariableTn[5] + (*iterCellFieldData_).scalarVariableTn[8] * deltaX + 3.0 * (*iterCellFieldData_).scalarVariableTn[9] * deltaY) * std::pow(ry, 2);
			}

			(*iterCellFieldData).scalarVariableTnLimited[3] = CfvMath::W12(dPhi3, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnLimited[4] = CfvMath::W12(dPhi4, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnLimited[5] = CfvMath::W12(dPhi5, (*iterCellFieldData).cellCellNumber - 1);
		}

		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			(*iterCellFieldData).scalarVariableTn[3] = (*iterCellFieldData).scalarVariableTnLimited[3];
			(*iterCellFieldData).scalarVariableTn[4] = (*iterCellFieldData).scalarVariableTnLimited[4];
			(*iterCellFieldData).scalarVariableTn[5] = (*iterCellFieldData).scalarVariableTnLimited[5];
		}
		// 2nd
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			dPhi1[0] = (*iterCellFieldData).scalarVariableTn[1];
			dPhi2[0] = (*iterCellFieldData).scalarVariableTn[2];
			for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber; ++ii)
			{
				int cr = (*iterCellFieldData).cellCellIndex[ii];
				iterCellFieldData_ = cellFieldData->begin() + cr - 1;
				real rx = (*iterCellFieldData).lengthReference.x / (*iterCellFieldData_).lengthReference.x;
				real ry = (*iterCellFieldData).lengthReference.y / (*iterCellFieldData_).lengthReference.y;
				real deltaX = ((*iterCellFieldData).baryCenter.x - ((*iterCellFieldData_).baryCenter.x + (*iterCellFieldData).cellFaceSideOff[ii].x)) / (*iterCellFieldData_).lengthReference.x;
				real deltaY = ((*iterCellFieldData).baryCenter.y - ((*iterCellFieldData_).baryCenter.y + (*iterCellFieldData).cellFaceSideOff[ii].y)) / (*iterCellFieldData_).lengthReference.y;
				dPhi1[ii] = ((*iterCellFieldData_).scalarVariableTn[1] + 2.0 * (*iterCellFieldData_).scalarVariableTn[3] * deltaX + (*iterCellFieldData_).scalarVariableTn[4] * deltaY + 3.0 * (*iterCellFieldData_).scalarVariableTn[6] * std::pow(deltaX, 2) + 2.0 * (*iterCellFieldData_).scalarVariableTn[7] * deltaX * deltaY + (*iterCellFieldData_).scalarVariableTn[8] * std::pow(deltaY, 2)) * rx;

				dPhi2[ii] = ((*iterCellFieldData_).scalarVariableTn[2] + 2.0 * (*iterCellFieldData_).scalarVariableTn[5] * deltaY + (*iterCellFieldData_).scalarVariableTn[4] * deltaX + (*iterCellFieldData_).scalarVariableTn[7] * std::pow(deltaX, 2) + 2.0 * (*iterCellFieldData_).scalarVariableTn[8] * deltaX * deltaY + 3.0 * (*iterCellFieldData_).scalarVariableTn[9] * std::pow(deltaY, 2)) * ry;
			}
			(*iterCellFieldData).scalarVariableTnLimited[1] = CfvMath::W12(dPhi1, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnLimited[2] = CfvMath::W12(dPhi2, (*iterCellFieldData).cellCellNumber - 1);
		}

		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			(*iterCellFieldData).scalarVariableTn[1] = (*iterCellFieldData).scalarVariableTnLimited[1];
			(*iterCellFieldData).scalarVariableTn[2] = (*iterCellFieldData).scalarVariableTnLimited[2];
		}

		return true;
	}

	template <unsigned int O>
	bool evolution<O>::excuteWBAPLimiter4thCR(
		parameter *parameter,
		cellFieldDataVector *cellFieldData)
	{

		real dPhi1[5]; // 4 + 1   0<->i
		real dPhi2[5];

		real dPhi3[5];
		real dPhi4[5];
		real dPhi5[5];

		real dPhi6[5];
		real dPhi7[5];
		real dPhi8[5];
		real dPhi9[5];

		cellFieldDataVector::iterator iterCellFieldData;
		cellFieldDataVector::iterator iterCellFieldData_;
		// 4th
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			dPhi6[0] = (*iterCellFieldData).scalarVariableTnCR[6];
			dPhi7[0] = (*iterCellFieldData).scalarVariableTnCR[7];
			dPhi8[0] = (*iterCellFieldData).scalarVariableTnCR[8];
			dPhi9[0] = (*iterCellFieldData).scalarVariableTnCR[9];
			for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber; ++ii)
			{
				int cr = (*iterCellFieldData).cellCellIndex[ii];
				iterCellFieldData_ = cellFieldData->begin() + cr - 1;
				real rx = (*iterCellFieldData).lengthReference.x / (*iterCellFieldData_).lengthReference.x;
				real ry = (*iterCellFieldData).lengthReference.y / (*iterCellFieldData_).lengthReference.y;
				dPhi6[ii] = (*iterCellFieldData_).scalarVariableTnCR[6] * std::pow(rx, 3);
				dPhi7[ii] = (*iterCellFieldData_).scalarVariableTnCR[7] * std::pow(rx, 2) * ry;
				dPhi8[ii] = (*iterCellFieldData_).scalarVariableTnCR[8] * rx * std::pow(ry, 2);
				dPhi9[ii] = (*iterCellFieldData_).scalarVariableTnCR[9] * std::pow(ry, 3);
			}

			(*iterCellFieldData).scalarVariableTnCRLimited[6] = CfvMath::W12(dPhi6, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnCRLimited[7] = CfvMath::W12(dPhi7, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnCRLimited[8] = CfvMath::W12(dPhi8, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnCRLimited[9] = CfvMath::W12(dPhi9, (*iterCellFieldData).cellCellNumber - 1);
		}

		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			(*iterCellFieldData).scalarVariableTnCR[6] = (*iterCellFieldData).scalarVariableTnCRLimited[6];
			(*iterCellFieldData).scalarVariableTnCR[7] = (*iterCellFieldData).scalarVariableTnCRLimited[7];
			(*iterCellFieldData).scalarVariableTnCR[8] = (*iterCellFieldData).scalarVariableTnCRLimited[8];
			(*iterCellFieldData).scalarVariableTnCR[9] = (*iterCellFieldData).scalarVariableTnCRLimited[9];
		}
		// 3rd
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			dPhi3[0] = (*iterCellFieldData).scalarVariableTnCR[3];
			dPhi4[0] = (*iterCellFieldData).scalarVariableTnCR[4];
			dPhi5[0] = (*iterCellFieldData).scalarVariableTnCR[5];
			for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber; ++ii)
			{
				int cr = (*iterCellFieldData).cellCellIndex[ii];
				iterCellFieldData_ = cellFieldData->begin() + cr - 1;
				real rx = (*iterCellFieldData).lengthReference.x / (*iterCellFieldData_).lengthReference.x;
				real ry = (*iterCellFieldData).lengthReference.y / (*iterCellFieldData_).lengthReference.y;
				real deltaX = ((*iterCellFieldData).baryCenter.x - ((*iterCellFieldData_).baryCenter.x + (*iterCellFieldData).cellFaceSideOff[ii].x)) / (*iterCellFieldData_).lengthReference.x;
				real deltaY = ((*iterCellFieldData).baryCenter.y - ((*iterCellFieldData_).baryCenter.y + (*iterCellFieldData).cellFaceSideOff[ii].y)) / (*iterCellFieldData_).lengthReference.y;
				dPhi3[ii] = ((*iterCellFieldData_).scalarVariableTnCR[3] + 3.0 * (*iterCellFieldData_).scalarVariableTnCR[6] * deltaX + (*iterCellFieldData_).scalarVariableTnCR[7] * deltaY) * std::pow(rx, 2);
				dPhi4[ii] = ((*iterCellFieldData_).scalarVariableTnCR[4] + 2.0 * (*iterCellFieldData_).scalarVariableTnCR[7] * deltaX + 2.0 * (*iterCellFieldData_).scalarVariableTnCR[8] * deltaY) * rx * ry;
				dPhi5[ii] = ((*iterCellFieldData_).scalarVariableTnCR[5] + (*iterCellFieldData_).scalarVariableTnCR[8] * deltaX + 3.0 * (*iterCellFieldData_).scalarVariableTnCR[9] * deltaY) * std::pow(ry, 2);
			}

			(*iterCellFieldData).scalarVariableTnCRLimited[3] = CfvMath::W12(dPhi3, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnCRLimited[4] = CfvMath::W12(dPhi4, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnCRLimited[5] = CfvMath::W12(dPhi5, (*iterCellFieldData).cellCellNumber - 1);
		}

		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			(*iterCellFieldData).scalarVariableTnCR[3] = (*iterCellFieldData).scalarVariableTnCRLimited[3];
			(*iterCellFieldData).scalarVariableTnCR[4] = (*iterCellFieldData).scalarVariableTnCRLimited[4];
			(*iterCellFieldData).scalarVariableTnCR[5] = (*iterCellFieldData).scalarVariableTnCRLimited[5];
		}
		// 2nd
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			dPhi1[0] = (*iterCellFieldData).scalarVariableTnCR[3];
			dPhi2[0] = (*iterCellFieldData).scalarVariableTnCR[4];
			for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber; ++ii)
			{
				int cr = (*iterCellFieldData).cellCellIndex[ii];
				iterCellFieldData_ = cellFieldData->begin() + cr - 1;
				real rx = (*iterCellFieldData).lengthReference.x / (*iterCellFieldData_).lengthReference.x;
				real ry = (*iterCellFieldData).lengthReference.y / (*iterCellFieldData_).lengthReference.y;
				real deltaX = ((*iterCellFieldData).baryCenter.x - ((*iterCellFieldData_).baryCenter.x + (*iterCellFieldData).cellFaceSideOff[ii].x)) / (*iterCellFieldData_).lengthReference.x;
				real deltaY = ((*iterCellFieldData).baryCenter.y - ((*iterCellFieldData_).baryCenter.y + (*iterCellFieldData).cellFaceSideOff[ii].y)) / (*iterCellFieldData_).lengthReference.y;
				dPhi1[ii] = ((*iterCellFieldData_).scalarVariableTnCR[1] + 2.0 * (*iterCellFieldData_).scalarVariableTnCR[3] * deltaX + (*iterCellFieldData_).scalarVariableTnCR[4] * deltaY + 3.0 * (*iterCellFieldData_).scalarVariableTnCR[6] * std::pow(deltaX, 2) + 2.0 * (*iterCellFieldData_).scalarVariableTnCR[7] * deltaX * deltaY + (*iterCellFieldData_).scalarVariableTnCR[8] * std::pow(deltaY, 2)) * rx;

				dPhi2[ii] = ((*iterCellFieldData_).scalarVariableTnCR[2] + 2.0 * (*iterCellFieldData_).scalarVariableTnCR[5] * deltaY + (*iterCellFieldData_).scalarVariableTnCR[4] * deltaX + (*iterCellFieldData_).scalarVariableTnCR[7] * std::pow(deltaX, 2) + 2.0 * (*iterCellFieldData_).scalarVariableTnCR[8] * deltaX * deltaY + 3.0 * (*iterCellFieldData_).scalarVariableTnCR[9] * std::pow(deltaY, 2)) * ry;
			}
			(*iterCellFieldData).scalarVariableTnCRLimited[1] = CfvMath::W12(dPhi1, (*iterCellFieldData).cellCellNumber - 1);
			(*iterCellFieldData).scalarVariableTnCRLimited[2] = CfvMath::W12(dPhi2, (*iterCellFieldData).cellCellNumber - 1);
		}

		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			(*iterCellFieldData).scalarVariableTnCR[1] = (*iterCellFieldData).scalarVariableTnCRLimited[1];
			(*iterCellFieldData).scalarVariableTnCR[2] = (*iterCellFieldData).scalarVariableTnCRLimited[2];
		}

		return true;
	}

	//------------------------------------------------------------------------------------------

	template <unsigned int O>
	bool evolution<O>::getSourceTerm(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		faceFieldDataVector *faceFieldData,
		GaussIntegralCellO1Grid<vO> *gaussIntegralCell)
	{

		//		std::cout << "computing source term..." << std::endl;
		cellFieldDataVector::iterator iterCellFieldData;
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			real *fI = new real[static_cast<int>((*iterCellFieldData).PG)];
			real *weight = new real[static_cast<int>((*iterCellFieldData).PG)];
			real *cofJacobi = new real[static_cast<int>((*iterCellFieldData).PG)];
			real parametricVolume = (*iterCellFieldData).parametricVolume;
			// get the values on Gauss points
			for (int gg = 0; gg < static_cast<int>((*iterCellFieldData).PG); ++gg)
			{
				// diff base matrix i
				tensor2D<real, NDOFS, NDOFS> matrixDiffBaseI;
				// base moment i
				tensor1D<real, NDOFS> momentI;
				point p = (*iterCellFieldData).gaussPairVector_[gg].p;
				point baryCenterI = (*iterCellFieldData).baryCenter;
				point scaleI = (*iterCellFieldData).lengthReference;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					momentI[kk] = (*iterCellFieldData).baseMoment[kk];
				}
				CfvMath::getDiffBaseValue(
					p,
					baryCenterI,
					scaleI,
					momentI,
					matrixDiffBaseI);
				real phi = (*iterCellFieldData).scalarVariableTn[0];
				real tempX, tempY;
				point dPhi, dPhiTilde;
				real d2PhiX, d2PhiY;
				// phi
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					phi += matrixDiffBaseI[kk][0] * (*iterCellFieldData).scalarVariableTn[kk];
				}
				// dphi
				tempX = 0.0;
				tempY = 0.0;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					tempX += matrixDiffBaseI[kk][1] * (*iterCellFieldData).scalarVariableTn[kk];
					tempY += matrixDiffBaseI[kk][2] * (*iterCellFieldData).scalarVariableTn[kk];
				}
				dPhi.x = tempX;
				dPhi.y = tempY;
				// dPhiTilde
				tempX = 0.0;
				tempY = 0.0;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					tempX += matrixDiffBaseI[kk][1] * (*iterCellFieldData).scalarVariableTnCR[kk];
					tempY += matrixDiffBaseI[kk][2] * (*iterCellFieldData).scalarVariableTnCR[kk];
				}
				dPhiTilde.x = tempX;
				dPhiTilde.y = tempY;
				// d2phi
				tempX = 0.0;
				tempY = 0.0;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					tempX += matrixDiffBaseI[kk][3] * (*iterCellFieldData).scalarVariableTn[kk];
					tempY += matrixDiffBaseI[kk][5] * (*iterCellFieldData).scalarVariableTn[kk];
				}
				d2PhiX = tempX;
				d2PhiY = tempY;
				// fI
				//	fI[gg] = getMax(0.0, 1.0 - std::fabs(getInnerProduct(dPhi, dPhiTilde)) + (*iterCellFieldData).artificialViscousParameter * phi * (d2PhiX + d2PhiY));
				//	fI[gg] = 1.0 - std::fabs(getInnerProduct(dPhi, dPhiTilde)) + (*iterCellFieldData).artificialViscousParameter * phi * (d2PhiX + d2PhiY);

				//----------------------
				// 20210407
				fI[gg] = 1.0 - std::fabs(getInnerProduct(dPhi, dPhiTilde));
				// fI[gg] = 1.0 - std::fabs(getInnerProduct(dPhi, dPhi));

				weight[gg] = (*iterCellFieldData).parametricValue[gg].second;
				cofJacobi[gg] = (*iterCellFieldData).gaussPairVector_[gg].JacobiCof;
			}
			// get the cell Gauss integrals
			real result;
			gaussIntegralCell->getIntegral(
				static_cast<int>((*iterCellFieldData).PG),
				fI,
				weight,
				cofJacobi,
				parametricVolume,
				result);
			(*iterCellFieldData).timeMarchingRHSTn = result / (*iterCellFieldData).volume;

			delete[] fI;
			delete[] weight;
			delete[] cofJacobi;
			fI = NULL;
			weight = NULL;
			cofJacobi = NULL;
		}
		//		std::cout << " ..source term has been computed." << std::endl;
		return true;
	}

	template <unsigned int O>
	bool evolution<O>::getSourceTerm(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		faceFieldDataVector *faceFieldData,
		GaussIntegralCellO2Grid<vO> *gaussIntegralCell)
	{

		//		std::cout << "computing source term..." << std::endl;
		cellFieldDataVector::iterator iterCellFieldData;
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			real *fI = new real[static_cast<int>((*iterCellFieldData).PG)];
			real *weight = new real[static_cast<int>((*iterCellFieldData).PG)];
			real *cofJacobi = new real[static_cast<int>((*iterCellFieldData).PG)];
			real parametricVolume = (*iterCellFieldData).parametricVolume;
			// get the values on Gauss points
			for (int gg = 0; gg < static_cast<int>((*iterCellFieldData).PG); ++gg)
			{
				// diff base matrix i
				tensor2D<real, NDOFS, NDOFS> matrixDiffBaseI;
				// base moment i
				tensor1D<real, NDOFS> momentI;
				point p = (*iterCellFieldData).gaussPairVector_[gg].p;
				point baryCenterI = (*iterCellFieldData).baryCenter;
				point scaleI = (*iterCellFieldData).lengthReference;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					momentI[kk] = (*iterCellFieldData).baseMoment[kk];
				}
				CfvMath::getDiffBaseValue(
					p,
					baryCenterI,
					scaleI,
					momentI,
					matrixDiffBaseI);
				real phi = (*iterCellFieldData).scalarVariableTn[0];
				real tempX, tempY;
				point dPhi, dPhiTilde;
				real d2PhiX, d2PhiY;
				// phi
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					phi += matrixDiffBaseI[kk][0] * (*iterCellFieldData).scalarVariableTn[kk];
				}
				// dphi
				tempX = 0.0;
				tempY = 0.0;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					tempX += matrixDiffBaseI[kk][1] * (*iterCellFieldData).scalarVariableTn[kk];
					tempY += matrixDiffBaseI[kk][2] * (*iterCellFieldData).scalarVariableTn[kk];
				}
				dPhi.x = tempX;
				dPhi.y = tempY;
				// dPhiTilde
				tempX = 0.0;
				tempY = 0.0;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					tempX += matrixDiffBaseI[kk][1] * (*iterCellFieldData).scalarVariableTnCR[kk];
					tempY += matrixDiffBaseI[kk][2] * (*iterCellFieldData).scalarVariableTnCR[kk];
				}
				dPhiTilde.x = tempX;
				dPhiTilde.y = tempY;
				// d2phi
				tempX = 0.0;
				tempY = 0.0;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					tempX += matrixDiffBaseI[kk][3] * (*iterCellFieldData).scalarVariableTn[kk];
					tempY += matrixDiffBaseI[kk][5] * (*iterCellFieldData).scalarVariableTn[kk];
				}
				d2PhiX = tempX;
				d2PhiY = tempY;
				// fI
				fI[gg] = 1.0 - std::fabs(getInnerProduct(dPhi, dPhiTilde)) + (*iterCellFieldData).artificialViscousParameter * phi * (d2PhiX + d2PhiY);
				weight[gg] = (*iterCellFieldData).parametricValue[gg].second;
				cofJacobi[gg] = (*iterCellFieldData).gaussPairVector_[gg].JacobiCof;
			}
			// get the cell Gauss integrals
			real result;
			gaussIntegralCell->getIntegral(
				static_cast<int>((*iterCellFieldData).PG),
				fI,
				weight,
				cofJacobi,
				parametricVolume,
				result);
			(*iterCellFieldData).timeMarchingRHSTn = result / (*iterCellFieldData).volume;

			delete[] fI;
			delete[] weight;
			delete[] cofJacobi;
			fI = NULL;
			weight = NULL;
			cofJacobi = NULL;
		}
		//		std::cout << " ..source term has been computed." << std::endl;
		return true;
	}

	//------------------------------------------------------------------------------------------

	template <unsigned int O>
	bool evolution<O>::excuteInnerLoopExplicitSSPRKStep3(
		const int iStep,
		cellFieldDataVector *cellFieldData)
	{
		// cellFieldDataVector::iterator iterCellFieldData;
		real cof3[4] = {0.0, 1.0, 0.25, 2.0 / 3.0};
		if (iStep == 1)
		{
#pragma omp parallel for schedule(static)
			for (int iCell = 0; iCell < cellFieldData->size(); iCell++)
			{
				(*cellFieldData)[iCell].scalarVariableTm[0] = (*cellFieldData)[iCell].scalarVariableTn[0];
			}
		}
#pragma omp parallel for schedule(static)
		for (int iCell = 0; iCell < cellFieldData->size(); iCell++)
		{
			(*cellFieldData)[iCell].scalarVariableTn[0] =
				(1.0 - cof3[iStep]) * (*cellFieldData)[iCell].scalarVariableTm[0] + cof3[iStep] * ((*cellFieldData)[iCell].scalarVariableTn[0] + (*cellFieldData)[iCell].dTPhysical * (*cellFieldData)[iCell].timeMarchingRHSTn);
			(*cellFieldData)[iCell].timeMarchingRHSRK[iStep] = (*cellFieldData)[iCell].timeMarchingRHSTn;
		}

		// for (int iCell = 0; iCell < cellFieldData->size(); iCell++)
		// {
			
		// }

		return true;
	}

	template <unsigned int O>
	bool evolution<O>::excuteInnerLoopExplicitSSPRKStep4(
		const int iStep,
		cellFieldDataVector *cellFieldData)
	{
		// cellFieldDataVector::iterator iterCellFieldData;
		real cof4[5] = {0.0, 0.5, 0.5, 1.0, 1.0 / 6.0};
		if (iStep == 1)
		{
#pragma omp parallel for schedule(static)
			for (int iCell = 0; iCell < cellFieldData->size(); iCell++)
			{
				(*cellFieldData)[iCell].scalarVariableTm[0] = (*cellFieldData)[iCell].scalarVariableTn[0];
			}
		}
		if (iStep < 4)
		{
#pragma omp parallel for schedule(static)
			for (int iCell = 0; iCell < cellFieldData->size(); iCell++)
			{
				(*cellFieldData)[iCell].scalarVariableTn[0] =
					(1.0 - cof4[iStep]) * (*cellFieldData)[iCell].scalarVariableTm[0] + cof4[iStep] * ((*cellFieldData)[iCell].scalarVariableTn[0] + (*cellFieldData)[iCell].dTPhysical * (*cellFieldData)[iCell].timeMarchingRHSTn);
				(*cellFieldData)[iCell].timeMarchingRHSRK[iStep] = (*cellFieldData)[iCell].timeMarchingRHSTn;
			}
			// for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
			// {
				
			// }
		}
		else if (iStep == 4)
		{
#pragma omp parallel for schedule(static)
			for (int iCell = 0; iCell < cellFieldData->size(); iCell++)
			{
				(*cellFieldData)[iCell].scalarVariableTn[0] =
					(*cellFieldData)[iCell].scalarVariableTm[0] + cof4[iStep] * (*cellFieldData)[iCell].dTPhysical * ((*cellFieldData)[iCell].timeMarchingRHSTn			   // current result,	<-> (4)
																													  + (*cellFieldData)[iCell].timeMarchingRHSRK[1]	   // Tm result,		<-> (1)
																													  + 2.0 * (*cellFieldData)[iCell].timeMarchingRHSRK[2] // Tm+(1) result,	<-> (2)
																													  + 2.0 * (*cellFieldData)[iCell].timeMarchingRHSRK[3] // Tm+(2) result,	<-> (3)
																													 );
			}
		}
		return true;
	}

	//------------------------------------------------------------------------------------------

	template <unsigned int O>
	bool evolution<O>::updateArtificalViscosity(
		const int iStep,
		const real refValue,
		cellFieldDataVector *cellFieldData)
	{

		cellFieldDataVector::iterator iterCellFieldData;
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			if (iStep < 0)
			{
				(*iterCellFieldData).artificialViscousParameter = refValue * 20.0;
			}
			else
			{
				// diff base matrix i
				tensor2D<real, NDOFS, NDOFS> matrixDiffBaseI;
				// base moment i
				tensor1D<real, NDOFS> momentI;
				point p = (*iterCellFieldData).baryCenter;
				point baryCenterI = (*iterCellFieldData).baryCenter;
				point scaleI = (*iterCellFieldData).lengthReference;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					momentI[kk] = (*iterCellFieldData).baseMoment[kk];
				}
				CfvMath::getDiffBaseValue(
					p,
					baryCenterI,
					scaleI,
					momentI,
					matrixDiffBaseI);
				// real phi = (*iterCellFieldData).scalarVariableTn[0];
				point dPhi;
				dPhi.setZero();
				real tempX, tempY;
				tempX = 0.0;
				tempY = 0.0;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					tempX += matrixDiffBaseI[kk][1] * (*iterCellFieldData).scalarVariableTn[kk];
					tempY += matrixDiffBaseI[kk][2] * (*iterCellFieldData).scalarVariableTn[kk];
				}
				dPhi.x = tempX;
				dPhi.y = tempY;
				//(x)<(y)?(x):(y)
				if (std::fabs(1.0 - std::pow(dPhi.length(), 2)) > 0.1 && std::fabs(1.0 - std::pow(dPhi.length(), 2)) < 1.0)
				{
					(*iterCellFieldData).artificialViscousParameter = refValue * 10.0;
				}
				else if (std::fabs(1.0 - std::pow(dPhi.length(), 2)) >= 1.0)
				{
					(*iterCellFieldData).artificialViscousParameter = refValue * 20.0;
				}
				else
				{
					(*iterCellFieldData).artificialViscousParameter =
						refValue < std::fabs(1.0 - std::pow(dPhi.length(), 2)) ? refValue : std::fabs(1.0 - std::pow(dPhi.length(), 2));
				}
			}
		}
		return true;
	}

	template <unsigned int O>
	bool evolution<O>::getArtificialViscosityTerm(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		faceFieldDataVector *faceFieldData,
		GaussIntegralFaceO1Grid<fO> *gaussIntegralFace)
	{

		cellFieldDataVector::iterator iterCell_;
		faceFieldDataVector::iterator iterFaceFieldData;

		// real nu = 0.008; un = 1
		real nu = 1.0; // 0.1;//0.02
		real un;
		real d;
		real SI;
		for (iterFaceFieldData = faceFieldData->begin(); iterFaceFieldData != faceFieldData->end(); ++iterFaceFieldData)
		{
			int cl = (*iterFaceFieldData).faceCellIndex[1];
			int cr = (*iterFaceFieldData).faceCellIndex[2];

			point *uNV = new point[static_cast<int>((*iterFaceFieldData).fPG)];
			real *fluxF = new real[static_cast<int>((*iterFaceFieldData).fPG)];
			real *weight = new real[static_cast<int>((*iterFaceFieldData).fPG)];
			real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData).fPG)];
			real parametricArea;
			// get the values on the Gauss points
			for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData).fPG); ++gg)
			{
				real qL, qnL, omegaL;
				real qR, qnR, omegaR;
				point dq;
				point dqL(0.0, 0.0);
				point dqR(0.0, 0.0);
				// cell L
				iterCell_ = cellFieldData->begin() + cl - 1;
				// base moment i
				tensor1D<real, NDOFS> momentI;
				// base i
				tensor2D<real, NDOFS, NDOFS> matrixDiffBaseI;
				point baryCenterI = (*iterCell_).baryCenter;
				point scaleI = (*iterCell_).lengthReference;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					momentI[kk] = (*iterCell_).baseMoment[kk];
				}
				point p = (*iterFaceFieldData).gaussPairVector_[gg].p;
				CfvMath::getDiffBaseValue(
					p,
					baryCenterI,
					scaleI,
					momentI,
					matrixDiffBaseI);
				qL = (*iterCell_).scalarVariableTn[0];
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					qL += matrixDiffBaseI[kk][0] * (*iterCell_).scalarVariableTn[kk];
					dqL.x += matrixDiffBaseI[kk][1] * (*iterCell_).scalarVariableTn[kk];
					dqL.y += matrixDiffBaseI[kk][2] * (*iterCell_).scalarVariableTn[kk];
				}

				omegaL = (*iterCell_).volume;
				//��λ��Ҫע��
				uNV[gg] = (1.0 / ((*iterFaceFieldData).gaussPairVector_[gg].normalVector.length())) * (*iterFaceFieldData).gaussPairVector_[gg].normalVector;

				qnL = getInnerProduct(dqL, uNV[gg]);
				SI = (*iterCell_).smoothIndicator;
				if (cr > 0)
				{
					// cell R
					iterCell_ = cellFieldData->begin() + cr - 1;
					// base moment j
					tensor1D<real, NDOFS> momentJ;
					// base j
					tensor2D<real, NDOFS, NDOFS> matrixDiffBaseJ;
					point baryCenterJ = (*iterCell_).baryCenter;
					point scaleJ = (*iterCell_).lengthReference;
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
					{
						momentJ[kk] = (*iterCell_).baseMoment[kk];
					}
					p = (*iterFaceFieldData).gaussPairVector_[gg].p;
					CfvMath::getDiffBaseValue(
						p + (*iterFaceFieldData).sideOff, //���ڱ߽�����
						baryCenterJ,
						scaleJ,
						momentJ,
						matrixDiffBaseJ);
					qR = (*iterCell_).scalarVariableTn[0];
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
					{
						qR += matrixDiffBaseJ[kk][0] * (*iterCell_).scalarVariableTn[kk];
						dqR.x += matrixDiffBaseJ[kk][1] * (*iterCell_).scalarVariableTn[kk];
						dqR.y += matrixDiffBaseJ[kk][2] * (*iterCell_).scalarVariableTn[kk];
						d = (baryCenterJ - baryCenterI).length();
					}
					omegaR = (*iterCell_).volume;
					dq = 0.5 * (dqL + dqR);
					qnR = getInnerProduct(dqR, uNV[gg]);
					SI += (*iterCell_).smoothIndicator;
					////------------------------------------

					un = d;
				}

				real qnL_abs = std::fabs(qnL);
				real qnR_abs = std::fabs(qnR);
				real ftemp;
				// 30P30N
				// if (SI < 0.1){
				//	ftemp = 0.05 * ScalarCfv::getMin(std::fabs(1.0 - dq.length()), 0.4 / std::pow((rO + 1), 2));
				// }
				// else{
				//	ftemp = 10.0 * 0.4 / std::pow((rO + 1), 2);
				// }

				ftemp = 10.0 * 0.4 / std::pow((rO + 1), 2);

				// 1 done
				// ftemp = 1.0 * 10.0 * 0.4 / std::pow((rO + 1), 2);

				// 3 done
				// ftemp = 1.0 / 4.0 * 10.0 * 0.4 / std::pow((rO + 1), 2);

				// 5 done
				// ftemp = 1.0 / 16.0 * 10.0 * 0.4 / std::pow((rO + 1), 2);

				// 7 done
				//  ftemp = 1.0 / 64.0 * 10.0 * 0.4 / std::pow((rO + 1), 2);

				// 9 done
				//  ftemp = 1.0 / 256.0 * 10.0 * 0.4 / std::pow((rO + 1), 2);

				fluxF[gg] = nu * un * ftemp * getInnerProduct(dq, uNV[gg]);

				weight[gg] = (*iterFaceFieldData).parametricValue[gg].second;
				cofJacobi[gg] = (*iterFaceFieldData).gaussPairVector_[gg].JacobiCof;
				parametricArea = (*iterFaceFieldData).parametricArea;
			}

			// get integral
			real *fI = new real[static_cast<int>((*iterFaceFieldData).fPG)];
			real result;
			for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData).fPG); ++gg)
			{
				fI[gg] = fluxF[gg];
			}
			gaussIntegralFace->getIntegral(
				static_cast<int>((*iterFaceFieldData).fPG),
				fI,
				weight,
				cofJacobi,
				parametricArea,
				result);

			// left cell
			iterCell_ = cellFieldData->begin() + cl - 1;
			//(*iterCell_).timeMarchingRHSTn += (result / (*iterCell_).volume);
			(*iterCell_).timeMarchingRHSTn += (result / (*iterCell_).volume);
			// right cell
			if (cr > 0 && (*iterFaceFieldData).faceProperty != Periodic)
			{
				iterCell_ = cellFieldData->begin() + cr - 1;
				//(*iterCell_).timeMarchingRHSTn -= (result / (*iterCell_).volume);
				(*iterCell_).timeMarchingRHSTn -= (result / (*iterCell_).volume);
			}
			delete[] fI;
			fI = NULL;

			delete[] uNV;
			delete[] cofJacobi;
			delete[] weight;
			delete[] fluxF;

			fluxF = NULL;
			uNV = NULL;
			cofJacobi = NULL;
			weight = NULL;
		}

		//--------------------------------------------------

		return true;
	}

	template <unsigned int O>
	bool evolution<O>::applyBoundaryLimiter( // time marching order 4
		parameter *parameter,
		cellFieldDataVector *cellFieldData)
	{
		real eps = 1e-6;
		cellFieldDataVector::iterator iterCell;
		for (iterCell = cellFieldData->begin(); iterCell != cellFieldData->end(); ++iterCell)
		{
			// if ((*iterCell).boundaryCellType_ == WallCell){
			if ((*iterCell).scalarVariableTn[0] < 1e-10)
			{
				// limit on avg
				(*iterCell).scalarVariableTn[0] = eps;
				(*iterCell).scalarVariableTm[0] = eps;
				// limit on slope
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					(*iterCell).scalarVariableTn[kk] = 0.0;
					(*iterCell).scalarVariableTm[kk] = 0.0;
					(*iterCell).scalarVariableTnCR[kk] = 0.0;
				}
			}
			//}
		}
		return true;
	}
}
#endif