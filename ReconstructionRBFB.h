#ifndef _RECONSTRUCTIONRBFB_H
#define _RECONSTRUCTIONRBFB_H
//#include "TypeDefine.h"
//#include "Point.h"

//#include "Geometry.h"
//#include "GaussIntegral.h"
#include "Variable.h"
#include "Math.h"
#include "Parameter.h"
#include "Reconstruction.h"

namespace ScalarCfv
{
	template <unsigned int O>
	class reconstructionRBFB1 : public reconstruction<O>
	{
	public:
		reconstructionRBFB1(){};
		~reconstructionRBFB1(){};
		const static unsigned int NDOFS = GLOBAL_NDOFS(O);
		const static unsigned int NDIFFS = GLOBAL_NDIFFS(O);
		const static unsigned int NDOFSCR = GLOBAL_NDOFSCR(O);
		const static unsigned int NDIFFSCR = GLOBAL_NDIFFSCR(O);

	public:
		inline bool initBaseMomentAndRelaxFactor(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			cellGaussDataVector *cellGaussData,
			GaussIntegralCellO1Grid<vO> *gaussIntegralCell) override;

		// bool initBaseMomentAndRelaxFactor(
		// 	parameter *parameter,
		// 	cellFieldDataVector *cellFieldData,
		// 	cellGaussDataVector *cellGaussData,
		// 	GaussIntegralCellO2Grid<vO> *gaussIntegralCell) override;

		inline bool initFaceWeight(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			faceFieldDataVector *faceFieldData) override;

		inline bool initReconstructionMatrixAndVector(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			cellGaussDataVector *cellGaussData,
			faceFieldDataVector *faceFieldData,
			faceGaussDataVector *faceGaussData,
			GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
			GaussIntegralFaceO1Grid<fO> *gaussIntegralFace) override;

		inline bool getBoundaryValueCheck(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			cellGaussDataVector *cellGaussData,
			faceFieldDataVector *faceFieldData,
			faceGaussDataVector *faceGaussData,
			GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
			GaussIntegralFaceO1Grid<fO> *gaussIntegralFace,
			std::ostream &out) override;

		// bool initReconstructionMatrixAndVector(
		// 	parameter *parameter,
		// 	cellFieldDataVector *cellFieldData,
		// 	cellGaussDataVector *cellGaussData,
		// 	faceFieldDataVector *faceFieldData,
		// 	faceGaussDataVector *faceGaussData,
		// 	GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
		// 	GaussIntegralFaceO2Grid<fO> *gaussIntegralFace) override;

		inline bool excuteReconstruction(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			cellGaussDataVector *cellGaussData,
			faceFieldDataVector *faceFieldData,
			faceGaussDataVector *faceGaussData,
			GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
			GaussIntegralFaceO1Grid<fO> *gaussIntegralFace) override;

		// bool excuteReconstruction(
		// 	parameter *parameter,
		// 	cellFieldDataVector *cellFieldData,
		// 	cellGaussDataVector *cellGaussData,
		// 	faceFieldDataVector *faceFieldData,
		// 	faceGaussDataVector *faceGaussData,
		// 	GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
		// 	GaussIntegralFaceO2Grid<fO> *gaussIntegralFace) override;

		// CR
		inline bool initReconstructionMatrixAndVectorCR(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			cellGaussDataVector *cellGaussData,
			faceFieldDataVector *faceFieldData,
			faceGaussDataVector *faceGaussData,
			GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
			GaussIntegralFaceO1Grid<fO> *gaussIntegralFace) override;

		// bool initReconstructionMatrixAndVectorCR(
		// 	parameter *parameter,
		// 	cellFieldDataVector *cellFieldData,
		// 	cellGaussDataVector *cellGaussData,
		// 	faceFieldDataVector *faceFieldData,
		// 	faceGaussDataVector *faceGaussData,
		// 	GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
		// 	GaussIntegralFaceO2Grid<fO> *gaussIntegralFace) override;

		inline bool excuteReconstructionCR(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			cellGaussDataVector *cellGaussData,
			faceFieldDataVector *faceFieldData,
			faceGaussDataVector *faceGaussData,
			GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
			GaussIntegralFaceO1Grid<fO> *gaussIntegralFace) override;

		// bool excuteReconstructionCR(
		// 	parameter *parameter,
		// 	cellFieldDataVector *cellFieldData,
		// 	cellGaussDataVector *cellGaussData,
		// 	faceFieldDataVector *faceFieldData,
		// 	faceGaussDataVector *faceGaussData,
		// 	GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
		// 	GaussIntegralFaceO2Grid<fO> *gaussIntegralFace) override;

		unsigned ReturnNDOFS() override { return NDOFS; };
		unsigned ReturnNDIFFS() override { return NDIFFS; };
		unsigned ReturnNDOFSCR() override { return NDOFSCR; };
		unsigned ReturnNDIFFSCR() override { return NDIFFSCR; };
	};

	//----------------------------------------------------------------------------------
	template <unsigned int O>
	inline bool reconstructionRBFB1<O>::initBaseMomentAndRelaxFactor(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		cellGaussDataVector *cellGaussData,
		GaussIntegralCellO1Grid<vO> *gaussIntegralCell)
	{
		std::cout << "initializing base moment and relax factor..." << std::endl;
		cellFieldDataVector::iterator iterCellFieldData;
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			// relax factor
			if ((*iterCellFieldData).boundaryCellType_ == InnerCell)
			{
				(*iterCellFieldData).relaxFactor = 1;
				(*iterCellFieldData).relaxFactorCR = 1;
			}
			else
			{
				(*iterCellFieldData).relaxFactor = 1;
				(*iterCellFieldData).relaxFactorCR = 1;
			}

			// base momment
			tensor1D<real, NDOFS> *phi = new tensor1D<real, NDOFS>[static_cast<int>((*iterCellFieldData).PG)];
			//------------------------------------------
			//|phi_1(GP1) phi_2(GP1) phi_3(GP1) ...
			//------------------------------------------
			//|phi_1(GP2) phi_2(GP2) phi_3(GP2) ...
			//------------------------------------------
			//|phi_1(GP3) phi_2(GP3) phi_3(GP3) ...
			//------------------------------------------

			iterCellFieldData->MeanIJ.setZero();
			for (int ii = 0; ii < static_cast<int>((*iterCellFieldData).PG); ++ii)
			{
				if (iterCellFieldData->cellType_ == Triangle)
					iterCellFieldData->MeanIJ +=
						CfvMath::getIJacobi((*iterCellFieldData).parametricValue[ii].first, *iterCellFieldData) *
						(*iterCellFieldData).parametricValue[ii].second;
				else if (iterCellFieldData->cellType_ == Quadrilateral)
					iterCellFieldData->MeanIJ +=
						CfvMath::getIJacobi((*iterCellFieldData).parametricValue[ii].first, *iterCellFieldData) *
						(*iterCellFieldData).parametricValue[ii].second;
				else
					assert(false);
			}

			// CfvMath::RegularizeJacobianSym(iterCellFieldData->MeanIJ);
			// CfvMath::RegularizeJacobian(InvMeanIJ);
#ifdef RBFB1_USE_UNITARY_MEANIJ
			Eigen::Matrix2d InvMeanIJ = iterCellFieldData->MeanIJ.inverse(); // = dxj/dxii
			real IJ0 = InvMeanIJ(0, Eigen::all).norm();
			real IJ1 = InvMeanIJ(1, Eigen::all).norm();
			InvMeanIJ(0, Eigen::all) /= IJ0;
			InvMeanIJ(1, Eigen::all) /= IJ1;
			// InvMeanIJ *= std::max(IJ0, IJ1);
			// InvMeanIJ *= std::min(IJ0, IJ1);
			// not multiplying any scale to

			// real IJ0 = std::sqrt(std::pow(iterCellFieldData->MeanIJ(0, 0), 2) + std::pow(iterCellFieldData->MeanIJ(1, 0), 2));
			// real IJ1 = std::sqrt(std::pow(iterCellFieldData->MeanIJ(0, 1), 2) + std::pow(iterCellFieldData->MeanIJ(1, 1), 2));
			// iterCellFieldData->MeanIJ(Eigen::all, 0) /= IJ0;
			// iterCellFieldData->MeanIJ(Eigen::all, 1) /= IJ1;
			iterCellFieldData->MeanIJ = InvMeanIJ.inverse();
#endif

			// std::cout << iterCellFieldData->MeanIJ << std::endl;
			// get the function values at Gaussian points
			for (int ii = 0; ii < static_cast<int>((*iterCellFieldData).PG); ++ii)
			{
				tensor1D<real, NDOFS> phiG;
				RBFB1GetMoment(
					(*iterCellFieldData).parametricValue[ii].first,
					(*iterCellFieldData).baryCenter,
					(*iterCellFieldData).lengthReference,
					phiG,
					*iterCellFieldData);

				for (int jj = 1; jj < NDOFS; ++jj)
				{							//ע�������ָ�귶Χ����20200315
					phi[ii][jj] = phiG[jj]; // phi: first->Gauss points; second->base function
				}
			}

			// get gauss integral
			for (int ii = 1; ii < NDOFS; ++ii)
			{
				real *f = new real[static_cast<int>((*iterCellFieldData).PG)];
				// tensor1D<real, static_cast<int>((*iterCellFieldData).PG)> f;
				for (int jj = 0; jj < static_cast<int>((*iterCellFieldData).PG); ++jj)
				{
					f[jj] = phi[jj][ii];
				}
				gaussIntegralCell->getIntegral(
					f,
					(*iterCellFieldData).index,
					cellGaussData,
					(*iterCellFieldData).baseMoment[ii]);
				(*iterCellFieldData).baseMoment[ii] /= (*iterCellFieldData).volume;
				delete[] f;
				f = NULL;
			}
			delete[] phi;
			phi = NULL;
			if ((std::fabs((*iterCellFieldData).baseMoment[1]) + std::fabs((*iterCellFieldData).baseMoment[2])) > parameter->EPS)
			{
				// std::cout << "	error: first order x- and y-moment should be small." << std::endl;
				// std::cout << (*iterCellFieldData).index << "\t"
				//	<< (*iterCellFieldData).baseMoment[1] << "\t"
				//	<< (*iterCellFieldData).baseMoment[2] << std::endl;
				// exit(1);
				// (*iterCellFieldData).baseMoment[1] = 0.0;
				// (*iterCellFieldData).baseMoment[2] = 0.0;
				// WARNING: DO NOT APPLY IN PARAMETRIC SCHEME
			}

			//
			//
			//
			// base momment CR
			tensor1D<real, NDOFSCR> *phiCR = new tensor1D<real, NDOFSCR>[static_cast<int>((*iterCellFieldData).PG)];
			for (int ii = 0; ii < static_cast<int>((*iterCellFieldData).PG); ++ii)
			{
				tensor1D<real, NDOFSCR> phiG;
				RBFB1GetMomentCR(
					(*iterCellFieldData).parametricValue[ii].first,
					(*iterCellFieldData).baryCenter,
					(*iterCellFieldData).lengthReference,
					phiG,
					*iterCellFieldData);
				for (int jj = 1; jj < NDOFSCR; ++jj)
				{							  //ע�������ָ�귶Χ����20200315
					phiCR[ii][jj] = phiG[jj]; // phi: first->Gauss points; second->base function
				}
			}

			// get gauss integral
			for (int ii = 1; ii < NDOFSCR; ++ii)
			{
				real *f = new real[static_cast<int>((*iterCellFieldData).PG)];
				// tensor1D<real, static_cast<int>((*iterCellFieldData).PG)> f;
				for (int jj = 0; jj < static_cast<int>((*iterCellFieldData).PG); ++jj)
				{
					f[jj] = phiCR[jj][ii];
				}
				gaussIntegralCell->getIntegral(
					f,
					(*iterCellFieldData).index,
					cellGaussData,
					(*iterCellFieldData).baseMomentCR[ii]);
				(*iterCellFieldData).baseMomentCR[ii] /= (*iterCellFieldData).volume;
				delete[] f;
				f = NULL;
			}
			delete[] phiCR;
			phiCR = NULL;
		}
		std::cout << " ..base moment and relax factor have been initialized." << std::endl;
		return true;
	}

	/*
		template <unsigned int O>
		bool reconstructionRBFB1<O>::initBaseMomentAndRelaxFactor(
			parameter *parameter,
			cellFieldDataVector *cellFieldData,
			cellGaussDataVector *cellGaussData,
			GaussIntegralCellO2Grid<vO> *gaussIntegralCell)
		{
			std::cout << "initializing base moment and relax factor..." << std::endl;

			// std::cout << "===================================" << std::endl;

			cellFieldDataVector::iterator iterCellFieldData;
			for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
			{
				// relax factor
				if ((*iterCellFieldData).boundaryCellType_ == InnerCell)
				{
					(*iterCellFieldData).relaxFactor = 1.0;
					(*iterCellFieldData).relaxFactorCR = 1.0;
				}
				else
				{
					(*iterCellFieldData).relaxFactor = 0.5;
					(*iterCellFieldData).relaxFactorCR = 0.5;
				}

				// base momment
				tensor1D<real, NDOFS> *phi = new tensor1D<real, NDOFS>[static_cast<int>((*iterCellFieldData).PG)];
				//------------------------------------------
				//|phi_1(GP1) phi_2(GP1) phi_3(GP1) ...
				//------------------------------------------
				//|phi_1(GP2) phi_2(GP2) phi_3(GP2) ...
				//------------------------------------------
				//|phi_1(GP3) phi_2(GP3) phi_3(GP3) ...
				//------------------------------------------

				// get the function values at Gaussian points
				for (int ii = 0; ii < static_cast<int>((*iterCellFieldData).PG); ++ii)
				{
					tensor1D<real, NDOFS> phiG;

					//	std::cout << "------------------------" << std::endl;
					//	std::cout << "PG: " << ii << std::endl;
					//	std::cout << (*iterCellFieldData).gaussPairVector_[1].p.x << "\t" << (*iterCellFieldData).gaussPairVector_[1].p.y
					//		<< "\t" << (*iterCellFieldData).gaussPairVector_[9].p.x << "\t" << (*iterCellFieldData).gaussPairVector_[9].p.y << std::endl;
					//	std::cout << (*iterCellFieldData).baryCenter.x << "\t" << (*iterCellFieldData).baryCenter.y << std::endl;
					//	std::cout << (*iterCellFieldData).lengthReference.x << "\t" << (*iterCellFieldData).lengthReference.y << std::endl;
					//	std::cout << phi[ii][4] << "\t" << phi[ii][5] << "\t" << phi[ii][6] << std::endl;
					//	std::cout << phi[ii][7] << "\t" << phi[ii][8] << "\t" << phi[ii][9] << std::endl;

					CfvMath::getMomentRBFB1(
						(*iterCellFieldData).parametricValue[ii].first,
						(*iterCellFieldData).baryCenter,
						(*iterCellFieldData).lengthReference,
						phiG,
						(*iterCellFieldData));
					for (int jj = 1; jj < NDOFS; ++jj)
					{							//ע�������ָ�귶Χ����20200315
						phi[ii][jj] = phiG[jj]; // phi: first->Gauss points; second->base function
					}

					//	std::cout << "------------------------" << std::endl;
					//	std::cout << "PG: " << ii << std::endl;
					//	std::cout << phi[ii][1] << "\t" << phi[ii][1] << "\t" << phi[ii][3] << std::endl;
					//	std::cout << phi[ii][4] << "\t" << phi[ii][5] << "\t" << phi[ii][6] << std::endl;
					//	std::cout << phi[ii][7] << "\t" << phi[ii][8] << "\t" << phi[ii][9] << std::endl;
				}

				// system("pause");

				// get gauss integral
				for (int ii = 1; ii < NDOFS; ++ii)
				{
					real *f = new real[static_cast<int>((*iterCellFieldData).PG)];
					// tensor1D<real, static_cast<int>((*iterCellFieldData).PG)> f;
					for (int jj = 0; jj < static_cast<int>((*iterCellFieldData).PG); ++jj)
					{
						f[jj] = phi[jj][ii];
					}
					gaussIntegralCell->getIntegral(
						f,
						(*iterCellFieldData).index,
						cellGaussData,
						(*iterCellFieldData).baseMoment[ii]);
					(*iterCellFieldData).baseMoment[ii] /= (*iterCellFieldData).volume;
					delete[] f;
					f = NULL;
				}
				delete[] phi;
				phi = NULL;
				if ((std::fabs((*iterCellFieldData).baseMoment[1]) + std::fabs((*iterCellFieldData).baseMoment[2])) > parameter->EPS)
				{
					// std::cout << "	error: first order x- and y-moment should be small." << std::endl;
					// std::cout << (*iterCellFieldData).index << "\t"
					//	<< (*iterCellFieldData).baseMoment[1] << "\t"
					//	<< (*iterCellFieldData).baseMoment[2] << std::endl;
					// exit(1);
					// (*iterCellFieldData).baseMoment[1] = 0.0;
					// (*iterCellFieldData).baseMoment[2] = 0.0;
					// WARNING: DO NOT APPLY IN PARAMETRIC SCHEME
				}
			}
			std::cout << " ..base moment and relax factor have been initialized." << std::endl;
			return true;
		}
	*/
	//----------------------------------------------------------------------------------

	template <unsigned int O>
	inline bool reconstructionRBFB1<O>::initFaceWeight(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		faceFieldDataVector *faceFieldData)
	{
		std::cout << "initializing face weight..." << std::endl;
		cellFieldDataVector::iterator iterCellFieldData;
		faceFieldDataVector::iterator iterFaceFieldData;
		// reconstruction coefficient
		for (iterFaceFieldData = faceFieldData->begin(); iterFaceFieldData != faceFieldData->end(); ++iterFaceFieldData)
		{
			int kf = iterFaceFieldData - faceFieldData->begin() + 1;
			int cl = (*iterFaceFieldData).faceCellIndex[1];
			int cr = (*iterFaceFieldData).faceCellIndex[2];
			int ff = 1;
			point faceL = iterFaceFieldData->faceNode[1].second - iterFaceFieldData->faceNode[2].second;
			for (; ff < (*cellFieldData)[cl - 1].cellFaceNumber + 1; ff++)
				if ((*cellFieldData)[cl - 1].cellFaceIndex[ff] == kf)
					break;
			assert(ff < (*cellFieldData)[cl - 1].cellFaceNumber + 1);
			iterFaceFieldData->interFacialJacobi.setZero();
			for (int gg = 0; gg < int(iterFaceFieldData->fPG); gg++)
			{
				point pparaml = CfvMath::GetFaceParam((*cellFieldData)[cl - 1].cellType_, ff, iterFaceFieldData->parametricValue[gg].first);
				iterFaceFieldData->interFacialJacobi += CfvMath::getFaceJacobi(pparaml, (*cellFieldData)[cl - 1], ff) * iterFaceFieldData->parametricValue[gg].second;
			}
			// iterFaceFieldData->interFacialJacobi = (*cellFieldData)[cl - 1].MeanIJ.inverse(); // dxj/dxii

			point unitNormalVector = (1.0 / (*iterFaceFieldData).area) * (*iterFaceFieldData).normalVector;
			point midPoint = 0.5 * ((*iterFaceFieldData).faceNode[1].second + (*iterFaceFieldData).faceNode[2].second);
			// point sideOff = (*iterFaceFieldData).sideOff;
			point delta;

			real w[NDIFFS] = {0.0}; // Warning Current Coincidentally FaceFieldData.NDOFS == self.NDIFFS

			real refLR = 0.0;
			real Omega = 0.0;
			real d = 0.0;
			if (cr > 0)
			{
				int ffr = 1;
				for (; ffr < (*cellFieldData)[cr - 1].cellFaceNumber + 1; ffr++)
					if ((*cellFieldData)[cr - 1].cellFaceIndex[ffr] == kf)
						break;
				assert(ffr < (*cellFieldData)[cr - 1].cellFaceNumber + 1);
				for (int gg = 0; gg < int(iterFaceFieldData->fPG); gg++)
				{
					point pparamr = CfvMath::GetFaceParam((*cellFieldData)[cr - 1].cellType_, ffr, iterFaceFieldData->parametricValue[gg].first);
					iterFaceFieldData->interFacialJacobi -= CfvMath::getFaceJacobi(pparamr, (*cellFieldData)[cr - 1], ffr) * iterFaceFieldData->parametricValue[gg].second;
				}
				// iterFaceFieldData->interFacialJacobi += (*cellFieldData)[cr - 1].MeanIJ.inverse();
				iterFaceFieldData->interFacialJacobi *= 0.5;

				(*iterFaceFieldData).weightVF = 1.0; // deflult value
				iterCellFieldData = cellFieldData->begin() + cl - 1;
				real omegaL = (*iterCellFieldData).volume;
				point scaleL = (*iterCellFieldData).lengthReference;
				point baryCenterL = (*iterCellFieldData).baryCenter;
				iterCellFieldData = cellFieldData->begin() + cr - 1;
				real omegaR = (*iterCellFieldData).volume;
				point scaleR = (*iterCellFieldData).lengthReference;
				point baryCenterR = (*iterCellFieldData).baryCenter;

				// correct
				delta.x = std::fabs((baryCenterL - baryCenterR + (*iterFaceFieldData).sideOff).x);
				delta.y = std::fabs((baryCenterL - baryCenterR + (*iterFaceFieldData).sideOff).y);

				// for parametric
				delta = -1.0 * (baryCenterL - baryCenterR + (*iterFaceFieldData).sideOff);

				refLR = 0.5 * (omegaL + omegaR) / (*iterFaceFieldData).area;
				Omega = 0.5 * (omegaL + omegaR);
				for (int ii = 0; ii < static_cast<int>((*iterFaceFieldData).NDOFS); ++ii)
				{
					w[ii] = 1.0;
				}
			}
			else if (cr == Wall)
			{
				if (parameter->isViscous == true)
				{
					// iterCellFieldData = cellFieldData->begin() + cl - 1;
					// real l1 = std::fabs((*iterCellFieldData).lengthReference.x / ((*iterCellFieldData).lengthReference.y + parameter->EPS));
					// real l2 = std::fabs((*iterCellFieldData).lengthReference.y / ((*iterCellFieldData).lengthReference.x + parameter->EPS));
					//(*iterFaceFieldData).weightVF = getMax(l1, l2)*getMax(l1, l2);//deflult value
					(*iterFaceFieldData).weightVF = 1.0; // deflult value
				}
				else if (parameter->isViscous == false)
				{
					(*iterFaceFieldData).weightVF = 0.0; // deflult value
				}
				iterCellFieldData = cellFieldData->begin() + cl - 1;
				real omegaL = (*iterCellFieldData).volume;
				point baryCenterL = (*iterCellFieldData).baryCenter;
				point scaleL = (*iterCellFieldData).lengthReference;
				Omega = omegaL;

				Omega = omegaL;
				delta.x = std::fabs((baryCenterL - midPoint).x);
				delta.y = std::fabs((baryCenterL - midPoint).y);

				// for parametric
				delta = midPoint - baryCenterL;
				delta = delta * 2.0;

				w[0] = 1.0;
				refLR = omegaL / (*iterFaceFieldData).area;
#ifdef RBFB1_USE_DELTA_INTERFACEJ
				// w[0] *= std::pow(faceL.length() / delta.length(), -1.0);
#endif
			}
			else if (cr == FarField)
			{
				(*iterFaceFieldData).weightVF = 1.0; // deflult value
				iterCellFieldData = cellFieldData->begin() + cl - 1;
				real omegaL = (*iterCellFieldData).volume;
				point baryCenterL = (*iterCellFieldData).baryCenter;
				point scaleL = (*iterCellFieldData).lengthReference;
				Omega = omegaL;

				Omega = omegaL;
				delta.x = std::fabs((baryCenterL - midPoint).x);
				delta.y = std::fabs((baryCenterL - midPoint).y);

				// for parametric
				delta = midPoint - baryCenterL;

				w[0] = 1.0;
				refLR = omegaL / (*iterFaceFieldData).area;
			}
			else if (cr == Symmetric)
			{
				(*iterFaceFieldData).weightVF = 1.0;
				iterCellFieldData = cellFieldData->begin() + cl - 1;
				real omegaL = (*iterCellFieldData).volume;
				point baryCenterL = (*iterCellFieldData).baryCenter;
				real normalProjection = std::fabs(getInnerProduct(unitNormalVector, (baryCenterL - midPoint)));

				point scaleL = (*iterCellFieldData).lengthReference;
				Omega = omegaL;
				// for (int ii = 0; ii < static_cast<int>((*iterFaceFieldData).NDOFS); ++ii){
				//	w[ii] = 1.0;
				// }

				//	Omega = omegaL;
				Omega = omegaL;
				delta.x = 2.0 * std::fabs((normalProjection * unitNormalVector).x);
				delta.y = 2.0 * std::fabs((normalProjection * unitNormalVector).y);

				// for parametric
				delta = 2 * (midPoint - baryCenterL);

				w[0] = 1.0;
				w[1] = 1.0;
				w[2] = 1.0;
				refLR = omegaL / (*iterFaceFieldData).area;
			}
			else if (cr == OutFlow)
			{
				std::cout << "	warning: OutFlow is not proposed currently." << std::endl;
			}
			else if (cr == InFlow)
			{
				std::cout << "	warning: InFlow is not proposed currently." << std::endl;
			}

			// real dx = delta.x;
			// real dy = delta.y;
			real dx = unitNormalVector.x * delta.length();
			real dy = unitNormalVector.y * delta.length();
			// real dx = unitNormalVector.x;
			// real dy = unitNormalVector.y;
#ifdef RBFB1_USE_DELTA_Jacobi_INTERFACEJ
			real J1norm = +unitNormalVector.x * iterFaceFieldData->interFacialJacobi(0, 0) + unitNormalVector.y * iterFaceFieldData->interFacialJacobi(0, 1);
			real J2tang = -unitNormalVector.y * iterFaceFieldData->interFacialJacobi(1, 0) + unitNormalVector.x * iterFaceFieldData->interFacialJacobi(1, 1);
			J2tang *= 0.01;
			iterFaceFieldData->interFacialJacobi(0, 0) = unitNormalVector.x * J1norm;
			iterFaceFieldData->interFacialJacobi(0, 1) = unitNormalVector.y * J1norm;
			iterFaceFieldData->interFacialJacobi(1, 0) = -unitNormalVector.y * J2tang;
			iterFaceFieldData->interFacialJacobi(1, 1) = unitNormalVector.x * J2tang;
			// std::cout << "FF" << unitNormalVector.length() << ' ' << J1norm / std::sqrt(std::pow(iterFaceFieldData->interFacialJacobi(0, 0), 2)+std::pow(iterFaceFieldData->interFacialJacobi(0, 1), 2));
			//  iterFaceFieldData->interFacialJacobi *= delta.length();
#endif
#ifdef RBFB1_USE_DELTA_NT_INTERFACEJ
			real wtgt = 1;
			iterFaceFieldData->interFacialJacobi(0, 0) = unitNormalVector.x * delta.length();
			iterFaceFieldData->interFacialJacobi(0, 1) = unitNormalVector.y * delta.length();
			iterFaceFieldData->interFacialJacobi(1, 0) = -unitNormalVector.y * faceL.length() * wtgt;
			iterFaceFieldData->interFacialJacobi(1, 1) = unitNormalVector.x * faceL.length() * wtgt;
#endif
#ifdef RBFB1_USE_DELTA_INTERFACEJ
			// real wtgt = std::pow(std::min(std::max(delta.length() / faceL.length(), 1e-0), 1e3), 1);
			// real wtgt = std::sqrt(delta.length() / faceL.length());
			real wtgt = 1.;
			iterFaceFieldData->interFacialJacobi(0, 0) = unitNormalVector.x * delta.length();
			iterFaceFieldData->interFacialJacobi(0, 1) = unitNormalVector.y * delta.length();
			iterFaceFieldData->interFacialJacobi(1, 0) = -unitNormalVector.y * delta.length() * wtgt;
			iterFaceFieldData->interFacialJacobi(1, 1) = unitNormalVector.x * delta.length() * wtgt;
#endif

			real alpha = 1.0;
			w[0] *= alpha;
			w[1] *= alpha;
			w[2] *= alpha;

			for (int ii = 0; ii < static_cast<int>((*iterFaceFieldData).NDOFS); ++ii)
			{ // 0: for b_vec
				int m = CfvMath::mMapping[ii];
				int n = CfvMath::nMapping[ii];
				real Cm = CfvMath::getCombination(m, n);
				real Pmn = CfvMath::getFactorial(m, n);

				(*iterFaceFieldData).faceWeightVF[ii] =
					w[ii] * (*iterFaceFieldData).weightVF * std::pow(dx, 2 * m) * std::pow(dy, 2 * n) * std::pow(Cm, 2) / std::pow(Pmn, 2);
				// for parametric
				(*iterFaceFieldData).faceWeightVF[ii] =
					w[ii] * (*iterFaceFieldData).weightVF * std::pow(dx, 1 * m) * std::pow(dy, 1 * n) * std::pow(Cm, 1) / std::pow(Pmn, 1);
			}
		}
		std::cout << " ..face weight has been initialized." << std::endl;
		return true;
	}

	//	here was commented code form Reconstruction.h

	//----------------------------------------------------------------------------------

	// Upgrade V1
	template <unsigned int O>
	inline bool reconstructionRBFB1<O>::initReconstructionMatrixAndVector(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		cellGaussDataVector *cellGaussData,
		faceFieldDataVector *faceFieldData,
		faceGaussDataVector *faceGaussData,
		GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO1Grid<fO> *gaussIntegralFace)
	{
		std::cout << "initializing reconstruction matrix and vector..." << std::endl;

		cellFieldDataVector::iterator iterCellFieldData;
		cellFieldDataVector::iterator iterCellFieldData_;
		faceFieldDataVector::iterator iterFaceFieldData_;
		// matrix A B b
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			// tensor2D<real, NDOFS, NDOFS> Aii;
			const int Aiisize = NDOFS;
			Eigen::MatrixXd Aii(NDOFS, NDOFS);
			Aii.setZero();

#ifdef TRIAL
			std::cout << "\n Some Cell: \n";
			for (auto &i : iterCellFieldData->parametricValue)
				std::cout << i.first << '\n';
#endif
			for (int ff = 1; ff < (*iterCellFieldData).cellFaceNumber + 1; ++ff)
			{ //ע��ѭ�����ޣ���20200314
#ifdef TRIAL
				std::cout << "SomeFace: \n";
#endif
				// for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
				// {
				// }
				int kf = (*iterCellFieldData).cellFaceIndex[ff];
				iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
				int cl = (*iterFaceFieldData_).faceCellIndex[1];
				int cr = (*iterFaceFieldData_).faceCellIndex[2];
				int cff = 0, cffr = 1;
				// alternate cell��ֻ���ڲ���Ԫ���ܻ���Ҫ
				if (cr == (*iterCellFieldData).index)
				{
					int ck = cr;
					cr = cl; //�߽粻�õ���
					cl = ck;
					cff = 1, cffr = 0;
					// std::cout << "		alternate left and right cell." << std::endl;
				}

				//�����
				iterCellFieldData_ = cellFieldData->begin() + cl - 1;
				// tensor2D<real, NDOFS, NDOFS> Bij;
				const int Siz = NDOFS;
				Eigen::MatrixXd Bij(NDOFS, NDOFS);
				Bij.setZero();
				// tensor1D<real, NDOFS> bi;
				Eigen::VectorXd bi(NDOFS);
				bi.setZero();
				if (cr > 0)
				{
					int ffr = 1;
					for (; ffr < (*cellFieldData)[cr - 1].cellFaceNumber + 1; ffr++)
						if ((*cellFieldData)[cr - 1].cellFaceIndex[ffr] == kf)
							break;
					assert(ffr < (*cellFieldData)[cr - 1].cellFaceNumber + 1);

					// find the neibhour cell index
					int cc;
					for (cc = 1; cc < (*iterCellFieldData).cellCellNumber + 1; ++cc)
					{
						if (cr == (*iterCellFieldData_).cellCellIndex[cc])
						{
							point match = (*iterCellFieldData).cellFaceSideOff[cc] - (*iterFaceFieldData_).sideOff;
							if (std::fabs(match.length()) < parameter->EPS)
							{
								break;
							}
						}
					}
					if (cc >= (*iterCellFieldData).cellCellNumber + 1)
					{
						std::cout << "	error: fail to find the neighbor cell in reconstruction." << std::endl;
						exit(1);
					}

					// integrated functions
					tensor2D<real, NDOFS, NDOFS> *fII = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					tensor2D<real, NDOFS, NDOFS> *fIJ = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					tensor1D<real, NDOFS> *fI = new tensor1D<real, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix and vector on each face Gauss point.

					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{

						// diff base matrix i
						tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
						tensor2D<real, NDOFSCR, NDIFFSCR> matrixDiffBaseICR;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						tensor1D<real, NDOFSCR> momentICR;
						// point p = (*iterFaceFieldData_).gaussPairVector_[gg].p + (*iterFaceFieldData_).sideOff;//20200318��������ƫ����
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						// right cell ����ƫ����
						//  ----------------       --------
						//  |     ||      |		  |     ||
						//  |     ||	     |   ...  |     ||
						//  |     ||	     |		  |     ||
						//-----------------		 ---------
						//  |     ||      |		  |     ||
						//  | R   ||  L   |   ...  | R   ||
						//  |     ||	     |		  |     ||
						//-----------------		 ---------
						//  |     ||      |		  |     ||
						//  |     ||	     |   ...  |     ||
						//  |     ||	     |		  |     ||
						//-----------------		 ---------
						// cell R���ǰ���ˮƽ�������ڴ����ģ�����ı���Ӧ�ü���ƫ����
						//	point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;//20200318����
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
#ifdef TRIAL
						std::cout << (*iterFaceFieldData_).parametricValue[gg].first << std::endl;
#endif
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
							momentICR[kk] = (*iterCellFieldData_).baseMomentCR[kk];
						point pparaml = CfvMath::GetFaceParam((*cellFieldData)[cl - 1].cellType_, ff, iterFaceFieldData_->parametricValue[gg].first);
						point pparamr = CfvMath::GetFaceParam((*cellFieldData)[cr - 1].cellType_, ffr, iterFaceFieldData_->parametricValue[gg].first, true);
#ifdef TRIAL
						std::cout << "LR" << CfvMath::getPoint(pparaml, (*cellFieldData)[cl - 1]) << CfvMath::getPoint(pparamr, (*cellFieldData)[cr - 1]) << std::endl;
#endif
						assert(CfvMath::getPoint(pparaml, (*cellFieldData)[cl - 1]).x == CfvMath::getPoint(pparamr, (*cellFieldData)[cr - 1]).x);
						RBFB1GetDiffBaseValue(
							pparaml,
							baryCenterI,
							scaleI,
							momentI,
							matrixDiffBaseI,
							*iterCellFieldData_);
						if (iterFaceFieldData_->diffBaseValueData[cff][gg][0][0] == UNINITReal)
							CfvMath::VVMatCopy(matrixDiffBaseI, iterFaceFieldData_->diffBaseValueData[cff][gg],
											   0, NDOFS, 0, NDIFFS);
						RBFB1GetDiffBaseValueCR(
							pparaml,
							baryCenterI,
							scaleI,
							momentICR,
							matrixDiffBaseICR,
							*iterCellFieldData_);
						if (iterFaceFieldData_->diffBaseValueDataCR[cff][gg][0][0] == UNINITReal)
							CfvMath::VVMatCopy(matrixDiffBaseICR, iterFaceFieldData_->diffBaseValueDataCR[cff][gg],
											   0, NDOFSCR, 0, NDIFFSCR);
						// diff base matrix j
						tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseJ;
						tensor2D<real, NDOFSCR, NDIFFSCR> matrixDiffBaseJCR;
						// base moment j
						tensor1D<real, NDOFS> momentJ;
						tensor1D<real, NDOFSCR> momentJCR;
						iterCellFieldData_ = cellFieldData->begin() + cr - 1;
						point baryCenterJ = (*iterCellFieldData_).baryCenter;
						point scaleJ = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							momentJ[kk] = (*iterCellFieldData_).baseMoment[kk];
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
							momentJCR[kk] = (*iterCellFieldData_).baseMomentCR[kk];
						RBFB1GetDiffBaseValue(
							pparamr, // 20200318����
							baryCenterJ,
							scaleJ,
							momentJ,
							matrixDiffBaseJ,
							*iterCellFieldData_);
						if (iterFaceFieldData_->diffBaseValueData[cffr][gg][0][0] == UNINITReal)
							CfvMath::VVMatCopy(matrixDiffBaseJ, iterFaceFieldData_->diffBaseValueData[cffr][gg],
											   0, NDOFS, 0, NDIFFS);
						RBFB1GetDiffBaseValueCR(
							pparamr, // 20200318����
							baryCenterJ,
							scaleJ,
							momentJCR,
							matrixDiffBaseJCR,
							*iterCellFieldData_);
						if (iterFaceFieldData_->diffBaseValueDataCR[cffr][gg][0][0] == UNINITReal)
							CfvMath::VVMatCopy(matrixDiffBaseJCR, iterFaceFieldData_->diffBaseValueDataCR[cffr][gg],
											   0, NDOFSCR, 0, NDIFFSCR);
						// ll: row    rr: col
						// get the integrated function values on each face Gauss point
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// matrix
						CfvMath::VVMatConjProd2DDiffCombine(matrixDiffBaseI, matrixDiffBaseI, fII[gg], iterFaceFieldData_->faceWeightVF,
															1, NDOFS, 0, NDIFFS, 1, NDOFS,
															*iterFaceFieldData_);
						CfvMath::VVMatConjProd2DDiffCombine(matrixDiffBaseI, matrixDiffBaseJ, fIJ[gg], iterFaceFieldData_->faceWeightVF,
															1, NDOFS, 0, NDIFFS, 1, NDOFS,
															*iterFaceFieldData_);
						// vector
						for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
						{
							// fI[gg][ll] = (*iterFaceFieldData_).faceWeightVF[0] * matrixDiffBaseI[ll][0];
							//  parametric
							fI[gg][ll] = std::pow((*iterFaceFieldData_).faceWeightVF[0], 2) * matrixDiffBaseI[ll][0];
						}
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii(ll, rr) += result; //+ means sum of each face

							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fIJ[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Bij(ll, rr) = result;

							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
					// get the face Gauss point for coefficient vector. integral along the interface between cell i and cell j
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real parametricArea;
						real result;
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							f[gg] = fI[gg][ll];
							weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
							cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
							parametricArea = (*iterFaceFieldData_).parametricArea;
						}
						gaussIntegralFace->getIntegral(
							static_cast<int>((*iterFaceFieldData_).fPG),
							f,
							weight,
							cofJacobi,
							parametricArea,
							result);
						bi[ll] = result;
						delete[] f;
						delete[] weight;
						delete[] cofJacobi;
						f = NULL;
						weight = NULL;
						cofJacobi = NULL;
					}
					// assign the matrix and vector with respect to each face
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					// CfvMath::VEVecMatCopy(bi, iterCellFieldData_->vectorbij, cc, 1, NDOFS);
					iterCellFieldData_->vectorbij.row(cc) = bi;
					// CfvMath::VVMatCopy(Bij, iterCellFieldData_->matrixBij[cc], 1, NDOFS, 1, NDOFS);
					iterCellFieldData_->matrixBij[cc] = Bij;
					delete[] fII;
					delete[] fIJ;
					delete[] fI;
					fII = NULL;
					fIJ = NULL;
					fI = NULL;
				}
				else if (cr == FarField || cr == Symmetric || cr == Wall)
				{
					tensor2D<real, NDOFS, NDOFS> *fII = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
						tensor2D<real, NDOFSCR, NDIFFSCR> matrixDiffBaseICR;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						tensor1D<real, NDOFSCR> momentICR;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						point pparaml;
						pparaml = CfvMath::GetFaceParam((*cellFieldData)[cl - 1].cellType_, ff, iterFaceFieldData_->parametricValue[gg].first);
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
							momentICR[kk] = (*iterCellFieldData_).baseMomentCR[kk];
						RBFB1GetDiffBaseValue(
							pparaml,
							baryCenterI,
							scaleI,
							momentI,
							matrixDiffBaseI,
							*iterCellFieldData_);
						if (iterFaceFieldData_->diffBaseValueData[cff][gg][0][0] == UNINITReal)
							CfvMath::VVMatCopy(matrixDiffBaseI, iterFaceFieldData_->diffBaseValueData[cff][gg],
											   0, NDOFS, 0, NDIFFS);
						RBFB1GetDiffBaseValueCR(
							pparaml,
							baryCenterI,
							scaleI,
							momentICR,
							matrixDiffBaseICR,
							*iterCellFieldData_);
						if (iterFaceFieldData_->diffBaseValueDataCR[cff][gg][0][0] == UNINITReal)
							CfvMath::VVMatCopy(matrixDiffBaseICR, iterFaceFieldData_->diffBaseValueDataCR[cff][gg],
											   0, NDOFSCR, 0, NDIFFSCR);
						// ll: row    rr: col
						// get the integrated function values on each face Gauss point
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// matrix
						CfvMath::VVMatConjProd2DDiffCombine(matrixDiffBaseI, matrixDiffBaseI, fII[gg], iterFaceFieldData_->faceWeightVF,
															1, NDOFS, 0, NDIFFS, 1, NDOFS,
															*iterFaceFieldData_);
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii(ll, rr) += result; //+ means sum of each face
							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
					delete[] fII;
					fII = NULL;
				}
				else if (cr == OutFlow)
				{
					std::cout << "	warning: OutFlow is not proposed currently." << std::endl;
				}
				else if (cr == InFlow)
				{
					std::cout << "	warning: InFlow is not proposed currently." << std::endl;
				}
			}
			// iterCellFieldData_ = cellFieldData->begin() + cl - 1;
			// assign Aii
			// CfvMath::VVMatCopy(Aii, iterCellFieldData->matrixAii, 1, NDOFS, 1, NDOFS);
			iterCellFieldData->matrixAii = Aii;
			// tensor2D<real, NDOFS, NDOFS> inverseAii;
			// CfvMath::getMatrixGeneralInverse(Aii, inverseAii);
			// assign inverse Aii
			// CfvMath::VVMatCopy(inverseAii, iterCellFieldData->matrixAiiInverse, 1, NDOFS, 1, NDOFS);
			// auto SVDResult = iterCellFieldData->matrixAii.bottomRightCorner(iterCellFieldData->matrixAii.rows() - 1, iterCellFieldData->matrixAii.cols() - 1).bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);

			// auto SVDResult = Aii.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
			// iterCellFieldData->matrixAiiInverse=
			// 	SVDResult.solve(Eigen::MatrixXd::Identity(iterCellFieldData->matrixAii.cols(), iterCellFieldData->matrixAii.rows()));
			CfvMath::EigenLeastSquareInverse(Aii, iterCellFieldData->matrixAiiInverse);

			// assign inverseAii*Bij & inverseAii*bi, for convenient calculation
			for (int cc = 1; cc < (*iterCellFieldData).cellCellNumber + 1; ++cc)
			{ // cellFaceNumber -> cellCellNumber
				// assign inverseAii*Bij
				// assign inverseAii*bi
				// CfvMath::VVMatMat(inverseAii, iterCellFieldData->matrixBij[cc], iterCellFieldData->matrixAiiInverseBij[cc],
				// 				  1, NDOFS, 1, NDOFS, 1, NDOFS);
				iterCellFieldData->matrixAiiInverseBij[cc] = iterCellFieldData->matrixAiiInverse * iterCellFieldData->matrixBij[cc];
				// CfvMath::VEMatVec(inverseAii, (iterCellFieldData->vectorbij.row(cc)), iterCellFieldData->vectorAiiInversebij[cc],
				// 				  1, NDOFS, 1, NDOFS);
				iterCellFieldData->vectorAiiInversebij.row(cc) =
					iterCellFieldData->vectorbij.row(cc) *
					iterCellFieldData->matrixAiiInverse.transpose();
			}
		}
		std::cout << " ..reconstruction matrix and vector has been initialized." << std::endl;
		return true;
	}

	// Upgrade V1
	template <unsigned int O>
	inline bool reconstructionRBFB1<O>::getBoundaryValueCheck(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		cellGaussDataVector *cellGaussData,
		faceFieldDataVector *faceFieldData,
		faceGaussDataVector *faceGaussData,
		GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO1Grid<fO> *gaussIntegralFace,
		std::ostream &out)
	{
		std::cout << "evaluating BC rec value..." << std::endl;

		cellFieldDataVector::iterator iterCellFieldData;
		cellFieldDataVector::iterator iterCellFieldData_;
		faceFieldDataVector::iterator iterFaceFieldData_;
		// matrix A B b
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			for (int ff = 1; ff < (*iterCellFieldData).cellFaceNumber + 1; ++ff)
			{ //?????????????20200314
				int kf = (*iterCellFieldData).cellFaceIndex[ff];
				iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
				int cl = (*iterFaceFieldData_).faceCellIndex[1];
				int cr = (*iterFaceFieldData_).faceCellIndex[2];
				// alternate cell???????????????????
				if (cr == (*iterCellFieldData).index)
				{
					int ck = cr;
					cr = cl; //??�^?????
					cl = ck;
					// std::cout << "		alternate left and right cell." << std::endl;
				}

				iterCellFieldData_ = cellFieldData->begin() + cl - 1;
				if (cr > 0)
				{
					// do nothing
				}
				else if (cr == Wall)
				{
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// base moment i
						tensor1D<real, NDOFS> momentI;
						tensor1D<real, NDOFSCR> momentICR;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						point pparaml;
						pparaml = CfvMath::GetFaceParam((*cellFieldData)[cl - 1].cellType_, ff, iterFaceFieldData_->parametricValue[gg].first);
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
							momentICR[kk] = (*iterCellFieldData_).baseMomentCR[kk];

						tensor1D<real, NDOFS> baseValueI;
						RBFB1GetBaseValue(
							pparaml,
							baryCenterI,
							scaleI,
							momentI,
							baseValueI,
							*iterCellFieldData_);
						tensor1D<real, NDOFSCR> baseValueICR;
						RBFB1GetBaseValueCR(
							pparaml,
							baryCenterI,
							scaleI,
							momentICR,
							baseValueICR,
							*iterCellFieldData_);

						real pointRecValue = (*iterCellFieldData).scalarVariableTn[0];
						real pointRecValueCR = (*iterCellFieldData).scalarVariableTnCR[0];
						for (unsigned kk = 1; kk < NDOFS; ++kk)
							pointRecValue += (*iterCellFieldData).scalarVariableTn[kk] * baseValueI[kk];
						for (unsigned kk = 1; kk < NDOFSCR; ++kk)
							pointRecValueCR += (*iterCellFieldData).scalarVariableTnCR[kk] * baseValueICR[kk];
						// TODO
						out << p.x << '\t' << p.y << '\t' << pointRecValue << '\t' << pointRecValueCR << '\n';
					}
				}
				else if (cr == FarField)
				{
				}
				else if (cr == Symmetric)
				{
				}
				else if (cr == OutFlow)
				{
					std::cout << "	warning: OutFlow is not proposed currently." << std::endl;
				}
				else if (cr == InFlow)
				{
					std::cout << "	warning: InFlow is not proposed currently." << std::endl;
				}
			}
			// iterCellFieldData_ = cellFieldData->begin() + cl - 1;
		}
		out << std::endl;
		std::cout << " ..evaluating BC rec value done." << std::endl;
		return true;
	}

	// Upgrade V1
	/*
	template <unsigned int O>
	bool reconstructionRBFB1<O>::initReconstructionMatrixAndVector__NotReady(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		cellGaussDataVector *cellGaussData,
		faceFieldDataVector *faceFieldData,
		faceGaussDataVector *faceGaussData,
		GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO2Grid<fO> *gaussIntegralFace)
	{
		std::cout << "initializing reconstruction matrix and vector..." << std::endl;

		cellFieldDataVector::iterator iterCellFieldData;
		cellFieldDataVector::iterator iterCellFieldData_;
		faceFieldDataVector::iterator iterFaceFieldData_;
		// matrix A B b
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			tensor2D<real, NDOFS, NDOFS> Aii;

			for (int ff = 1; ff < (*iterCellFieldData).cellFaceNumber + 1; ++ff)
			{ //ע��ѭ�����ޣ���20200314
				int kf = (*iterCellFieldData).cellFaceIndex[ff];
				iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
				int cl = (*iterFaceFieldData_).faceCellIndex[1];
				int cr = (*iterFaceFieldData_).faceCellIndex[2];
				// alternate cell��ֻ���ڲ���Ԫ���ܻ���Ҫ
				if (cr == (*iterCellFieldData).index)
				{
					int ck = cr;
					cr = cl; //�߽粻�õ���
					cl = ck;
					// std::cout << "		alternate left and right cell." << std::endl;
				}

				//�����
				iterCellFieldData_ = cellFieldData->begin() + cl - 1;
				tensor2D<real, NDOFS, NDOFS> Bij;
				tensor1D<real, NDOFS> bi;
				if (cr > 0)
				{
					// find the neibhour cell index
					int cc;
					for (cc = 1; cc < (*iterCellFieldData).cellCellNumber + 1; ++cc)
					{
						if (cr == (*iterCellFieldData_).cellCellIndex[cc])
						{
							point match = (*iterCellFieldData).cellFaceSideOff[cc] - (*iterFaceFieldData_).sideOff;
							if (std::fabs(match.length()) < parameter->EPS)
							{
								break;
							}
						}
					}
					if (cc >= (*iterCellFieldData).cellCellNumber + 1)
					{
						std::cout << "	error: fail to find the neighbor cell in reconstruction." << std::endl;
						exit(1);
					}

					// integrated functions
					tensor2D<real, NDOFS, NDOFS> *fII = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					tensor2D<real, NDOFS, NDOFS> *fIJ = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					tensor1D<real, NDOFS> *fI = new tensor1D<real, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix and vector on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						// point p = (*iterFaceFieldData_).gaussPairVector_[gg].p + (*iterFaceFieldData_).sideOff;//20200318��������ƫ����
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						// right cell ����ƫ����
						//  ----------------       --------
						//  |     ||      |		  |     ||
						//  |     ||	     |   ...  |     ||
						//  |     ||	     |		  |     ||
						//-----------------		 ---------
						//  |     ||      |		  |     ||
						//  | R   ||  L   |   ...  | R   ||
						//  |     ||	     |		  |     ||
						//-----------------		 ---------
						//  |     ||      |		  |     ||
						//  |     ||	     |   ...  |     ||
						//  |     ||	     |		  |     ||
						//-----------------		 ---------
						// cell R���ǰ���ˮƽ�������ڴ����ģ�����ı���Ӧ�ü���ƫ����
						//	point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;//20200318����
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getDiffBaseValueRBFB1(
							p,
							baryCenterI,
							scaleI,
							momentI,
							matrixDiffBaseI);
						// diff base matrix j
						tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseJ;
						// base moment j
						tensor1D<real, NDOFS> momentJ;
						iterCellFieldData_ = cellFieldData->begin() + cr - 1;
						point baryCenterJ = (*iterCellFieldData_).baryCenter;
						point scaleJ = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentJ[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getDiffBaseValueRBFB1(
							p + (*iterFaceFieldData_).sideOff, // 20200318����
							baryCenterJ,
							scaleJ,
							momentJ,
							matrixDiffBaseJ);
						// ll: row    rr: col
						// get the integrated function values on each face Gauss point
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// matrix
						for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
						{
							for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
							{
								fII[gg][ll][rr] = 0.0;
								fIJ[gg][ll][rr] = 0.0;
								// mapping col can be 0, but row cannot be 0
								for (int kk = 0; kk < static_cast<int>(NDIFFS); ++kk)
								{
									fII[gg][ll][rr] += (*iterFaceFieldData_).faceWeightVF[kk] * matrixDiffBaseI[ll][kk] * matrixDiffBaseI[rr][kk];
									fIJ[gg][ll][rr] += (*iterFaceFieldData_).faceWeightVF[kk] * matrixDiffBaseI[ll][kk] * matrixDiffBaseJ[rr][kk];
								}
							}
						}
						// vector
						for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
						{
							fI[gg][ll] = (*iterFaceFieldData_).faceWeightVF[0] * matrixDiffBaseI[ll][0];
						}
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii[ll][rr] += result; //+ means sum of each face

							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fIJ[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Bij[ll][rr] = result;

							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
					// get the face Gauss point for coefficient vector. integral along the interface between cell i and cell j
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real parametricArea;
						real result;
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							f[gg] = fI[gg][ll];
							weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
							cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
							parametricArea = (*iterFaceFieldData_).parametricArea;
						}
						gaussIntegralFace->getIntegral(
							static_cast<int>((*iterFaceFieldData_).fPG),
							f,
							weight,
							cofJacobi,
							parametricArea,
							result);
						bi[ll] = result;
						delete[] f;
						delete[] weight;
						delete[] cofJacobi;
						f = NULL;
						weight = NULL;
						cofJacobi = NULL;
					}
					// assign the matrix and vector with respect to each face
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							(*iterCellFieldData_).matrixBij[cc][ll][rr] = Bij[ll][rr];
						}
						(*iterCellFieldData_).vectorbij[cc][ll] = bi[ll];
					}
					delete[] fII;
					delete[] fIJ;
					delete[] fI;
					fII = NULL;
					fIJ = NULL;
					fI = NULL;
				}
				else if (cr == Wall)
				{
					tensor2D<real, NDOFS, NDOFS> *fII = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getDiffBaseValueRBFB1(
							p,
							baryCenterI,
							scaleI,
							momentI,
							matrixDiffBaseI);
						// ll: row    rr: col
						// get the integrated function values on each face Gauss point
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// matrix
						for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
						{
							for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
							{
								fII[gg][ll][rr] = 0.0;
								// mapping
								for (int kk = 0; kk < static_cast<int>(NDIFFS); ++kk)
								{
									fII[gg][ll][rr] += (*iterFaceFieldData_).faceWeightVF[kk] * matrixDiffBaseI[ll][kk] * matrixDiffBaseI[rr][kk];
								}
							}
						}
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii[ll][rr] += result; //+ means sum of each face
							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
					delete[] fII;
					fII = NULL;
				}
				else if (cr == FarField)
				{
					tensor2D<real, NDOFS, NDOFS> *fII = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getDiffBaseValueRBFB1(
							p,
							baryCenterI,
							scaleI,
							momentI,
							matrixDiffBaseI);
						// ll: row    rr: col
						// get the integrated function values on each face Gauss point
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// matrix
						for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
						{
							for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
							{
								fII[gg][ll][rr] = 0.0;
								// mapping
								for (int kk = 0; kk < static_cast<int>(NDIFFS); ++kk)
								{
									fII[gg][ll][rr] += (*iterFaceFieldData_).faceWeightVF[kk] * matrixDiffBaseI[ll][kk] * matrixDiffBaseI[rr][kk];
								}
							}
						}
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii[ll][rr] += result; //+ means sum of each face
							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
					delete[] fII;
					fII = NULL;
				}
				else if (cr == Symmetric)
				{
					tensor2D<real, NDOFS, NDOFS> *fII = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getDiffBaseValueRBFB1(
							p,
							baryCenterI,
							scaleI,
							momentI,
							matrixDiffBaseI);
						// ll: row    rr: col
						// get the integrated function values on each face Gauss point
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// matrix
						for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
						{
							for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
							{
								fII[gg][ll][rr] = 0.0;
								// mapping
								for (int kk = 0; kk < static_cast<int>(NDIFFS); ++kk)
								{
									fII[gg][ll][rr] += (*iterFaceFieldData_).faceWeightVF[kk] * matrixDiffBaseI[ll][kk] * matrixDiffBaseI[rr][kk];
								}
							}
						}
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii[ll][rr] += result; //+ means sum of each face
							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
					delete[] fII;
					fII = NULL;
				}
				else if (cr == OutFlow)
				{
					std::cout << "	warning: OutFlow is not proposed currently." << std::endl;
				}
				else if (cr == InFlow)
				{
					std::cout << "	warning: InFlow is not proposed currently." << std::endl;
				}
			}
			// iterCellFieldData_ = cellFieldData->begin() + cl - 1;
			// assign Aii
			for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
			{
				for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
				{
					(*iterCellFieldData).matrixAii[ll][rr] = Aii[ll][rr];
				}
			}
			tensor2D<real, NDOFS, NDOFS> inverseAii;
			CfvMath::getMatrixGeneralInverse(Aii, inverseAii);
			// assign inverse Aii
			for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
			{
				for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
				{
					(*iterCellFieldData).matrixAiiInverse[ll][rr] = inverseAii[ll][rr];
				}
			}
			// assign inverseAii*Bij & inverseAii*bi, for convenient calculation
			for (int cc = 1; cc < (*iterCellFieldData).cellCellNumber + 1; ++cc)
			{ // cellFaceNumber -> cellCellNumber
				// assign inverseAii*Bij
				for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
				{
					for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
					{
						(*iterCellFieldData).matrixAiiInverseBij[cc][ll][rr] = 0.0;
						for (int pp = 1; pp < static_cast<int>(NDOFS); ++pp)
						{
							(*iterCellFieldData).matrixAiiInverseBij[cc][ll][rr] +=
								inverseAii[ll][pp] * (*iterCellFieldData).matrixBij[cc][pp][rr];
						}
					}
				}
				// assign inverseAii*bi
				for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
				{
					(*iterCellFieldData).vectorAiiInversebij[cc][ll] = 0.0;
					for (int pp = 1; pp < static_cast<int>(NDOFS); ++pp)
					{
						(*iterCellFieldData).vectorAiiInversebij[cc][ll] +=
							inverseAii[ll][pp] * (*iterCellFieldData).vectorbij[cc][pp];
					}
				}
			}
		}
		std::cout << " ..reconstruction matrix and vector has been initialized." << std::endl;
		return true;
	}
	*/

	//----------------------------------------------------------------------------------
	// Upgrade V1
	template <unsigned int O>
	inline bool reconstructionRBFB1<O>::excuteReconstruction(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		cellGaussDataVector *cellGaussData,
		faceFieldDataVector *faceFieldData,
		faceGaussDataVector *faceGaussData,
		GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO1Grid<fO> *gaussIntegralFace)
	{
		//		std::cout << "excuting reconstruction ..." << std::endl;

		// #pragma omp parallel for schedule(guided)
		for (int iCell = 0; iCell < cellFieldData->size(); iCell++)
		{
			cellFieldDataVector::iterator iterCellFieldData;
			cellFieldDataVector::iterator iterCellFieldData_;
			faceFieldDataVector::iterator iterFaceFieldData_;
			iterCellFieldData = cellFieldData->begin() + iCell;
			if ((*iterCellFieldData).boundaryCellType_ == InnerCell)
			{
				// assign gradient values
				for (int ii = 1; ii < static_cast<int>(NDOFS); ++ii)
				{
					(*iterCellFieldData).scalarVariableTn[ii] =
						(1.0 - (*iterCellFieldData).relaxFactor) * (*iterCellFieldData).scalarVariableTn[ii];
				}
				// deal with Bij*uj
				for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber + 1; ++ii)
				{
					int cr = (*iterCellFieldData).cellCellIndex[ii];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					// CfvMath::VVMatVec(iterCellFieldData->matrixAiiInverseBij[ii], iterCellFieldData_->scalarVariableTn, iterCellFieldData->scalarVariableTn, iterCellFieldData->relaxFactor,
					// 				  1, NDOFS, 1, NDOFS, false);
					iterCellFieldData->scalarVariableTn += iterCellFieldData->relaxFactor *
														   (iterCellFieldData->matrixAiiInverseBij[ii] * iterCellFieldData_->scalarVariableTn);
				}
				// deal with bi
				for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber + 1; ++ii)
				{
					int cr = (*iterCellFieldData).cellCellIndex[ii];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					// CfvMath::VVVecAdd(iterCellFieldData->vectorAiiInversebij[ii], iterCellFieldData->scalarVariableTn,
					// 				  iterCellFieldData->relaxFactor * ((*iterCellFieldData_).scalarVariableTn[0] - (*iterCellFieldData).scalarVariableTn[0]),
					// 				  1, NDOFS);
					iterCellFieldData->scalarVariableTn += iterCellFieldData->vectorAiiInversebij.row(ii) *
														   iterCellFieldData->relaxFactor * ((*iterCellFieldData_).scalarVariableTn[0] - (*iterCellFieldData).scalarVariableTn[0]);
				}
			}
			else
			{
				// tensor1D<real, NDOFS> vectorBoundaryCorrection;
				Eigen::VectorXd vectorBoundaryCorrection(NDOFS);
				vectorBoundaryCorrection.setZero();
				for (int ii = 1; ii < (*iterCellFieldData).cellFaceNumber + 1; ++ii)
				{
					int ff = (*iterCellFieldData).cellFaceIndex[ii];
					iterFaceFieldData_ = faceFieldData->begin() + ff - 1;
					int cl = (*iterFaceFieldData_).faceCellIndex[1];
					int cr = (*iterFaceFieldData_).faceCellIndex[2];	  // for boundary face, 2 always correspond to R cell
					iterCellFieldData_ = cellFieldData->begin() + cl - 1; // 20200315!!

					tensor1D<real, NDOFS> *fI = new tensor1D<real, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					if (cr == Wall)
					{
						// get the vector on each face Gauss point.
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							auto &matrixDiffBaseII = iterFaceFieldData_->diffBaseValueData[0][gg];
							real uI = (*iterCellFieldData_).scalarVariableTn[0];
							real uBV = 0.0;
							// get the function values at Gauss points
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								// parametric
								fI[gg][ll] = std::pow((*iterFaceFieldData_).faceWeightVF[0], 2) * matrixDiffBaseII[ll][0] * (uBV - uI);
							}

							//	//--20200513
							//	//1st order derivative correction
							//	real duBVdx = 0.0;
							//	for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll){
							//		duBVdx += matrixDiffBaseI[ll][1] * (*iterCellFieldData_).scalarVariableTn[ll];
							//	}
							//	real duBVdy = 0.0;
							//	for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll){
							//		duBVdy += matrixDiffBaseI[ll][2] * (*iterCellFieldData_).scalarVariableTn[ll];
							//	}
							//	point uNV = (1.0 / (*iterFaceFieldData_).area) * (*iterFaceFieldData_).normalVector;
							//	real duBVdn = duBVdx * uNV.x + duBVdy * uNV.y;
							//	//	duBVdx = duBVdx - 2.0*duBVdn*uNV.x;
							//	//	duBVdy = duBVdy - 2.0*duBVdn*uNV.y;
							//	duBVdx = duBVdx;
							//	duBVdy = duBVdy;
							//	for (int ll = 1; ll < 3; ++ll){//static_cast<int>(NDOFS)
							//		fI[gg][ll] += (*iterFaceFieldData_).faceWeightVF[1] * matrixDiffBaseI[ll][1] * duBVdx;
							//	}
							//	for (int ll = 1; ll < 3; ++ll){//static_cast<int>(NDOFS)
							//		fI[gg][ll] += (*iterFaceFieldData_).faceWeightVF[2] * matrixDiffBaseI[ll][2] * duBVdy;
							//	}

							//
						}
					}
					else if (cr == FarField)
					{
						// get the vector on each face Gauss point.
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							auto &matrixDiffBaseII = iterFaceFieldData_->diffBaseValueData[0][gg];

							real uI = (*iterCellFieldData_).scalarVariableTn[0];
							real uBV = uI;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								uBV += matrixDiffBaseII[ll][0] * (*iterCellFieldData_).scalarVariableTn[ll];
							}
							// get the function values at Gauss points
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								// parametric
								fI[gg][ll] = std::pow((*iterFaceFieldData_).faceWeightVF[0], 2) * matrixDiffBaseII[ll][0] * (uBV - uI);
							}

							////--20200513
							////1st order derivative correction
							real duBVdx = 0.0;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								duBVdx += matrixDiffBaseII[ll][1] * (*iterCellFieldData_).scalarVariableTn[ll];
							}
							real duBVdy = 0.0;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								duBVdy += matrixDiffBaseII[ll][2] * (*iterCellFieldData_).scalarVariableTn[ll];
							}
							point uNV = (1.0 / (*iterFaceFieldData_).area) * (*iterFaceFieldData_).normalVector;
							real duBVdn = duBVdx * uNV.x + duBVdy * uNV.y;
							//	duBVdx = duBVdx - 2.0*duBVdn*uNV.x;
							//	duBVdy = duBVdy - 2.0*duBVdn*uNV.y;
							duBVdx = duBVdx;
							duBVdy = duBVdy;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								// parametric
								fI[gg][ll] += std::pow((*iterFaceFieldData_).faceWeightVF[1], 2) * matrixDiffBaseII[ll][1] * duBVdx;
							}
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								// parametric
								fI[gg][ll] += std::pow((*iterFaceFieldData_).faceWeightVF[2], 2) * matrixDiffBaseII[ll][2] * duBVdy;
							}

							//
						}
					}
					else if (cr == Symmetric)
					{
						// get the vector on each face Gauss point.
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							auto &matrixDiffBaseII = iterFaceFieldData_->diffBaseValueData[0][gg];

							real uI = (*iterCellFieldData_).scalarVariableTn[0];
							real uBV = uI;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								uBV += matrixDiffBaseII[ll][0] * (*iterCellFieldData_).scalarVariableTn[ll];
							}
							// get the function values at Gauss points
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								// parametric
								fI[gg][ll] = std::pow((*iterFaceFieldData_).faceWeightVF[0], 2) * matrixDiffBaseII[ll][0] * (uBV - uI);
							}
							// 1st order derivative correction
							real duBVdx = 0.0;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								duBVdx += matrixDiffBaseII[ll][1] * (*iterCellFieldData_).scalarVariableTn[ll];
							}
							real duBVdy = 0.0;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								duBVdy += matrixDiffBaseII[ll][2] * (*iterCellFieldData_).scalarVariableTn[ll];
							}
							point uNV = (1.0 / (*iterFaceFieldData_).area) * (*iterFaceFieldData_).normalVector;
							real duBVdn = duBVdx * uNV.x + duBVdy * uNV.y;
							duBVdx = duBVdx - 2.0 * duBVdn * uNV.x;
							duBVdy = duBVdy - 2.0 * duBVdn * uNV.y;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								// parametric
								fI[gg][ll] += std::pow((*iterFaceFieldData_).faceWeightVF[1], 2) * matrixDiffBaseII[ll][1] * duBVdx;
							}
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								// parametric
								fI[gg][ll] += std::pow((*iterFaceFieldData_).faceWeightVF[2], 2) * matrixDiffBaseII[ll][2] * duBVdy;
							}
						}
					}
					else if (cr == OutFlow)
					{
						std::cout << "	warning: OutFlow is not proposed currently." << std::endl;
					}
					else if (cr == InFlow)
					{
						std::cout << "	warning: InFlow is not proposed currently." << std::endl;
					}
					// get integral
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real parametricArea;
						real result;
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							f[gg] = fI[gg][ll];
							weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
							cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
							parametricArea = (*iterFaceFieldData_).parametricArea;
						}
						gaussIntegralFace->getIntegral(
							static_cast<int>((*iterFaceFieldData_).fPG),
							f,
							weight,
							cofJacobi,
							parametricArea,
							result);
						vectorBoundaryCorrection[ll] += result; //��һ����һ����
						delete[] f;
						delete[] weight;
						delete[] cofJacobi;
						f = NULL;
						weight = NULL;
						cofJacobi = NULL;
					}
					delete[] fI;
					fI = NULL;
				}
				// proceed with reconstruction
				// assign gradient values
				for (int ii = 1; ii < static_cast<int>(NDOFS); ++ii)
				{
					(*iterCellFieldData).scalarVariableTn[ii] =
						(1.0 - (*iterCellFieldData).relaxFactor) * (*iterCellFieldData).scalarVariableTn[ii];
				}
				// deal with Bij*uj
				for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber + 1; ++ii)
				{
					int cr = (*iterCellFieldData).cellCellIndex[ii];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;

					// CfvMath::VVMatVec(iterCellFieldData->matrixAiiInverseBij[ii], iterCellFieldData_->scalarVariableTn, iterCellFieldData->scalarVariableTn, iterCellFieldData->relaxFactor,
					// 				  1, NDOFS, 1, NDOFS, false);
					iterCellFieldData->scalarVariableTn += iterCellFieldData->relaxFactor *
														   (iterCellFieldData->matrixAiiInverseBij[ii] * iterCellFieldData_->scalarVariableTn);
				}
				// deal with bi
				for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber + 1; ++ii)
				{
					int cr = (*iterCellFieldData).cellCellIndex[ii];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					// CfvMath::VVVecAdd(iterCellFieldData->vectorAiiInversebij[ii], iterCellFieldData->scalarVariableTn,
					// 				  iterCellFieldData->relaxFactor * ((*iterCellFieldData_).scalarVariableTn[0] - (*iterCellFieldData).scalarVariableTn[0]),
					// 				  1, NDOFS);
					iterCellFieldData->scalarVariableTn +=
						iterCellFieldData->vectorAiiInversebij.row(ii) *
						iterCellFieldData->relaxFactor * ((*iterCellFieldData_).scalarVariableTn[0] - (*iterCellFieldData).scalarVariableTn[0]);
				}
				// deal with boundary correction
				// CfvMath::VVMatVec(iterCellFieldData->matrixAiiInverse, vectorBoundaryCorrection, iterCellFieldData->scalarVariableTn, iterCellFieldData->relaxFactor,
				// 				  1, NDOFS, 1, NDOFS, false);
				iterCellFieldData->scalarVariableTn += iterCellFieldData->relaxFactor *
													   iterCellFieldData->matrixAiiInverse *
													   vectorBoundaryCorrection;

				// // RBFB1:
				// //	attempt a restriction
				// switch (NDOFS)
				// {
				// case 7:
				// 	for (int i = 0; i < 4; i++)
				// 		iterCellFieldData->scalarVariableTn[4 + i] *= 0.2;
				// 	// real Mvalue{0};
				// 	// for (int i = 0; i < 4; i++)
				// 	// 	Mvalue += iterCellFieldData->scalarVariableTn[4 + i];
				// 	// for (int i = 0; i < 4; i++)
				// 	// 	iterCellFieldData->scalarVariableTn[4 + i] = iterCellFieldData->scalarVariableTn[4 + i] * Mvalue > 0 ? iterCellFieldData->scalarVariableTn[4 + i] : 0;
				// 	// break;

				// case 3:
				// case 6:
				// 	break;

				// default:
				// 	assert(false);
				// 	break;
				// }
			}
		}
		//		std::cout << " ..reconstruction has been completed." << std::endl;
		return true;
	}

	// Upgrade V1
	/*
	template <unsigned int O>
	bool reconstructionRBFB1<O>::excuteReconstruction_NotReady(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		cellGaussDataVector *`+3,
		faceFieldDataVector *faceFieldData,
		faceGaussDataVector *faceGaussData,
		GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO2Grid<fO> *gaussIntegralFace)
	{
		//		std::cout << "excuting reconstruction ..." << std::end

		cellFieldDataVector::iterator iterCellFieldData;
		cellFieldDataVector::iterator iterCellFieldData_;
		faceFieldDataVector::iterator iterFaceFieldData_;
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			if ((*iterCellFieldData).boundaryCellType_ == InnerCell)
			{
				// assign gradient values
				for (int ii = 1; ii < static_cast<int>(NDOFS); ++ii)
				{
					(*iterCellFieldData).scalarVariableTn[ii] =
						(1.0 - (*iterCellFieldData).relaxFactor) * (*iterCellFieldData).scalarVariableTn[ii];
				}
				// deal with Bij*uj
				for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber + 1; ++ii)
				{
					int cr = (*iterCellFieldData).cellCellIndex[ii];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					for (int jj = 1; jj < static_cast<int>(NDOFS); ++jj)
					{
						real sumLocal = 0.0;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							sumLocal += (*iterCellFieldData).matrixAiiInverseBij[ii][jj][kk] * (*iterCellFieldData_).scalarVariableTn[kk];
						}
						(*iterCellFieldData).scalarVariableTn[jj] += (*iterCellFieldData).relaxFactor * sumLocal;
						//+ means sum of neighbour cells
						//---------------------------------
						//
						//	real sumLocal2 = (*iterCellFieldData).vectorAiiInversebij[ii][jj] *
						//		((*iterCellFieldData_).scalarVariableTn[0] - (*iterCellFieldData).scalarVariableTn[0]);
						//	(*iterCellFieldData).scalarVariableTn[jj] += (*iterCellFieldData).relaxFactor * sumLocal2;
						//
						//---------------------------------
					}
				}
				// deal with bi
				for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber + 1; ++ii)
				{
					int cr = (*iterCellFieldData).cellCellIndex[ii];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					for (int jj = 1; jj < static_cast<int>(NDOFS); ++jj)
					{
						real sumLocal = (*iterCellFieldData).vectorAiiInversebij[ii][jj] *
										((*iterCellFieldData_).scalarVariableTn[0] - (*iterCellFieldData).scalarVariableTn[0]);
						(*iterCellFieldData).scalarVariableTn[jj] += (*iterCellFieldData).relaxFactor * sumLocal; //+ means sum of neighbour cells
					}
				}
			}
			else
			{
				tensor1D<real, NDOFS> vectorBoundaryCorrection;
				for (int ii = 1; ii < (*iterCellFieldData).cellFaceNumber + 1; ++ii)
				{
					int ff = (*iterCellFieldData).cellFaceIndex[ii];
					iterFaceFieldData_ = faceFieldData->begin() + ff - 1;
					int cl = (*iterFaceFieldData_).faceCellIndex[1];
					int cr = (*iterFaceFieldData_).faceCellIndex[2];	  // for boundary face, 2 always correspond to R cell
					iterCellFieldData_ = cellFieldData->begin() + cl - 1; // 20200315!!

					tensor1D<real, NDOFS> *fI = new tensor1D<real, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					if (cr == Wall)
					{
						// get the vector on each face Gauss point.
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							// base moment i
							tensor1D<real, NDOFS> momentI;
							// base i
							tensor1D<real, NDOFS> vectorBaseI;
							point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
							point baryCenterI = (*iterCellFieldData_).baryCenter;
							point scaleI = (*iterCellFieldData_).lengthReference;
							for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							{
								momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
							}
							CfvMath::getBaseValueRBFB1(
								p,
								baryCenterI,
								scaleI,
								momentI,
								vectorBaseI);
							real uI = (*iterCellFieldData_).scalarVariableTn[0];
							real uBV = 0.0;
							// get the function values at Gauss points
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								fI[gg][ll] = (*iterFaceFieldData_).faceWeightVF[0] * vectorBaseI[ll] * (uBV - uI);
							}
						}
					}
					else if (cr == FarField)
					{
						// get the vector on each face Gauss point.
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							// base moment i
							tensor1D<real, NDOFS> momentI;
							// base i
							tensor1D<real, NDOFS> vectorBaseI;
							point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
							point baryCenterI = (*iterCellFieldData_).baryCenter;
							point scaleI = (*iterCellFieldData_).lengthReference;
							for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							{
								momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
							}
							CfvMath::getBaseValueRBFB1(
								p,
								baryCenterI,
								scaleI,
								momentI,
								vectorBaseI);
							real uI = (*iterCellFieldData_).scalarVariableTn[0];
							real uBV = uI;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								uBV += vectorBaseI[ll] * (*iterCellFieldData_).scalarVariableTn[ll];
							}
							// get the function values at Gauss points
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								fI[gg][ll] = (*iterFaceFieldData_).faceWeightVF[0] * vectorBaseI[ll] * (uBV - uI);
							}
						}
					}
					else if (cr == Symmetric)
					{
						// get the vector on each face Gauss point.
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							// base moment i
							tensor1D<real, NDOFS> momentI;
							// base i
							tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
							point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;

							point baryCenterI = (*iterCellFieldData_).baryCenter;
							point scaleI = (*iterCellFieldData_).lengthReference;
							for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							{
								momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
							}
							CfvMath::getDiffBaseValueRBFB1(
								p,
								baryCenterI,
								scaleI,
								momentI,
								matrixDiffBaseI);
							real uI = (*iterCellFieldData_).scalarVariableTn[0];
							real uBV = uI;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								uBV += matrixDiffBaseI[ll][0] * (*iterCellFieldData_).scalarVariableTn[ll];
							}
							// get the function values at Gauss points
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								fI[gg][ll] = (*iterFaceFieldData_).faceWeightVF[0] * matrixDiffBaseI[ll][0] * (uBV - uI);
							}
							// 1st order derivative correction
							real duBVdx = 0.0;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								duBVdx += matrixDiffBaseI[ll][1] * (*iterCellFieldData_).scalarVariableTn[ll];
							}
							real duBVdy = 0.0;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								duBVdy += matrixDiffBaseI[ll][2] * (*iterCellFieldData_).scalarVariableTn[ll];
							}
							point uNV = (1.0 / (*iterFaceFieldData_).area) * (*iterFaceFieldData_).normalVector;
							real duBVdn = duBVdx * uNV.x + duBVdy * uNV.y;
							duBVdx = duBVdx - 2.0 * duBVdn * uNV.x;
							duBVdy = duBVdy - 2.0 * duBVdn * uNV.y;
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								fI[gg][ll] += (*iterFaceFieldData_).faceWeightVF[1] * matrixDiffBaseI[ll][1] * duBVdx;
							}
							for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
							{
								fI[gg][ll] += (*iterFaceFieldData_).faceWeightVF[2] * matrixDiffBaseI[ll][2] * duBVdy;
							}
						}
					}
					else if (cr == OutFlow)
					{
						std::cout << "	warning: OutFlow is not proposed currently." << std::endl;
					}
					else if (cr == InFlow)
					{
						std::cout << "	warning: InFlow is not proposed currently." << std::endl;
					}
					// get integral
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real parametricArea;
						real result;
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							f[gg] = fI[gg][ll];
							weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
							cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
							parametricArea = (*iterFaceFieldData_).parametricArea;
						}
						gaussIntegralFace->getIntegral(
							static_cast<int>((*iterFaceFieldData_).fPG),
							f,
							weight,
							cofJacobi,
							parametricArea,
							result);
						vectorBoundaryCorrection[ll] += result; //��һ����һ����
						delete[] f;
						delete[] weight;
						delete[] cofJacobi;
						f = NULL;
						weight = NULL;
						cofJacobi = NULL;
					}
					delete[] fI;
					fI = NULL;
				}
				// proceed with reconstruction
				// assign gradient values
				for (int ii = 1; ii < static_cast<int>(NDOFS); ++ii)
				{
					(*iterCellFieldData).scalarVariableTn[ii] =
						(1.0 - (*iterCellFieldData).relaxFactor) * (*iterCellFieldData).scalarVariableTn[ii];
				}
				// deal with Bij*uj
				for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber + 1; ++ii)
				{
					int cr = (*iterCellFieldData).cellCellIndex[ii];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					for (int jj = 1; jj < static_cast<int>(NDOFS); ++jj)
					{
						real sumLocal = 0.0;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							sumLocal += (*iterCellFieldData).matrixAiiInverseBij[ii][jj][kk] * (*iterCellFieldData_).scalarVariableTn[kk];
						}
						(*iterCellFieldData).scalarVariableTn[jj] += (*iterCellFieldData).relaxFactor * sumLocal;
					}
				}
				// deal with bi
				for (int ii = 1; ii < (*iterCellFieldData).cellCellNumber + 1; ++ii)
				{
					int cr = (*iterCellFieldData).cellCellIndex[ii];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					for (int jj = 1; jj < static_cast<int>(NDOFS); ++jj)
					{
						real sumLocal = (*iterCellFieldData).vectorAiiInversebij[ii][jj] *
										((*iterCellFieldData_).scalarVariableTn[0] - (*iterCellFieldData).scalarVariableTn[0]);
						(*iterCellFieldData).scalarVariableTn[jj] += (*iterCellFieldData).relaxFactor * sumLocal;
					}
				}
				// deal with boundary correction
				for (int jj = 1; jj < static_cast<int>(NDOFS); ++jj)
				{
					real sumLocal = 0.0;
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
					{
						sumLocal += (*iterCellFieldData).matrixAiiInverse[jj][kk] * vectorBoundaryCorrection[kk];
					}
					(*iterCellFieldData).scalarVariableTn[jj] += (*iterCellFieldData).relaxFactor * sumLocal;
				}
			}
		}
		//		std::cout << " ..reconstruction has been completed." << std::endl;
		return true;
	}
	*/

	//----------------------------------------------------------------------------------
	template <unsigned int O>
	inline bool reconstructionRBFB1<O>::initReconstructionMatrixAndVectorCR(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		cellGaussDataVector *cellGaussData,
		faceFieldDataVector *faceFieldData,
		faceGaussDataVector *faceGaussData,
		GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO1Grid<fO> *gaussIntegralFace)
	{
		//		std::cout << "initializing CR reconstruction matrix and vector..." << std::endl;

		cellFieldDataVector::iterator iterCellFieldData;
		cellFieldDataVector::iterator iterCellFieldData_;
		faceFieldDataVector::iterator iterFaceFieldData_;
		// matrix A
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			// tensor2D<real, NDOFSCR, NDOFSCR> Aii;
			Eigen::MatrixXd Aii(NDOFSCR, NDOFSCR);
			Aii.setZero();
			for (int ff = 1; ff < (*iterCellFieldData).cellFaceNumber + 1; ++ff)
			{ //ע��ѭ�����ޣ���20200314
				int kf = (*iterCellFieldData).cellFaceIndex[ff];
				iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
				int cl = (*iterFaceFieldData_).faceCellIndex[1];
				int cr = (*iterFaceFieldData_).faceCellIndex[2];
				// alternate cell��ֻ���ڲ���Ԫ���ܻ���Ҫ
				if (cr == (*iterCellFieldData).index)
				{
					int ck = cr;
					cr = cl; //�߽粻�õ���
					cl = ck;
				}

				//�����
				iterCellFieldData_ = cellFieldData->begin() + cl - 1;
				if (cr > 0 || cr == Wall || cr == FarField || Symmetric)
				{
					tensor2D<real, NDOFSCR, NDOFSCR> *fII = new tensor2D<real, NDOFSCR, NDOFSCR>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix and vector on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
#ifdef RBFB1_CR_DIFFBOUND
						tensor2D<real, NDOFSCR, NDIFFSCR> vectorBaseI;
#else
						tensor1D<real, NDOFSCR> vectorBaseI;
#endif
						// base moment i
						tensor1D<real, NDOFSCR> momentI;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						point pparaml;
						pparaml = CfvMath::GetFaceParam((*cellFieldData)[cl - 1].cellType_, ff, iterFaceFieldData_->parametricValue[gg].first);
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						// CfvMath::getBaseValue(
						// 	p,
						// 	baryCenterI,
						// 	scaleI,
						// 	momentI,
						// 	vectorBaseI);
#ifdef RBFB1_CR_DIFFBOUND
						RBFB1GetDiffBaseValueCR(
							pparaml,
							baryCenterI,
							scaleI,
							momentI,
							vectorBaseI,
							*iterCellFieldData_
#ifdef RBFB1_CR_INTERPOLATE
							,
							*faceFieldData
#endif
						);
#else
						RBFB1GetBaseValueCR(
							pparaml,
							baryCenterI,
							scaleI,
							momentI,
							vectorBaseI,
							*iterCellFieldData_
#ifdef RBFB1_CR_INTERPOLATE
							,
							*faceFieldData
#endif
						);
#endif
						// ll: row    rr: col
						// get the integrated function values on each face Gauss point
						for (int ll = 1; ll < static_cast<int>(NDOFSCR); ++ll)
						{
							for (int rr = 1; rr < static_cast<int>(NDOFSCR); ++rr)
							{
#ifdef RBFB1_CR_DIFFBOUND
								fII[gg][ll][rr] = std::pow((*iterFaceFieldData_).faceWeightVF[0], 2) * vectorBaseI[ll][0] * vectorBaseI[rr][0] +
												  (vectorBaseI[ll][1] * (*iterFaceFieldData_).faceWeightVF[1] + vectorBaseI[ll][2] * (*iterFaceFieldData_).faceWeightVF[2]) *
													  (vectorBaseI[rr][1] * (*iterFaceFieldData_).faceWeightVF[1] + vectorBaseI[rr][2] * (*iterFaceFieldData_).faceWeightVF[2]) * RBFB1_CR_DIFFDIM;
#else
								// parametric
								fII[gg][ll][rr] = std::pow((*iterFaceFieldData_).faceWeightVF[0], 2) * vectorBaseI[ll] * vectorBaseI[rr];
#endif
							}
						}
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFSCR); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFSCR); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii(ll, rr) += result; //+ means sum of each face
							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
					delete[] fII;
					fII = NULL;
				}
				else if (cr == OutFlow)
				{
					std::cout << "	warning: OutFlow is not proposed currently." << std::endl;
				}
				else if (cr == InFlow)
				{
					std::cout << "	warning: InFlow is not proposed currently." << std::endl;
				}
			}
			// iterCellFieldData_ = cellFieldData->begin() + cl - 1;
			// assign Aii
			// redundancy //iterCellFieldData->matrixAiiCR = Aii;
			// CfvMath::VVMatCopy(Aii, iterCellFieldData->matrixAiiCR, 1, NDOFSCR, 1, NDOFSCR);
			CfvMath::EigenLeastSquareInverse(Aii, iterCellFieldData->matrixAiiInverseCR);
			// tensor2D<real, NDOFSCR, NDOFSCR> inverseAii;
			// CfvMath::getMatrixGeneralInverse(Aii, inverseAii);
			// assign inverse Aii
			// CfvMath::VVMatCopy(inverseAii, iterCellFieldData->matrixAiiInverseCR, 1, NDOFSCR, 1, NDOFSCR);
		}
		//		std::cout << " ..CR reconstruction matrix and vector has been initialized." << std::endl;
		return true;
	}

	/*
	template <unsigned int O>
	bool reconstructionRBFB1<O>::initReconstructionMatrixAndVectorCR_NotReady(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		cellGaussDataVector *cellGaussData,
		faceFieldDataVector *faceFieldData,
		faceGaussDataVector *faceGaussData,
		GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO2Grid<fO> *gaussIntegralFace)
	{
		//		std::cout << "initializing CR reconstruction matrix and vector..." << std::endl;

		cellFieldDataVector::iterator iterCellFieldData;
		cellFieldDataVector::iterator iterCellFieldData_;
		faceFieldDataVector::iterator iterFaceFieldData_;
		// matrix A
		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			tensor2D<real, NDOFS, NDOFS> Aii;
			for (int ff = 1; ff < (*iterCellFieldData).cellFaceNumber + 1; ++ff)
			{ //ע��ѭ�����ޣ���20200314
				int kf = (*iterCellFieldData).cellFaceIndex[ff];
				iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
				int cl = (*iterFaceFieldData_).faceCellIndex[1];
				int cr = (*iterFaceFieldData_).faceCellIndex[2];
				// alternate cell��ֻ���ڲ���Ԫ���ܻ���Ҫ
				if (cr == (*iterCellFieldData).index)
				{
					int ck = cr;
					cr = cl; //�߽粻�õ���
					cl = ck;
				}

				//�����
				iterCellFieldData_ = cellFieldData->begin() + cl - 1;
				if (cr > 0)
				{
					// find the neibhour cell index
					int cc;
					for (cc = 1; cc < (*iterCellFieldData).cellCellNumber + 1; ++cc)
					{
						if (cr == (*iterCellFieldData_).cellCellIndex[cc])
						{
							point match = (*iterCellFieldData).cellFaceSideOff[cc] - (*iterFaceFieldData_).sideOff;
							if (std::fabs(match.length()) < parameter->EPS)
							{
								break;
							}
						}
					}
					if (cc >= (*iterCellFieldData).cellCellNumber + 1)
					{
						std::cout << "	error: fail to find the neighbor cell in reconstruction." << std::endl;
						exit(1);
					}

					// integrated functions
					tensor2D<real, NDOFS, NDOFS> *fII = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix and vector on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor1D<real, NDOFS> vectorBaseI;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getBaseValueRBFB1(
							p,
							baryCenterI,
							scaleI,
							momentI,
							vectorBaseI);
						// matrix
						for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
						{
							for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
							{
								fII[gg][ll][rr] = (*iterFaceFieldData_).faceWeightVF[0] * vectorBaseI[ll] * vectorBaseI[rr];
							}
						}
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii[ll][rr] += result; //+ means sum of each face
							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
				}
				else if (cr == Wall)
				{
					tensor2D<real, NDOFS, NDOFS> *fII = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor1D<real, NDOFS> vectorBaseI;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getBaseValueRBFB1(
							p,
							baryCenterI,
							scaleI,
							momentI,
							vectorBaseI);
						// ll: row    rr: col
						// get the integrated function values on each face Gauss point
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// matrix
						for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
						{
							for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
							{
								fII[gg][ll][rr] = (*iterFaceFieldData_).faceWeightVF[0] * vectorBaseI[ll] * vectorBaseI[rr];
							}
						}
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii[ll][rr] += result; //+ means sum of each face
							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
					delete[] fII;
					fII = NULL;
				}
				else if (cr == FarField)
				{
					tensor2D<real, NDOFS, NDOFS> *fII = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor1D<real, NDOFS> vectorBaseI;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getBaseValueRBFB1(
							p,
							baryCenterI,
							scaleI,
							momentI,
							vectorBaseI);
						// ll: row    rr: col
						// get the integrated function values on each face Gauss point
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// matrix
						for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
						{
							for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
							{
								fII[gg][ll][rr] = (*iterFaceFieldData_).faceWeightVF[0] * vectorBaseI[ll] * vectorBaseI[rr];
							}
						}
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii[ll][rr] += result; //+ means sum of each face
							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
					delete[] fII;
					fII = NULL;
				}
				else if (cr == Symmetric)
				{
					tensor2D<real, NDOFS, NDOFS> *fII = new tensor2D<real, NDOFS, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the matrix on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor1D<real, NDOFS> vectorBaseI;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getBaseValueRBFB1(
							p,
							baryCenterI,
							scaleI,
							momentI,
							vectorBaseI);
						// ll: row    rr: col
						// get the integrated function values on each face Gauss point
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// matrix
						for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
						{
							for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
							{
								fII[gg][ll][rr] = (*iterFaceFieldData_).faceWeightVF[0] * vectorBaseI[ll] * vectorBaseI[rr];
							}
						}
					}
					// get the face Gauss point for coefficient matrix. integral along the interface between cell i and cell j
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
						{
							real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
							real parametricArea;
							real result;
							for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
							{
								f[gg] = fII[gg][ll][rr];
								weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
								cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
								parametricArea = (*iterFaceFieldData_).parametricArea;
							}
							gaussIntegralFace->getIntegral(
								static_cast<int>((*iterFaceFieldData_).fPG),
								f,
								weight,
								cofJacobi,
								parametricArea,
								result);
							Aii[ll][rr] += result; //+ means sum of each face
							delete[] f;
							delete[] weight;
							delete[] cofJacobi;
							f = NULL;
							weight = NULL;
							cofJacobi = NULL;
						}
					}
					delete[] fII;
					fII = NULL;
				}
				else if (cr == OutFlow)
				{
					std::cout << "	warning: OutFlow is not proposed currently." << std::endl;
				}
				else if (cr == InFlow)
				{
					std::cout << "	warning: InFlow is not proposed currently." << std::endl;
				}
			}
			// iterCellFieldData_ = cellFieldData->begin() + cl - 1;
			// assign Aii
			for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
			{
				for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
				{
					(*iterCellFieldData).matrixAiiCR[ll][rr] = Aii[ll][rr];
				}
			}
			tensor2D<real, NDOFS, NDOFS> inverseAii;
			CfvMath::getMatrixGeneralInverse(Aii, inverseAii);
			// assign inverse Aii
			for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
			{
				for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
				{
					(*iterCellFieldData).matrixAiiInverseCR[ll][rr] = inverseAii[ll][rr];
				}
			}
		}
		//		std::cout << " ..CR reconstruction matrix and vector has been initialized." << std::endl;
		return true;
	}
	*/

	//----------------------------------------------------------------------------------
	// Update V1
	template <unsigned int O>
	inline bool reconstructionRBFB1<O>::excuteReconstructionCR(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		cellGaussDataVector *cellGaussData,
		faceFieldDataVector *faceFieldData,
		faceGaussDataVector *faceGaussData,
		GaussIntegralCellO1Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO1Grid<fO> *gaussIntegralFace)
	{

#pragma omp parallel for schedule(guided)
		for (int iCell = 0; iCell < cellFieldData->size(); iCell++)
		{
			cellFieldDataVector::iterator iterCellFieldData;
			cellFieldDataVector::iterator iterCellFieldData_;
			faceFieldDataVector::iterator iterFaceFieldData_;
			iterCellFieldData = cellFieldData->begin() + iCell;
			// tensor1D<real, NDOFSCR> bi;
			Eigen::VectorXd bi(NDOFSCR);
			bi.setZero();
			// deal with boundary
			for (int ff = 1; ff < (*iterCellFieldData).cellFaceNumber + 1; ++ff)
			{ //ע��ѭ�����ޣ���20200314
				int kf = (*iterCellFieldData).cellFaceIndex[ff];
				iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
				int cl = (*iterFaceFieldData_).faceCellIndex[1];
				int cr = (*iterFaceFieldData_).faceCellIndex[2];
				// alternate cell��ֻ���ڲ���Ԫ���ܻ���Ҫ
				if (cr == (*iterCellFieldData).index)
				{
					int ck = cr;
					cr = cl; //�߽粻�õ���
					cl = ck;
					// std::cout << "		alternate left and right cell." << std::endl;
				}
				real *fIL = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
				tensor1D<real, NDOFSCR> *fIBase = new tensor1D<real, NDOFSCR>[static_cast<int>((*iterFaceFieldData_).fPG)];
#ifdef RBFB1_CR_DIFFBOUND_NO
				real *fILDx = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
				real *fILDy = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
				tensor1D<real, NDOFSCR> *fIBaseDx = new tensor1D<real, NDOFSCR>[static_cast<int>((*iterFaceFieldData_).fPG)];
				tensor1D<real, NDOFSCR> *fIBaseDy = new tensor1D<real, NDOFSCR>[static_cast<int>((*iterFaceFieldData_).fPG)];
#endif
				if (cr > 0)
				{
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						fIL[gg] = 0.0;
				}
				else if (cr == Wall)
				{
					// get the function values on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// // diff base matrix i
						// tensor1D<real, NDOFS> vectorBaseI;
						// base moment i
						// tensor1D<real, NDOFS> momentI;
						// point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						// point pparaml = CfvMath::GetFaceParam((*cellFieldData)[cl - 1].cellType_, ff, iterFaceFieldData_->parametricValue[gg].first);
						// iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// point baryCenterI = (*iterCellFieldData_).baryCenter;
						// point scaleI = (*iterCellFieldData_).lengthReference;
						// for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						// {
						// 	momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						// }
						// CfvMath::getBaseValueRBFA1(
						// 	p,
						// 	baryCenterI,
						// 	scaleI,
						// 	momentI,
						// 	vectorBaseI);
						// CfvMath::getBaseValueRBFB1(
						// 	pparaml,
						// 	baryCenterI,
						// 	scaleI,
						// 	momentI,
						// 	vectorBaseI,
						// 	*iterCellFieldData_);

						auto &vectorBaseIICR = iterFaceFieldData_->diffBaseValueDataCR[0][gg];

						fIL[gg] = 0.0;
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
							fIBase[gg][kk] = vectorBaseIICR[kk][0];

#ifdef RBFB1_CR_DIFFBOUND_NO
						auto &vectorBaseII = iterFaceFieldData_->diffBaseValueData[0][gg];
						fILDx[gg] = 0., fILDy[gg] = 0.;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							fILDx[gg] += vectorBaseII[kk][1] * (*iterCellFieldData_).scalarVariableTn[kk],
								fILDy[gg] += vectorBaseII[kk][2] * (*iterCellFieldData_).scalarVariableTn[kk];
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
							fIBaseDx[gg][kk] = vectorBaseIICR[kk][1], fIBaseDy[gg][kk] = vectorBaseIICR[kk][2];
#endif
					}
				}
				else if (cr == FarField)
				{
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						// tensor1D<real, NDOFS> vectorBaseI;
						// base moment i
						// tensor1D<real, NDOFS> momentI;
						// point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						// point pparaml = CfvMath::GetFaceParam((*cellFieldData)[cl - 1].cellType_, ff, iterFaceFieldData_->parametricValue[gg].first);
						// iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// point baryCenterI = (*iterCellFieldData_).baryCenter;
						// point scaleI = (*iterCellFieldData_).lengthReference;
						// for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						// {
						// 	momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						// }
						// CfvMath::getBaseValue(
						// 	p,
						// 	baryCenterI,
						// 	scaleI,
						// 	momentI,
						// 	vectorBaseI);
						// CfvMath::getBaseValueRBFB1(
						// 	pparaml,
						// 	baryCenterI,
						// 	scaleI,
						// 	momentI,
						// 	vectorBaseI,
						// 	*iterCellFieldData_);
						auto &vectorBaseIICR = iterFaceFieldData_->diffBaseValueDataCR[0][gg];
						fIL[gg] = (*iterCellFieldData_).scalarVariableTn[0];
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
							fIBase[gg][kk] = vectorBaseIICR[kk][0];

#ifdef RBFB1_CR_DIFFBOUND_NO
						auto &vectorBaseII = iterFaceFieldData_->diffBaseValueData[0][gg];
						fILDx[gg] = 0., fILDy[gg] = 0.;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							fILDx[gg] += vectorBaseII[kk][1] * (*iterCellFieldData_).scalarVariableTn[kk],
								fILDy[gg] += vectorBaseII[kk][2] * (*iterCellFieldData_).scalarVariableTn[kk];
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
							fIBaseDx[gg][kk] = vectorBaseIICR[kk][1], fIBaseDy[gg][kk] = vectorBaseIICR[kk][2];
#endif
					}
				}
				else if (cr == Symmetric)
				{
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						// tensor1D<real, NDOFS> vectorBaseI;
						// base moment i
						// tensor1D<real, NDOFS> momentI;
						// point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						// point pparaml = CfvMath::GetFaceParam((*cellFieldData)[cl - 1].cellType_, ff, iterFaceFieldData_->parametricValue[gg].first);
						// iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						// point baryCenterI = (*iterCellFieldData_).baryCenter;
						// point scaleI = (*iterCellFieldData_).lengthReference;
						// for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						// {
						// 	momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						// }
						// CfvMath::getBaseValue(
						// 	p,
						// 	baryCenterI,
						// 	scaleI,
						// 	momentI,
						// 	vectorBaseI);
						// CfvMath::getBaseValueRBFB1(
						// 	pparaml,
						// 	baryCenterI,
						// 	scaleI,
						// 	momentI,
						// 	vectorBaseI,
						// 	*iterCellFieldData_);
						auto &vectorBaseII = iterFaceFieldData_->diffBaseValueData[0][gg];
						auto &vectorBaseIICR = iterFaceFieldData_->diffBaseValueDataCR[0][gg];

						fIL[gg] = (*iterCellFieldData_).scalarVariableTn[0];
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							fIL[gg] += vectorBaseII[kk][0] * (*iterCellFieldData_).scalarVariableTn[kk];
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
							fIBase[gg][kk] = vectorBaseIICR[kk][0];

#ifdef RBFB1_CR_DIFFBOUND_NO
						fILDx[gg] = 0., fILDy[gg] = 0.;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							fILDx[gg] += vectorBaseII[kk][1] * (*iterCellFieldData_).scalarVariableTn[kk],
								fILDy[gg] += vectorBaseII[kk][2] * (*iterCellFieldData_).scalarVariableTn[kk];
						for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
							fIBaseDx[gg][kk] = vectorBaseIICR[kk][1], fIBaseDy[gg][kk] = vectorBaseIICR[kk][2];
#endif
					}
				}
				else if (cr == OutFlow)
				{
					std::cout << "	warning: OutFlow is not proposed currently." << std::endl;
				}
				else if (cr == InFlow)
				{
					std::cout << "	warning: InFlow is not proposed currently." << std::endl;
				}

				// get the RHS term
				iterCellFieldData_ = cellFieldData->begin() + cl - 1;
				real valueL = (*iterCellFieldData_).scalarVariableTn[0];
				for (int ll = 1; ll < static_cast<int>(NDOFSCR); ++ll)
				{
					real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					real parametricArea;
					real result;
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
#ifdef RBFB1_INCREMENT_CR
						f[gg] = 0.;
#else
						f[gg] = (fIL[gg] - valueL) * fIBase[gg][ll] * std::pow((*iterFaceFieldData_).faceWeightVF[0], 2);
#endif
						// integrated function value
						weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
						cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
						parametricArea = (*iterFaceFieldData_).parametricArea;
					}
					gaussIntegralFace->getIntegral(
						static_cast<int>((*iterFaceFieldData_).fPG),
						f,
						weight,
						cofJacobi,
						parametricArea,
						result);
					bi[ll] += result; //+ means sum of each face
					delete[] f;
					delete[] weight;
					delete[] cofJacobi;
					f = NULL;
					weight = NULL;
					cofJacobi = NULL;
				}
				delete[] fIL;
				delete[] fIBase;
#ifdef RBFB1_CR_DIFFBOUND_NO
				delete[] fILDx;
				delete[] fILDy;
				delete[] fIBaseDx;
				delete[] fIBaseDy;
#endif
				fIBase = NULL;
				fIL = NULL;
			}
			// deal with inner
			for (int cc = 1; cc < (*iterCellFieldData).cellCellNumber + 1; ++cc)
			{
				int cl = (*iterCellFieldData).cellCellIndex[0];
				int cr = (*iterCellFieldData).cellCellIndex[cc];
				int kf;
				//ȷ����
				int ff, ffr;
				for (ff = 1; ff < (*iterCellFieldData).cellFaceNumber + 1; ++ff)
				{
					kf = (*iterCellFieldData).cellFaceIndex[ff];
					iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
					int clMatch = (*iterFaceFieldData_).faceCellIndex[1];
					int crMatch = (*iterFaceFieldData_).faceCellIndex[2];
					if (cr == clMatch || cr == crMatch)
						break;
				}
				for (ffr = 1; ffr < (*cellFieldData)[cr - 1].cellFaceNumber + 1; ffr++)
					if ((*cellFieldData)[cr - 1].cellFaceIndex[ffr] == kf)
						break;
				assert(ffr < (*cellFieldData)[cr - 1].cellFaceNumber + 1);

				iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
				int cff = 0, cffr = 1;
				if (iterFaceFieldData_->faceCellIndex[1] == cr)
					cff = 1, cffr = 0;

				point p1 = (*iterFaceFieldData_).faceNode[1].second;
				point p2 = (*iterFaceFieldData_).faceNode[2].second;
				point pMid = 0.5 * (p1 + p2);
				point pparaml = CfvMath::GetFaceParam((*cellFieldData)[cl - 1].cellType_, ff, point(0.5, 0));
				point pparamr = CfvMath::GetFaceParam((*cellFieldData)[cr - 1].cellType_, ffr, point(0.5, 0), true);
				point outDirection = pMid - (*iterCellFieldData).baryCenter;
				point normalVector = getInnerProduct(outDirection, (*iterFaceFieldData_).normalVector) * (*iterFaceFieldData_).normalVector;
				point uNV = (1.0 / (normalVector.length() + parameter->EPS)) * normalVector; // cell L -> cell R
				//ʹ��ӭ�緽ʽ�ж���ֵ����ֵ
				// cell L
				iterCellFieldData_ = cellFieldData->begin() + cl - 1;
				// base moment i
				tensor1D<real, NDOFS> momentI;
				// base i
				tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
				point baryCenterI = (*iterCellFieldData_).baryCenter;
				point scaleI = (*iterCellFieldData_).lengthReference;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
					momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
				// CfvMath::getDiffBaseValue(
				// 	pMid,
				// 	baryCenterI,
				// 	scaleI,
				// 	momentI,
				// 	matrixDiffBaseI);

				if (iterFaceFieldData_->diffBaseValueDataMid[cff][0][0] == UNINITReal)
				{
					RBFB1GetDiffBaseValue(
						pparaml,
						baryCenterI,
						scaleI,
						momentI,
						matrixDiffBaseI,
						*iterCellFieldData_);
					CfvMath::VVMatCopy(matrixDiffBaseI, iterFaceFieldData_->diffBaseValueDataMid[cff],
									   0, NDOFS, 0, NDIFFS);
				}
				auto &matrixDiffBaseII = iterFaceFieldData_->diffBaseValueDataMid[cff];

				real dxL = 0.0;
				real dyL = 0.0;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					dxL += matrixDiffBaseII[kk][1] * (*iterCellFieldData_).scalarVariableTn[kk]; // based on VR
					dyL += matrixDiffBaseII[kk][2] * (*iterCellFieldData_).scalarVariableTn[kk];
				}
				// cell R
				iterCellFieldData_ = cellFieldData->begin() + cr - 1;
				// base moment i
				tensor1D<real, NDOFS> momentJ;
				// base i
				tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseJ;
				point baryCenterJ = (*iterCellFieldData_).baryCenter;
				point scaleJ = (*iterCellFieldData_).lengthReference;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
					momentJ[kk] = (*iterCellFieldData_).baseMoment[kk];
				// CfvMath::getDiffBaseValue(
				// 	pMid + (*iterFaceFieldData_).sideOff, //���ڱ߽�����
				// 	baryCenterJ,
				// 	scaleJ,
				// 	momentJ,
				// 	matrixDiffBaseJ);
				if (iterFaceFieldData_->diffBaseValueDataMid[cffr][0][0] == UNINITReal)
				{
					RBFB1GetDiffBaseValue(
						pparamr, //���ڱ߽�����
						baryCenterJ,
						scaleJ,
						momentJ,
						matrixDiffBaseJ,
						*iterCellFieldData_);
					CfvMath::VVMatCopy(matrixDiffBaseJ, iterFaceFieldData_->diffBaseValueDataMid[cffr],
									   0, NDOFS, 0, NDIFFS);
				}
				auto &matrixDiffBaseJJ = iterFaceFieldData_->diffBaseValueDataMid[cffr];
				real dxR = 0.0;
				real dyR = 0.0;
				for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
				{
					dxR += matrixDiffBaseJJ[kk][1] * (*iterCellFieldData_).scalarVariableTn[kk]; // based on VR
					dyR += matrixDiffBaseJJ[kk][2] * (*iterCellFieldData_).scalarVariableTn[kk];
				}
				point windDirection;
				windDirection.x = 0.5 * (dxL + dxR);
				windDirection.y = 0.5 * (dyL + dyR);
				real signLR = CfvMath::getSign(getInnerProduct(windDirection, uNV));
				real *fIL = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
				real *fIR = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
				tensor1D<real, NDOFSCR> *fIBase = new tensor1D<real, NDOFSCR>[static_cast<int>((*iterFaceFieldData_).fPG)];
#ifdef RBFB1_CR_DIFFBOUND
				real *fILDx = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
				real *fIRDx = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
				tensor1D<real, NDOFSCR> *fIBaseDx = new tensor1D<real, NDOFSCR>[static_cast<int>((*iterFaceFieldData_).fPG)];
				real *fILDy = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
				real *fIRDy = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
				tensor1D<real, NDOFSCR> *fIBaseDy = new tensor1D<real, NDOFSCR>[static_cast<int>((*iterFaceFieldData_).fPG)];
#endif
				// get the function values on each face Gauss point.
				for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
				{
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					auto &vectorBaseII = iterFaceFieldData_->diffBaseValueData[cff][gg];
					auto &vectorBaseIICR = iterFaceFieldData_->diffBaseValueDataCR[cff][gg];
					fIL[gg] = (*iterCellFieldData_).scalarVariableTn[0];
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						fIL[gg] += vectorBaseII[kk][0] * (*iterCellFieldData_).scalarVariableTn[kk];
					for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
						fIBase[gg][kk] = vectorBaseIICR[kk][0];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					auto &vectorBaseJJ = iterFaceFieldData_->diffBaseValueData[cffr][gg];
					fIR[gg] = (*iterCellFieldData_).scalarVariableTn[0];
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						fIR[gg] += vectorBaseJJ[kk][0] * (*iterCellFieldData_).scalarVariableTn[kk];

#ifdef RBFB1_CR_DIFFBOUND
					fILDx[gg] = 0., fILDy[gg] = 0.;
					fIRDx[gg] = 0., fIRDy[gg] = 0.;
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						fILDx[gg] += vectorBaseII[kk][1] * (*iterCellFieldData_).scalarVariableTn[kk],
							fILDy[gg] += vectorBaseII[kk][2] * (*iterCellFieldData_).scalarVariableTn[kk];
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						fIRDx[gg] += vectorBaseJJ[kk][1] * (*iterCellFieldData_).scalarVariableTn[kk],
							fIRDy[gg] += vectorBaseJJ[kk][2] * (*iterCellFieldData_).scalarVariableTn[kk];
					for (int kk = 1; kk < static_cast<int>(NDOFSCR); ++kk)
						fIBaseDx[gg][kk] = vectorBaseIICR[kk][1], fIBaseDy[gg][kk] = vectorBaseIICR[kk][2];
#endif
				}
				// get the RHS term
				iterCellFieldData_ = cellFieldData->begin() + cl - 1;
				real valueL = (*iterCellFieldData_).scalarVariableTn[0];
				for (int ll = 1; ll < static_cast<int>(NDOFSCR); ++ll)
				{
					real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					real parametricArea;
					real result;
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
#ifdef RBFB1_INCREMENT_CR
						f[gg] = (0.5 * (1.0 + signLR) * fIL[gg] + 0.5 * (1.0 - signLR) * fIR[gg] - fIL[gg]) *
								fIBase[gg][ll] * std::pow((*iterFaceFieldData_).faceWeightVF[0], 2);
#else
						f[gg] = (0.5 * (1.0 + signLR) * fIL[gg] + 0.5 * (1.0 - signLR) * fIR[gg] - valueL) *
								fIBase[gg][ll] * std::pow((*iterFaceFieldData_).faceWeightVF[0], 2);
#endif
#ifdef RBFB1_CR_DIFFBOUND
						real fbdx = (0.5 * (1.0 + signLR) * fILDx[gg] + 0.5 * (1.0 - signLR) * fIRDx[gg]);
						real fbdy = (0.5 * (1.0 + signLR) * fILDy[gg] + 0.5 * (1.0 - signLR) * fIRDy[gg]);
#ifdef RBFB1_INCREMENT_CR
						fbdx -= fILDx[gg];
						fbdy -= fILDy[gg];
#endif
						f[gg] += ((*iterFaceFieldData_).faceWeightVF[1] * fbdx + (*iterFaceFieldData_).faceWeightVF[2] * fbdy) *
								 ((*iterFaceFieldData_).faceWeightVF[1] * fIBaseDx[gg][ll] + (*iterFaceFieldData_).faceWeightVF[2] * fIBaseDy[gg][ll]) * RBFB1_CR_DIFFDIM;
#endif
						// integrated function value
						weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
						cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
						parametricArea = (*iterFaceFieldData_).parametricArea;
					}
					gaussIntegralFace->getIntegral(
						static_cast<int>((*iterFaceFieldData_).fPG),
						f,
						weight,
						cofJacobi,
						parametricArea,
						result);
					bi[ll] += result; //+ means sum of each face
					delete[] f;
					delete[] weight;
					delete[] cofJacobi;
					f = NULL;
					weight = NULL;
					cofJacobi = NULL;
				}
				delete[] fIL;
				delete[] fIR;
				delete[] fIBase;
#ifdef RBFB1_CR_DIFFBOUND
				delete[] fILDx;
				delete[] fIRDx;
				delete[] fIBaseDx;
				delete[] fILDy;
				delete[] fIRDy;
				delete[] fIBaseDy;
#endif
				fIL = NULL;
				fIR = NULL;
				fIBase = NULL;
			}

			// assign gradient values : relax
			// for (int ii = 1; ii < static_cast<int>(NDOFSCR); ++ii)
			// 	(*iterCellFieldData).scalarVariableTnCR[ii] =
			// 		(1.0 - (*iterCellFieldData).relaxFactorCR) * (*iterCellFieldData).scalarVariableTnCR[ii];
			(*iterCellFieldData).scalarVariableTnCR = (1.0 - (*iterCellFieldData).relaxFactorCR) * (*iterCellFieldData).scalarVariableTnCR;
			// update gradient values
			// for (int ll = 1; ll < static_cast<int>(NDOFSCR); ++ll)
			// 	for (int rr = 1; rr < static_cast<int>(NDOFSCR); ++rr)
			// 		(*iterCellFieldData).scalarVariableTnCR[ll] += (*iterCellFieldData).matrixAiiInverseCR[ll][rr] * bi[rr];
			(*iterCellFieldData).scalarVariableTnCR += (*iterCellFieldData).matrixAiiInverseCR * bi * (*iterCellFieldData).relaxFactorCR;
		}
		//		std::cout << " ..CR reconstruction has been completed." << std::endl;
		return true;
	}

	/*
	template <unsigned int O>
	bool reconstructionRBFB1<O>::excuteReconstructionCR_NotReady(
		parameter *parameter,
		cellFieldDataVector *cellFieldData,
		cellGaussDataVector *cellGaussData,
		faceFieldDataVector *faceFieldData,
		faceGaussDataVector *faceGaussData,
		GaussIntegralCellO2Grid<vO> *gaussIntegralCell,
		GaussIntegralFaceO2Grid<fO> *gaussIntegralFace)
	{

		cellFieldDataVector::iterator iterCellFieldData;
		cellFieldDataVector::iterator iterCellFieldData_;
		faceFieldDataVector::iterator iterFaceFieldData_;

		for (iterCellFieldData = cellFieldData->begin(); iterCellFieldData != cellFieldData->end(); ++iterCellFieldData)
		{
			tensor1D<real, NDOFS> bi;
			if ((*iterCellFieldData).boundaryCellType_ == InnerCell)
			{
				for (int cc = 1; cc < (*iterCellFieldData).cellCellNumber + 1; ++cc)
				{
					int cl = (*iterCellFieldData).cellCellIndex[0];
					int cr = (*iterCellFieldData).cellCellIndex[cc];
					int kf;
					//ȷ����
					for (int ff = 1; ff < (*iterCellFieldData).cellFaceNumber + 1; ++ff)
					{
						kf = (*iterCellFieldData).cellFaceIndex[ff];
						iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
						int clMatch = (*iterFaceFieldData_).faceCellIndex[1];
						int crMatch = (*iterFaceFieldData_).faceCellIndex[2];
						if (cr == clMatch || cr == crMatch)
						{
							break;
						}
					}
					iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
					point p1 = (*iterFaceFieldData_).faceNode[1].second;
					point p2 = (*iterFaceFieldData_).faceNode[2].second;
					point pMid = 0.5 * (p1 + p2);
					point outDirection = pMid - (*iterCellFieldData).baryCenter;
					point normalVector = getInnerProduct(outDirection, (*iterFaceFieldData_).normalVector) * (*iterFaceFieldData_).normalVector;
					point uNV = (1.0 / (normalVector.length() + parameter->EPS)) * normalVector; // cell L -> cell R
					//ʹ��ӭ�緽ʽ�ж���ֵ����ֵ
					// cell L
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					// base moment i
					tensor1D<real, NDOFS> momentI;
					// base i
					tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
					point baryCenterI = (*iterCellFieldData_).baryCenter;
					point scaleI = (*iterCellFieldData_).lengthReference;
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
					{
						momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
					}
					CfvMath::getDiffBaseValueRBFB1(
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
					// cell R
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					// base moment i
					tensor1D<real, NDOFS> momentJ;
					// base i
					tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseJ;
					point baryCenterJ = (*iterCellFieldData_).baryCenter;
					point scaleJ = (*iterCellFieldData_).lengthReference;
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
					{
						momentJ[kk] = (*iterCellFieldData_).baseMoment[kk];
					}
					CfvMath::getDiffBaseValueRBFB1(
						pMid + (*iterFaceFieldData_).sideOff, //���ڱ߽�����
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
					point windDirection;
					windDirection.x = 0.5 * (dxL + dxR);
					windDirection.y = 0.5 * (dyL + dyR);
					real signLR = CfvMath::getSign(getInnerProduct(windDirection, uNV));
					real *fIL = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					real *fIR = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					tensor1D<real, NDOFS> *fIBase = new tensor1D<real, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the function values on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor1D<real, NDOFS> vectorBaseI;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getBaseValueRBFB1(
							p,
							baryCenterI,
							scaleI,
							momentI,
							vectorBaseI);
						fIL[gg] = (*iterCellFieldData_).scalarVariableTn[0];
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							fIL[gg] += vectorBaseI[kk] * (*iterCellFieldData_).scalarVariableTn[kk];
							fIBase[gg][kk] = vectorBaseI[kk];
						}
						// diff base matrix j
						tensor1D<real, NDOFS> vectorBaseJ;
						// base moment j
						tensor1D<real, NDOFS> momentJ;
						iterCellFieldData_ = cellFieldData->begin() + cr - 1;
						point baryCenterJ = (*iterCellFieldData_).baryCenter;
						point scaleJ = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentJ[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getBaseValueRBFB1(
							p + (*iterFaceFieldData_).sideOff, //���ڱ߽�����
							baryCenterJ,
							scaleJ,
							momentJ,
							vectorBaseJ);
						fIR[gg] = (*iterCellFieldData_).scalarVariableTn[0];
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							fIR[gg] += vectorBaseJ[kk] * (*iterCellFieldData_).scalarVariableTn[kk];
						}
					}
					// get the RHS term
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					real valueL = (*iterCellFieldData_).scalarVariableTn[0];
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real parametricArea;
						real result;
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							f[gg] = (0.5 * (1.0 + signLR) * fIL[gg] + 0.5 * (1.0 - signLR) * fIR[gg] - valueL) * fIBase[gg][ll] * (*iterFaceFieldData_).faceWeightVF[0];
							// integrated function value
							weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
							cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
							parametricArea = (*iterFaceFieldData_).parametricArea;
						}
						gaussIntegralFace->getIntegral(
							static_cast<int>((*iterFaceFieldData_).fPG),
							f,
							weight,
							cofJacobi,
							parametricArea,
							result);
						bi[ll] += result; //+ means sum of each face
						delete[] f;
						delete[] weight;
						delete[] cofJacobi;
						f = NULL;
						weight = NULL;
						cofJacobi = NULL;
					}
					delete[] fIL;
					delete[] fIR;
					delete[] fIBase;
					fIL = NULL;
					fIR = NULL;
					fIBase = NULL;
				}
			}
			else
			{
				// deal with boundary
				for (int ff = 1; ff < (*iterCellFieldData).cellFaceNumber + 1; ++ff)
				{ //ע��ѭ�����ޣ���20200314
					int kf = (*iterCellFieldData).cellFaceIndex[ff];
					iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
					int cl = (*iterFaceFieldData_).faceCellIndex[1];
					int cr = (*iterFaceFieldData_).faceCellIndex[2];
					// alternate cell��ֻ���ڲ���Ԫ���ܻ���Ҫ
					if (cr == (*iterCellFieldData).index)
					{
						int ck = cr;
						cr = cl; //�߽粻�õ���
						cl = ck;
						// std::cout << "		alternate left and right cell." << std::endl;
					}
					real *fIL = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					tensor1D<real, NDOFS> *fIBase = new tensor1D<real, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					if (cr == Wall)
					{
						// get the function values on each face Gauss point.
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							// diff base matrix i
							tensor1D<real, NDOFS> vectorBaseI;
							// base moment i
							tensor1D<real, NDOFS> momentI;
							point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
							iterCellFieldData_ = cellFieldData->begin() + cl - 1;
							point baryCenterI = (*iterCellFieldData_).baryCenter;
							point scaleI = (*iterCellFieldData_).lengthReference;
							for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							{
								momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
							}
							CfvMath::getBaseValueRBFB1(
								p,
								baryCenterI,
								scaleI,
								momentI,
								vectorBaseI);

							fIL[gg] = 0.0;
							for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							{
								// fIL[gg] += vectorBaseI[kk] * (*iterCellFieldData_).scalarVariableTn[kk];
								fIBase[gg][kk] = vectorBaseI[kk];
							}
						}
					}
					else if (cr == FarField)
					{
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							// diff base matrix i
							tensor1D<real, NDOFS> vectorBaseI;
							// base moment i
							tensor1D<real, NDOFS> momentI;
							point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
							iterCellFieldData_ = cellFieldData->begin() + cl - 1;
							point baryCenterI = (*iterCellFieldData_).baryCenter;
							point scaleI = (*iterCellFieldData_).lengthReference;
							for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							{
								momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
							}
							CfvMath::getBaseValueRBFB1(
								p,
								baryCenterI,
								scaleI,
								momentI,
								vectorBaseI);

							fIL[gg] = (*iterCellFieldData_).scalarVariableTn[0];
							for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							{
								// fIL[gg] += vectorBaseI[kk] * (*iterCellFieldData_).scalarVariableTn[kk];
								fIBase[gg][kk] = vectorBaseI[kk];
							}
						}
					}
					else if (cr == Symmetric)
					{
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							// diff base matrix i
							tensor1D<real, NDOFS> vectorBaseI;
							// base moment i
							tensor1D<real, NDOFS> momentI;
							point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
							iterCellFieldData_ = cellFieldData->begin() + cl - 1;
							point baryCenterI = (*iterCellFieldData_).baryCenter;
							point scaleI = (*iterCellFieldData_).lengthReference;
							for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							{
								momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
							}
							CfvMath::getBaseValueRBFB1(
								p,
								baryCenterI,
								scaleI,
								momentI,
								vectorBaseI);
							fIL[gg] = (*iterCellFieldData_).scalarVariableTn[0];
							for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
							{
								fIL[gg] += vectorBaseI[kk] * (*iterCellFieldData_).scalarVariableTn[kk];
								fIBase[gg][kk] = vectorBaseI[kk];
							}
						}
					}
					else if (cr == OutFlow)
					{
						std::cout << "	warning: OutFlow is not proposed currently." << std::endl;
					}
					else if (cr == InFlow)
					{
						std::cout << "	warning: InFlow is not proposed currently." << std::endl;
					}
					// get the RHS term
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					real valueL = (*iterCellFieldData_).scalarVariableTn[0];
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real parametricArea;
						real result;
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							f[gg] = (fIL[gg] - valueL) * fIBase[gg][ll] * (*iterFaceFieldData_).faceWeightVF[0];
							// integrated function value
							weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
							cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
							parametricArea = (*iterFaceFieldData_).parametricArea;
						}
						gaussIntegralFace->getIntegral(
							static_cast<int>((*iterFaceFieldData_).fPG),
							f,
							weight,
							cofJacobi,
							parametricArea,
							result);
						bi[ll] += result; //+ means sum of each face
						delete[] f;
						delete[] weight;
						delete[] cofJacobi;
						f = NULL;
						weight = NULL;
						cofJacobi = NULL;
					}
					delete[] fIL;
					delete[] fIBase;
					fIBase = NULL;
					fIL = NULL;
				}
				// deal with inner
				for (int cc = 1; cc < (*iterCellFieldData).cellCellNumber + 1; ++cc)
				{
					int cl = (*iterCellFieldData).cellCellIndex[0];
					int cr = (*iterCellFieldData).cellCellIndex[cc];
					int kf;
					//ȷ����
					for (int ff = 1; ff < (*iterCellFieldData).cellFaceNumber + 1; ++ff)
					{
						kf = (*iterCellFieldData).cellFaceIndex[ff];
						iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
						int clMatch = (*iterFaceFieldData_).faceCellIndex[1];
						int crMatch = (*iterFaceFieldData_).faceCellIndex[2];
						if (cr == clMatch || cr == crMatch)
						{
							break;
						}
					}
					iterFaceFieldData_ = faceFieldData->begin() + kf - 1;
					point p1 = (*iterFaceFieldData_).faceNode[1].second;
					point p2 = (*iterFaceFieldData_).faceNode[2].second;
					point pMid = 0.5 * (p1 + p2);
					point outDirection = pMid - (*iterCellFieldData).baryCenter;
					point normalVector = getInnerProduct(outDirection, (*iterFaceFieldData_).normalVector) * (*iterFaceFieldData_).normalVector;
					point uNV = (1.0 / (normalVector.length() + parameter->EPS)) * normalVector; // cell L -> cell R
					//ʹ��ӭ�緽ʽ�ж���ֵ����ֵ
					// cell L
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					// base moment i
					tensor1D<real, NDOFS> momentI;
					// base i
					tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseI;
					point baryCenterI = (*iterCellFieldData_).baryCenter;
					point scaleI = (*iterCellFieldData_).lengthReference;
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
					{
						momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
					}
					CfvMath::getDiffBaseValueRBFB1(
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
					// cell R
					iterCellFieldData_ = cellFieldData->begin() + cr - 1;
					// base moment i
					tensor1D<real, NDOFS> momentJ;
					// base i
					tensor2D<real, NDOFS, NDIFFS> matrixDiffBaseJ;
					point baryCenterJ = (*iterCellFieldData_).baryCenter;
					point scaleJ = (*iterCellFieldData_).lengthReference;
					for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
					{
						momentJ[kk] = (*iterCellFieldData_).baseMoment[kk];
					}
					CfvMath::getDiffBaseValueRBFB1(
						pMid + (*iterFaceFieldData_).sideOff, //���ڱ߽�����
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
					point windDirection;
					windDirection.x = 0.5 * (dxL + dxR);
					windDirection.y = 0.5 * (dyL + dyR);
					real signLR = CfvMath::getSign(getInnerProduct(windDirection, uNV));
					real *fIL = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					real *fIR = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
					tensor1D<real, NDOFS> *fIBase = new tensor1D<real, NDOFS>[static_cast<int>((*iterFaceFieldData_).fPG)];
					// get the function values on each face Gauss point.
					for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
					{
						// diff base matrix i
						tensor1D<real, NDOFS> vectorBaseI;
						// base moment i
						tensor1D<real, NDOFS> momentI;
						point p = (*iterFaceFieldData_).gaussPairVector_[gg].p;
						iterCellFieldData_ = cellFieldData->begin() + cl - 1;
						point baryCenterI = (*iterCellFieldData_).baryCenter;
						point scaleI = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentI[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getBaseValueRBFB1(
							p,
							baryCenterI,
							scaleI,
							momentI,
							vectorBaseI);

						fIL[gg] = (*iterCellFieldData_).scalarVariableTn[0];
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							fIL[gg] += vectorBaseI[kk] * (*iterCellFieldData_).scalarVariableTn[kk];
							fIBase[gg][kk] = vectorBaseI[kk];
						}
						// diff base matrix j
						tensor1D<real, NDOFS> vectorBaseJ;
						// base moment j
						tensor1D<real, NDOFS> momentJ;
						iterCellFieldData_ = cellFieldData->begin() + cr - 1;
						point baryCenterJ = (*iterCellFieldData_).baryCenter;
						point scaleJ = (*iterCellFieldData_).lengthReference;
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							momentJ[kk] = (*iterCellFieldData_).baseMoment[kk];
						}
						CfvMath::getBaseValueRBFB1(
							p + (*iterFaceFieldData_).sideOff, //���ڱ߽�����
							baryCenterJ,
							scaleJ,
							momentJ,
							vectorBaseJ);
						fIR[gg] = (*iterCellFieldData_).scalarVariableTn[0];
						for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk)
						{
							fIR[gg] += vectorBaseJ[kk] * (*iterCellFieldData_).scalarVariableTn[kk];
						}
					}
					// get the RHS term
					iterCellFieldData_ = cellFieldData->begin() + cl - 1;
					real valueL = (*iterCellFieldData_).scalarVariableTn[0];
					for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
					{
						real *f = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *weight = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real *cofJacobi = new real[static_cast<int>((*iterFaceFieldData_).fPG)];
						real parametricArea;
						real result;
						for (int gg = 0; gg < static_cast<int>((*iterFaceFieldData_).fPG); ++gg)
						{
							f[gg] = (0.5 * (1.0 + signLR) * fIL[gg] + 0.5 * (1.0 - signLR) * fIR[gg] - valueL) * fIBase[gg][ll] * (*iterFaceFieldData_).faceWeightVF[0];
							// integrated function value
							weight[gg] = (*iterFaceFieldData_).parametricValue[gg].second;
							cofJacobi[gg] = (*iterFaceFieldData_).gaussPairVector_[gg].JacobiCof;
							parametricArea = (*iterFaceFieldData_).parametricArea;
						}
						gaussIntegralFace->getIntegral(
							static_cast<int>((*iterFaceFieldData_).fPG),
							f,
							weight,
							cofJacobi,
							parametricArea,
							result);
						bi[ll] += result; //+ means sum of each face
						delete[] f;
						delete[] weight;
						delete[] cofJacobi;
						f = NULL;
						weight = NULL;
						cofJacobi = NULL;
					}
					delete[] fIL;
					delete[] fIR;
					delete[] fIBase;
					fIL = NULL;
					fIR = NULL;
					fIBase = NULL;
				}
			}
			// assign gradient values : relax
			for (int ii = 1; ii < static_cast<int>(NDOFS); ++ii)
			{
				(*iterCellFieldData).scalarVariableTnCR[ii] =
					(1.0 - (*iterCellFieldData).relaxFactorCR) * (*iterCellFieldData).scalarVariableTnCR[ii];
			}
			// update gradient values
			for (int ll = 1; ll < static_cast<int>(NDOFS); ++ll)
			{
				for (int rr = 1; rr < static_cast<int>(NDOFS); ++rr)
				{
					(*iterCellFieldData).scalarVariableTnCR[ll] += (*iterCellFieldData).matrixAiiInverseCR[ll][rr] * bi[rr];
				}
			}
		}
		//		std::cout << " ..CR reconstruction has been completed." << std::endl;
		return true;
	}
	*/

}
#endif