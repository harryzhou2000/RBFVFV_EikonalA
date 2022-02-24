#include "problemCircle.h"

namespace ScalarCfv
{
	bool problemCIRCLE::timeMarchingExplicitSSPRK(
		vertexVector* node,
		GaussIntegralCellO1Grid<vO>* gaussIntegralCell,
		GaussIntegralFaceO1Grid<fO>* gaussIntegralFace,
		const std::string & filenameInputBackup,
		const std::string & filenameOutputBackup,
		const std::string & filenameOutputSln,
		const std::string & filenameOutputResidual,
		const std::string & filenameOutputError
		){

		real residual;
		int IOflag = 0;

		real L1, L2, L8;

		allocateArray();
		if (parameter_->isContinue == true){
			readBackUp(filenameInputBackup);
		}
		else{
			initializeArray();
		}

		reconstruction_->initBaseMomentAndRelaxFactor(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			gaussIntegralCell
			);
		reconstruction_->initFaceWeight(
			parameter_,
			cellFieldData_,
			faceFieldData_
			);
		reconstruction_->initReconstructionMatrixAndVector(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace
			);
		reconstruction_->initReconstructionMatrixAndVectorCR(
			parameter_,
			cellFieldData_,
			cellGaussData_,
			faceFieldData_,
			faceGaussData_,
			gaussIntegralCell,
			gaussIntegralFace
			);
		exportCurrentSln(node, 0, filenameOutputSln);
		//main loop
		parameter_->refTn = 0.0;//global time marching
		for (int iStep = parameter_->nStart; iStep <= parameter_->nEnd; ++iStep){
			evolution_->getTimeStep(
				parameter_,
				cellFieldData_,
				faceFieldData_);
			if (parameter_->isLocalTimeMarching == false){
				parameter_->refTn += parameter_->refDeltaTn;
			}

			for (int irk = 1; irk <= parameter_->nStepTimeMarching; ++irk){
				for (int isubrk = 1; isubrk <= parameter_->nInnerStepTimeMarching; ++isubrk){
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
					//deal with AV?
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


					if (parameter_->nStepTimeMarching == 3){
						evolution_->excuteInnerLoopExplicitSSPRKStep3(
							irk,
							cellFieldData_);
					}
					else if (parameter_->nStepTimeMarching == 4){
						evolution_->excuteInnerLoopExplicitSSPRKStep4(
							irk,
							cellFieldData_);
					}
					//evolution_->applyBoundaryLimiter(
					//	parameter_,
					//	cellFieldData_);
				}

				if (iStep%parameter_->nScreenOutput == 0){
					std::cout << "	RK step:" << "\t" << "\t" << irk << "has been done." << std::endl;

				}
			}
			//post -> residual
			getResidual(residual);
			exportResidual(
				IOflag,
				iStep,
				residual,
				filenameOutputResidual);

			//post -> GRE
			//getGRE(L1, L2, L8);

			getGRE(gaussIntegralCell, L1, L2, L8);
			exportGRE(
				IOflag,
				iStep,
				L1, L2, L8,
				filenameOutputError);

			//evolution_->excuteFilter(parameter_,
			//	cellFieldData_);

			if (iStep%parameter_->nScreenOutput == 0){
				std::cout << "--------------------------------------" << std::endl;
				std::cout << "Current step:" << "\t" << "\t" << iStep << std::endl;

			}
			if (iStep%parameter_->nFileOutput == 0){
				exportCurrentSln(node, iStep, filenameOutputSln);
			}

			if (iStep%parameter_->nBackUp == 0){
				//exportBackUp(filenameOutputBackup);
				exportBackUp(iStep, filenameOutputBackup);
			}
		}
		return true;
	}

//	bool problemCIRCLE::getGRE(real & L1, real & L2, real & L8){
//		real L1t = 0.0;
//		real L2t = 0.0;
//		real L8t = 0.0;
//		real volt = 0.0;
//		cellFieldDataVector::iterator iterCell;
//		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell){
//			tensor1D<real, NDOFS> moment;
//			//base j
//			tensor2D<real, NDOFS, NDOFS> matrixDiffBase;
//			point p = (*iterCell).baryCenter;
//			point baryCenter = (*iterCell).baryCenter;
//			point scale = (*iterCell).lengthReference;
//			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk){
//				moment[kk] = (*iterCell).baseMoment[kk];
//			}
//			CfvMath::getDiffBaseValue(
//				p,
//				baryCenter,
//				scale,
//				moment,
//				matrixDiffBase);
//			real dist = (*iterCell).scalarVariableTn[0];
//			for (int kk = 1; kk < static_cast<int>(NDOFS); ++kk){
//				dist += matrixDiffBase[kk][0] * (*iterCell).scalarVariableTn[kk];
//			}
//			
//			real diff = std::fabs(dist - ((*iterCell).baryCenter.length() - 0.5));
//
//			L1t += diff * (*iterCell).volume;
//			L2t += diff*diff* (*iterCell).volume;
//			L8t = getMax(L8t, diff);
//			volt += (*iterCell).volume;
//		}
//		L1 = L1t / volt;
//		L2 = std::sqrt(L2t / volt);
//		L8 = L8t;
//		return true;
//	}

	bool problemCIRCLE::getGRE(GaussIntegralCellO1Grid<vO>* gaussIntegralCell, real & L1, real & L2, real & L8){
		real L1t = 0.0;
		real L2t = 0.0;
		real L8t = 0.0;
		real volt = 0.0;
		cellFieldDataVector::iterator iterCell;
		for (iterCell = cellFieldData_->begin(); iterCell != cellFieldData_->end(); ++iterCell){

			real* fI
				= new real[static_cast<int>((*iterCell).PG)];
			real* weight = new real[static_cast<int>((*iterCell).PG)];
			real* cofJacobi = new real[static_cast<int>((*iterCell).PG)];
			real parametricVolume = (*iterCell).parametricVolume;
			for (int gg = 0; gg < static_cast<int>((*iterCell).PG); ++gg){
				point p = (*iterCell).gaussPairVector_[gg].p;
				fI[gg] = p.length() - 0.5;
				weight[gg] = (*iterCell).parametricValue[gg].second;
				cofJacobi[gg] = (*iterCell).gaussPairVector_[gg].JacobiCof;
			}
			
			real result;
			real* f = new real[static_cast<int>((*iterCell).PG)];
			for (int gg = 0; gg < static_cast<int>((*iterCell).PG); ++gg){
				f[gg] = fI[gg];
			}
			gaussIntegralCell->getIntegral(
				static_cast<int>((*iterCell).PG),
				f,
				weight,
				cofJacobi,
				parametricVolume,
				result);
			result /= (*iterCell).volume;

			real diff = std::fabs(result - (*iterCell).scalarVariableTn[0]);

			L1t += diff * (*iterCell).volume;
			L2t += diff*diff* (*iterCell).volume;
			L8t = getMax(L8t, diff);
			volt += (*iterCell).volume;
		}
		L1 = L1t / volt;
		L2 = std::sqrt(L2t / volt);
		L8 = L8t;
		return true;
	}

	bool problemCIRCLE::exportGRE(
		int & flag,
		const int step,
		real L1, real L2, real L8,
		const std::string & filename){
	
		std::ofstream fileOut;
		fileOut.open(filename.c_str(), std::ios_base::app);
		if (!fileOut){
			std::cerr << "		error: unable to open \"" << filename << "\"\n";
			exit(1);
		}
		if (flag == 0){
			fileOut << "VARIABLES = \"Steps\", \"L1(Err)\", \"L2(Err)\", \"L8(Err)\"\n";
		}
		fileOut << step << "\t" << L1 << "\t" << L2 << "\t" << L8 << std::endl;
		flag += 1;
		return true;
	
	}


}