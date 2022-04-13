#ifndef _PARAMETER_H
#define _PARAMETER_H
#include "TypeDefine.h"
#include "Point.h"
#include <fstream>
#include <iomanip> //for setprecision
#include <string>
#include <iostream>
#include <sstream> //for istringstream

namespace ScalarCfv
{
	//extern std::string fileIn_Mesh("D:/calculations/new_test/grid.in");
	//extern std::string fileOut_Mesh("D:/calculations/new_test/grid.plt");
	//extern std::string fileIn_Parameter("D:/calculations/new_test/inputParameters.txt");
	extern std::string fineOut_BndCheck;
	extern std::string fileIn_Mesh;		 // ("D:/calculations/new_test/grid.in");
	extern std::string fileOut_Mesh;	 // ("D:/calculations/new_test/grid.plt");
	extern std::string fileIn_Parameter; // ("D:/calculations/new_test/inputParameters.txt");
	//extern std::string fileIn_OldSln;// ("D:/calculations/new_test/inputParameters.txt");
	extern std::string fileOut_Sln; // ("D:/calculations/new_test/inputParameters.txt");
	extern std::string fileOut_Residual;
	extern std::string fileOut_BackUp;
	extern std::string fileIn_BackUp;
	//debug
	extern std::string fileOut_debug;	   // ("D:/calculations/new_test/inputParameters.txt");
	extern std::string fileOut_debug_sln;  // ("D:/calculations/new_test/inputParameters.txt");
	extern std::string fileOut_debug_dsln; // ("D:/calculations/new_test/inputParameters.txt");
	
	class parameter
	{
	public:
		parameter() : ndim(2),
					  nodeNumber(0),
					  nodeNumberGhost(0),
					  faceNumber(0),
					  faceNumberGhost(0),
					  cellNumber(0),
					  cellNumberGhost(0),
					  innerCellNumber(0),
					  boundaryCellNumber(0),
					  boundaryMixedCellNumber(0),
					  boundaryWallCellNumber(0),
					  boundaryFarfieldCellNumber(0),
					  boundarySymmetricCellNumber(0),
					  boundaryOutflowCellNumber(0),
					  boundaryInflowCellNumber(0),
					  boundaryNumber(0),
					  boundaryPeriodicCellNumber(0),
					  boundaryFaceNumber(0),
					  reconstructionOrder(0),
					  nStart(0),
					  nEnd(0),
					  nScreenOutput(0),
					  nFileOutput(0),
					  nBackUp(0),
					  nStepTimeMarching(0),
					  nInnerStepTimeMarching(0),
					  isImplicitTimeMarching(false),
					  isContinue(false),
					  isPeriodicBoundary(false),
					  isViscous(false),
					  EPS(1e-12),
					  unit(1.0),
					  refTn(0.0),
					  refDeltaTn(0.0),
					  residualTn(0.0),
					  cofAV(0.0){};
		~parameter(){};

	public:
		//topology parameter
		int inputMeshOrder;

		int nodeNumber;
		int nodeNumberGhost;
		int faceNumber;
		int faceNumberGhost;
		int cellNumber;
		int cellNumberGhost;

		int innerCellNumber;
		int boundaryCellNumber;

		int boundaryWallCellNumber;
		int boundaryFarfieldCellNumber;
		int boundarySymmetricCellNumber;
		int boundaryOutflowCellNumber;
		int boundaryInflowCellNumber;
		int boundaryPeriodicCellNumber;
		int boundaryMixedCellNumber;

		int boundaryNumber;
		//int boundaryNumberPeriodic;
		int boundaryFaceNumber;
		point offSidePeriodic;
		//geomertic parameter
		real unit;

		//gas parameter
		real Re;
		real PrLam;
		real PrTur;
		real muReference;
		real tmpWallReference;

		//flow field parameter
		real rhoReference;
		real uxReference;
		real uyReference;
		real preReference;
		real tmpReference;
		real maReference;

		//time parameter
		real CFL;
		real totalTimePhysical;
		real dTReferencePhysical;

		//control parameter
		int ndim;
		int reconstructionOrder;
		int nStart;		   //����
		int nEnd;		   //����
		int nScreenOutput; //����
		int nFileOutput;   //����
		int nBackUp;	   //����
		//int nRKLoop;
		int nStepTimeMarching;		//���� S3P3 & S4P4 &S1P1
		int nInnerStepTimeMarching; //����
		bool isContinue;
		bool isPeriodicBoundary;
		bool isViscous;
		real EPS;
		bool isImplicitTimeMarching; //����
		bool isLocalTimeMarching;	 //����
		//int timeMarchingStep;
		real cofAV;

		real refTn;
		real refDeltaTn;
		real residualTn;
		//integral control parameter
		//unsigned int pgFace;
		//unsigned int pnFace;
		//unsigned int pgVolTri;
		//unsigned int pnVolTri;
		//unsigned int pgVolQuad;
		//unsigned int pnVolQuad;

	public:
		bool readParameterFile(const std::string &filename);
	};
}
#endif