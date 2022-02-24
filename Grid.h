#ifndef _GRID_H
#define _GRID_H
//#include "TypeDefine.h"
//#include "Point.h"
//#include "Geometry.h"
#include "GaussIntegral.h"
#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>

namespace ScalarCfv
{
	class grid
	{
	public:
		grid(vertexVector* node,
			faceVector* face,
			cellVector* cell,
			boundaryVector* boundary,
			parameter* parameter
			) {
			this->node_ = node;
			this->face_ = face;
			this->cell_ = cell;
			this->boundary_ = boundary;
			this->parameter_ = parameter;
		};
		~grid() {};

	public:
		vertexVector* node_;
		faceVector* face_;
		cellVector* cell_;
		boundaryVector* boundary_;
		parameter* parameter_;
		std::vector<std::pair<int, int>> boundaryPair;		//boundaryNum->boundaryType vs. boundaryFaceNum
		
		int boundaryNumberPeriodic;
		//boundaryNumPeriodic->boundaryTag vs. boundaryOffside
		std::vector<std::pair<int, int>> boundaryPairPeriodicTag;
		std::vector<point> boundaryPeriodicOffside;
	public:
		virtual bool readGridFile(const std::string & filename);
		//virtual bool allocateNode();
		virtual bool allocateFace();
		virtual bool allocateCell();
		//virtual bool allocateBoundary();

		virtual bool setNodeTopology();
		virtual bool setFaceTopology();
		virtual int checkFace(const int n1, const int n2);
		virtual bool setCellTopology();
		//virtual bool setNodeFaceTopology();//需要等setFaceTopology执行完后才能执行
		virtual bool setBoundaryTopology();
		virtual bool setFaceProperty();
		virtual bool saveGrid(const std::string & filename);

		//virtual bool getGeometryInfo();

		virtual bool getGeometryInfo(
			cellGaussDataVector* cellGaussData,
			GaussIntegralFaceO1Grid<fO>* gaussFaceIntegral_,
			GaussIntegralCellO1Grid<vO>* gaussCellIntegral_);

		virtual bool getGeometryInfo(
			cellGaussDataVector* cellGaussData,
			GaussIntegralFaceO2Grid<fO>* gaussFaceIntegral_,
			GaussIntegralCellO2Grid<vO>* gaussCellIntegral_);

		virtual bool freeGrid();

		//debug
		virtual bool exportTest(const std::string & filename);
		virtual bool exportTestGeoInfo(const std::string & filename);
		//debug 20200304
		virtual bool debugTest(
			faceGaussDataVector* faceGaussData,
			GaussIntegralFaceO1Grid<fO>* gaussFaceIntegral_);

		virtual bool debugTest(
			faceGaussDataVector* faceGaussData,
			GaussIntegralFaceO2Grid<fO>* gaussFaceIntegral_);

		virtual bool adjustNV();
	};

}
#endif