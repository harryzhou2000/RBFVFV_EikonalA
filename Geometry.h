#ifndef _GEOMETRY_H
#define _GEOMETRY_H
#include "TypeDefine.h"
#include "Point.h"
#include "Parameter.h"
#include <cmath>

namespace ScalarCfv
{
	class vertex
	{
	public:
		vertex() :
			index(0),
			//nodeEdgeNumber(0),
			nodeFaceNumber(0),
			nodeCellNumber(0),
			nodeCellfaceNumber(0),
			isO1Node(true) {};
		~vertex() {};
	
		//vertex parameters
		int index;
		point nodePhysical;
	
		//vertex topology
		//int nodeEdgeNumber;
		int nodeFaceNumber;
		int nodeCellNumber;
		int nodeCellfaceNumber;
		//std::vector<int> nodeEdgeIndex;
		std::vector<int> nodeFaceIndex;
		std::vector<int> nodeCellIndex;
		std::vector<int> nodeCellfaceIndex;
	
		//high order mesh parameters
		bool isO1Node;
	};
	typedef std::vector<vertex> vertexVector;

	//class edge
	//{};
	
	//boundary!
	class face
	{
	public:
		face() :
			index(0),
			faceProperty(0),
			faceNodeNumber(0),
			faceCellNumber(0),
			spectralRadius(0.0),
			area(0.0),
			areaReference(0.0),
			o2PointNumber(0) {};
		~face() {};
	
		//face parameters
		int index;
		int faceProperty;
		int faceNodeNumber;
		int faceCellNumber;
		
		real spectralRadius;
		real area;
		real areaReference;
		point normalVector;
		point sideOff;
	
		//face topology
		std::vector<std::pair<int, point> > faceNode;  //face2NodeIndex, face2NodePoint
		//std::vector<std::pair<int, int>> face2Cell;

		//std::vector<int> face2NodeIndex;
		std::vector<int> faceCellIndex;

		//o2mesh
		int o2PointNumber;
		std::vector<std::pair<int, point> > faceNodeO2Point; //.first is not used!
	};
	typedef std::vector<face> faceVector;
	
	class cell
	{
	public:
		//cell(parameter* parameter) :
		//	index(0),
		//	volume(0.0){
		//	this->parameter_ = parameter;
		//};
		cell() :
			index(0),
			volume(0.0),
			cellNodeNumber(0),
			cellFaceNumber(0),
			cellCellNumber(0),
			isBoundaryCell(false),
			isCountByboundaryCellType(false),
			boundaryCellType_(InnerCell),
			cellInnerFaceNumber(0),
			cellBoundaryFaceNumber(0),
			o2PointNumber(0) {};
		~cell() {};
	
		//cell parameters
		//parameter* parameter_;

		int index;
		bool isBoundaryCell;
		bool isCountByboundaryCellType;
		cellType cellType_;
		boundaryCellType boundaryCellType_;
		real volume;
		point baryCenter;
		point cellSideOff;
		point lengthReference;
	
		//cell topology
		int cellNodeNumber;
		int cellFaceNumber;
		int cellCellNumber;
		//std::vector<int> cellNodeIndex;

		std::vector<std::pair<int, point> > cellNode; //cellNodeIndex, cellNodePoint
		std::vector<int> cellFaceIndex;
		std::vector<int> cellCellIndex;

		std::vector<point> cellFaceSideOff;
		std::vector<real> cellFaceWeight;

		int cellInnerFaceNumber;
		int cellBoundaryFaceNumber;

		//o2mesh
		int o2PointNumber;
		std::vector<std::pair<int, point>> cellNodeO2Point; //.first is not used!

	};
	typedef std::vector<cell> cellVector;

	class boundary : public face{
	
	public:
		boundary():face(){};
		~boundary() {};
	public:
		boundaryType boundaryType_;
		int faceIndex;
		int boundaryTag;
		int matchIndex;
		//std::vector<std::pair<int, point>> faceNode;
	};
	typedef std::vector<boundary> boundaryVector;

}
#endif