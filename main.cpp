#include "main.h"

//���Խڵ�1����ɼ������˼�Gauss���ֺ�Ԥ��02.20����

int main(int argc, char **argv)
{

	ScalarCfv::vertexVector node;
	ScalarCfv::faceVector face;
	ScalarCfv::cellVector cell;
	ScalarCfv::boundaryVector boundary;
	ScalarCfv::parameter parameter;

	ScalarCfv::faceGaussDataVector faceGaussData;
	ScalarCfv::cellGaussDataVector cellGaussData;

#ifdef O2mesh
	ScalarCfv::GaussIntegralFaceO2Grid<fO> gaussIntegralFace;
	ScalarCfv::GaussIntegralCellO2Grid<vO> gaussIntegralCell;
#endif
#ifdef O1mesh
	ScalarCfv::GaussIntegralFaceO1Grid<fO> gaussIntegralFace;
	ScalarCfv::GaussIntegralCellO1Grid<vO> gaussIntegralCell;
#endif

	ScalarCfv::cellFieldDataVector cellFieldData;
	ScalarCfv::faceFieldDataVector faceFieldData;
#ifdef USE_RBFA1
	ScalarCfv::reconstructionRBFA1<rO> reconstruction;
#else
	ScalarCfv::reconstruction<rO> reconstruction;
#endif
	ScalarCfv::evolution<rO> evolution;

	ScalarCfv::fieldsolver fieldSolver(
		&parameter,
		&cellFieldData,
		&faceFieldData,
		&cellGaussData,
		&faceGaussData,
		&evolution,
		&reconstruction);

	//ScalarCfv::problemCIRCLE fieldSolver(
	//	&parameter,
	//	&cellFieldData,
	//	&faceFieldData,
	//	&cellGaussData,
	//	&faceGaussData,
	//	&evolution,
	//	&reconstruction);

	//ScalarCfv::problemAT fieldSolver(
	//	&parameter,
	//	&cellFieldData,
	//	&faceFieldData,
	//	&cellGaussData,
	//	&faceGaussData,
	//	&evolution,
	//	&reconstruction);

	ScalarCfv::grid grid(&node,
						 &face,
						 &cell,
						 &boundary,
						 &parameter);

	parameter.readParameterFile(ScalarCfv::fileIn_Parameter);
	grid.readGridFile(ScalarCfv::fileIn_Mesh);

	grid.allocateCell();
	grid.allocateFace();
	grid.setNodeTopology();
	grid.setFaceTopology();
	grid.setBoundaryTopology();
	grid.setFaceProperty();
	grid.setCellTopology();
	grid.saveGrid(ScalarCfv::fileOut_Mesh);
	grid.adjustNV();

	gaussIntegralCell.allocateArray(&parameter, &cell, &cellGaussData);
	gaussIntegralCell.initializeArray(&cell, &cellGaussData);
	grid.getGeometryInfo(&cellGaussData, &gaussIntegralFace, &gaussIntegralCell);
	gaussIntegralCell.syncData(&cell, &cellGaussData);
	grid.exportTestGeoInfo(ScalarCfv::fileOut_debug);

	gaussIntegralFace.allocateArray(&parameter, &face, &faceGaussData);
	gaussIntegralFace.initializeArray(&face, &faceGaussData);
	//	grid.debugTest(&faceGaussData);

	//#define TEST

#ifdef TEST

	fieldSolver.allocateArray();
	fieldSolver.initializeArray();

	//face info

	//ScalarCfv::faceFieldDataVector::iterator iterFaceFieldData;
	//for (iterFaceFieldData = faceFieldData.begin(); iterFaceFieldData != faceFieldData.end(); ++iterFaceFieldData){
	//	std::cout << "------------------------" << std::endl;
	//	std::cout << "node: " << (*iterFaceFieldData).faceNode[1].first << "\t" << (*iterFaceFieldData).faceNode[2].first << std::endl;
	//	std::cout.precision(12);
	//	for (int ii = 0; ii < (*iterFaceFieldData).fPG; ++ii){
	//		std::cout << "\t" << (*iterFaceFieldData).gaussPairVector_[ii].p.x << "\t" << (*iterFaceFieldData).gaussPairVector_[ii].p.y << std::endl;
	//	}
	//	std::cout << "--------------" << std::endl;
	//	for (int ii = 0; ii < (*iterFaceFieldData).fPG; ++ii){
	//		std::cout << "\t" << (*iterFaceFieldData).gaussPairVector_[ii].JacobiCof << std::endl;
	//	}
	//	std::cout << "--------------" << std::endl;
	//	for (int ii = 0; ii < (*iterFaceFieldData).fPG; ++ii){
	//		std::cout << "\t" << (*iterFaceFieldData).gaussPairVector_[ii].normalVector.x << "\t" << (*iterFaceFieldData).gaussPairVector_[ii].normalVector.y << std::endl;
	//	}
	//}

	//cell info

	ScalarCfv::cellFieldDataVector::iterator iterCellFieldData;
	for (iterCellFieldData = cellFieldData.begin(); iterCellFieldData != cellFieldData.end(); ++iterCellFieldData)
	{
		std::cout << "------------------------" << std::endl;
		std::cout << "cell: " << (*iterCellFieldData).index << std::endl;
		std::cout << "PG: " << (*iterCellFieldData).PG << std::endl;
		std::cout.precision(12);
		for (int ii = 0; ii < (*iterCellFieldData).PG; ++ii)
		{
			std::cout << "\t" << (*iterCellFieldData).gaussPairVector_[ii].p.x << "\t" << (*iterCellFieldData).gaussPairVector_[ii].p.y << std::endl;
		}
		std::cout << "--------------" << std::endl;
		for (int ii = 0; ii < (*iterCellFieldData).PG; ++ii)
		{
			std::cout << "\t" << (*iterCellFieldData).gaussPairVector_[ii].JacobiCof << std::endl;
		}
	}

#else

	fieldSolver.timeMarchingExplicitSSPRK(
		&node,
		&gaussIntegralCell,
		&gaussIntegralFace,
		ScalarCfv::fileIn_BackUp,
		ScalarCfv::fileOut_BackUp,
		ScalarCfv::fileOut_Sln,
		ScalarCfv::fileOut_Residual);

	//fieldSolver.timeMarchingExplicitSSPRK(
	//	&node,
	//	&gaussIntegralCell,
	//	&gaussIntegralFace,
	//	ScalarCfv::fileIn_BackUp,
	//	ScalarCfv::fileOut_BackUp,
	//	ScalarCfv::fileOut_Sln,
	//	ScalarCfv::fileOut_Residual,
	//	ScalarCfv::fileOut_Error
	//	);

#endif

	//fieldSolver.timeMarchingExplicit_DEBUG(
	//	&node,
	//  &gaussIntegralCell,
	//	&gaussIntegralFace,
	//	ScalarCfv::fileIn_BackUp,
	//	ScalarCfv::fileOut_BackUp,
	//	ScalarCfv::fileOut_Sln,
	//	ScalarCfv::fileOut_Residual
	//	);

	std::cout << "end" << std::endl;
	return 0;
}