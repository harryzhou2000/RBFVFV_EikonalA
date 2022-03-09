#include "Grid.h"
#include "EigenTensor.hpp"
// 20200220 ���������ļ����ڵ���-1����0�ſ�ʼ
namespace ScalarCfv
{

	bool grid::readGridFile(const std::string &filename)
	{
		std::ifstream fileIn;
		fileIn.open(filename.c_str());
		if (!fileIn)
		{
			std::cerr << "	error: fail to open grid file: unable to open \"" << filename << "\"\n";
			exit(1);
		}
		std::cout << "reading grid file..." << std::endl;

		fileIn >> parameter_->inputMeshOrder;
		if (parameter_->inputMeshOrder != mO)
		{
			std::cout << "	error: input mesh order incorrect." << std::endl;
			exit(1);
		}

		fileIn >> parameter_->cellNumber;
		fileIn >> parameter_->cellNumberGhost;
		fileIn >> parameter_->nodeNumber;
		fileIn >> parameter_->nodeNumberGhost;

		std::cout << "------------------------------------------------------------------" << std::endl;
		std::cout << "	cellNumber:"
				  << "\t" << parameter_->cellNumber << std::endl;
		std::cout << "	nodeNumber:"
				  << "\t" << parameter_->nodeNumber << std::endl;

		cell_->resize(parameter_->cellNumber);
		node_->resize(parameter_->nodeNumber);

		//������Ԫ-�ڵ�
		cellVector::iterator iterCell;
		int ii = 0;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			(*iterCell).index = ii + 1;			  // index
			fileIn >> (*iterCell).cellFaceNumber; // cellFaceNumber

			// std::cout << (*iterCell).cellFaceNumber << std::endl;
			if ((*iterCell).cellFaceNumber == 3)
			{											 // cellNodeNumber
				(*iterCell).cellType_ = Triangle;		 // cellType_
				(*iterCell).cellNodeNumber = 3;			 // cellNode.first
				(*iterCell).cellFaceIndex.resize(3 + 1); //����cellFaceIndex
				(*iterCell).cellNode.resize(3 + 1);		 //����cellNode
			}
			else if ((*iterCell).cellFaceNumber == 4)
			{

				// std::cout << "true" << std::endl;

				(*iterCell).cellType_ = Quadrilateral;
				(*iterCell).cellNodeNumber = 4;
				(*iterCell).cellFaceIndex.resize(4 + 1);
				(*iterCell).cellNode.resize(4 + 1);
			}
			(*iterCell).cellNode.resize((*iterCell).cellNodeNumber + 1);
			for (int jj = 1; jj < (*iterCell).cellNodeNumber + 1; ++jj)
			{
				fileIn >> (*iterCell).cellNode[jj].first; // cellNode.second is the node position
														  //(*iterCell).cellNode[jj].first -= 1;			//����Ƿ���ȷ����20200220
			}
			++ii;
		}
		std::cout << "	node info has been read." << std::endl;

		//�����ڵ�����
		vertexVector::iterator iterNode;
		ii = 0;
		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{								// index
			(*iterNode).index = ii + 1; // nodePhysical
			fileIn >> (*iterNode).nodePhysical.x;
			fileIn >> (*iterNode).nodePhysical.y;
			(*iterNode).nodePhysical = (1.0 / parameter_->unit) * (*iterNode).nodePhysical;
			auto x = (*iterNode).nodePhysical.x, y = (*iterNode).nodePhysical.y;

			// DEBUG: rotate the grid 
			real theta = 1.0 * std::acos(-1) / 4.0;
			(*iterNode).nodePhysical.x = std::cos(theta) * x + std::sin(theta) * y;
			(*iterNode).nodePhysical.y = -std::sin(theta) * x + std::cos(theta) * y;
			++ii;
		}
		std::cout << "	node coordinate has been read." << std::endl;

		//��Ҫ�޸������ļ���20200218
		fileIn >> parameter_->boundaryNumber;
		fileIn >> parameter_->boundaryFaceNumber; // total boundary Face Number

		boundaryPair.resize(parameter_->boundaryNumber + 1);
		for (ii = 1; ii < parameter_->boundaryNumber + 1; ++ii)
		{
			fileIn >> boundaryPair[ii].first;  // tag
			fileIn >> boundaryPair[ii].second; // face number
		}
		boundary_->resize(parameter_->boundaryFaceNumber);
		boundaryVector::iterator iterBoundary;
		int boundaryTypeIn;
		ii = 0;
		for (iterBoundary = boundary_->begin(); iterBoundary != boundary_->end(); ++iterBoundary)
		{
			(*iterBoundary).index = ii + 1;
			fileIn >> (*iterBoundary).boundaryTag;
			fileIn >> boundaryTypeIn;
			if (boundaryTypeIn == 1)
			{											// boundaryTag
				(*iterBoundary).boundaryType_ = Wall;	// boundaryType_
				(*iterBoundary).faceNode.resize(2 + 1); // faceNode.first  i.e., node tag
				fileIn >> (*iterBoundary).faceNode[1].first;
				fileIn >> (*iterBoundary).faceNode[2].first;
			}
			else if (boundaryTypeIn == 2)
			{
				(*iterBoundary).boundaryType_ = FarField;
				(*iterBoundary).faceNode.resize(2 + 1);
				fileIn >> (*iterBoundary).faceNode[1].first;
				fileIn >> (*iterBoundary).faceNode[2].first;
			}
			else if (boundaryTypeIn == 3)
			{
				(*iterBoundary).boundaryType_ = Symmetric;
				(*iterBoundary).faceNode.resize(2 + 1);
				fileIn >> (*iterBoundary).faceNode[1].first;
				fileIn >> (*iterBoundary).faceNode[2].first;
			}
			else if (boundaryTypeIn == 4)
			{
				(*iterBoundary).boundaryType_ = OutFlow;
				(*iterBoundary).faceNode.resize(2 + 1);
				fileIn >> (*iterBoundary).faceNode[1].first;
				fileIn >> (*iterBoundary).faceNode[2].first;
			}
			else if (boundaryTypeIn == 5)
			{
				(*iterBoundary).boundaryType_ = InFlow;
				(*iterBoundary).faceNode.resize(2 + 1);
				fileIn >> (*iterBoundary).faceNode[1].first;
				fileIn >> (*iterBoundary).faceNode[2].first;
			}
			++ii;
		}
		std::cout << "	boundary info has been read." << std::endl;

		if (parameter_->isPeriodicBoundary)
		{
			fileIn >> boundaryNumberPeriodic;
			// boundaryPairPeriodic.resize(parameter_->boundaryNumber);
			boundaryPairPeriodicTag.resize(boundaryNumberPeriodic + 1);
			boundaryPeriodicOffside.resize(boundaryNumberPeriodic + 1);
			for (ii = 1; ii < boundaryNumberPeriodic + 1; ++ii)
			{
				fileIn >> boundaryPairPeriodicTag[ii].first;  // tag1 ��1��ʼ
				fileIn >> boundaryPairPeriodicTag[ii].second; // tag2
				fileIn >> boundaryPeriodicOffside[ii].x;	  // offside_x
				fileIn >> boundaryPeriodicOffside[ii].y;	  // offside_y

				for (iterBoundary = boundary_->begin(); iterBoundary != boundary_->end(); ++iterBoundary)
				{
					if ((*iterBoundary).boundaryTag == boundaryPairPeriodicTag[ii].first || (*iterBoundary).boundaryTag == boundaryPairPeriodicTag[ii].second)
					{
						(*iterBoundary).boundaryType_ = Periodic;
					}
				}
			}
			std::cout << "	*** periodic boundary info has been read." << std::endl;
		}

		if (parameter_->inputMeshOrder > 1)
		{
			for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
			{
				int np;
				fileIn >> np; // tag1 ��1��ʼ
				(*iterCell).o2PointNumber = np;
				(*iterCell).cellNodeO2Point.resize(np + 1);
				for (int ii = 1; ii < np + 1; ++ii)
				{
					(*iterCell).cellNodeO2Point[ii].first = ii;
					fileIn >> (*iterCell).cellNodeO2Point[ii].second.x; //
				}
				fileIn >> np;
				for (int ii = 1; ii < np + 1; ++ii)
				{
					fileIn >> (*iterCell).cellNodeO2Point[ii].second.y;
				}
			}
		}

		fileIn.close();
		std::cout << " ..grid file has been closed." << std::endl;

		return true;
	}

	bool grid::allocateCell()
	{
		//
		cellVector::iterator iterCell;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			(*iterCell).cellCellNumber = (*iterCell).cellFaceNumber;		  //�˴���ʱ������ʼ�� 20200219
			(*iterCell).cellCellIndex.resize((*iterCell).cellCellNumber + 1); //ע��˴���(*iterCell).cellCellNumber->current cell
			(*iterCell).cellFaceSideOff.resize((*iterCell).cellCellNumber + 1);
			(*iterCell).cellFaceWeight.resize((*iterCell).cellCellNumber + 1);
			for (int jj = 0; jj < (*iterCell).cellCellNumber + 1; ++jj)
			{
				(*iterCell).cellCellIndex[jj] = 0;
				(*iterCell).cellFaceWeight[jj] = 0.0;
				(*iterCell).cellFaceSideOff[jj].setZero();
			}
			//(*iterCell).cellCellIndex[(*iterCell).cellCellNumber] = 0;//20200225����
		}
		// cellCellNumber
		//����cellCellIndex
		//����&��ʼ��cellFaceSideOff
		//����&��ʼ��cellFaceWeight
		return true;
	}

	bool grid::allocateFace()
	{
		int totalFaceNumber = 0;
		cellVector::iterator iterCell;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			totalFaceNumber += (*iterCell).cellFaceNumber;
		}
		parameter_->faceNumber = (totalFaceNumber + parameter_->boundaryFaceNumber) / 2;
		parameter_->faceNumberGhost = parameter_->faceNumber;
		face_->resize(parameter_->faceNumber);

		faceVector::iterator iterFace;
		int ii = 0;
		for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace)
		{
			(*iterFace).faceNodeNumber = 2;
			(*iterFace).index = ii + 1;
			(*iterFace).faceCellIndex.resize(2 + 1);
			(*iterFace).faceNode.resize(2 + 1);
			for (int jj = 0; jj < 2 + 1; ++jj)
			{
				(*iterFace).faceNode[jj].first = 0;
				(*iterFace).faceNode[jj].second.setZero();
			}
			++ii;
		}
		//����faceCellIndex
		//
		return true;
	}

	bool grid::setNodeTopology()
	{

		std::cout << "seting node topology..." << std::endl;

		vertexVector::iterator iterNode;
		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{
			(*iterNode).nodeCellNumber = 0;
		}

		cellVector::iterator iterCell;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			for (int jj = 1; jj < (*iterCell).cellFaceNumber + 1; ++jj)
			{
				iterNode = node_->begin(); //ͨ���������ҵ�Ԫ��

				iterNode += (*iterCell).cellNode[jj].first - 1; //ע��ָ���ƶ��Ƿ���ȷ����20200220
				(*iterNode).nodeCellNumber += 1;
			}
		}

		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{
			//(*iterNode).nodeCellNumber = 0;
			(*iterNode).nodeCellIndex.resize((*iterNode).nodeCellNumber + 1);
			for (int jj = 1; jj < (*iterNode).nodeCellNumber + 1; ++jj)
			{
				(*iterNode).nodeCellIndex[jj] = 0;
			}
			(*iterNode).nodeCellNumber = 0;
		}
		// int ii = 0;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			for (int jj = 1; jj < (*iterCell).cellFaceNumber + 1; ++jj)
			{
				iterNode = node_->begin();						//ͨ���������ҵ�Ԫ��
				iterNode += (*iterCell).cellNode[jj].first - 1; //ע��ָ���ƶ��Ƿ���ȷ����20200220
				(*iterNode).nodeCellNumber += 1;				//ע���Լ�λ�ã���20200225
				(*iterNode).nodeCellIndex[(*iterNode).nodeCellNumber] = (*iterCell).index;
			}
			//++ii;
		}
		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{
			// std::cout << (*iterNode).nodeCellNumber << std::endl;
			//(*iterNode).nodeFaceNumber = 0;
			(*iterNode).nodeCellfaceNumber = 0;
			(*iterNode).nodeFaceIndex.resize((*iterNode).nodeCellNumber + 1);
			(*iterNode).nodeCellfaceIndex.resize(2 * (*iterNode).nodeCellNumber + 1);
			// for (int jj = 0; jj < (*iterNode).nodeCellNumber + 1; ++jj){
			//	(*iterNode).nodeFaceIndex[jj] = 0;
			// }
			for (int jj = 0; jj < 2 * (*iterNode).nodeCellNumber + 1; ++jj)
			{ //ע��ȷ�����Ϊɶ2��
				(*iterNode).nodeCellfaceIndex[jj] = 0;
			}
			// 20200225
			//����Ϊ����Ҫ2��
			//(*iterNode).nodeFaceIndex.resize((*iterNode).nodeCellNumber);
			// for (int jj = 0; jj < (*iterNode).nodeCellNumber; ++jj){
			//	(*iterNode).nodeFaceIndex[jj] = 0;
			// }

			// std::cout << (*iterNode).index << "\t" << (*iterNode).nodeFaceIndex.size() << std::endl;
		}
		std::cout << " ..node topology has been set." << std::endl;

		// system("pause");

		return true;
	}

	int grid::checkFace(int n1, int n2)
	{
		vertexVector::iterator iterNode;
		faceVector::iterator iterFace;
		iterNode = node_->begin() + n1 - 1; //ע��ָ���ƶ��Ƿ���ȷ����20200220//20200225ȷ�ϣ�����-1����Ϊ��0��ʼ
		for (int ii = 1; ii < (*iterNode).nodeCellfaceNumber + 1; ++ii)
		{
			int f1 = (*iterNode).nodeCellfaceIndex[ii];
			iterFace = face_->begin() + f1 - 1; //ע��ָ���ƶ��Ƿ���ȷ����20200220
			bool isTrue1 = ((*iterFace).faceNode[1].first == n1) && ((*iterFace).faceNode[2].first == n2);
			bool isTrue2 = ((*iterFace).faceNode[1].first == n2) && ((*iterFace).faceNode[2].first == n1);
			if (isTrue1 || isTrue2)
			{
				return f1;
			}
		}
		return 0;
	}

	// debug
	bool grid::exportTest(const std::string &filename)
	{
		std::ofstream fileOut;
		fileOut.open(filename.c_str());
		if (!fileOut)
		{
			std::cerr << "	error: unable to open \"" << filename << "\"\n";
			// return false;
			exit(1);
		}

		fileOut << parameter_->inputMeshOrder << std::endl;

		fileOut << parameter_->cellNumber << " " << parameter_->cellNumberGhost << " " << parameter_->nodeNumber << " " << parameter_->nodeNumberGhost << std::endl;

		cellVector::iterator iterCell;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			fileOut << (*iterCell).cellFaceNumber << " ";
			for (int jj = 1; jj < (*iterCell).cellFaceNumber + 1; ++jj)
			{
				fileOut << (*iterCell).cellNode[jj].first << " ";
			}
			fileOut << std::endl;
		}

		vertexVector::iterator iterNode;
		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{
			fileOut << (*iterNode).nodePhysical.x << " "
					<< (*iterNode).nodePhysical.y << std::endl;
		}

		fileOut << parameter_->boundaryNumber << std::endl;
		fileOut << parameter_->boundaryFaceNumber << std::endl;
		for (int ii = 1; ii < parameter_->boundaryNumber + 1; ++ii)
		{
			fileOut << boundaryPair[ii].first << " ";		 // tag
			fileOut << boundaryPair[ii].second << std::endl; // face number
		}

		boundaryVector::iterator iterBoundary;
		for (iterBoundary = boundary_->begin(); iterBoundary != boundary_->end(); ++iterBoundary)
		{
			fileOut << (*iterBoundary).boundaryTag << " ";
			fileOut << -(*iterBoundary).boundaryType_ << " ";
			fileOut << (*iterBoundary).faceNode[1].first << " ";
			fileOut << (*iterBoundary).faceNode[2].first << std::endl;
		}
		fileOut.close();
		return true;
	}

	bool grid::setFaceTopology()
	{

		std::cout << "seting face topology..." << std::endl;

		faceVector::iterator iterFace;
		for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace)
		{ // index
			(*iterFace).faceCellIndex[1] = -100;
			(*iterFace).faceCellIndex[2] = -100;
		}

		int ii = 0;
		int faceNumberCheck = 0;
		cellVector::iterator iterCell;
		vertexVector::iterator iterNode;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			for (int jj = 1; jj < (*iterCell).cellFaceNumber + 1; ++jj)
			{
				int n1 = (*iterCell).cellNode[jj].first;								  //�Ǳ�ţ�
				int n2 = (*iterCell).cellNode[jj % (*iterCell).cellFaceNumber + 1].first; // 20200227 -1!!
				int faceIndex = checkFace(n2, n1);
				// std::cout << (*iterCell).index << "\t" << faceIndex << std::endl;

				if (faceIndex == 0)
				{
					++ii;
					iterNode = node_->begin() + n1 - 1;
					point p1 = (*iterNode).nodePhysical;
					if ((*iterNode).nodeCellfaceNumber > 50)
					{ // nodeFaceNumber �˴���Ϊ�ֲ��ģ���Ϊ�˽���face��ϵ�õģ���Ϊ����Ԫѭ������
						std::cout << "	error: node is shared by too many edges/faces!" << std::endl;
						exit(1);
					}
					(*iterNode).nodeCellfaceNumber += 1; //ע���Լ�λ�ã���20200225
					(*iterNode).nodeCellfaceIndex[(*iterNode).nodeCellfaceNumber] = ii;

					iterNode = node_->begin() + n2 - 1;
					point p2 = (*iterNode).nodePhysical;

					if ((*iterNode).nodeCellfaceNumber > 50)
					{
						std::cout << "	error: node is shared by too many edges/faces!" << std::endl;
						exit(1);
					}
					(*iterNode).nodeCellfaceNumber += 1;
					(*iterNode).nodeCellfaceIndex[(*iterNode).nodeCellfaceNumber] = ii;

					iterFace = face_->begin();
					iterFace += (ii - 1);
					(*iterFace).faceCellIndex[1] = (*iterCell).index;
					(*iterFace).normalVector.x = (p2.y - p1.y);
					(*iterFace).normalVector.y = -(p2.x - p1.x);
					(*iterFace).faceNode[1].first = n1;
					(*iterFace).faceNode[2].first = n2;

					(*iterFace).faceNode[1].second = p1; // 20200226����
					(*iterFace).faceNode[2].second = p2;

					(*iterFace).area = (*iterFace).normalVector.length();
				}
				else
				{
					iterFace = face_->begin();
					iterFace += (faceIndex - 1);
					(*iterFace).faceCellIndex[2] = (*iterCell).index;
				}
			}
			if ((*iterCell).index == parameter_->cellNumber)
			{						  // [���+1=��Ŀ]//20200225ע������ı���Ƿ�������
				faceNumberCheck = ii; //��ţ�
			}
		}

		// std::cout << faceNumberCheck << "\t" << parameter_->faceNumber << std::endl;

		if (faceNumberCheck != parameter_->faceNumber)
		{
			std::cout << "	error: face number mismatching!" << std::endl;
			std::cout << faceNumberCheck << "\t" << parameter_->faceNumber << std::endl;
			exit(1);
		}

		std::cout << "	node to cell-face has been processed." << std::endl;

		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{
			(*iterNode).nodeFaceNumber = 0;
		}
		for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace)
		{
			for (int jj = 1; jj < 3; ++jj)
			{
				int n1 = (*iterFace).faceNode[jj].first;
				iterNode = node_->begin() + n1 - 1;
				(*iterNode).nodeFaceNumber += 1;
			}
		}
		// std::cout <<"1---------------------------------"  << std::endl;
		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{
			(*iterNode).nodeFaceIndex.resize((*iterNode).nodeFaceNumber + 1);
			for (int jj = 1; jj < (*iterNode).nodeFaceNumber + 1; ++jj)
			{
				(*iterNode).nodeFaceIndex[jj] = 0;
			}
			(*iterNode).nodeFaceNumber = 0;
		}

		// std::cout <<"2---------------------------------"  << std::endl;
		for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace)
		{

			for (int jj = 1; jj < 3; ++jj)
			{
				int n1 = (*iterFace).faceNode[jj].first;

				iterNode = node_->begin() + n1 - 1; //ͨ���������ҵ�Ԫ��
				// iterNode += (n1 - 1);
				(*iterNode).nodeFaceNumber += 1;

				// std::cout << "-----------------------" << std::endl;
				// std::cout << (*iterNode).index << std::endl;
				// std::cout << (*iterNode).nodeFaceNumber << std::endl;
				// std::cout << (*iterNode).nodeFaceIndex.size() << std::endl;

				(*iterNode).nodeFaceIndex[(*iterNode).nodeFaceNumber] = (*iterFace).index;

				// std::cout << (*iterNode).nodeFaceNumber << std::endl;

				// std::cout << (*iterNode).nodeFaceIndex[(*iterNode).nodeFaceNumber] << std::endl;
			}
		}

		// std::cout <<"3---------------------------------"  << std::endl;

		std::cout << " ..face topology has been set." << std::endl;
		// system("pause");
		return true;
	}

	bool grid::setBoundaryTopology()
	{

		std::cout << "creating boundary topology..." << std::endl;

		faceVector::iterator iterFace;
		vertexVector::iterator iterNode;
		boundaryVector::iterator iterBoundary;

		for (iterBoundary = boundary_->begin(); iterBoundary != boundary_->end(); ++iterBoundary)
		{
			int n1 = (*iterBoundary).faceNode[1].first;
			int n2 = (*iterBoundary).faceNode[2].first;
			// std::cout << n1 << "\t" << n2  << std::endl;
			int jj;
			iterNode = node_->begin();
			iterNode += (n1 - 1);
			for (jj = 1; jj < (*iterNode).nodeFaceNumber + 1; ++jj)
			{ //�����⣡��20200228
				int f1 = (*iterNode).nodeFaceIndex[jj];
				iterFace = face_->begin();
				iterFace += (f1 - 1);
				int m1 = (*iterFace).faceNode[1].first;
				int m2 = (*iterFace).faceNode[2].first;
				bool isMatch = false;
				if (n2 == m1 || n2 == m2)
				{
					isMatch = true;
					(*iterBoundary).faceIndex = f1;
				}
				// std::cout << m1 << "\t" << m2 << "\t" << isMatch << std::endl;
			}
		}
		std::cout << " ..boundary topology has been created." << std::endl;

		if (parameter_->isPeriodicBoundary)
		{
			std::cout << "creating periodic boundary topology..." << std::endl;
			for (int ii = 1; ii < boundaryNumberPeriodic + 1; ++ii)
			{
				int itag1 = boundaryPairPeriodicTag[ii].first; // itag1��itag2��Ϊ���ڱ�
				int itag2 = boundaryPairPeriodicTag[ii].second;
				// std::cout << itag1 << "\t" << itag2 << std::endl;
				for (iterBoundary = boundary_->begin(); iterBoundary != boundary_->end(); ++iterBoundary)
				{
					if (itag1 == (*iterBoundary).boundaryTag)
					{
						int n1 = (*iterBoundary).faceNode[1].first;
						int n2 = (*iterBoundary).faceNode[2].first;

						iterNode = node_->begin() + n1 - 1;
						point p1 = (*iterNode).nodePhysical;
						iterNode = node_->begin() + n2 - 1;
						point p2 = (*iterNode).nodePhysical;
						point p = 0.5 * (p1 + p2);
						int jj = 0;

						// std::cout << p.x << "\t" << p.y << std::endl;
						int flag = 1;
						boundaryVector::iterator iterBoundaryMatch;
						for (iterBoundaryMatch = boundary_->begin(); iterBoundaryMatch != boundary_->end(); ++iterBoundaryMatch)
						{
							if (itag2 == (*iterBoundaryMatch).boundaryTag)
							{
								++jj;

								int nMatch1 = (*iterBoundaryMatch).faceNode[1].first;
								int nMatch2 = (*iterBoundaryMatch).faceNode[2].first;
								iterNode = node_->begin() + nMatch1 - 1;
								point pMatch1 = (*iterNode).nodePhysical;
								iterNode = node_->begin() + nMatch2 - 1;
								point pMatch2 = (*iterNode).nodePhysical;
								point pMatch = 0.5 * (pMatch1 + pMatch2);
								point pOffside = pMatch - p;
								pOffside -= boundaryPeriodicOffside[ii]; // 20200318
								//�����Ǽ�ȥƫ��������Ϊƫ�����ǰ���������ʽ���ġ���ʹ��ʱ���Ǽ���ƫ��������Ҫע�⣡��
								real pOffsideValue = pOffside.length();
								if (pOffsideValue < parameter_->EPS)
								{
									flag = 0;
									//�˴�����flag�жϣ���������ָ�ꡣ��Ϊԭ���ĺ����ǰ�ָ��ѭ���ģ�
									//���û���ҵ���ָ��ᳬ����Χ�������ǰ��յ�����ѭ���ģ�û�ҵ��Ļ�ָ�겻�ᳬ����Χ��
									break;
								}
							}
						}
						// iterBoundaryMatch -= 1;//20200317 ע���Լ�˳�򣡣�

						// if (jj > boundaryPair[itag2].second){//20200225ע�������ָ�꣡��
						if (flag == 1)
						{ // 20200225ע�������ָ�꣡��
							std::cout << "	error: fail to match periodic boundary!" << std::endl;
							exit(1);
						}
						(*iterBoundary).matchIndex = (*iterBoundaryMatch).faceIndex; // indexΪ�ֲ�ָ��
						int f1 = (*iterBoundary).faceIndex;
						int fMatch1 = (*iterBoundaryMatch).faceIndex;
						iterFace = face_->begin() + f1 - 1;
						faceVector::iterator iterFaceMatch;
						iterFaceMatch = face_->begin() + fMatch1 - 1;
						(*iterFace).faceCellIndex[2] = (*iterFaceMatch).faceCellIndex[1];
						(*iterFace).faceProperty = Periodic;
						(*iterFace).sideOff = boundaryPeriodicOffside[ii];
					}
				}
			}
			std::cout << " ..periodic boundary topology has been created." << std::endl;
		}
		return true;
	}

	bool grid::setFaceProperty()
	{

		std::cout << "setting face property topology..." << std::endl;

		faceVector::iterator iterFace;
		boundaryVector::iterator iterBoundary;
		for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace)
		{
			int n1 = (*iterFace).faceNode[1].first;
			int n2 = (*iterFace).faceNode[2].first;
			if ((*iterFace).faceCellIndex[2] < 0)
			{
				for (int ii = 1; ii < parameter_->boundaryNumber + 1; ++ii)
				{
					int jj = 1;
					for (iterBoundary = boundary_->begin(); iterBoundary != boundary_->end(); ++iterBoundary)
					{
						if ((*iterBoundary).boundaryTag == ii)
						{ //ע��߽��ǩ�ţ���20200220//20200225 boundaryTag��1��ʼ
							bool isTrue1 = ((*iterBoundary).faceNode[1].first == n1 && (*iterBoundary).faceNode[2].first == n2);
							bool isTrue2 = ((*iterBoundary).faceNode[1].first == n2 && (*iterBoundary).faceNode[2].first == n1);
							if (isTrue1 || isTrue2)
							{
								(*iterFace).faceCellIndex[2] = (*iterBoundary).boundaryType_; //ע�������20200226 -(*iterBoundary).boundaryType_
								// std::cout << (*iterFace).faceCellIndex[2] << std::endl;
								break; //�����ȷ��break��������forѭ������20200221
							}
							++jj;
						}
					}
					if (jj <= boundaryPair[ii].second)
					{
						break;
					}
					if (ii > parameter_->boundaryNumber)
					{
						std::cout << "	error: fail to find boundary type!" << std::endl;
						exit(1);
					}
				}
			}
		}

		std::cout << " ..face property has been set." << std::endl;
		//	system("pause");
		return true;
	}

	bool grid::saveGrid(const std::string &filename)
	{
		std::cout << "saving grid..." << std::endl;

		std::ofstream fileOut;
		fileOut.open(filename.c_str());
		if (!fileOut)
		{
			std::cerr << "	error: unable to open \"" << filename << "\"\n";
			// return false;
			exit(1);
		}
		std::cout << "	in Tecplot format, to file : " << filename << "..." << std::endl;
		fileOut << "VARIABLES = \"x\", \"y\", \n";
		fileOut << "ZONE N =" << parameter_->nodeNumber << ","
				<< "E=" << parameter_->cellNumber << ","
				<< "F=FEPOINT,ET=QUADRILATERAL" << std::endl;
		// fileOut << "ZONE N =" << parameter_->nodeNumber << "," << "E=" << parameter_->cellNumber << "," << "F=FEPOINT,ET=TRIANGLE" << std::endl;;
		fileOut << std::setprecision(15);

		vertexVector::iterator iterNode;
		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{
			fileOut << (*iterNode).nodePhysical.x << " "
					<< (*iterNode).nodePhysical.y << std::endl;
		}
		cellVector::iterator iterCell;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			if ((*iterCell).cellType_ == Triangle)
			{
				// std::cout << "Triangle" << "\t" << (*iterCell).cellType_ << "\t" << (*iterCell).cellNodeNumber << std::endl;

				for (int jj = 1; jj < (*iterCell).cellNodeNumber + 1; ++jj)
				{
					fileOut << (*iterCell).cellNode[jj].first << " ";
				}
				fileOut << (*iterCell).cellNode[(*iterCell).cellNodeNumber].first;
				fileOut << std::endl;
			}
			else if ((*iterCell).cellType_ == Quadrilateral)
			{
				// std::cout << "Quadrilateral" << "\t" << (*iterCell).cellType_ << "\t" << (*iterCell).cellNodeNumber << std::endl;
				for (int jj = 1; jj < (*iterCell).cellNodeNumber + 1; ++jj)
				{
					fileOut << (*iterCell).cellNode[jj].first << " ";
				}
				fileOut << std::endl;
			}
		}
		fileOut.close();
		std::cout << " ..grid has been saved." << std::endl;
		return true;
	}

	bool grid::setCellTopology()
	{
		std::cout << "creating cell topology..." << std::endl;
		vertexVector::iterator iterNode;
		faceVector::iterator iterFace;
		cellVector::iterator iterCell;

		// 20200226����
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			for (int ii = 1; ii < (*iterCell).cellNodeNumber + 1; ++ii)
			{
				int n1 = (*iterCell).cellNode[ii].first;
				iterNode = node_->begin() + n1 - 1;
				point p1 = (*iterNode).nodePhysical;
				(*iterCell).cellNode[ii].second = p1;
			}
		}

		// std::cout << "1---------------------------------" << std::endl;

		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			(*iterCell).cellFaceNumber = 0;
			(*iterCell).cellCellNumber = 0;
		}

		for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace)
		{
			int cl = (*iterFace).faceCellIndex[1];
			int cr = (*iterFace).faceCellIndex[2];
			if (cl < 0)
			{
				std::cout << "	error: left cell should be inner cell!" << std::endl;
				exit(1);
			}
			iterCell = cell_->begin() + cl - 1;
			(*iterCell).cellFaceNumber += 1; //ע���Լ�λ�ã���20200225
			(*iterCell).cellFaceIndex[(*iterCell).cellFaceNumber] = (*iterFace).index;

			if (cr > 0 && (*iterFace).faceProperty != Periodic)
			{
				iterCell = cell_->begin() + cr - 1;
				(*iterCell).cellFaceNumber += 1;
				(*iterCell).cellFaceIndex[(*iterCell).cellFaceNumber] = (*iterFace).index;
			}
		}

		// std::cout << "2---------------------------------" << std::endl;

		parameter_->innerCellNumber = 0;
		parameter_->boundaryCellNumber = 0;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			for (int jj = 1; jj < (*iterCell).cellFaceNumber + 1; ++jj)
			{
				int f1 = (*iterCell).cellFaceIndex[jj]; //���
				for (int kk = 1; kk < 3; ++kk)
				{
					iterFace = face_->begin() + f1 - 1; //��������λ��
					int c1 = (*iterFace).faceCellIndex[kk];
					if (c1 != (*iterCell).index && c1 > 0)
					{
						(*iterCell).cellCellNumber += 1; //ע��˴��Լ�λ�ã���20200225

						// std::cout << "--------------------" << std::endl;
						// std::cout << (*iterCell).cellCellNumber << std::endl;
						// std::cout << (*iterCell).cellCellIndex.size() << std::endl;

						(*iterCell).cellCellIndex[(*iterCell).cellCellNumber] = c1;
						(*iterCell).cellFaceSideOff[(*iterCell).cellCellNumber] = (*iterFace).sideOff;
					}
				}
			}

			if ((*iterCell).cellCellNumber == (*iterCell).cellFaceNumber)
			{
				parameter_->innerCellNumber += 1;
				(*iterCell).isBoundaryCell = false;
				//
			}
			else if ((*iterCell).cellCellNumber > 0 && (*iterCell).cellCellNumber < (*iterCell).cellFaceNumber)
			{
				parameter_->boundaryCellNumber += 1;
				(*iterCell).isBoundaryCell = true;
			}
			else
			{
				std::cout << "	error: incorrect cell-cell vs. cell-face!" << std::endl;
				exit(1);
			}
		}

		// std::cout << "3---------------------------------" << std::endl;

		if (parameter_->innerCellNumber + parameter_->boundaryCellNumber != parameter_->cellNumber)
		{
			std::cout << "	error: cell count mistake!" << std::endl;
			exit(1);
		}

		// std::cout << "---------------------------------------" << std::endl;
		std::cout << "	inner cell number:"
				  << "\t" << parameter_->innerCellNumber << std::endl;
		std::cout << "	boundary cell number:"
				  << "\t" << parameter_->boundaryCellNumber << std::endl;

		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			(*iterCell).cellCellIndex[0] = (*iterCell).index; //ע��20200225����
			(*iterCell).cellSideOff.setZero();
		}

		// std::cout << "---------------------------------------" << std::endl;

		faceVector::iterator iterFace_;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			for (int jj = 1; jj < (*iterCell).cellFaceNumber + 1; ++jj)
			{
				int ff = (*iterCell).cellFaceIndex[jj];
				iterFace_ = face_->begin() + ff - 1;
				int cl = (*iterFace_).faceCellIndex[1];
				int cr = (*iterFace_).faceCellIndex[2];
				if (cl > 0 && cr > 0)
				{
					(*iterCell).cellInnerFaceNumber += 1;
				}
				else
				{
					(*iterCell).cellBoundaryFaceNumber += 1;
				}
			}
		}
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			if ((*iterCell).cellInnerFaceNumber == (*iterCell).cellFaceNumber)
			{
				(*iterCell).boundaryCellType_ = InnerCell;
			}
			else if ((*iterCell).cellBoundaryFaceNumber > 1)
			{
				(*iterCell).boundaryCellType_ = MixedCell;
				parameter_->boundaryMixedCellNumber += 1;
			}
			else if ((*iterCell).cellBoundaryFaceNumber == 1)
			{
				for (int jj = 1; jj < (*iterCell).cellFaceNumber + 1; ++jj)
				{
					int ff = (*iterCell).cellFaceIndex[jj];
					iterFace_ = face_->begin() + ff - 1;
					int cl = (*iterFace_).faceCellIndex[1];
					int cr = (*iterFace_).faceCellIndex[2];
					if (cr == Wall)
					{
						(*iterCell).boundaryCellType_ = WallCell;
						parameter_->boundaryWallCellNumber += 1;
					}
					else if (cr == FarField)
					{
						(*iterCell).boundaryCellType_ = FarFieldCell;
						parameter_->boundaryFarfieldCellNumber += 1;
					}
					else if (cr == Symmetric)
					{
						(*iterCell).boundaryCellType_ = SymmetricCell;
						parameter_->boundarySymmetricCellNumber += 1;
					}
					else if (cr == OutFlow)
					{
						(*iterCell).boundaryCellType_ = OutFlowCell;
						parameter_->boundaryOutflowCellNumber += 1;
					}
					else if (cr == InFlow)
					{
						(*iterCell).boundaryCellType_ = InFlowCell;
						parameter_->boundaryInflowCellNumber += 1;
					}
				}
			}
		}

		//	for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace){
		//		int cl = (*iterFace).faceCellIndex[1];
		//		int cr = (*iterFace).faceCellIndex[2];
		//
		//		iterCell = cell_->begin() + cl - 1;
		//		//std::cout << cl << "\t" << cr << std::endl;
		//		//std::cout << cr << "\t" << (*iterCell).index << std::endl;
		//		//std::cout << (*iterCell).isCountByboundaryCellType << std::endl;
		//		//std::cout << "---------------------------------------" << std::endl;
		//		if (cr == Wall && (*iterCell).isCountByboundaryCellType == false){
		//			(*iterCell).boundaryCellType_ = WallCell;
		//			parameter_->boundaryWallCellNumber += 1;
		//			(*iterCell).isCountByboundaryCellType = true;
		//		}
		//		else if (cr == FarField && (*iterCell).isCountByboundaryCellType == false){
		//			(*iterCell).boundaryCellType_ = FarFieldCell;
		//			parameter_->boundaryFarfieldCellNumber += 1;
		//			(*iterCell).isCountByboundaryCellType = true;
		//		}
		//		else if (cr == Symmetric && (*iterCell).isCountByboundaryCellType == false){
		//			(*iterCell).boundaryCellType_ = SymmetricCell;
		//			parameter_->boundarySymmetricCellNumber += 1;
		//			(*iterCell).isCountByboundaryCellType = true;
		//		}
		//		else if (cr == OutFlow && (*iterCell).isCountByboundaryCellType == false){
		//			(*iterCell).boundaryCellType_ = OutFlowCell;
		//			parameter_->boundaryOutflowCellNumber += 1;
		//			(*iterCell).isCountByboundaryCellType = true;
		//		}
		//		else if (cr == InFlow && (*iterCell).isCountByboundaryCellType == false){
		//			(*iterCell).boundaryCellType_ = InFlowCell;
		//			parameter_->boundaryInflowCellNumber += 1;
		//			(*iterCell).isCountByboundaryCellType = true;
		//		}
		//		else if ((*iterFace).faceProperty == Periodic && (*iterCell).isCountByboundaryCellType == false){
		//			(*iterCell).boundaryCellType_ = InnerCell;//ע��������ع����ֵ�Ӱ�죡��20200318
		//			parameter_->boundaryPeriodicCellNumber += 1;
		//			(*iterCell).isCountByboundaryCellType = true;
		//		}
		//	}

		std::cout << "------------------------------------------------------------------" << std::endl;
		std::cout << "	wall cell number:"
				  << "\t" << parameter_->boundaryWallCellNumber << std::endl;
		std::cout << "	far field cell number:"
				  << "\t" << parameter_->boundaryFarfieldCellNumber << std::endl;
		std::cout << "	symmetric cell number:"
				  << "\t" << parameter_->boundarySymmetricCellNumber << std::endl;
		std::cout << "	outflow cell number:"
				  << "\t" << parameter_->boundaryOutflowCellNumber << std::endl;
		std::cout << "	inflow cell number:"
				  << "\t" << parameter_->boundaryInflowCellNumber << std::endl;
		std::cout << "	mixed cell number:"
				  << "\t" << parameter_->boundaryMixedCellNumber << std::endl;

		if ((parameter_->boundaryWallCellNumber + parameter_->boundaryFarfieldCellNumber + parameter_->boundarySymmetricCellNumber + parameter_->boundaryOutflowCellNumber + parameter_->boundaryInflowCellNumber + parameter_->boundaryMixedCellNumber) != parameter_->boundaryCellNumber)
		{
			std::cout << "	warning: boundary cell count mistake!" << std::endl;
			//��һ����Ԫ�������߽��ʱ������ֲ�ƥ�����������������ο�
			// exit(1);
		}

		if (parameter_->inputMeshOrder > 1)
		{
			for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
			{
				for (int jj = 1; jj < (*iterCell).cellFaceNumber + 1; ++jj)
				{
					int n1 = (*iterCell).cellNode[jj].first;								  //�Ǳ�ţ�
					int n2 = (*iterCell).cellNode[jj % (*iterCell).cellFaceNumber + 1].first; // 20200227 -1!!
					int faceIndex = checkFace(n2, n1);
					iterFace = face_->begin() + faceIndex - 1;
					(*iterFace).o2PointNumber = 1;
					(*iterFace).faceNodeO2Point.resize((*iterFace).o2PointNumber + 1);
					(*iterFace).faceNodeO2Point[1].second = (*iterCell).cellNodeO2Point[jj].second;
				}
				// the following codes indicate that the cell-face order is not correspond to
				//	the face order that dertermined by cell-node order. so we need to match them.
				//	std::cout << "------------------------" << std::endl;
				//	std::cout << "cell:" << "\t" << (*iterCell).index << "\t" << (*iterCell).cellFaceNumber << std::endl;
				//	for (int ii = 1; ii < (*iterCell).cellNodeNumber + 1; ++ii){
				//		std::cout << (*iterCell).cellNode[ii].first << std::endl;
				//
				//	}
				//	std::cout << "-------" << std::endl;
				//	for (int ii = 1; ii < (*iterCell).cellFaceNumber + 1; ++ii){
				//		int ff = (*iterCell).cellFaceIndex[ii];
				//		iterFace = face_->begin() + ff - 1;
				//		std::cout << (*iterFace).faceNode[1].first << "\t" << (*iterFace).faceNode[2].first << std::endl;
				//
				//	}
			}

			//	//test
			//	for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace){
			//		int n1 = (*iterFace).faceNode[1].first;
			//		int n2 = (*iterFace).faceNode[2].first;
			//		point p3 = (*iterFace).faceNodeO2Point[1].second;
			//
			//		std::cout << "----" << std::endl;
			//		std::cout << n1 << "\t" << n2 << "\t" << p3.x << "\t" << p3.y << std::endl;
			//	}
		}

		std::cout << " ..cell topology has been created." << std::endl;

		//	system("pause");

		return true;
	}

	bool grid::freeGrid()
	{
		vertexVector::iterator iterNode;
		faceVector::iterator iterFace;
		cellVector::iterator iterCell;

		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{
			std::vector<int>().swap((*iterNode).nodeFaceIndex);
			std::vector<int>().swap((*iterNode).nodeCellIndex);
		}

		for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace)
		{
			std::vector<std::pair<int, point>>().swap((*iterFace).faceNode);
			std::vector<int>().swap((*iterFace).faceCellIndex);
		}

		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			std::vector<std::pair<int, point>>().swap((*iterCell).cellNode);
			std::vector<int>().swap((*iterCell).cellFaceIndex);
			std::vector<int>().swap((*iterCell).cellCellIndex);
			std::vector<point>().swap((*iterCell).cellFaceSideOff);
			std::vector<real>().swap((*iterCell).cellFaceWeight);
		}

		//ȷ�������ͷŷ�ʽ�Ƿ���ȷ����20200224
		std::vector<vertex>().swap(*node_);
		std::vector<face>().swap(*face_);
		std::vector<cell>().swap(*cell_);

		return true;
	}

	bool grid::getGeometryInfo(
		cellGaussDataVector *cellGaussData,
		GaussIntegralFaceO1Grid<fO> *gaussFaceIntegral_,
		GaussIntegralCellO1Grid<vO> *gaussCellIntegral_)
	{

		std::cout << "computing geometric info..." << std::endl;
		cellVector::iterator iterCell;
		cellGaussDataVector::iterator iterGaussData;

		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			iterGaussData = cellGaussData->begin() + (*iterCell).index - 1;
			// get the function values at Gaussian points
			int PG;
			if ((*iterGaussData).cellType_ == Triangle)
			{
				PG = static_cast<int>((*iterGaussData).PGTri[vO]);
			}
			else if ((*iterGaussData).cellType_ == Quadrilateral)
			{
				PG = static_cast<int>((*iterGaussData).PGQuad[vO]);
			}

			// std::cout << "PG" << "\t" << PG << std::endl;

			real *f = new real[PG];
			real *fx = new real[PG];
			real *fy = new real[PG];
			for (int ii = 0; ii < PG; ++ii)
			{
				f[ii] = 1.0;										// (*iterCell).cellNode[ii].second;
				fx[ii] = (*iterGaussData).gaussPairVector_[ii].p.x; // (*iterCell).cellNode[ii].second;
				fy[ii] = (*iterGaussData).gaussPairVector_[ii].p.y;

				// std::cout << "--------------------------------------" << std::endl;
				// std::cout << f[ii - 1] << "\t" << fx[ii - 1] << "\t" << fy[ii - 1] << std::endl;
			}
			gaussCellIntegral_->getIntegral(
				f,
				(*iterCell).index,
				cellGaussData,
				(*iterCell).volume);
			gaussCellIntegral_->getIntegral(
				fx,
				(*iterCell).index,
				cellGaussData,
				(*iterCell).baryCenter.x);
			gaussCellIntegral_->getIntegral(
				fy,
				(*iterCell).index,
				cellGaussData,
				(*iterCell).baryCenter.y);
			(*iterCell).baryCenter = (1.0 / (*iterCell).volume) * (*iterCell).baryCenter;

			real *px = new real[(*iterCell).cellNodeNumber];
			real *py = new real[(*iterCell).cellNodeNumber];
			for (int ii = 1; ii < (*iterCell).cellNodeNumber + 1; ++ii)
			{
				px[ii - 1] = (*iterCell).cellNode[ii].second.x;
				py[ii - 1] = (*iterCell).cellNode[ii].second.y;
			}

			real xMax = getMax(px, (*iterCell).cellNodeNumber);
			real xMin = getMin(px, (*iterCell).cellNodeNumber);
			real yMax = getMax(py, (*iterCell).cellNodeNumber);
			real yMin = getMin(py, (*iterCell).cellNodeNumber);

			(*iterCell).lengthReference.x = (xMax - xMin) / 2.0;
			(*iterCell).lengthReference.y = (yMax - yMin) / 2.0;

			delete[] px;
			delete[] py;
			px = NULL;
			py = NULL;

			delete[] f;
			delete[] fx;
			delete[] fy;
			f = NULL;
			fx = NULL;
			fy = NULL;
		}

		// baryCenter??
		std::cout << " ..geometric info has been computed." << std::endl;
		return true;
	}

	bool grid::getGeometryInfo(
		cellGaussDataVector *cellGaussData,
		GaussIntegralFaceO2Grid<fO> *gaussFaceIntegral_,
		GaussIntegralCellO2Grid<vO> *gaussCellIntegral_)
	{

		std::cout << "computing geometric info..." << std::endl;
		cellVector::iterator iterCell;
		cellGaussDataVector::iterator iterGaussData;

		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			iterGaussData = cellGaussData->begin() + (*iterCell).index - 1;
			// get the function values at Gaussian points
			int PG;
			if ((*iterGaussData).cellType_ == Triangle)
			{
				PG = static_cast<int>((*iterGaussData).PGTri[vO]);
			}
			else if ((*iterGaussData).cellType_ == Quadrilateral)
			{
				PG = static_cast<int>((*iterGaussData).PGQuad[vO]);
			}

			// std::cout << "PG" << "\t" << PG << std::endl;

			real *f = new real[PG];
			real *fx = new real[PG];
			real *fy = new real[PG];
			for (int ii = 0; ii < PG; ++ii)
			{
				f[ii] = 1.0;										// (*iterCell).cellNode[ii].second;
				fx[ii] = (*iterGaussData).gaussPairVector_[ii].p.x; // (*iterCell).cellNode[ii].second;
				fy[ii] = (*iterGaussData).gaussPairVector_[ii].p.y;

				// std::cout << "--------------------------------------" << std::endl;
				// std::cout << f[ii - 1] << "\t" << fx[ii - 1] << "\t" << fy[ii - 1] << std::endl;
			}
			gaussCellIntegral_->getIntegral(
				f,
				(*iterCell).index,
				cellGaussData,
				(*iterCell).volume);
			gaussCellIntegral_->getIntegral(
				fx,
				(*iterCell).index,
				cellGaussData,
				(*iterCell).baryCenter.x);
			gaussCellIntegral_->getIntegral(
				fy,
				(*iterCell).index,
				cellGaussData,
				(*iterCell).baryCenter.y);
			(*iterCell).baryCenter = (1.0 / (*iterCell).volume) * (*iterCell).baryCenter;

			real *px = new real[(*iterCell).cellNodeNumber];
			real *py = new real[(*iterCell).cellNodeNumber];
			for (int ii = 1; ii < (*iterCell).cellNodeNumber + 1; ++ii)
			{
				px[ii - 1] = (*iterCell).cellNode[ii].second.x;
				py[ii - 1] = (*iterCell).cellNode[ii].second.y;
			}

			real xMax = getMax(px, (*iterCell).cellNodeNumber);
			real xMin = getMin(px, (*iterCell).cellNodeNumber);
			real yMax = getMax(py, (*iterCell).cellNodeNumber);
			real yMin = getMin(py, (*iterCell).cellNodeNumber);

			(*iterCell).lengthReference.x = (xMax - xMin) / 2.0;
			(*iterCell).lengthReference.y = (yMax - yMin) / 2.0;

			delete[] px;
			delete[] py;
			px = NULL;
			py = NULL;

			delete[] f;
			delete[] fx;
			delete[] fy;
			f = NULL;
			fx = NULL;
			fy = NULL;
		}

		// baryCenter??
		std::cout << " ..geometric info has been computed." << std::endl;
		return true;
	}

	bool grid::debugTest(
		faceGaussDataVector *faceGaussData,
		GaussIntegralFaceO1Grid<fO> *gaussFaceIntegral_)
	{

		std::cout << "computing face geometric info..." << std::endl;
		faceVector::iterator iterFace;
		faceGaussDataVector::iterator iterGaussData;

		for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace)
		{
			iterGaussData = faceGaussData->begin() + (*iterFace).index - 1;
			int PG = static_cast<int>((*iterGaussData).PG[fO]);

			point *flux = new point[PG];
			real *fI = new real[PG];
			real *weight = new real[PG];
			real *cofJacobi = new real[PG];
			real parametricArea;
			for (int ii = 0; ii < PG; ++ii)
			{
				flux[ii] = (*iterGaussData).normalVector;
				flux[ii] = (1.0 / (*iterGaussData).area) * flux[ii];

				fI[ii] = 1.0;
				weight[ii] = (*iterGaussData).parametricValue[ii].second;
				cofJacobi[ii] = (*iterGaussData).gaussPairVector_[ii].JacobiCof;
				parametricArea = (*iterGaussData).parametricArea;
			}

			//ʹ�õ�λ��������Ϊͨ�����в��ԣ���20200304

			real result;
			gaussFaceIntegral_->getIntegral(
				flux,
				(*iterFace).index,
				faceGaussData,
				result);

			std::cout << " ------------------------" << std::endl;
			std::cout << "face node: " << (*iterGaussData).faceNode[1].first << "\t" << (*iterGaussData).faceNode[2].first << std::endl;
			std::cout << "integral result: " << result << std::endl;
			std::cout << "reference: " << (*iterFace).area << std::endl;
			std::cout << "difference: " << result - (*iterFace).area << std::endl;
			std::cout << " --------" << std::endl;

			gaussFaceIntegral_->getIntegral(
				PG,
				fI,
				weight,
				cofJacobi,
				parametricArea,
				result);
			std::cout << "integral result: " << result << std::endl;
			std::cout << "reference: " << (*iterFace).area << std::endl;
			std::cout << "difference: " << result - (*iterFace).area << std::endl;

			// std::cout << PG << std::endl;
			// for (int ii = 0; ii < PG; ++ii){
			//	std::cout << "gauss point "<< ii << std::endl;
			//	std::cout << (*iterGaussData).gaussPairVector_[ii].p.x << "\t" << (*iterGaussData).gaussPairVector_[ii].p.y << std::endl;
			//	std::cout << (*iterGaussData).gaussPairVector_[ii].normalVector.x << "\t" << (*iterGaussData).gaussPairVector_[ii].normalVector.y << std::endl;
			//	std::cout << (*iterGaussData).gaussPairVector_[ii].JacobiCof << std::endl;
			// }

			delete[] fI;
			delete[] flux;
			delete[] weight;
			delete[] cofJacobi;
			fI = NULL;
			flux = NULL;
			weight = NULL;
			cofJacobi = NULL;
		}

		for (iterGaussData = faceGaussData->begin(); iterGaussData != faceGaussData->end(); ++iterGaussData)
		{
			std::cout << " --------" << std::endl;
			std::cout << (*iterGaussData).faceNodeNumber << std::endl;
			for (int ii = 1; ii < (*iterGaussData).faceNodeNumber + 1; ++ii)
			{
				std::cout << (*iterGaussData).faceNode[ii].second.x << "\t" << (*iterGaussData).faceNode[ii].second.y << std::endl;
			}
		}

		std::cout << " ..face geometric info has been computed." << std::endl;
		return true;
	}

	bool grid::debugTest(
		faceGaussDataVector *faceGaussData,
		GaussIntegralFaceO2Grid<fO> *gaussFaceIntegral_)
	{

		std::cout << "computing face geometric info..." << std::endl;
		faceVector::iterator iterFace;
		faceGaussDataVector::iterator iterGaussData;

		for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace)
		{
			iterGaussData = faceGaussData->begin() + (*iterFace).index - 1;
			int PG = static_cast<int>((*iterGaussData).PG[fO]);

			point *flux = new point[PG];
			real *fI = new real[PG];
			real *weight = new real[PG];
			real *cofJacobi = new real[PG];
			real parametricArea;
			for (int ii = 0; ii < PG; ++ii)
			{
				flux[ii] = (*iterGaussData).normalVector;
				flux[ii] = (1.0 / (*iterGaussData).area) * flux[ii];

				fI[ii] = 1.0;
				weight[ii] = (*iterGaussData).parametricValue[ii].second;
				cofJacobi[ii] = (*iterGaussData).gaussPairVector_[ii].JacobiCof;
				parametricArea = (*iterGaussData).parametricArea;
			}

			//ʹ�õ�λ��������Ϊͨ�����в��ԣ���20200304

			real result;
			gaussFaceIntegral_->getIntegral(
				flux,
				(*iterFace).index,
				faceGaussData,
				result);

			std::cout << " ------------------------" << std::endl;
			std::cout << "face node: " << (*iterGaussData).faceNode[1].first << "\t" << (*iterGaussData).faceNode[2].first << std::endl;
			std::cout << "integral result: " << result << std::endl;
			std::cout << "reference: " << (*iterFace).area << std::endl;
			std::cout << "difference: " << result - (*iterFace).area << std::endl;
			std::cout << " --------" << std::endl;

			gaussFaceIntegral_->getIntegral(
				PG,
				fI,
				weight,
				cofJacobi,
				parametricArea,
				result);
			std::cout << "integral result: " << result << std::endl;
			std::cout << "reference: " << (*iterFace).area << std::endl;
			std::cout << "difference: " << result - (*iterFace).area << std::endl;

			// std::cout << PG << std::endl;
			// for (int ii = 0; ii < PG; ++ii){
			//	std::cout << "gauss point "<< ii << std::endl;
			//	std::cout << (*iterGaussData).gaussPairVector_[ii].p.x << "\t" << (*iterGaussData).gaussPairVector_[ii].p.y << std::endl;
			//	std::cout << (*iterGaussData).gaussPairVector_[ii].normalVector.x << "\t" << (*iterGaussData).gaussPairVector_[ii].normalVector.y << std::endl;
			//	std::cout << (*iterGaussData).gaussPairVector_[ii].JacobiCof << std::endl;
			// }

			delete[] fI;
			delete[] flux;
			delete[] weight;
			delete[] cofJacobi;
			fI = NULL;
			flux = NULL;
			weight = NULL;
			cofJacobi = NULL;
		}

		for (iterGaussData = faceGaussData->begin(); iterGaussData != faceGaussData->end(); ++iterGaussData)
		{
			std::cout << " --------" << std::endl;
			std::cout << (*iterGaussData).faceNodeNumber << std::endl;
			for (int ii = 1; ii < (*iterGaussData).faceNodeNumber + 1; ++ii)
			{
				std::cout << (*iterGaussData).faceNode[ii].second.x << "\t" << (*iterGaussData).faceNode[ii].second.y << std::endl;
			}
		}

		std::cout << " ..face geometric info has been computed." << std::endl;
		return true;
	}

	// debug
	bool grid::exportTestGeoInfo(const std::string &filename)
	{
		std::cout << "saving grid..." << std::endl;

		int lineControl = 5;

		std::ofstream fileOut;
		fileOut.open(filename.c_str());
		if (!fileOut)
		{
			std::cerr << "	error: unable to open \"" << filename << "\"\n";
			// return false;
			exit(1);
		}
		std::cout << "	in Tecplot format, to file : " << filename << "..." << std::endl;
		fileOut << "VARIABLES = \"x\", \"y\", \"vol\", \"b.x\", \"b.y\"\n";
		fileOut << "ZONE N =" << parameter_->nodeNumber << ","
				<< "E=" << parameter_->cellNumber << ","
				<< "VARLOCATION=([1-2]=NODAL,[3-5]=CELLCENTERED)" << std::endl;
		fileOut << ", "
				<< "DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL" << std::endl;
		// fileOut << "ZONE N =" << parameter_->nodeNumber << "," << "E=" << parameter_->cellNumber << "," << "F=FEPOINT,ET=TRIANGLE" << std::endl;;
		fileOut << std::setprecision(15);

		vertexVector::iterator iterNode;
		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{
			fileOut << (*iterNode).nodePhysical.x << "\t";
			if (((*iterNode).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << std::endl;
		for (iterNode = node_->begin(); iterNode != node_->end(); ++iterNode)
		{
			fileOut << (*iterNode).nodePhysical.y << "\t";
			if (((*iterNode).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << std::endl;
		cellVector::iterator iterCell;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			fileOut << (*iterCell).volume << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << std::endl;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			fileOut << (*iterCell).baryCenter.x << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << std::endl;
		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			fileOut << (*iterCell).baryCenter.y << "\t";
			if (((*iterCell).index % lineControl) == 0)
				fileOut << '\n';
		}
		fileOut << std::endl;

		for (iterCell = cell_->begin(); iterCell != cell_->end(); ++iterCell)
		{
			if ((*iterCell).cellType_ == Triangle)
			{
				// std::cout << "Triangle" << "\t" << (*iterCell).cellType_ << "\t" << (*iterCell).cellNodeNumber << std::endl;

				for (int jj = 1; jj < (*iterCell).cellNodeNumber + 1; ++jj)
				{
					fileOut << (*iterCell).cellNode[jj].first << " ";
				}
				fileOut << (*iterCell).cellNode[(*iterCell).cellNodeNumber].first;
				fileOut << std::endl;
			}
			else if ((*iterCell).cellType_ == Quadrilateral)
			{
				// std::cout << "Quadrilateral" << "\t" << (*iterCell).cellType_ << "\t" << (*iterCell).cellNodeNumber << std::endl;
				for (int jj = 1; jj < (*iterCell).cellNodeNumber + 1; ++jj)
				{
					fileOut << (*iterCell).cellNode[jj].first << " ";
				}
				fileOut << std::endl;
			}
		}
		fileOut << std::endl;
		fileOut.close();
		std::cout << "grid has been saved." << std::endl;
		return true;
	}

	bool grid::adjustNV()
	{

		cellVector::iterator iterCell;
		faceVector::iterator iterFace;

		// left -> right cell

		for (iterFace = face_->begin(); iterFace != face_->end(); ++iterFace)
		{ // index
			int cl = (*iterFace).faceCellIndex[1];
			int cr = (*iterFace).faceCellIndex[2];
			point fSideOff = (*iterFace).sideOff;
			if (cr < 0)
			{
				iterCell = cell_->begin() + cl - 1;
				point bc;
				bc.setZero();
				for (int ii = 1; ii < (*iterCell).cellNodeNumber + 1; ++ii)
				{
					bc += (*iterCell).cellNode[ii].second;
				}
				bc /= static_cast<real>((*iterCell).cellNodeNumber);
				point pMid = 0.5 * ((*iterFace).faceNode[1].second + (*iterFace).faceNode[2].second);
				point outDirection = pMid - bc;
				if (getInnerProduct(outDirection, (*iterFace).normalVector) < 0.0)
				{
					(*iterFace).normalVector = -1.0 * (*iterFace).normalVector;
				}
			}
			else
			{
				iterCell = cell_->begin() + cl - 1;
				point bcL;
				bcL.setZero();
				for (int ii = 1; ii < (*iterCell).cellNodeNumber + 1; ++ii)
				{
					bcL += (*iterCell).cellNode[ii].second;
				}
				bcL /= static_cast<real>((*iterCell).cellNodeNumber);
				iterCell = cell_->begin() + cr - 1;
				point bcR;
				bcR.setZero();
				for (int ii = 1; ii < (*iterCell).cellNodeNumber + 1; ++ii)
				{
					bcR += (*iterCell).cellNode[ii].second;
				}
				bcR /= static_cast<real>((*iterCell).cellNodeNumber);
				point outDirection = bcR - bcL - fSideOff; // 20201117
				if (getInnerProduct(outDirection, (*iterFace).normalVector) < 0.0)
				{
					(*iterFace).normalVector = -1.0 * (*iterFace).normalVector;
				}
			}
		}

		// system("pause");

		return true;
	}

	bool grid::adjustNodeFaceOrder()
	{
		for (auto &c : *cell_)
		{
			if (c.cellType_ == Quadrilateral)
			{
				int a[4][2] = {{1, 2}, {2, 3}, {3, 4}, {4, 1}};
				int f[4][2], n[4][2];
				for (int ifacecell = 1; ifacecell <= c.cellFaceNumber; ifacecell++)
				{
					auto iface = c.cellFaceIndex[ifacecell];
					f[ifacecell - 1][0] = (*face_)[iface - 1].faceNode[1].first;
					f[ifacecell - 1][1] = (*face_)[iface - 1].faceNode[2].first;
					n[ifacecell - 1][0] = c.cellNode[a[ifacecell - 1][0]].first;
					n[ifacecell - 1][1] = c.cellNode[a[ifacecell - 1][1]].first;
				}
				int perm[4];
				for (int i = 0; i < 4; i++)
				{
					int j;
					for (j = 0; j < 4; j++)
						if ((f[j][0] == n[i][0] && f[j][1] == n[i][1]) || (f[j][0] == n[i][1] && f[j][1] == n[i][0]))
							break;
					assert(j < 4);
					perm[i] = j;
				}
				auto oldcellFaceIndex = c.cellFaceIndex;
				auto oldcellFaceWeight = c.cellFaceWeight;
				auto oldcellFaceSideoff = c.cellFaceSideOff;
				for (int i = 0; i < 4; i++)
				{
					c.cellFaceIndex[i + 1] = oldcellFaceIndex[perm[i] + 1];
					c.cellFaceWeight[i + 1] = oldcellFaceWeight[perm[i] + 1];
					c.cellFaceSideOff[i + 1] = oldcellFaceSideoff[perm[i] + 1];
				}
			}
		}
	}
}