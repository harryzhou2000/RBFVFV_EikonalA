#ifndef _GAUSSINTEGRAL_H
#define _GAUSSINTEGRAL_H
#include "TypeDefine.h"
#include "Point.h"
#include "Parameter.h"
#include "Geometry.h"
//#include "Variable.h"

//PG=3
//20180218������3����������

namespace ScalarCfv
{
	
	//-------------------------------------------------------------------------
	class gaussPair
	{
	public:
		point p;
		point normalVector;
		real  JacobiCof;
		gaussPair(real x_, real y_, real nx_, real ny_, real jc_) :p(x_, y_), normalVector(nx_, ny_), JacobiCof(jc_){};
		gaussPair(): p(0, 0), normalVector(0, 0), JacobiCof(0){};
	};
	typedef std::vector<gaussPair> gaussPairVector;
	typedef std::vector<gaussPairVector> gaussPairVectorVector;

	//-------------------------------------------------------------------------
	//O1 grid 

	class faceGaussData : public face{

	public:
		faceGaussData() :face(){};
		~faceGaussData() {
			gaussPairVector().swap(gaussPairVector_);
			std::vector<unsigned int>().swap(PG);
			std::vector<std::pair<point, real>>().swap(parametricValue);
			//point.x <-> kesi;  point.y <->eta;
		};
	public:
		gaussPairVector gaussPairVector_;
		std::vector<unsigned int> PG;
		std::vector<std::pair<point, real>> parametricValue;
		real parametricArea;
	};
	typedef std::vector<faceGaussData> faceGaussDataVector;

	template<unsigned int OG>
	class GaussIntegralFaceO1Grid
	{
	public:
		~GaussIntegralFaceO1Grid() {};
		GaussIntegralFaceO1Grid() {};

		//unsigned int PG[OG + 1];

		//typedef std::vector<std::pair<int, point>> pairVectorPhysical;
		//typedef std::vector<std::pair<real, real>> pairVectorParameter;
		//pairVectorParameter parametricValue;
		//gaussPairVectorVector faceData;

	public:
		virtual bool allocateArray(
			parameter* parameter,
			faceVector* faceVector,
			faceGaussDataVector* faceGaussData);
		virtual bool initializeArray(
			faceVector* faceVector,
			faceGaussDataVector* faceGaussData);

		//flux
		virtual bool getIntegral(
			point f[],	//Gauss���ϵĺ���ֵ���������ͨ�������������
			int idxFace,		//���λ��ָ�꣬��ΪҪʹ��faceData
			faceGaussDataVector* faceGaussData, 
			real& result);

		//overload scaler
		virtual bool getIntegral(
			unsigned int PG,
			real f[],		
			real weight[],
			real cofJacobi[],
			real parametricArea,
			real& result);


		virtual bool freeArray(void);
		virtual point getGaussPointFace(
			int idxFace,
			std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
			int idx,
			faceGaussDataVector* faceGaussData);									//���λ��ָ�꣬��ΪҪʹ��faceData
		virtual point getNormalVector(
			int idxFace,
			std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
			int idx,
			faceGaussDataVector* faceGaussData);									//���λ��ָ�꣬��ΪҪʹ��faceData

		virtual real getJacobiCof(
			int idxFace,
			std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
			int idx,
			faceGaussDataVector* faceGaussData);									//���λ��ָ�꣬��ΪҪʹ��faceData

	};

	template<unsigned int OG>
	bool GaussIntegralFaceO1Grid<OG>::allocateArray(
		parameter* parameter,
		faceVector* faceVector,
		faceGaussDataVector* faceGaussData){

		std::cout << "allocating face Gauss integral array..." << std::endl;

		faceGaussData->resize(parameter->faceNumber);

		faceVector::iterator iterFace;
		faceGaussDataVector::iterator iterGaussData;

		for (iterFace = faceVector->begin(); iterFace != faceVector->end(); ++iterFace){
			int ii = (*iterFace).index;
			iterGaussData = faceGaussData->begin() + ii - 1;
			(*iterGaussData).index = (*iterFace).index;
			(*iterGaussData).faceProperty = (*iterFace).faceProperty;
			(*iterGaussData).faceNodeNumber = (*iterFace).faceNodeNumber;
			(*iterGaussData).faceCellNumber = (*iterFace).faceCellNumber;
			(*iterGaussData).spectralRadius = (*iterFace).spectralRadius;
			(*iterGaussData).area = (*iterFace).area;
			(*iterGaussData).areaReference = (*iterFace).areaReference;
			(*iterGaussData).normalVector = orientFlagF*(*iterFace).normalVector;
			(*iterGaussData).sideOff = (*iterFace).sideOff;


			// switch (OG){
			// case 1:
			// 	(*iterGaussData).PG.resize(1 + 1);
			// 	(*iterGaussData).PG[1] = (OG + 1) / 2;
			// 	break;
			// case 3:
			// 	(*iterGaussData).PG.resize(3 + 1);
			// 	(*iterGaussData).PG[3] = (OG + 1) / 2;
			// 	break;
			// case 5:
			// 	(*iterGaussData).PG.resize(5 + 1);
			// 	(*iterGaussData).PG[5] = (OG + 1) / 2;
			// 	break;
			// default:
			// 	std::cout << "	error: currently only support PG from 1-3!" << std::endl;
			// 	exit(1);
			// }
			(*iterGaussData).PG.resize(OG + 1);
			(*iterGaussData).PG[OG] = (OG + 1) / 2;

			(*iterGaussData).parametricValue.resize((*iterGaussData).PG[OG]);

			(*iterGaussData).gaussPairVector_.resize((*iterGaussData).PG[OG]);

			(*iterGaussData).faceNodeNumber = (*iterFace).faceNodeNumber;
			(*iterGaussData).faceNode.resize((*iterFace).faceNodeNumber + 1);

			(*iterGaussData).faceCellIndex.resize(2 + 1);
			(*iterGaussData).faceCellIndex[1] = (*iterFace).faceCellIndex[1];
			(*iterGaussData).faceCellIndex[2] = (*iterFace).faceCellIndex[2];
			for (int jj = 1; jj < (*iterFace).faceNodeNumber + 1; ++jj){
				(*iterGaussData).faceNode[jj].first = (*iterFace).faceNode[jj].first;
				(*iterGaussData).faceNode[jj].second = (*iterFace).faceNode[jj].second;
			}

		}

		std::cout << " ..face Gauss integral array has been allocated." << std::endl;
		return true;
	}

	template<unsigned int OG>
	point GaussIntegralFaceO1Grid<OG>::getGaussPointFace(
		int idxFace,
		std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
		int idx,
		faceGaussDataVector* faceGaussData){

		faceGaussDataVector::iterator iterGaussData;
		iterGaussData = faceGaussData->begin() + idxFace - 1;
		real kesi = (*iterGaussData).parametricValue[idx].first.x;
		real baseValue1 = 1.0 - kesi;
		real baseValue2 = kesi;
		point p = baseValue1*physicalPoint[0].second + baseValue2*physicalPoint[1].second;
		return p;
	}

	template<unsigned int OG>
	point GaussIntegralFaceO1Grid<OG>::getNormalVector(
		int idxFace,
		std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
		int idx,
		faceGaussDataVector* faceGaussData){

		faceGaussDataVector::iterator iterGaussData;
		iterGaussData = faceGaussData->begin() + idxFace - 1;
		real kesi = (*iterGaussData).parametricValue[idx].first.x;
		real dBaseValue1dKesi = -1.0;
		real dBaseValue2dKesi = 1.0;
		point normalVector;

		normalVector.x = -(physicalPoint[0].second.y*dBaseValue1dKesi + physicalPoint[1].second.y*dBaseValue2dKesi);
		normalVector.y = physicalPoint[0].second.x*dBaseValue1dKesi + physicalPoint[1].second.x*dBaseValue2dKesi;

		//ע��ȷ������20200303!!
		if (getInnerProduct(normalVector, (*iterGaussData).normalVector) < 0){
			normalVector = orientFlagG*normalVector;
		}


		//debug
		//std::cout << " ------------------------" << std::endl;
		//std::cout << "face node: " << (*iterGaussData).faceNode[1].first << "\t" << (*iterGaussData).faceNode[2].first << std::endl;
		//std::cout << (*iterGaussData).normalVector.x << "\t" << (*iterGaussData).normalVector.y << std::endl;
		//std::cout << normalVector.x << "\t" << normalVector.y << std::endl;
		//std::cout << "difference: " << "\t" << normalVector.x - (*iterGaussData).normalVector.x
		//	<< "\t" << normalVector.y - (*iterGaussData).normalVector.y << std::endl;

		return normalVector;
	}

	template<unsigned int OG>
	real GaussIntegralFaceO1Grid<OG>::getJacobiCof(
		int idxFace,
		std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
		int idx,
		faceGaussDataVector* faceGaussData){

		faceGaussDataVector::iterator iterGaussData;
		iterGaussData = faceGaussData->begin() + idxFace - 1;
		real kesi = (*iterGaussData).parametricValue[idx].first.x;
		real dBaseValue1dKesi = -1.0;
		real dBaseValue2dKesi = 1.0;
		real dxdKesi, dydKesi, cofJacobi;
		point origin = physicalPoint[0].second;

		//point origin;
		//origin.setZero();

		dxdKesi = -(physicalPoint[0].second.x - origin.x)*dBaseValue1dKesi 
			+ (physicalPoint[1].second.x - origin.x)*dBaseValue2dKesi;
		dydKesi = -(physicalPoint[0].second.y - origin.y)*dBaseValue1dKesi 
			+ (physicalPoint[1].second.y - origin.y)*dBaseValue2dKesi;

		cofJacobi = std::fabs(std::sqrt(std::pow(dxdKesi, 2) + std::pow(dydKesi, 2)));

		return cofJacobi;
	}

	template<unsigned int OG>
	bool GaussIntegralFaceO1Grid<OG>::initializeArray(
		faceVector* faceVector,
		faceGaussDataVector* faceGaussData){
		//PG = 3

		std::cout << "initializing face Gauss integral array..." << std::endl;


		faceGaussDataVector::iterator iterGaussData;
		for (iterGaussData = faceGaussData->begin(); iterGaussData != faceGaussData->end(); ++iterGaussData){
			switch (OG){
			case 1:
				(*iterGaussData).parametricValue[0].first.x = (0 + 1)*0.5;
				(*iterGaussData).parametricValue[0].second = 2.0*0.5;
				break;
			case 3:
				(*iterGaussData).parametricValue[0].first.x = (-0.577350269189626 + 1)*0.5;
				(*iterGaussData).parametricValue[0].second = 1.0*0.5;
				(*iterGaussData).parametricValue[1].first.x = (0.577350269189626 + 1)*0.5;
				(*iterGaussData).parametricValue[1].second = 1.0*0.5;
				break;
			case 5:
				(*iterGaussData).parametricValue[0].first.x = (-0.774596669241483 + 1)*0.5;
				(*iterGaussData).parametricValue[0].second = 0.555555555555555*0.5;
				(*iterGaussData).parametricValue[1].first.x = (0 + 1)*0.5;
				(*iterGaussData).parametricValue[1].second = 0.888888888888888*0.5;
				(*iterGaussData).parametricValue[2].first.x = (0.774596669241483 + 1)*0.5;
				(*iterGaussData).parametricValue[2].second = 0.555555555555555*0.5;
				break;
			case 7:
				(*iterGaussData).parametricValue[0].first.x = (-0.861136311594054 + 1) * 0.5;
				(*iterGaussData).parametricValue[0].second = 0.347854845137452 * 0.5;
				(*iterGaussData).parametricValue[1].first.x = (-0.339981043584857 + 1) * 0.5;
				(*iterGaussData).parametricValue[1].second = 0.652145154862546 * 0.5;
				(*iterGaussData).parametricValue[2].first.x = (0.339981043584857 + 1) * 0.5;
				(*iterGaussData).parametricValue[2].second = 0.652145154862546 * 0.5;
				(*iterGaussData).parametricValue[3].first.x = (0.861136311594054 + 1) * 0.5;
				(*iterGaussData).parametricValue[3].second = 0.347854845137452 * 0.5;
				break;
			default:
				std::cout << "	error: currently only support PG from 1-3!" << std::endl;
				exit(1);
			}
			(*iterGaussData).parametricArea = 1;

			std::pair<int, point>* node = new std::pair<int, point>[(*iterGaussData).faceNodeNumber];
			for (int kk = 0; kk < (*iterGaussData).faceNodeNumber; ++kk){
				node[kk].first = (*iterGaussData).faceNode[kk + 1].first;
				node[kk].second = (*iterGaussData).faceNode[kk + 1].second;
			}
			for (int jj = 0; jj < static_cast<int>((*iterGaussData).PG[OG]); ++jj){
				point physicalPoint;
				point nromalVector;
				real cofJacobi;
				physicalPoint = getGaussPointFace((*iterGaussData).index, node, jj, faceGaussData);
				nromalVector = getNormalVector((*iterGaussData).index, node, jj, faceGaussData);
				cofJacobi = getJacobiCof((*iterGaussData).index, node, jj, faceGaussData);

				(*iterGaussData).gaussPairVector_[jj].p = physicalPoint;
				(*iterGaussData).gaussPairVector_[jj].normalVector = nromalVector;
				(*iterGaussData).gaussPairVector_[jj].JacobiCof = cofJacobi;
			}
			delete[] node;
			node = NULL;
		}
		std::cout << " ..face Gauss integral array initialized." << std::endl;
		return true;
		
	}


	template<unsigned int OG>
	bool GaussIntegralFaceO1Grid<OG>::getIntegral(
		point f[],	//Gauss���ϵĺ���ֵ���������ͨ�������������
		int idxFace,		//���λ��ָ�꣬��ΪҪʹ��faceData
		faceGaussDataVector* faceGaussData,
		real& result){
		result = 0;
		faceGaussDataVector::iterator iterGaussData;
		iterGaussData = faceGaussData->begin() + idxFace - 1;
		for (unsigned int ii = 0; ii < static_cast<int>((*iterGaussData).PG[OG]); ++ii){
			//parametric space |s|=1
			//ע��ȷ�ϻ��ֹ�ʽ�Ƿ���ȷ����20200218
			result += getInnerProduct(f[ii], (*iterGaussData).gaussPairVector_[ii].normalVector) 
				* (*iterGaussData).parametricValue[ii].second * (*iterGaussData).parametricArea;
		}
		return true;
	}
 
	//overload
	template<unsigned int OG>
	bool GaussIntegralFaceO1Grid<OG>::getIntegral(
		unsigned int PG,
		real f[],
		real weight[],
		real cofJacobi[],
		real parametricArea,
		real& result){
		result = 0;
		for (int ii = 0; ii < static_cast<int>(PG); ++ii){
			//parametric space |s|=1
			//ע��ȷ�ϻ��ֹ�ʽ�Ƿ���ȷ����20200218
			result += f[ii] * weight[ii] * cofJacobi[ii] * parametricArea;
		}
		return true;
	}

	template<unsigned int OG>
	bool GaussIntegralFaceO1Grid<OG>::freeArray(void){

		//std::vector<std::pair<real, real>>().swap(parametricValue);
		//gaussPairVectorVector().swap(faceData);

		return true;
	}

	//O2 grid
	template<unsigned int OG>
	class GaussIntegralFaceO2Grid : public GaussIntegralFaceO1Grid<OG>
	{
	public:
		~GaussIntegralFaceO2Grid() {};
		GaussIntegralFaceO2Grid() {};

	public:
		//bool allocateArray(const parameter* parameter) override;
		//bool initializeArray(const faceVector* faceVector) override;
		//bool getIntegral(const point f[], const int idx, const point vecRef) override;
		////bool freeArray(void) override;
		//bool allocateArray(parameter* parameter) override;
		//bool initializeArray(faceVector* faceVector) override;
		//bool getIntegral(
		//	point f[],	//Gauss���ϵĺ���ֵ���������ͨ�������������
		//	int idx,		//���λ��ָ�꣬��ΪҪʹ��faceData
		//	point vecRef, //�ο�����
		//	real& result) override;
		//bool freeArray(void) override;
		//point getGaussPointFace(
		//	std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O2������������Ϊ3
		//	int idx) override;						//���λ��ָ�꣬��ΪҪʹ��faceData
		//point getNormalVector(
		//	std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O2������������Ϊ3
		//	int idx) override;						//���λ��ָ�꣬��ΪҪʹ��faceData
		bool allocateArray(
			parameter* parameter,
			faceVector* faceVector,
			faceGaussDataVector* faceGaussData) override;
		bool initializeArray(
			faceVector* faceVector,
			faceGaussDataVector* faceGaussData) override;

		//flux
		virtual bool getIntegral(
			point f[],	//Gauss���ϵĺ���ֵ���������ͨ�������������
			int idxFace,		//���λ��ָ�꣬��ΪҪʹ��faceData
			faceGaussDataVector* faceGaussData,
			real& result) override;

		//overload scaler
		virtual bool getIntegral(
			unsigned int PG,
			real f[],
			real weight[],
			real cofJacobi[],
			real parametricArea,
			real& result) override;

		bool freeArray(void) override;
		point getGaussPointFace(
			int idxFace,
			std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
			int idx,
			faceGaussDataVector* faceGaussData) override;									//���λ��ָ�꣬��ΪҪʹ��faceData
		point getNormalVector(
			int idxFace,
			std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
			int idx,
			faceGaussDataVector* faceGaussData) override;

		real getJacobiCof(
			int idxFace,
			std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
			int idx,
			faceGaussDataVector* faceGaussData) override;									//���λ��ָ�꣬��ΪҪʹ��faceData
	};


	template<unsigned int OG>
	bool GaussIntegralFaceO2Grid<OG>::allocateArray(
		parameter* parameter,
		faceVector* faceVector,
		faceGaussDataVector* faceGaussData){

		std::cout << "allocating face Gauss integral array..." << std::endl;

		faceGaussData->resize(parameter->faceNumber);

		faceVector::iterator iterFace;
		faceGaussDataVector::iterator iterGaussData;

		for (iterFace = faceVector->begin(); iterFace != faceVector->end(); ++iterFace){
			int ii = (*iterFace).index;
			iterGaussData = faceGaussData->begin() + ii - 1;
			(*iterGaussData).index = (*iterFace).index;
			(*iterGaussData).faceProperty = (*iterFace).faceProperty;
			(*iterGaussData).faceNodeNumber = (*iterFace).faceNodeNumber;
			(*iterGaussData).faceCellNumber = (*iterFace).faceCellNumber;
			(*iterGaussData).spectralRadius = (*iterFace).spectralRadius;
			(*iterGaussData).area = (*iterFace).area;
			(*iterGaussData).areaReference = (*iterFace).areaReference;
			(*iterGaussData).normalVector = orientFlagF*(*iterFace).normalVector;
			(*iterGaussData).sideOff = (*iterFace).sideOff;


			switch (OG){
			case 1:
				(*iterGaussData).PG.resize(1 + 1);
				(*iterGaussData).PG[1] = (OG + 1) / 2;
				break;
			case 3:
				(*iterGaussData).PG.resize(3 + 1);
				(*iterGaussData).PG[3] = (OG + 1) / 2;
				break;
			case 5:
				(*iterGaussData).PG.resize(5 + 1);
				(*iterGaussData).PG[5] = (OG + 1) / 2;
				break;
			default:
				std::cout << "	error: currently only support PG from 1-3!" << std::endl;
				exit(1);
			}

			(*iterGaussData).parametricValue.resize((*iterGaussData).PG[OG]);

			(*iterGaussData).gaussPairVector_.resize((*iterGaussData).PG[OG]);

			(*iterGaussData).faceNodeNumber = (*iterFace).faceNodeNumber;
			(*iterGaussData).faceNode.resize((*iterFace).faceNodeNumber + 1);

			(*iterGaussData).faceCellIndex.resize(2 + 1);
			(*iterGaussData).faceCellIndex[1] = (*iterFace).faceCellIndex[1];
			(*iterGaussData).faceCellIndex[2] = (*iterFace).faceCellIndex[2];
			for (int jj = 1; jj < (*iterFace).faceNodeNumber + 1; ++jj){
				(*iterGaussData).faceNode[jj].first = (*iterFace).faceNode[jj].first;
				(*iterGaussData).faceNode[jj].second = (*iterFace).faceNode[jj].second;
			}

			if ((*iterFace).o2PointNumber != 0){
				(*iterGaussData).o2PointNumber = (*iterFace).o2PointNumber;
				std::pair<int, point> o2Point;
				o2Point.first = (*iterFace).faceNodeO2Point[1].first;
				o2Point.second = (*iterFace).faceNodeO2Point[1].second;
				(*iterGaussData).faceNodeNumber += 1;
				(*iterGaussData).faceNode.push_back(o2Point);
			}


		//	std::cout << "----" << std::endl;
		//	std::cout << (*iterGaussData).index << "\t" << (*iterGaussData).faceNodeNumber << std::endl;
		//	std::cout << (*iterGaussData).faceNode[1].first << "\t" << (*iterGaussData).faceNode[2].first 
		//		<< "\t" << (*iterGaussData).faceNode[3].second.x << "\t" << (*iterGaussData).faceNode[3].second.y << std::endl;
		}

		//system("pause");

		std::cout << " ..face Gauss integral array has been allocated." << std::endl;
		return true;
	}


	template<unsigned int OG>
	point GaussIntegralFaceO2Grid<OG>::getGaussPointFace(
		int idxFace,
		std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
		int idx,
		faceGaussDataVector* faceGaussData){

		faceGaussDataVector::iterator iterGaussData;
		iterGaussData = faceGaussData->begin() + idxFace - 1;
		real kesi = (*iterGaussData).parametricValue[idx].first.x;
		real baseValue1 = (-1.0 + kesi) * (-1.0 + 2.0*kesi);
		real baseValue2 = kesi * (-1.0 + 2.0*kesi);
		real baseValue3 = -4 * (-1.0 + kesi) * kesi;
		point p = baseValue1*physicalPoint[0].second 
			+ baseValue2*physicalPoint[1].second
			+ baseValue3*physicalPoint[2].second;
		return p;
	}


	template<unsigned int OG>
	point GaussIntegralFaceO2Grid<OG>::getNormalVector(
		int idxFace,
		std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
		int idx,
		faceGaussDataVector* faceGaussData){

		faceGaussDataVector::iterator iterGaussData;
		iterGaussData = faceGaussData->begin() + idxFace - 1;
		real kesi = (*iterGaussData).parametricValue[idx].first.x;
		real dBaseValue1dKesi = -3.0 + 4.0*kesi;
		real dBaseValue2dKesi = -1.0 + 4.0*kesi;
		real dBaseValue3dKesi =  4.0 - 8.0*kesi;

		point normalVector;

		normalVector.x = -(physicalPoint[0].second.y*dBaseValue1dKesi 
			+ physicalPoint[1].second.y*dBaseValue2dKesi
			+ physicalPoint[2].second.y*dBaseValue3dKesi);
		normalVector.y = physicalPoint[0].second.x*dBaseValue1dKesi 
			+ physicalPoint[1].second.x*dBaseValue2dKesi
			+ physicalPoint[2].second.x*dBaseValue3dKesi;

		//ע��ȷ������20200303!!
		if (getInnerProduct(normalVector, (*iterGaussData).normalVector) < 0){
			normalVector = orientFlagG*normalVector;
		}


		//debug
		//std::cout << " ------------------------" << std::endl;
		//std::cout << "face node: " << (*iterGaussData).faceNode[1].first << "\t" << (*iterGaussData).faceNode[2].first << std::endl;
		//std::cout << (*iterGaussData).normalVector.x << "\t" << (*iterGaussData).normalVector.y << std::endl;
		//std::cout << normalVector.x << "\t" << normalVector.y << std::endl;
		//std::cout << "difference: " << "\t" << normalVector.x - (*iterGaussData).normalVector.x
		//	<< "\t" << normalVector.y - (*iterGaussData).normalVector.y << std::endl;

		return normalVector;
	}

	template<unsigned int OG>
	real GaussIntegralFaceO2Grid<OG>::getJacobiCof(
		int idxFace,
		std::pair<int, point> physicalPoint[],	//������Ľڵ㣬�����ռ��е����꣬O1������������Ϊ2
		int idx,
		faceGaussDataVector* faceGaussData){

		faceGaussDataVector::iterator iterGaussData;
		iterGaussData = faceGaussData->begin() + idxFace - 1;
		real kesi = (*iterGaussData).parametricValue[idx].first.x;
		real dBaseValue1dKesi = -3.0 + 4.0*kesi;
		real dBaseValue2dKesi = -1.0 + 4.0*kesi;
		real dBaseValue3dKesi = 4.0 - 8.0*kesi;
		real dxdKesi, dydKesi, cofJacobi;
		point origin = physicalPoint[0].second;

		//point origin;
		//origin.setZero();

		dxdKesi = (physicalPoint[0].second.x - origin.x)*dBaseValue1dKesi
			+ (physicalPoint[1].second.x - origin.x)*dBaseValue2dKesi
			+ (physicalPoint[2].second.x - origin.x)*dBaseValue3dKesi;
		dydKesi = (physicalPoint[0].second.y - origin.y)*dBaseValue1dKesi
			+ (physicalPoint[1].second.y - origin.y)*dBaseValue2dKesi
			+ (physicalPoint[2].second.y - origin.y)*dBaseValue3dKesi;

		cofJacobi = std::fabs(std::sqrt(std::pow(dxdKesi, 2) + std::pow(dydKesi, 2)));

		return cofJacobi;
	}

	template<unsigned int OG>
	bool GaussIntegralFaceO2Grid<OG>::initializeArray(
		faceVector* faceVector,
		faceGaussDataVector* faceGaussData){
		//PG = 3

		std::cout << "initializing face Gauss integral array..." << std::endl;


		faceGaussDataVector::iterator iterGaussData;
		for (iterGaussData = faceGaussData->begin(); iterGaussData != faceGaussData->end(); ++iterGaussData){
			switch (OG){
			case 1:
				(*iterGaussData).parametricValue[0].first.x = (0 + 1)*0.5;
				(*iterGaussData).parametricValue[0].second = 2.0*0.5;
				break;
			case 3:
				(*iterGaussData).parametricValue[0].first.x = (-0.577350269189626 + 1)*0.5;
				(*iterGaussData).parametricValue[0].second = 1.0*0.5;
				(*iterGaussData).parametricValue[1].first.x = (0.577350269189626 + 1)*0.5;
				(*iterGaussData).parametricValue[1].second = 1.0*0.5;
				break;
			case 5:
				(*iterGaussData).parametricValue[0].first.x = (-0.774596669241483 + 1)*0.5;
				(*iterGaussData).parametricValue[0].second = 0.555555555555555*0.5;
				(*iterGaussData).parametricValue[1].first.x = (0 + 1)*0.5;
				(*iterGaussData).parametricValue[1].second = 0.888888888888888*0.5;
				(*iterGaussData).parametricValue[2].first.x = (0.774596669241483 + 1)*0.5;
				(*iterGaussData).parametricValue[2].second = 0.555555555555555*0.5;
				break;
			default:
				std::cout << "	error: currently only support PG from 1-3!" << std::endl;
				exit(1);
			}
			(*iterGaussData).parametricArea = 1;

			std::pair<int, point>* node = new std::pair<int, point>[(*iterGaussData).faceNodeNumber];
			for (int kk = 0; kk < (*iterGaussData).faceNodeNumber; ++kk){
				node[kk].first = (*iterGaussData).faceNode[kk + 1].first;
				node[kk].second = (*iterGaussData).faceNode[kk + 1].second;
			}
			for (int jj = 0; jj < static_cast<int>((*iterGaussData).PG[OG]); ++jj){
				point physicalPoint;
				point nromalVector;
				real cofJacobi;
				physicalPoint = getGaussPointFace((*iterGaussData).index, node, jj, faceGaussData);
				nromalVector = getNormalVector((*iterGaussData).index, node, jj, faceGaussData);
				cofJacobi = getJacobiCof((*iterGaussData).index, node, jj, faceGaussData);

				(*iterGaussData).gaussPairVector_[jj].p = physicalPoint;
				(*iterGaussData).gaussPairVector_[jj].normalVector = nromalVector;
				(*iterGaussData).gaussPairVector_[jj].JacobiCof = cofJacobi;
			}
			delete[] node;
			node = NULL;
		}
		std::cout << " ..face Gauss integral array initialized." << std::endl;
		return true;

	}


	template<unsigned int OG>
	bool GaussIntegralFaceO2Grid<OG>::getIntegral(
		point f[],	//Gauss���ϵĺ���ֵ���������ͨ�������������
		int idxFace,		//���λ��ָ�꣬��ΪҪʹ��faceData
		faceGaussDataVector* faceGaussData,
		real& result){
		result = 0;
		faceGaussDataVector::iterator iterGaussData;
		iterGaussData = faceGaussData->begin() + idxFace - 1;
		for (unsigned int ii = 0; ii < static_cast<int>((*iterGaussData).PG[OG]); ++ii){
			//parametric space |s|=1
			//ע��ȷ�ϻ��ֹ�ʽ�Ƿ���ȷ����20200218
			result += getInnerProduct(f[ii], (*iterGaussData).gaussPairVector_[ii].normalVector)
				* (*iterGaussData).parametricValue[ii].second * (*iterGaussData).parametricArea;
		}
		return true;
	}

	//overload
	template<unsigned int OG>
	bool GaussIntegralFaceO2Grid<OG>::getIntegral(
		unsigned int PG,
		real f[],
		real weight[],
		real cofJacobi[],
		real parametricArea,
		real& result){
		result = 0;
		for (int ii = 0; ii < static_cast<int>(PG); ++ii){
			//parametric space |s|=1
			//ע��ȷ�ϻ��ֹ�ʽ�Ƿ���ȷ����20200218
			result += f[ii] * weight[ii] * cofJacobi[ii] * parametricArea;
		}
		return true;
	}

	template<unsigned int OG>
	bool GaussIntegralFaceO2Grid<OG>::freeArray(void){

		//std::vector<std::pair<real, real>>().swap(parametricValue);
		//gaussPairVectorVector().swap(faceData);

		return true;
	}

	//-------------------------------------------------------------------------

	class cellGaussData : public cell{

	public:
		cellGaussData() :cell(){};
		~cellGaussData() {
			gaussPairVector().swap(gaussPairVector_);
			std::vector<int>().swap(PGTri);
			std::vector<int>().swap(PGQuad);
			std::vector<std::pair<point, real>>().swap(parametricValueTri);
			std::vector<std::pair<point, real>>().swap(parametricValueQuad);
		};
	public:
		gaussPairVector gaussPairVector_;
		std::vector<int> PGTri;
		std::vector<int> PGQuad;
		std::vector<std::pair<point, real>> parametricValueTri;
		std::vector<std::pair<point, real>> parametricValueQuad;
		real parametricVolume;
	};
	typedef std::vector<cellGaussData> cellGaussDataVector;


	template<unsigned int OG>
	class GaussIntegralCellO1Grid
	{
	public:

		GaussIntegralCellO1Grid() {};
		~GaussIntegralCellO1Grid() {};

	public:
		virtual bool allocateArray(
			parameter* parameter, 
			cellVector* cellVector, 
			cellGaussDataVector* cellGaussData);
		virtual bool initializeArray(
			cellVector* cellVector,
			cellGaussDataVector* cellGaussData);

		virtual bool getIntegral(
			real f[],			//Gauss���ϵĺ���ֵ
			int idxCell,									//��Ԫλ��ָ�꣬��ΪҪʹ��cellData
			cellGaussDataVector* cellGaussData,
			real& result);									//���ֽ��

		virtual bool getIntegral(
			unsigned int PG,
			real f[],
			real JacobiCof[],
			real parametricValue[],
			real parametricVolume,
			real& result);									//���ֽ��

		virtual bool freeArray(void);
		virtual point getGaussPointCell(
			int idxCell,								//��Ԫ����ָ�꣬��Ҫ�ж��������ε�Ԫ�����ı��ε�Ԫ
			std::pair<int, point> physicalPoint[],	//���ɵ�Ԫ�Ľڵ㣬�����ռ��е�����
			int idx,								//��Ԫλ��ָ�꣬��ΪҪʹ��cellData
			cellGaussDataVector* cellGaussData);									
		virtual real getJacobiCof(
			int idxCell,							//��Ԫ����ָ�꣬��Ҫ�ж��������ε�Ԫ�����ı��ε�Ԫ
			std::pair<int, point> physicalPoint[],	//���ɵ�Ԫ�Ľڵ㣬�����ռ��е�����
			int idx,
			cellGaussDataVector* cellGaussData);									//��Ԫλ��ָ�꣬��ΪҪʹ��cellData

		virtual bool syncData(
			cellVector* cellVector,
			cellGaussDataVector* cellGaussData);

		//debug 20200302
		//virtual bool exportTest(cellVector* cellVector, const std::string & filename);
	};

	template<unsigned int OG>
	bool GaussIntegralCellO1Grid<OG>::allocateArray(
		parameter* parameter,
		cellVector* cellVector,
		cellGaussDataVector* cellGaussData){

		std::cout << "allocating cell Gauss integral array..." << std::endl;

		cellGaussData->resize(parameter->cellNumber);

		cellVector::iterator iterCell;
		cellGaussDataVector::iterator iterGaussData;
		for (iterCell = cellVector->begin(); iterCell != cellVector->end(); ++iterCell){
			int ii = (*iterCell).index;
			iterGaussData = cellGaussData->begin() + ii - 1;
			(*iterGaussData).index = (*iterCell).index;
			(*iterGaussData).isBoundaryCell = (*iterCell).isBoundaryCell;
			(*iterGaussData).isCountByboundaryCellType = (*iterCell).isCountByboundaryCellType;
			(*iterGaussData).cellType_ = (*iterCell).cellType_;
			(*iterGaussData).boundaryCellType_ = (*iterCell).boundaryCellType_;
			(*iterGaussData).volume = (*iterCell).volume;
			(*iterGaussData).baryCenter = (*iterCell).baryCenter;
			(*iterGaussData).cellSideOff = (*iterCell).cellSideOff;
			(*iterGaussData).lengthReference = (*iterCell).lengthReference;
			(*iterGaussData).cellNodeNumber = (*iterCell).cellNodeNumber;
			(*iterGaussData).cellFaceNumber = (*iterCell).cellFaceNumber;
			(*iterGaussData).cellCellNumber = (*iterCell).cellCellNumber;

			switch (OG){
			case 1:
				(*iterGaussData).PGTri.resize(1 + 1);
				(*iterGaussData).PGTri[1] = 1;
				(*iterGaussData).PGQuad.resize(1 + 1);
				(*iterGaussData).PGQuad[1] = 1;
				break;
			case 2:
				(*iterGaussData).PGTri.resize(2 + 1);
				(*iterGaussData).PGTri[2] = 3;
				(*iterGaussData).PGQuad.resize(2 + 1);
				(*iterGaussData).PGQuad[2] = 4;
				break;
			case 3:
				(*iterGaussData).PGTri.resize(3 + 1);
				(*iterGaussData).PGTri[3] = 4;
				(*iterGaussData).PGQuad.resize(3 + 1);
				(*iterGaussData).PGQuad[3] = 9;
				break;
			default:
				std::cout << "	error: Currently only support PG from 1-3!" << std::endl;
				exit(1);
			}

			(*iterGaussData).parametricValueTri.resize((*iterGaussData).PGTri[OG]);
			(*iterGaussData).parametricValueQuad.resize((*iterGaussData).PGQuad[OG]);

			if ((*iterGaussData).cellType_ == Triangle){
				(*iterGaussData).gaussPairVector_.resize((*iterGaussData).PGTri[OG]);
				(*iterGaussData).cellNodeNumber = 3;
				(*iterGaussData).cellNode.resize(3 + 1);
			}
			else if ((*iterGaussData).cellType_ == Quadrilateral){
				(*iterGaussData).gaussPairVector_.resize((*iterGaussData).PGQuad[OG]);
				(*iterGaussData).cellNodeNumber = 4;
				(*iterGaussData).cellNode.resize(4 + 1);
			}
			for (int jj = 1; jj < (*iterGaussData).cellNodeNumber + 1; ++jj){
				(*iterGaussData).cellNode[jj].first = (*iterCell).cellNode[jj].first;
				(*iterGaussData).cellNode[jj].second = (*iterCell).cellNode[jj].second;
			}
			

			(*iterGaussData).cellFaceIndex.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceIndex[jj] = (*iterCell).cellFaceIndex[jj];
				(*iterGaussData).cellFaceIndex[jj] = (*iterCell).cellFaceIndex[jj];
			}

			(*iterGaussData).cellCellIndex.resize((*iterGaussData).cellCellNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellCellNumber + 1; ++jj){
				(*iterGaussData).cellCellIndex[jj] = (*iterCell).cellCellIndex[jj];
				(*iterGaussData).cellCellIndex[jj] = (*iterCell).cellCellIndex[jj];
			}

			(*iterGaussData).cellFaceSideOff.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceSideOff[jj] = (*iterCell).cellFaceSideOff[jj];
				(*iterGaussData).cellFaceSideOff[jj] = (*iterCell).cellFaceSideOff[jj];
			}

			(*iterGaussData).cellFaceWeight.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceWeight[jj] = (*iterCell).cellFaceWeight[jj];
				(*iterGaussData).cellFaceWeight[jj] = (*iterCell).cellFaceWeight[jj];
			}
		}

		std::cout << " ..cell Gauss integral array has been allocated." << std::endl;

		return true;
	}

	template<unsigned int OG>
	point GaussIntegralCellO1Grid<OG>::getGaussPointCell(
		int idxCell,
		std::pair<int, point> physicalPoint[],
		int idx,
		cellGaussDataVector* cellGaussData){
		//tri
		//quad
		point p;
		cellGaussDataVector::iterator iterGaussData;
		iterGaussData = cellGaussData->begin() + idxCell - 1;
		if ((*iterGaussData).cellType_ == Triangle){
			real kesi = (*iterGaussData).parametricValueTri[idx].first.x;
			real eta = (*iterGaussData).parametricValueTri[idx].first.y;
			real baseValue1 = 1.0 - kesi - eta;
			real baseValue2 = kesi;
			real baseValue3 = eta;
			p = baseValue1*physicalPoint[0].second
				+ baseValue2*physicalPoint[1].second
				+ baseValue3*physicalPoint[2].second;
		}
		else if ((*iterGaussData).cellType_ == Quadrilateral){
			real kesi = (*iterGaussData).parametricValueQuad[idx].first.x;
			real eta = (*iterGaussData).parametricValueQuad[idx].first.y;
			real baseValue1 = (-1.0 + kesi)*(-1.0 + eta);
			real baseValue2 = kesi*(1.0 - eta);
			real baseValue3 = kesi*eta;
			real baseValue4 = eta*(1.0 - kesi);
			p = baseValue1*physicalPoint[0].second
				+ baseValue2*physicalPoint[1].second
				+ baseValue3*physicalPoint[2].second
				+ baseValue4*physicalPoint[3].second;
		}
		return p;
	}

	template<unsigned int OG>
	real GaussIntegralCellO1Grid<OG>::getJacobiCof(
		int idxCell,
		std::pair<int, point> physicalPoint[],
		int idx,
		cellGaussDataVector* cellGaussData){
		//tri
		//quad
		real J11, J12, J21, J22;
		cellGaussDataVector::iterator iterGaussData;
		iterGaussData = cellGaussData->begin() + idxCell - 1;
		if ((*iterGaussData).cellType_ == Triangle){
			real kesi = (*iterGaussData).parametricValueTri[idx].first.x;
			real eta = (*iterGaussData).parametricValueTri[idx].first.y;
			real dBaseValue1dKesi = -1.0, dBaseValue1dEta = -1.0;
			real dBaseValue2dKesi = 1.0, dBaseValue2dEta = 0.0;
			real dBaseValue3dKesi = 0.0, dBaseValue3dEta = 1.0;

			//
			point origin = physicalPoint[0].second;
			
			//J11 = physicalPoint[0].second.x*dBaseValue1dKesi
			//	+ physicalPoint[1].second.x*dBaseValue2dKesi
			//	+ physicalPoint[2].second.x*dBaseValue3dKesi;
			//J12 = physicalPoint[0].second.x*dBaseValue1dEta
			//	+ physicalPoint[1].second.x*dBaseValue2dEta
			//	+ physicalPoint[2].second.x*dBaseValue3dEta;
			//J21 = physicalPoint[0].second.y*dBaseValue1dKesi
			//	+ physicalPoint[1].second.y*dBaseValue2dKesi
			//	+ physicalPoint[2].second.y*dBaseValue3dKesi;
			//J22 = physicalPoint[0].second.y*dBaseValue1dEta
			//	+ physicalPoint[1].second.y*dBaseValue2dEta
			//	+ physicalPoint[2].second.y*dBaseValue3dEta;

			
			J11 = (physicalPoint[0].second.x - origin.x)*dBaseValue1dKesi
				+ (physicalPoint[1].second.x - origin.x)*dBaseValue2dKesi
				+ (physicalPoint[2].second.x - origin.x)*dBaseValue3dKesi;

			J12 = (physicalPoint[0].second.x - origin.x)*dBaseValue1dEta
				+ (physicalPoint[1].second.x - origin.x)*dBaseValue2dEta
				+ (physicalPoint[2].second.x - origin.x)*dBaseValue3dEta;

			J21 = (physicalPoint[0].second.y - origin.y)*dBaseValue1dKesi
				+ (physicalPoint[1].second.y - origin.y)*dBaseValue2dKesi
				+ (physicalPoint[2].second.y - origin.y)*dBaseValue3dKesi;

			J22 = (physicalPoint[0].second.y - origin.y)*dBaseValue1dEta
				+ (physicalPoint[1].second.y - origin.y)*dBaseValue2dEta
				+ (physicalPoint[2].second.y - origin.y)*dBaseValue3dEta;
			
		}
		else if ((*iterGaussData).cellType_ == Quadrilateral){
			real kesi = (*iterGaussData).parametricValueQuad[idx].first.x;
			real eta = (*iterGaussData).parametricValueQuad[idx].first.y;
			real dBaseValue1dKesi = -1.0 + eta, dBaseValue1dEta = -1.0 + kesi;
			real dBaseValue2dKesi = 1.0 - eta, dBaseValue2dEta = -kesi;
			real dBaseValue3dKesi = eta, dBaseValue3dEta = kesi;
			real dBaseValue4dKesi = -eta, dBaseValue4dEta = 1.0 - kesi;
			point origin = physicalPoint[0].second;
			//J11 = physicalPoint[0].second.x*dBaseValue1dKesi
			//	+ physicalPoint[1].second.x*dBaseValue2dKesi
			//	+ physicalPoint[2].second.x*dBaseValue3dKesi
			//	+ physicalPoint[3].second.x*dBaseValue4dKesi;
			//J12 = physicalPoint[0].second.x*dBaseValue1dEta
			//	+ physicalPoint[1].second.x*dBaseValue2dEta
			//	+ physicalPoint[2].second.x*dBaseValue3dEta
			//	+ physicalPoint[3].second.x*dBaseValue4dEta;
			//J21 = physicalPoint[0].second.y*dBaseValue1dKesi
			//	+ physicalPoint[1].second.y*dBaseValue2dKesi
			//	+ physicalPoint[2].second.y*dBaseValue3dKesi
			//	+ physicalPoint[3].second.y*dBaseValue4dKesi;
			//J22 = physicalPoint[0].second.y*dBaseValue1dEta
			//	+ physicalPoint[1].second.y*dBaseValue2dEta
			//	+ physicalPoint[2].second.y*dBaseValue3dEta
			//	+ physicalPoint[3].second.y*dBaseValue4dEta;

			J11 = (physicalPoint[0].second.x - origin.x)*dBaseValue1dKesi
				+ (physicalPoint[1].second.x - origin.x)*dBaseValue2dKesi
				+ (physicalPoint[2].second.x - origin.x)*dBaseValue3dKesi
				+ (physicalPoint[3].second.x - origin.x)*dBaseValue4dKesi;
			J12 = (physicalPoint[0].second.x - origin.x)*dBaseValue1dEta
				+ (physicalPoint[1].second.x - origin.x)*dBaseValue2dEta
				+ (physicalPoint[2].second.x - origin.x)*dBaseValue3dEta
				+ (physicalPoint[3].second.x - origin.x)*dBaseValue4dEta;
			J21 = (physicalPoint[0].second.y - origin.y)*dBaseValue1dKesi
				+ (physicalPoint[1].second.y - origin.y)*dBaseValue2dKesi
				+ (physicalPoint[2].second.y - origin.y)*dBaseValue3dKesi
				+ (physicalPoint[3].second.y - origin.y)*dBaseValue4dKesi;
			J22 = (physicalPoint[0].second.y - origin.y)*dBaseValue1dEta
				+ (physicalPoint[1].second.y - origin.y)*dBaseValue2dEta
				+ (physicalPoint[2].second.y - origin.y)*dBaseValue3dEta
				+ (physicalPoint[3].second.y - origin.y)*dBaseValue4dEta;
		}
		return std::fabs(J11*J22 - J12*J21);//ע��ȷ��Jacobi��ʽ�Ƿ���ȷ����20200218
	}

	template<unsigned int OG>
	bool GaussIntegralCellO1Grid<OG>::initializeArray(
		cellVector* cellVector,
		cellGaussDataVector* cellGaussData){
		//PG = 3
		std::cout << "initializing cell Gauss integral array..." << std::endl;

		cellGaussDataVector::iterator iterGaussData;
		for (iterGaussData = cellGaussData->begin(); iterGaussData != cellGaussData->end(); ++iterGaussData){
			switch (OG){
			case 1:
				
				(*iterGaussData).parametricValueTri[0].first.x = 1.0 / 3.0; //kesi
				(*iterGaussData).parametricValueTri[0].first.y = 1.0 / 3.0; //eta
				(*iterGaussData).parametricValueTri[0].second = 1.0;

				(*iterGaussData).parametricValueQuad[0].first.x = 0.5; //kesi
				(*iterGaussData).parametricValueQuad[0].first.y = 0.5; //eta
				(*iterGaussData).parametricValueQuad[0].second = 1.0;
				break;
			case 2:
				(*iterGaussData).parametricValueTri[0].first.x = 2.0 / 3.0; //kesi
				(*iterGaussData).parametricValueTri[0].first.y = 1.0 / 6.0; //eta
				(*iterGaussData).parametricValueTri[0].second = 0.333333333333333;

				(*iterGaussData).parametricValueTri[1].first.x = 1.0 / 6.0; //kesi
				(*iterGaussData).parametricValueTri[1].first.y = 2.0 / 3.0; //eta
				(*iterGaussData).parametricValueTri[1].second = 0.333333333333333;

				(*iterGaussData).parametricValueTri[2].first.x = 1.0 / 6.0; //kesi
				(*iterGaussData).parametricValueTri[2].first.y = 1.0 / 6.0; //eta
				(*iterGaussData).parametricValueTri[2].second = 0.333333333333333;

				(*iterGaussData).parametricValueQuad[0].first.x = (-std::sqrt(3.0) + 3.0) / 6.0; //kesi
				(*iterGaussData).parametricValueQuad[0].first.y = (-std::sqrt(3.0) + 3.0) / 6.0; //eta
				(*iterGaussData).parametricValueQuad[0].second = 0.25;

				(*iterGaussData).parametricValueQuad[1].first.x = (std::sqrt(3.0) + 3.0) / 6.0; //kesi
				(*iterGaussData).parametricValueQuad[1].first.y = (-std::sqrt(3.0) + 3.0) / 6.0; //eta
				(*iterGaussData).parametricValueQuad[1].second = 0.25;

				(*iterGaussData).parametricValueQuad[2].first.x = (std::sqrt(3.0) + 3.0) / 6.0; //kesi
				(*iterGaussData).parametricValueQuad[2].first.y = (std::sqrt(3.0) + 3.0) / 6.0; //eta
				(*iterGaussData).parametricValueQuad[2].second = 0.25;

				(*iterGaussData).parametricValueQuad[3].first.x = (-std::sqrt(3.0) + 3.0) / 6.0; //kesi
				(*iterGaussData).parametricValueQuad[3].first.y = (std::sqrt(3.0) + 3.0) / 6.0; //eta
				(*iterGaussData).parametricValueQuad[3].second = 0.25;
				break;
			case 3:
				(*iterGaussData).parametricValueTri[0].first.x = 1.0 / 3.0; //kesi
				(*iterGaussData).parametricValueTri[0].first.y = 1.0 / 3.0; //eta
				(*iterGaussData).parametricValueTri[0].second = -27.0 / 48.0;

				(*iterGaussData).parametricValueTri[1].first.x = 0.2; //kesi
				(*iterGaussData).parametricValueTri[1].first.y = 0.6; //eta
				(*iterGaussData).parametricValueTri[1].second = 25.0 / 48.0;

				(*iterGaussData).parametricValueTri[2].first.x = 0.2; //kesi
				(*iterGaussData).parametricValueTri[2].first.y = 0.2; //eta
				(*iterGaussData).parametricValueTri[2].second = 25.0 / 48.0;

				(*iterGaussData).parametricValueTri[3].first.x = 0.6; //kesi
				(*iterGaussData).parametricValueTri[3].first.y = 0.2; //eta
				(*iterGaussData).parametricValueTri[3].second = 25.0 / 48.0;

				(*iterGaussData).parametricValueQuad[0].first.x = (-0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[0].first.y = (-0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[0].second = 0.555555555555555*0.555555555555555 / 4.0;

				(*iterGaussData).parametricValueQuad[1].first.x = (0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[1].first.y = (-0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[1].second = 0.555555555555555*0.555555555555555 / 4.0;

				(*iterGaussData).parametricValueQuad[2].first.x = (0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[2].first.y = (0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[2].second = 0.555555555555555*0.555555555555555 / 4.0;

				(*iterGaussData).parametricValueQuad[3].first.x = (-0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[3].first.y = (0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[3].second = 0.555555555555555*0.555555555555555 / 4.0;

				(*iterGaussData).parametricValueQuad[4].first.x = 0.5; //kesi
				(*iterGaussData).parametricValueQuad[4].first.y = (-0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[4].second = 0.555555555555555*0.888888888888888 / 4.0;

				(*iterGaussData).parametricValueQuad[5].first.x = (0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[5].first.y = 0.5; //eta
				(*iterGaussData).parametricValueQuad[5].second = 0.555555555555555*0.888888888888888 / 4.0;

				(*iterGaussData).parametricValueQuad[6].first.x = 0.5; //kesi
				(*iterGaussData).parametricValueQuad[6].first.y = (0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[6].second = 0.555555555555555*0.888888888888888 / 4.0;

				(*iterGaussData).parametricValueQuad[7].first.x = (-0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[7].first.y = 0.5; //eta
				(*iterGaussData).parametricValueQuad[7].second = 0.555555555555555*0.888888888888888 / 4.0;

				(*iterGaussData).parametricValueQuad[8].first.x = 0.5; //kesi
				(*iterGaussData).parametricValueQuad[8].first.y = 0.5; //eta
				(*iterGaussData).parametricValueQuad[8].second = 0.888888888888888*0.888888888888888 / 4.0;
				break;
			default:
				std::cout << "	error: Currently only support PG from 1-3!" << std::endl;
				exit(1);
			}

			unsigned int PG;
			unsigned int P;
			if ((*iterGaussData).cellType_ == Triangle){
				P = 3;
				PG = static_cast<int>((*iterGaussData).PGTri[OG]);

				(*iterGaussData).parametricVolume = 0.5;
			}
			else if ((*iterGaussData).cellType_ == Quadrilateral){
				P = 4;
				PG = static_cast<int>((*iterGaussData).PGQuad[OG]);

				(*iterGaussData).parametricVolume = 1.0;
			}

			std::pair<int, point>* node = new std::pair<int, point>[P];
			for (int kk = 0; kk < (*iterGaussData).cellNodeNumber; ++kk){
				node[kk].first = (*iterGaussData).cellNode[kk + 1].first;
				node[kk].second = (*iterGaussData).cellNode[kk + 1].second;
			}

			for (int jj = 0; jj < static_cast<int>(PG); ++jj){
				point physicalPoint;
				real cofJacobi;
				//if ((*iterGaussData).cellType_ == Triangle){
				//	physicalPoint = getGaussPointCell((*iterGaussData).cellType_, node, jj);
				//	cofJacobi = getJacobiCof((*iterGaussData).cellType_, node, jj);//confirm wether the parameter node[] is proper?
				//}
				//else if ((*iterGaussData).cellType_ == Quadrilateral){
				//	physicalPoint = getGaussPointCell((*iterGaussData).cellType_, node, jj);
				//	cofJacobi = getJacobiCof((*iterGaussData).cellType_, node, jj);
				//}
				physicalPoint = getGaussPointCell((*iterGaussData).index, node, jj, cellGaussData);
				cofJacobi = getJacobiCof((*iterGaussData).index, node, jj, cellGaussData);
				(*iterGaussData).gaussPairVector_[jj].p = physicalPoint;
				(*iterGaussData).gaussPairVector_[jj].JacobiCof = cofJacobi;

				//cellData[(*iter).index - 1][jj].p = physicalPoint;
				//cellData[(*iter).index - 1][jj].JacobiCof = cofJacobi;

			//	std::cout << "----------------------" << std::endl;
			//	std::cout << (*iterGaussData).gaussPairVector_[jj].p.x << "\t" << (*iterGaussData).gaussPairVector_[jj].p.y << std::endl;
			//	std::cout << (*iterGaussData).gaussPairVector_[jj].JacobiCof << std::endl;

			}
			delete[] node;
			node = NULL;
		}

		//system("pause");

		std::cout << " ..cell Gauss integral array has been initialized." << std::endl;
		return true;
	}

	template<unsigned int OG>
	bool GaussIntegralCellO1Grid<OG>::syncData(
		cellVector* cellVector,
		cellGaussDataVector* cellGaussData){
	
	
		cellVector::iterator iterCell;
		cellGaussDataVector::iterator iterGaussData;
		for (iterCell = cellVector->begin(); iterCell != cellVector->end(); ++iterCell){
			int ii = (*iterCell).index;
			iterGaussData = cellGaussData->begin() + ii - 1;
			(*iterGaussData).index = (*iterCell).index;
			(*iterGaussData).isBoundaryCell = (*iterCell).isBoundaryCell;
			(*iterGaussData).isCountByboundaryCellType = (*iterCell).isCountByboundaryCellType;
			(*iterGaussData).cellType_ = (*iterCell).cellType_;
			(*iterGaussData).volume = (*iterCell).volume;
			(*iterGaussData).baryCenter = (*iterCell).baryCenter;
			(*iterGaussData).cellSideOff = (*iterCell).cellSideOff;
			(*iterGaussData).lengthReference = (*iterCell).lengthReference;
			(*iterGaussData).cellNodeNumber = (*iterCell).cellNodeNumber;
			(*iterGaussData).cellFaceNumber = (*iterCell).cellFaceNumber;
			(*iterGaussData).cellCellNumber = (*iterCell).cellCellNumber;

			//std::cout << "========" << std::endl;
			//std::cout << (*iterGaussData).index << "\t" << (*iterCell).index << std::endl;
			//std::cout << (*iterGaussData).volume << "\t" << (*iterCell).volume << std::endl;
			//std::cout << (*iterGaussData).cellNodeNumber << "\t" << (*iterCell).cellNodeNumber << std::endl;

			switch (OG){
			case 1:
				(*iterGaussData).PGTri.resize(1 + 1);
				(*iterGaussData).PGTri[1] = 1;
				(*iterGaussData).PGQuad.resize(1 + 1);
				(*iterGaussData).PGQuad[1] = 1;
				break;
			case 2:
				(*iterGaussData).PGTri.resize(2 + 1);
				(*iterGaussData).PGTri[2] = 3;
				(*iterGaussData).PGQuad.resize(2 + 1);
				(*iterGaussData).PGQuad[2] = 4;
				break;
			case 3:
				(*iterGaussData).PGTri.resize(3 + 1);
				(*iterGaussData).PGTri[3] = 4;
				(*iterGaussData).PGQuad.resize(3 + 1);
				(*iterGaussData).PGQuad[3] = 9;
				break;
			default:
				std::cout << "	error: Currently only support PG from 1-3!" << std::endl;
				exit(1);
			}

			(*iterGaussData).parametricValueTri.resize((*iterGaussData).PGTri[OG]);
			(*iterGaussData).parametricValueQuad.resize((*iterGaussData).PGQuad[OG]);

			if ((*iterGaussData).cellType_ == Triangle){
				(*iterGaussData).gaussPairVector_.resize((*iterGaussData).PGTri[OG]);
				(*iterGaussData).cellNodeNumber = 3;
				(*iterGaussData).cellNode.resize(3 + 1);
			}
			else if ((*iterGaussData).cellType_ == Quadrilateral){
				(*iterGaussData).gaussPairVector_.resize((*iterGaussData).PGQuad[OG]);
				(*iterGaussData).cellNodeNumber = 4;
				(*iterGaussData).cellNode.resize(4 + 1);
			}
			for (int jj = 1; jj < (*iterGaussData).cellNodeNumber + 1; ++jj){
				(*iterGaussData).cellNode[jj].first = (*iterCell).cellNode[jj].first;
				(*iterGaussData).cellNode[jj].second = (*iterCell).cellNode[jj].second;
			}


			(*iterGaussData).cellFaceIndex.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceIndex[jj] = (*iterCell).cellFaceIndex[jj];
				(*iterGaussData).cellFaceIndex[jj] = (*iterCell).cellFaceIndex[jj];
			}

			(*iterGaussData).cellCellIndex.resize((*iterGaussData).cellCellNumber + 1);
			for (int jj = 0; jj < (*iterGaussData).cellCellNumber + 1; ++jj){
				(*iterGaussData).cellCellIndex[jj] = (*iterCell).cellCellIndex[jj];
				(*iterGaussData).cellCellIndex[jj] = (*iterCell).cellCellIndex[jj];
			}

			(*iterGaussData).cellFaceSideOff.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceSideOff[jj] = (*iterCell).cellFaceSideOff[jj];
				(*iterGaussData).cellFaceSideOff[jj] = (*iterCell).cellFaceSideOff[jj];
			}

			(*iterGaussData).cellFaceWeight.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceWeight[jj] = (*iterCell).cellFaceWeight[jj];
				(*iterGaussData).cellFaceWeight[jj] = (*iterCell).cellFaceWeight[jj];
			}
		}

		std::cout << " ..cell Gauss integral array has been allocated." << std::endl;
		return true;
	
	}

	template<unsigned int OG>
	bool GaussIntegralCellO1Grid<OG>::getIntegral(
		real f[],
		int idxCell,
		cellGaussDataVector* cellGaussData,
		real& result){
		result = 0;
		cellGaussDataVector::iterator iterGaussData;
		iterGaussData = cellGaussData->begin() + idxCell - 1;
		if ((*iterGaussData).cellType_ == Triangle){
			for (int ii = 0; ii < static_cast<int>((*iterGaussData).PGTri[OG]); ++ii){
				result += f[ii] * (*iterGaussData).gaussPairVector_[ii].JacobiCof * (*iterGaussData).parametricValueTri[ii].second * (*iterGaussData).parametricVolume;
				//ע��ȷ�ϻ��ֹ�ʽ�Ƿ���ȷ����20200218

				//std::cout << "--------------------------------------" << std::endl;
				//std::cout << f[ii] << std::endl;
				//std::cout << (*iterGaussData).gaussPairVector_[ii].JacobiCof << std::endl;
				//std::cout << (*iterGaussData).parametricValueTri[ii].second << std::endl;
				//std::cout << result << std::endl;
			}
		}
		else if ((*iterGaussData).cellType_ == Quadrilateral){
			
			for (int ii = 0; ii < static_cast<int>((*iterGaussData).PGQuad[OG]); ++ii){
				result += f[ii] * (*iterGaussData).gaussPairVector_[ii].JacobiCof * (*iterGaussData).parametricValueQuad[ii].second * (*iterGaussData).parametricVolume;

				//std::cout << "--------------------------------------" << std::endl;
				//std::cout << f[ii] << std::endl;
				//std::cout << (*iterGaussData).gaussPairVector_[ii].JacobiCof << std::endl;
				//std::cout << (*iterGaussData).parametricValueQuad[ii].second << std::endl;
				//std::cout << result << std::endl;

			}
		}
		//std::cout << "--------------------------------------" << std::endl;
		//std::cout << f[0] << std::endl;
		//std::cout << (*iterGaussData).gaussPairVector_[0].JacobiCof << std::endl;
		//std::cout << (*iterGaussData).parametricValueTri[0].second << std::endl;
		//std::cout << result << std::endl;

		return true;
	}

	template<unsigned int OG>
	bool GaussIntegralCellO1Grid<OG>::getIntegral(	
		unsigned int PG,
		real f[],
		real JacobiCof[],
		real parametricValue[],
		real parametricVolume,
		real& result){
		result = 0;

		for (int ii = 0; ii < static_cast<int>(PG); ++ii){
			result += f[ii] * JacobiCof[ii] * parametricValue[ii] * parametricVolume;
		}

		return true;
	}

	template<unsigned int OG>
	bool GaussIntegralCellO1Grid<OG>::freeArray(void){

		//pairVectorParameter().swap(parametricValueTri);
		//pairVectorParameter().swap(parametricValueQuad);
		//gaussPairVectorVector().swap(cellData);

		return true;
	}


	//O2 grid
	template<unsigned int OG>
	class GaussIntegralCellO2Grid : public GaussIntegralCellO1Grid<OG>
	{
	public:
		~GaussIntegralCellO2Grid() {};
		GaussIntegralCellO2Grid() {};

	public:

		bool allocateArray(
			parameter* parameter,
			cellVector* cellVector,
			cellGaussDataVector* cellGaussData) override;
		bool initializeArray(
			cellVector* cellVector,
			cellGaussDataVector* cellGaussData) override;
		bool getIntegral(real f[],			//Gauss���ϵĺ���ֵ
			int idx,									//��Ԫλ��ָ�꣬��ΪҪʹ��cellData
			cellGaussDataVector* cellGaussData,
			real& result) override;									//���ֽ��

		bool getIntegral(
			unsigned int PG,
			real f[],
			real JacobiCof[],
			real parametricValue[],
			real parametricVolume,
			real& result) override;									//���ֽ��

		bool freeArray(void) override;
		point getGaussPointCell(
			int idxCell,								//��Ԫ����ָ�꣬��Ҫ�ж��������ε�Ԫ�����ı��ε�Ԫ
			std::pair<int, point> physicalPoint[],	//���ɵ�Ԫ�Ľڵ㣬�����ռ��е�����
			int idx,								//��Ԫλ��ָ�꣬��ΪҪʹ��cellData
			cellGaussDataVector* cellGaussData) override;
		real getJacobiCof(
			int idxCell,							//��Ԫ����ָ�꣬��Ҫ�ж��������ε�Ԫ�����ı��ε�Ԫ
			std::pair<int, point> physicalPoint[],	//���ɵ�Ԫ�Ľڵ㣬�����ռ��е�����
			int idx,
			cellGaussDataVector* cellGaussData) override;

		bool syncData(
			cellVector* cellVector,
			cellGaussDataVector* cellGaussData) override;
	};

	template<unsigned int OG>
	bool GaussIntegralCellO2Grid<OG>::allocateArray(
		parameter* parameter,
		cellVector* cellVector,
		cellGaussDataVector* cellGaussData){

		std::cout << "allocating cell Gauss integral array..." << std::endl;

		cellGaussData->resize(parameter->cellNumber);

		cellVector::iterator iterCell;
		cellGaussDataVector::iterator iterGaussData;
		for (iterCell = cellVector->begin(); iterCell != cellVector->end(); ++iterCell){
			int ii = (*iterCell).index;
			iterGaussData = cellGaussData->begin() + ii - 1;
			(*iterGaussData).index = (*iterCell).index;
			(*iterGaussData).isBoundaryCell = (*iterCell).isBoundaryCell;
			(*iterGaussData).isCountByboundaryCellType = (*iterCell).isCountByboundaryCellType;
			(*iterGaussData).cellType_ = (*iterCell).cellType_;
			(*iterGaussData).boundaryCellType_ = (*iterCell).boundaryCellType_;
			(*iterGaussData).volume = (*iterCell).volume;
			(*iterGaussData).baryCenter = (*iterCell).baryCenter;
			(*iterGaussData).cellSideOff = (*iterCell).cellSideOff;
			(*iterGaussData).lengthReference = (*iterCell).lengthReference;
			(*iterGaussData).cellNodeNumber = (*iterCell).cellNodeNumber;
			(*iterGaussData).cellFaceNumber = (*iterCell).cellFaceNumber;
			(*iterGaussData).cellCellNumber = (*iterCell).cellCellNumber;

			switch (OG){
			case 1:
				(*iterGaussData).PGTri.resize(1 + 1);
				(*iterGaussData).PGTri[1] = 1;
				(*iterGaussData).PGQuad.resize(1 + 1);
				(*iterGaussData).PGQuad[1] = 1;
				break;
			case 2:
				(*iterGaussData).PGTri.resize(2 + 1);
				(*iterGaussData).PGTri[2] = 3;
				(*iterGaussData).PGQuad.resize(2 + 1);
				(*iterGaussData).PGQuad[2] = 4;
				break;
			case 3:
				(*iterGaussData).PGTri.resize(3 + 1);
				(*iterGaussData).PGTri[3] = 4;
				(*iterGaussData).PGQuad.resize(3 + 1);
				(*iterGaussData).PGQuad[3] = 9;
				break;
			default:
				std::cout << "	error: Currently only support PG from 1-3!" << std::endl;
				exit(1);
			}

			(*iterGaussData).parametricValueTri.resize((*iterGaussData).PGTri[OG]);
			(*iterGaussData).parametricValueQuad.resize((*iterGaussData).PGQuad[OG]);

			if ((*iterGaussData).cellType_ == Triangle){
				(*iterGaussData).gaussPairVector_.resize((*iterGaussData).PGTri[OG]);
				(*iterGaussData).cellNodeNumber = 3;
				(*iterGaussData).cellNode.resize(3 + 1);
			}
			else if ((*iterGaussData).cellType_ == Quadrilateral){
				(*iterGaussData).gaussPairVector_.resize((*iterGaussData).PGQuad[OG]);
				(*iterGaussData).cellNodeNumber = 4;
				(*iterGaussData).cellNode.resize(4 + 1);
			}
			for (int jj = 1; jj < (*iterGaussData).cellNodeNumber + 1; ++jj){
				(*iterGaussData).cellNode[jj].first = (*iterCell).cellNode[jj].first;
				(*iterGaussData).cellNode[jj].second = (*iterCell).cellNode[jj].second;
			}


			(*iterGaussData).cellFaceIndex.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceIndex[jj] = (*iterCell).cellFaceIndex[jj];
				(*iterGaussData).cellFaceIndex[jj] = (*iterCell).cellFaceIndex[jj];
			}

			(*iterGaussData).cellCellIndex.resize((*iterGaussData).cellCellNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellCellNumber + 1; ++jj){
				(*iterGaussData).cellCellIndex[jj] = (*iterCell).cellCellIndex[jj];
				(*iterGaussData).cellCellIndex[jj] = (*iterCell).cellCellIndex[jj];
			}

			(*iterGaussData).cellFaceSideOff.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceSideOff[jj] = (*iterCell).cellFaceSideOff[jj];
				(*iterGaussData).cellFaceSideOff[jj] = (*iterCell).cellFaceSideOff[jj];
			}

			(*iterGaussData).cellFaceWeight.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceWeight[jj] = (*iterCell).cellFaceWeight[jj];
				(*iterGaussData).cellFaceWeight[jj] = (*iterCell).cellFaceWeight[jj];
			}

			if ((*iterCell).o2PointNumber != 0){
				(*iterGaussData).o2PointNumber = (*iterCell).o2PointNumber;
				std::vector<std::pair<int, point> > o2Point;
				o2Point.resize((*iterCell).o2PointNumber + 1);
				for (int jj = 1; jj < (*iterCell).o2PointNumber + 1; ++jj){
					o2Point[jj].first = (*iterCell).cellNodeO2Point[jj].first;
					o2Point[jj].second = (*iterCell).cellNodeO2Point[jj].second;
				}		
				(*iterGaussData).cellNodeNumber += (*iterCell).o2PointNumber;
				for (int jj = 1; jj < (*iterCell).o2PointNumber + 1; ++jj){
					(*iterGaussData).cellNode.push_back(o2Point[jj]);
				//	std::cout << (*iterGaussData).cellNode[4 + jj].second.x << "\t" << (*iterGaussData).cellNode[4 + jj].second.y << std::endl;
				}
			}

		//	std::cout << "----" << std::endl;
		//	std::cout << (*iterGaussData).index << "\t" << (*iterGaussData).cellNodeNumber << std::endl;
		//	for (int jj = 1; jj < (*iterGaussData).cellNodeNumber + 1; ++jj){
		//		std::cout << (*iterGaussData).cellNode[jj].second.x << "\t" << (*iterGaussData).cellNode[jj].second.y << std::endl;
		//	}
		}	
		std::cout << " ..cell Gauss integral array has been allocated." << std::endl;
		//system("pause");
		return true;
	}

	template<unsigned int OG>
	point GaussIntegralCellO2Grid<OG>::getGaussPointCell(
		int idxCell,
		std::pair<int, point> physicalPoint[],
		int idx,
		cellGaussDataVector* cellGaussData){
		//tri
		//quad
		point p;
		cellGaussDataVector::iterator iterGaussData;
		iterGaussData = cellGaussData->begin() + idxCell - 1;
		if ((*iterGaussData).cellType_ == Triangle){
			real kesi = (*iterGaussData).parametricValueTri[idx].first.x;
			real eta = (*iterGaussData).parametricValueTri[idx].first.y;
			real baseValue1 = (-1.0 + kesi + eta)*(-1.0 + 2.0*kesi + 2.0*eta);
			real baseValue2 = kesi*(-1.0 + 2.0*kesi);
			real baseValue3 = eta*(-1.0 + 2.0*eta);
			real baseValue4 = -4.0*(-1.0 + kesi + eta)*kesi;
			real baseValue5 = 4.0*kesi*eta;
			real baseValue6 = -4.0*(-1.0 + kesi + eta)*eta;

			p = baseValue1*physicalPoint[0].second
				+ baseValue2*physicalPoint[1].second
				+ baseValue3*physicalPoint[2].second
				+ baseValue4*physicalPoint[3].second
				+ baseValue5*physicalPoint[4].second
				+ baseValue6*physicalPoint[5].second;
		}
		else if ((*iterGaussData).cellType_ == Quadrilateral){
			real kesi = (*iterGaussData).parametricValueQuad[idx].first.x;
			real eta = (*iterGaussData).parametricValueQuad[idx].first.y;
			real baseValue1 = (-1.0 + eta)*(-1.0 + 2.0*eta)*(-1.0 + kesi)*(-1.0 + 2.0*kesi);
			real baseValue2 = (-1.0 + eta)*(-1.0 + 2.0*eta)*kesi*(-1.0 + 2.0*kesi);
			real baseValue3 = eta*(-1.0 + 2.0*eta)*kesi*(-1.0 + 2.0*kesi);
			real baseValue4 = eta*(-1.0 + 2.0*eta)*(-1.0 + kesi)*(-1.0 + 2.0*kesi);
			real baseValue5 = -4.0*(-1.0 + eta)*(-1.0 + 2.0*eta)*(-1.0 + kesi)*kesi;
			real baseValue6 = -4.0*kesi*eta*(-1.0 + eta)*(-1.0 + 2.0*kesi);
			real baseValue7 = -4.0*kesi*eta*(-1.0 + 2.0*eta)*(-1.0 + kesi);
			real baseValue8 = -4.0*(-1.0 + eta)*eta*(-1.0 + kesi)*(-1.0 + 2.0*kesi);
			real baseValue9 = 16.0*(-1.0 + eta)*eta*(-1.0 + kesi)*kesi;

			p = baseValue1*physicalPoint[0].second
				+ baseValue2*physicalPoint[1].second
				+ baseValue3*physicalPoint[2].second
				+ baseValue4*physicalPoint[3].second
				+ baseValue5*physicalPoint[4].second
				+ baseValue6*physicalPoint[5].second
				+ baseValue7*physicalPoint[6].second
				+ baseValue8*physicalPoint[7].second
				+ baseValue9*physicalPoint[8].second;
		}
		return p;
	}

	template<unsigned int OG>
	real GaussIntegralCellO2Grid<OG>::getJacobiCof(
		int idxCell,
		std::pair<int, point> physicalPoint[],
		int idx,
		cellGaussDataVector* cellGaussData){
		//tri
		//quad
		real J11, J12, J21, J22;
		cellGaussDataVector::iterator iterGaussData;
		iterGaussData = cellGaussData->begin() + idxCell - 1;
		if ((*iterGaussData).cellType_ == Triangle){
			real kesi = (*iterGaussData).parametricValueTri[idx].first.x;
			real eta = (*iterGaussData).parametricValueTri[idx].first.y;
			real dBaseValue1dKesi = -3.0 + 4.0*kesi + 4.0*eta,		dBaseValue1dEta = -3.0 + 4.0*kesi + 4.0*eta;
			real dBaseValue2dKesi = -1.0 + 4.0*kesi,				dBaseValue2dEta = 0.0;
			real dBaseValue3dKesi = 0.0,							dBaseValue3dEta = -1.0 + 4.0*eta;
			real dBaseValue4dKesi = -4.0*(-1.0 + 2.0*kesi + eta),	dBaseValue4dEta = -4.0*kesi;
			real dBaseValue5dKesi = 4.0*eta,						dBaseValue5dEta = 4.0*kesi;
			real dBaseValue6dKesi = -4.0*eta,						dBaseValue6dEta = -4.0*(-1.0 + kesi + 2.0*eta);

			point origin = physicalPoint[0].second;

			J11 = (physicalPoint[0].second.x - origin.x)*dBaseValue1dKesi
				+ (physicalPoint[1].second.x - origin.x)*dBaseValue2dKesi
				+ (physicalPoint[2].second.x - origin.x)*dBaseValue3dKesi
				+ (physicalPoint[3].second.x - origin.x)*dBaseValue4dKesi
				+ (physicalPoint[4].second.x - origin.x)*dBaseValue5dKesi
				+ (physicalPoint[5].second.x - origin.x)*dBaseValue6dKesi;

			J12 = (physicalPoint[0].second.x - origin.x)*dBaseValue1dEta
				+ (physicalPoint[1].second.x - origin.x)*dBaseValue2dEta
				+ (physicalPoint[2].second.x - origin.x)*dBaseValue3dEta
				+ (physicalPoint[3].second.x - origin.x)*dBaseValue4dEta
				+ (physicalPoint[4].second.x - origin.x)*dBaseValue5dEta
				+ (physicalPoint[5].second.x - origin.x)*dBaseValue6dEta;

			J21 = (physicalPoint[0].second.y - origin.y)*dBaseValue1dKesi
				+ (physicalPoint[1].second.y - origin.y)*dBaseValue2dKesi
				+ (physicalPoint[2].second.y - origin.y)*dBaseValue3dKesi
				+ (physicalPoint[3].second.y - origin.y)*dBaseValue4dKesi
				+ (physicalPoint[4].second.y - origin.y)*dBaseValue5dKesi
				+ (physicalPoint[5].second.y - origin.y)*dBaseValue6dKesi;

			J22 = (physicalPoint[0].second.y - origin.y)*dBaseValue1dEta
				+ (physicalPoint[1].second.y - origin.y)*dBaseValue2dEta
				+ (physicalPoint[2].second.y - origin.y)*dBaseValue3dEta
				+ (physicalPoint[3].second.y - origin.y)*dBaseValue4dEta
				+ (physicalPoint[4].second.y - origin.y)*dBaseValue5dEta
				+ (physicalPoint[5].second.y - origin.y)*dBaseValue6dEta;

		}
		else if ((*iterGaussData).cellType_ == Quadrilateral){
			real kesi = (*iterGaussData).parametricValueQuad[idx].first.x;
			real eta = (*iterGaussData).parametricValueQuad[idx].first.y;
			real dBaseValue1dKesi = (-3.0 + 4.0*kesi)*(1.0 - 3.0*eta + 2.0*eta*eta),		dBaseValue1dEta = (-3.0 + 4.0*eta)*(1.0 - 3.0*kesi + 2.0*kesi*kesi);
			real dBaseValue2dKesi = (-1.0 + 4.0*kesi)*(1.0 - 3.0*eta + 2.0*eta*eta),		dBaseValue2dEta = kesi*(-1.0 + 2.0*kesi)*(-3.0 + 4.0*eta);
			real dBaseValue3dKesi = (-1.0 + 4.0*kesi)*eta*(-1.0 + 2.0*eta),					dBaseValue3dEta = kesi*(-1.0 + 2.0*kesi)*(-1.0 + 4.0*eta);
			real dBaseValue4dKesi = (-3.0 + 4.0*kesi)*eta*(-1.0 + 2.0*eta),					dBaseValue4dEta = (1.0 - 3.0*kesi + 2.0*kesi*kesi)*(-1.0 + 4.0*eta);
			real dBaseValue5dKesi = -4.0*(-1.0 + 2.0*kesi)*(1.0 - 3.0*eta + 2.0*eta*eta),	dBaseValue5dEta = -4.0*(-1.0 + kesi)*kesi*(-3.0 + 4.0*eta);
			real dBaseValue6dKesi = -4.0*(-1.0 + 4.0*kesi)*eta*(-1.0 + eta),				dBaseValue6dEta = -4.0*(-1.0 + 2.0*kesi)*kesi*(-1.0 + 2.0*eta);
			real dBaseValue7dKesi = -4.0*(-1.0 + 2.0*kesi)*eta*(-1.0 + 2.0*eta),			dBaseValue7dEta = -4.0*(-1.0 + 4.0*eta)*kesi*(-1.0 + kesi);
			real dBaseValue8dKesi = -4.0*(-3.0 + 4.0*kesi)*eta*(-1.0 + eta),				dBaseValue8dEta = -4.0*(1.0 - 3.0*kesi + 2.0*kesi*kesi)*(-1.0 + 2.0*eta);
			real dBaseValue9dKesi = 16.0*(-1.0 + 2.0*kesi)*eta*(-1.0 + eta),				dBaseValue9dEta = 16.0*(-1.0 + 2.0*eta)*kesi*(-1.0 + kesi);

			point origin = physicalPoint[0].second;

			J11 = (physicalPoint[0].second.x - origin.x)*dBaseValue1dKesi
				+ (physicalPoint[1].second.x - origin.x)*dBaseValue2dKesi
				+ (physicalPoint[2].second.x - origin.x)*dBaseValue3dKesi
				+ (physicalPoint[3].second.x - origin.x)*dBaseValue4dKesi
				+ (physicalPoint[4].second.x - origin.x)*dBaseValue5dKesi
				+ (physicalPoint[5].second.x - origin.x)*dBaseValue6dKesi
				+ (physicalPoint[6].second.x - origin.x)*dBaseValue7dKesi
				+ (physicalPoint[7].second.x - origin.x)*dBaseValue8dKesi
				+ (physicalPoint[8].second.x - origin.x)*dBaseValue9dKesi;

			J12 = (physicalPoint[0].second.x - origin.x)*dBaseValue1dEta
				+ (physicalPoint[1].second.x - origin.x)*dBaseValue2dEta
				+ (physicalPoint[2].second.x - origin.x)*dBaseValue3dEta
				+ (physicalPoint[3].second.x - origin.x)*dBaseValue4dEta
				+ (physicalPoint[4].second.x - origin.x)*dBaseValue5dEta
				+ (physicalPoint[5].second.x - origin.x)*dBaseValue6dEta
				+ (physicalPoint[6].second.x - origin.x)*dBaseValue7dEta
				+ (physicalPoint[7].second.x - origin.x)*dBaseValue8dEta
				+ (physicalPoint[8].second.x - origin.x)*dBaseValue9dEta;

			J21 = (physicalPoint[0].second.y - origin.y)*dBaseValue1dKesi
				+ (physicalPoint[1].second.y - origin.y)*dBaseValue2dKesi
				+ (physicalPoint[2].second.y - origin.y)*dBaseValue3dKesi
				+ (physicalPoint[3].second.y - origin.y)*dBaseValue4dKesi
				+ (physicalPoint[4].second.y - origin.y)*dBaseValue5dKesi
				+ (physicalPoint[5].second.y - origin.y)*dBaseValue6dKesi
				+ (physicalPoint[6].second.y - origin.y)*dBaseValue7dKesi
				+ (physicalPoint[7].second.y - origin.y)*dBaseValue8dKesi
				+ (physicalPoint[8].second.y - origin.y)*dBaseValue9dKesi;

			J22 = (physicalPoint[0].second.y - origin.y)*dBaseValue1dEta
				+ (physicalPoint[1].second.y - origin.y)*dBaseValue2dEta
				+ (physicalPoint[2].second.y - origin.y)*dBaseValue3dEta
				+ (physicalPoint[3].second.y - origin.y)*dBaseValue4dEta
				+ (physicalPoint[4].second.y - origin.y)*dBaseValue5dEta
				+ (physicalPoint[5].second.y - origin.y)*dBaseValue6dEta
				+ (physicalPoint[6].second.y - origin.y)*dBaseValue7dEta
				+ (physicalPoint[7].second.y - origin.y)*dBaseValue8dEta
				+ (physicalPoint[8].second.y - origin.y)*dBaseValue9dEta;
		}
		return std::fabs(J11*J22 - J12*J21);//ע��ȷ��Jacobi��ʽ�Ƿ���ȷ����20200218
	}

	template<unsigned int OG>
	bool GaussIntegralCellO2Grid<OG>::initializeArray(
		cellVector* cellVector,
		cellGaussDataVector* cellGaussData){
		//PG = 3
		std::cout << "initializing cell Gauss integral array..." << std::endl;

		cellGaussDataVector::iterator iterGaussData;
		for (iterGaussData = cellGaussData->begin(); iterGaussData != cellGaussData->end(); ++iterGaussData){
			switch (OG){
			case 1:

				(*iterGaussData).parametricValueTri[0].first.x = 1.0 / 3.0; //kesi
				(*iterGaussData).parametricValueTri[0].first.y = 1.0 / 3.0; //eta
				(*iterGaussData).parametricValueTri[0].second = 1.0;

				(*iterGaussData).parametricValueQuad[0].first.x = 0.5; //kesi
				(*iterGaussData).parametricValueQuad[0].first.y = 0.5; //eta
				(*iterGaussData).parametricValueQuad[0].second = 1.0;
				break;
			case 2:
				(*iterGaussData).parametricValueTri[0].first.x = 2.0 / 3.0; //kesi
				(*iterGaussData).parametricValueTri[0].first.y = 1.0 / 6.0; //eta
				(*iterGaussData).parametricValueTri[0].second = 0.333333333333333;

				(*iterGaussData).parametricValueTri[1].first.x = 1.0 / 6.0; //kesi
				(*iterGaussData).parametricValueTri[1].first.y = 2.0 / 3.0; //eta
				(*iterGaussData).parametricValueTri[1].second = 0.333333333333333;

				(*iterGaussData).parametricValueTri[2].first.x = 1.0 / 6.0; //kesi
				(*iterGaussData).parametricValueTri[2].first.y = 1.0 / 6.0; //eta
				(*iterGaussData).parametricValueTri[2].second = 0.333333333333333;

				(*iterGaussData).parametricValueQuad[0].first.x = (-std::sqrt(3.0) + 3.0) / 6.0; //kesi
				(*iterGaussData).parametricValueQuad[0].first.y = (-std::sqrt(3.0) + 3.0) / 6.0; //eta
				(*iterGaussData).parametricValueQuad[0].second = 0.25;

				(*iterGaussData).parametricValueQuad[1].first.x = (std::sqrt(3.0) + 3.0) / 6.0; //kesi
				(*iterGaussData).parametricValueQuad[1].first.y = (-std::sqrt(3.0) + 3.0) / 6.0; //eta
				(*iterGaussData).parametricValueQuad[1].second = 0.25;

				(*iterGaussData).parametricValueQuad[2].first.x = (std::sqrt(3.0) + 3.0) / 6.0; //kesi
				(*iterGaussData).parametricValueQuad[2].first.y = (std::sqrt(3.0) + 3.0) / 6.0; //eta
				(*iterGaussData).parametricValueQuad[2].second = 0.25;

				(*iterGaussData).parametricValueQuad[3].first.x = (-std::sqrt(3.0) + 3.0) / 6.0; //kesi
				(*iterGaussData).parametricValueQuad[3].first.y = (std::sqrt(3.0) + 3.0) / 6.0; //eta
				(*iterGaussData).parametricValueQuad[3].second = 0.25;
				break;
			case 3:
				(*iterGaussData).parametricValueTri[0].first.x = 1.0 / 3.0; //kesi
				(*iterGaussData).parametricValueTri[0].first.y = 1.0 / 3.0; //eta
				(*iterGaussData).parametricValueTri[0].second = -27.0 / 48.0;

				(*iterGaussData).parametricValueTri[1].first.x = 0.2; //kesi
				(*iterGaussData).parametricValueTri[1].first.y = 0.6; //eta
				(*iterGaussData).parametricValueTri[1].second = 25.0 / 48.0;

				(*iterGaussData).parametricValueTri[2].first.x = 0.2; //kesi
				(*iterGaussData).parametricValueTri[2].first.y = 0.2; //eta
				(*iterGaussData).parametricValueTri[2].second = 25.0 / 48.0;

				(*iterGaussData).parametricValueTri[3].first.x = 0.6; //kesi
				(*iterGaussData).parametricValueTri[3].first.y = 0.2; //eta
				(*iterGaussData).parametricValueTri[3].second = 25.0 / 48.0;

				(*iterGaussData).parametricValueQuad[0].first.x = (-0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[0].first.y = (-0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[0].second = 0.555555555555555*0.555555555555555 / 4.0;

				(*iterGaussData).parametricValueQuad[1].first.x = (0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[1].first.y = (-0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[1].second = 0.555555555555555*0.555555555555555 / 4.0;

				(*iterGaussData).parametricValueQuad[2].first.x = (0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[2].first.y = (0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[2].second = 0.555555555555555*0.555555555555555 / 4.0;

				(*iterGaussData).parametricValueQuad[3].first.x = (-0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[3].first.y = (0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[3].second = 0.555555555555555*0.555555555555555 / 4.0;

				(*iterGaussData).parametricValueQuad[4].first.x = 0.5; //kesi
				(*iterGaussData).parametricValueQuad[4].first.y = (-0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[4].second = 0.555555555555555*0.888888888888888 / 4.0;

				(*iterGaussData).parametricValueQuad[5].first.x = (0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[5].first.y = 0.5; //eta
				(*iterGaussData).parametricValueQuad[5].second = 0.555555555555555*0.888888888888888 / 4.0;

				(*iterGaussData).parametricValueQuad[6].first.x = 0.5; //kesi
				(*iterGaussData).parametricValueQuad[6].first.y = (0.774596669241 + 1.0) / 2.0; //eta
				(*iterGaussData).parametricValueQuad[6].second = 0.555555555555555*0.888888888888888 / 4.0;

				(*iterGaussData).parametricValueQuad[7].first.x = (-0.774596669241 + 1.0) / 2.0; //kesi
				(*iterGaussData).parametricValueQuad[7].first.y = 0.5; //eta
				(*iterGaussData).parametricValueQuad[7].second = 0.555555555555555*0.888888888888888 / 4.0;

				(*iterGaussData).parametricValueQuad[8].first.x = 0.5; //kesi
				(*iterGaussData).parametricValueQuad[8].first.y = 0.5; //eta
				(*iterGaussData).parametricValueQuad[8].second = 0.888888888888888*0.888888888888888 / 4.0;
				break;
			default:
				std::cout << "	error: Currently only support PG from 1-3!" << std::endl;
				exit(1);
			}

			unsigned int PG;
			unsigned int P;
			if ((*iterGaussData).cellType_ == Triangle){
				P = 3 + (*iterGaussData).o2PointNumber;
				PG = static_cast<int>((*iterGaussData).PGTri[OG]);

				(*iterGaussData).parametricVolume = 0.5;
			}
			else if ((*iterGaussData).cellType_ == Quadrilateral){
				P = 4 + (*iterGaussData).o2PointNumber;
				PG = static_cast<int>((*iterGaussData).PGQuad[OG]);

				(*iterGaussData).parametricVolume = 1.0;
			}

			std::pair<int, point>* node = new std::pair<int, point>[P];
			for (int kk = 0; kk < (*iterGaussData).cellNodeNumber; ++kk){
				node[kk].first = (*iterGaussData).cellNode[kk + 1].first;
				node[kk].second = (*iterGaussData).cellNode[kk + 1].second;
			}

		//	std::cout << "========================" << std::endl;
		//	std::cout << (*iterGaussData).index << std::endl;
		//	std::cout << P << "\t" <<static_cast<int>(PG) << std::endl;
			for (int jj = 0; jj < static_cast<int>(PG); ++jj){
				point physicalPoint;
				real cofJacobi;

				physicalPoint = getGaussPointCell((*iterGaussData).index, node, jj, cellGaussData);
				cofJacobi = getJacobiCof((*iterGaussData).index, node, jj, cellGaussData);
				(*iterGaussData).gaussPairVector_[jj].p = physicalPoint;
				(*iterGaussData).gaussPairVector_[jj].JacobiCof = cofJacobi;

			//	std::cout << "----------------------" << std::endl;
			//	std::cout.precision(20);
			//	std::cout << (*iterGaussData).gaussPairVector_[jj].p.x << "\t" << (*iterGaussData).gaussPairVector_[jj].p.y << std::endl;
			//	std::cout << (*iterGaussData).gaussPairVector_[jj].JacobiCof << std::endl;
			}
			delete[] node;
			node = NULL;
		}
		std::cout << " ..cell Gauss integral array has been initialized." << std::endl;
		return true;
	}

	template<unsigned int OG>
	bool GaussIntegralCellO2Grid<OG>::syncData(
		cellVector* cellVector,
		cellGaussDataVector* cellGaussData){

		cellVector::iterator iterCell;
		cellGaussDataVector::iterator iterGaussData;
		for (iterCell = cellVector->begin(); iterCell != cellVector->end(); ++iterCell){
			int ii = (*iterCell).index;
			iterGaussData = cellGaussData->begin() + ii - 1;
			(*iterGaussData).index = (*iterCell).index;
			(*iterGaussData).isBoundaryCell = (*iterCell).isBoundaryCell;
			(*iterGaussData).isCountByboundaryCellType = (*iterCell).isCountByboundaryCellType;
			(*iterGaussData).cellType_ = (*iterCell).cellType_;
			(*iterGaussData).volume = (*iterCell).volume;
			(*iterGaussData).baryCenter = (*iterCell).baryCenter;
			(*iterGaussData).cellSideOff = (*iterCell).cellSideOff;
			(*iterGaussData).lengthReference = (*iterCell).lengthReference;
		//	(*iterGaussData).cellNodeNumber = (*iterCell).cellNodeNumber;
			(*iterGaussData).cellFaceNumber = (*iterCell).cellFaceNumber;
			(*iterGaussData).cellCellNumber = (*iterCell).cellCellNumber;

			//std::cout << "========" << std::endl;
			//std::cout << (*iterGaussData).index << "\t" << (*iterCell).index << std::endl;
			//std::cout << (*iterGaussData).volume << "\t" << (*iterCell).volume << std::endl;
			//std::cout << (*iterGaussData).cellNodeNumber << "\t" << (*iterCell).cellNodeNumber << std::endl;

			switch (OG){
			case 1:
				(*iterGaussData).PGTri.resize(1 + 1);
				(*iterGaussData).PGTri[1] = 1;
				(*iterGaussData).PGQuad.resize(1 + 1);
				(*iterGaussData).PGQuad[1] = 1;
				break;
			case 2:
				(*iterGaussData).PGTri.resize(2 + 1);
				(*iterGaussData).PGTri[2] = 3;
				(*iterGaussData).PGQuad.resize(2 + 1);
				(*iterGaussData).PGQuad[2] = 4;
				break;
			case 3:
				(*iterGaussData).PGTri.resize(3 + 1);
				(*iterGaussData).PGTri[3] = 4;
				(*iterGaussData).PGQuad.resize(3 + 1);
				(*iterGaussData).PGQuad[3] = 9;
				break;
			default:
				std::cout << "	error: Currently only support PG from 1-3!" << std::endl;
				exit(1);
			}

			(*iterGaussData).parametricValueTri.resize((*iterGaussData).PGTri[OG]);
			(*iterGaussData).parametricValueQuad.resize((*iterGaussData).PGQuad[OG]);

			if ((*iterGaussData).cellType_ == Triangle){
				(*iterGaussData).gaussPairVector_.resize((*iterGaussData).PGTri[OG]);
			//	(*iterGaussData).cellNodeNumber = (*iterCell).cellNodeNumber;
			//	(*iterGaussData).cellNode.resize((*iterCell).cellNodeNumber + 1);
			}
			else if ((*iterGaussData).cellType_ == Quadrilateral){
				(*iterGaussData).gaussPairVector_.resize((*iterGaussData).PGQuad[OG]);
			//	(*iterGaussData).cellNodeNumber = (*iterCell).cellNodeNumber;
			//	(*iterGaussData).cellNode.resize((*iterCell).cellNodeNumber + 1);
			}
			//for (int jj = 1; jj < (*iterGaussData).cellNodeNumber + 1; ++jj){
			//	(*iterGaussData).cellNode[jj].first = (*iterCell).cellNode[jj].first;
			//	(*iterGaussData).cellNode[jj].second = (*iterCell).cellNode[jj].second;
			//}


			(*iterGaussData).cellFaceIndex.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceIndex[jj] = (*iterCell).cellFaceIndex[jj];
				(*iterGaussData).cellFaceIndex[jj] = (*iterCell).cellFaceIndex[jj];
			}

			(*iterGaussData).cellCellIndex.resize((*iterGaussData).cellCellNumber + 1);
			for (int jj = 0; jj < (*iterGaussData).cellCellNumber + 1; ++jj){
				(*iterGaussData).cellCellIndex[jj] = (*iterCell).cellCellIndex[jj];
				(*iterGaussData).cellCellIndex[jj] = (*iterCell).cellCellIndex[jj];
			}

			(*iterGaussData).cellFaceSideOff.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceSideOff[jj] = (*iterCell).cellFaceSideOff[jj];
				(*iterGaussData).cellFaceSideOff[jj] = (*iterCell).cellFaceSideOff[jj];
			}

			(*iterGaussData).cellFaceWeight.resize((*iterGaussData).cellFaceNumber + 1);
			for (int jj = 1; jj < (*iterGaussData).cellFaceNumber + 1; ++jj){
				(*iterGaussData).cellFaceWeight[jj] = (*iterCell).cellFaceWeight[jj];
				(*iterGaussData).cellFaceWeight[jj] = (*iterCell).cellFaceWeight[jj];
			}
		}

		std::cout << " ..cell Gauss integral array has been allocated." << std::endl;
		return true;

	}

	template<unsigned int OG>
	bool GaussIntegralCellO2Grid<OG>::getIntegral(
		real f[],
		int idxCell,
		cellGaussDataVector* cellGaussData,
		real& result){
		result = 0;
		cellGaussDataVector::iterator iterGaussData;
		iterGaussData = cellGaussData->begin() + idxCell - 1;
		if ((*iterGaussData).cellType_ == Triangle){
			for (int ii = 0; ii < static_cast<int>((*iterGaussData).PGTri[OG]); ++ii){
				result += f[ii] * (*iterGaussData).gaussPairVector_[ii].JacobiCof * (*iterGaussData).parametricValueTri[ii].second * (*iterGaussData).parametricVolume;
			}
		}
		else if ((*iterGaussData).cellType_ == Quadrilateral){
			for (int ii = 0; ii < static_cast<int>((*iterGaussData).PGQuad[OG]); ++ii){
				result += f[ii] * (*iterGaussData).gaussPairVector_[ii].JacobiCof * (*iterGaussData).parametricValueQuad[ii].second * (*iterGaussData).parametricVolume;
			}
		}
		return true;
	}

	template<unsigned int OG>
	bool GaussIntegralCellO2Grid<OG>::getIntegral(
		unsigned int PG,
		real f[],
		real JacobiCof[],
		real parametricValue[],
		real parametricVolume,
		real& result){
		result = 0;

		for (int ii = 0; ii < static_cast<int>(PG); ++ii){
			result += f[ii] * JacobiCof[ii] * parametricValue[ii] * parametricVolume;
		}

		return true;
	}

	template<unsigned int OG>
	bool GaussIntegralCellO2Grid<OG>::freeArray(void){

		//pairVectorParameter().swap(parametricValueTri);
		//pairVectorParameter().swap(parametricValueQuad);
		//gaussPairVectorVector().swap(cellData);

		return true;
	}



}
#endif