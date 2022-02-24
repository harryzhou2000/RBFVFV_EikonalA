#include "Parameter.h"

namespace ScalarCfv
{
	const std::string BASE_PATH = "./Cases/test";
	extern std::string fineOut_BndCheck(BASE_PATH + "/bndCheck.txt");
	extern std::string fileIn_Mesh(BASE_PATH + "/grid.in");
	//		extern std::string fileIn_Mesh("D:/calculations/new_test/grid_o2.in");
	extern std::string fileOut_Mesh(BASE_PATH + "/grid.plt");
	//extern std::string fileIn_Parameter("D:/calculations/new_test/inputParameters.txt");
	extern std::string fileIn_Parameter(BASE_PATH + "/inputParameters_notes.txt");
	extern std::string fileOut_debug_sln(BASE_PATH + "/sln.plt");

	extern std::string fileIn_BackUp(BASE_PATH + "/sln.sav");
	extern std::string fileOut_BackUp(BASE_PATH + "/sln.sav");
	extern std::string fileOut_Sln(BASE_PATH + "/sln");
	extern std::string fileOut_Residual(BASE_PATH + "/residual.txt");
	//debug
	//extern std::string fileOut_debug("D:/calculations/new_test/debug.txt");
	extern std::string fileOut_debug(BASE_PATH + "/debug.plt");
	extern std::string fileOut_debug_dsln(BASE_PATH + "/dsln.plt");

	bool parameter::readParameterFile(const std::string &filename)
	{
#define IF_NOTES

#ifdef IF_NOTES
		real inputParameters[50]; //����������50������
		std::ifstream fileIn;
		fileIn.open(filename.c_str());
		if (!fileIn)
		{
			std::cerr << "error: fail to open grid file: unable to open \"" << filename << "\"\n";
			exit(1);
		}
		std::cout << "reading parameter file..." << std::endl;
		std::string DateLine;
		int ii = 1;
		try
		{
			while (!fileIn.eof())
			{
				std::getline(fileIn, DateLine);
				std::cout << DateLine << std::endl;
				std::string::size_type isNoteLine = DateLine.find("#");
				if (!isNoteLine)
				{
					continue;
				}
				std::istringstream istrStream(DateLine);
				real str2real;
				istrStream >> str2real;
				inputParameters[ii] = str2real;
				++ii;
			}
		}
		catch (std::exception e)
		{
			// (1) it will catch the error show show
			std::cerr << e.what() << std::endl;
		}
		catch (std::out_of_range e)
		{
			// (2) this will also catch the same error if (1) was not there but could
			// get you more details if you wanted since its more specific but i have
			// not digged into it further
			std::cerr << e.what() << std::endl;
		}
		catch (...)
		{
			// (3) just for sanity check if first two didn't catch it
			std::cerr << "something went wrong";
		}
		fileIn.close();

		unit = inputParameters[1];
		Re = inputParameters[2];
		PrLam = inputParameters[3];
		PrTur = inputParameters[4];
		muReference = inputParameters[5];

		tmpWallReference = inputParameters[6];
		rhoReference = inputParameters[7];
		uxReference = inputParameters[8];
		uyReference = inputParameters[9];
		preReference = inputParameters[10];

		tmpReference = inputParameters[11];
		maReference = inputParameters[12];
		CFL = inputParameters[13];
		totalTimePhysical = inputParameters[14];
		dTReferencePhysical = inputParameters[15];

		ndim = static_cast<int>(inputParameters[16]);
		reconstructionOrder = static_cast<int>(inputParameters[17]);
		nStart = static_cast<int>(inputParameters[18]);
		nEnd = static_cast<int>(inputParameters[19]);
		nScreenOutput = static_cast<int>(inputParameters[20]);

		nFileOutput = static_cast<int>(inputParameters[21]);
		nBackUp = static_cast<int>(inputParameters[22]);
		nStepTimeMarching = static_cast<int>(inputParameters[23]);
		nInnerStepTimeMarching = static_cast<int>(inputParameters[24]);
		isContinue = static_cast<bool>(inputParameters[25]);

		isPeriodicBoundary = static_cast<bool>(inputParameters[26]);
		isViscous = static_cast<bool>(inputParameters[27]);
		EPS = inputParameters[28];
		isImplicitTimeMarching = static_cast<bool>(inputParameters[29]);
		isLocalTimeMarching = static_cast<bool>(inputParameters[30]);

		cofAV = inputParameters[31];

		//timeMarchingStep		= static_cast<int>(inputParameters[31]);
		std::cout << "------------------------------------------------------------------" << std::endl;
		std::cout << "	testing..." << std::endl;
		std::cout << "\t"
				  << "reconstructionOrder"
				  << "\t" << reconstructionOrder << std::endl;
		std::cout << "\t"
				  << "isContinue"
				  << "\t"
				  << "\t" << isContinue << std::endl;
		std::cout << "\t"
				  << "isPeriodicBoundary"
				  << "\t" << isPeriodicBoundary << std::endl;
		std::cout << "\t"
				  << "isViscous"
				  << "\t"
				  << "\t" << isViscous << std::endl;
		std::cout << "\t"
				  << "EPS"
				  << "\t"
				  << "\t" << EPS << std::endl;
		std::cout << "\t"
				  << "isImplicitTimeMarching"
				  << "\t" << isImplicitTimeMarching << std::endl;
		std::cout << "\t"
				  << "nStepTimeMarching"
				  << "\t" << nStepTimeMarching << std::endl;
		std::cout << "------------------------------------------------------------------" << std::endl;
		std::cout << " ..parameter file has been read." << std::endl;
		//system("pause");
#else
		std::ifstream fileIn;
		fileIn.open(filename.c_str());
		if (!fileIn)
		{
			std::cerr << "error: fail to open grid file: unable to open \"" << filename << "\"\n";
			exit(1);
		}
		std::cout << "reading parameter file..." << std::endl;
		//std::string skipCharacters;
		//std::getline(fileIn, skipCharacters);
		fileIn >> unit;
		//std::getline(fileIn, skipCharacters);
		fileIn >> Re;
		fileIn >> PrLam;
		fileIn >> PrTur;
		fileIn >> muReference;

		fileIn >> tmpWallReference;
		fileIn >> rhoReference;
		fileIn >> uxReference;
		fileIn >> uyReference;
		fileIn >> preReference;

		fileIn >> tmpReference;
		fileIn >> maReference;
		fileIn >> CFL;
		fileIn >> totalTimePhysical;
		fileIn >> dTReferencePhysical;

		fileIn >> ndim;
		fileIn >> reconstructionOrder;
		fileIn >> nStart;
		fileIn >> nEnd;
		fileIn >> nScreenOutput;

		fileIn >> nFileOutput;
		fileIn >> nBackUp;
		fileIn >> nStepTimeMarching;
		fileIn >> nInnerStepTimeMarching;
		fileIn >> isContinue;

		fileIn >> isPeriodicBoundary;
		fileIn >> isViscous;
		fileIn >> EPS;
		fileIn >> isImplicitTimeMarching;
		fileIn >> isLocalTimeMarching;

		//fileIn >> timeMarchingStep;
		fileIn.close();

		std::cout << "------------------------------------------------------------------" << std::endl;
		std::cout << "	testing..." << std::endl;
		std::cout << "\t"
				  << "reconstructionOrder"
				  << "\t" << reconstructionOrder << std::endl;
		std::cout << "\t"
				  << "isContinue"
				  << "\t"
				  << "\t" << isContinue << std::endl;
		std::cout << "\t"
				  << "isPeriodicBoundary"
				  << "\t" << isPeriodicBoundary << std::endl;
		std::cout << "\t"
				  << "isViscous"
				  << "\t" << isViscous << std::endl;
		std::cout << "------------------------------------------------------------------" << std::endl;
		std::cout << " ..parameter file has been read." << std::endl;

#endif
		return true;
	}
}