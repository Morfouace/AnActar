///////////////////////////////////////
//                                   //
// P. Morfouace - GANIL 2017         //
// email: morfouace@ganil.fr         //
//                                   //
// MDetectorConfig class file:                //
//     - DetectorConfig                   //
//                                   //
///////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "MDetectorConfig.h"

using namespace std;


/////////////////////////////////////////////////
MDetectorConfig::MDetectorConfig()
{
}


//////////////////////////////////////////////////
MDetectorConfig::~MDetectorConfig()
{
}

//////////////////////////////////////////////////
void MDetectorConfig::ReadConfigurationFile(string filename)
{
	std::ifstream ifile(filename.c_str());
  if(!ifile.is_open()){
    std::cout << "file: " << filename << " not found!!! " << std::endl;
	}

	string line;
	string DataBuffer;
	bool ReadingStatus=false;
	bool check_padsX=false;
	bool check_padsY=false;
	bool check_sizeX=false;
	bool check_sizeY=false;
	bool check_gas=false;
	bool check_pressure=false;
	bool check_driftvelocity=false;

	while(!ifile.eof()){
		getline(ifile,line);
		if(line.compare(0,12,"config_actar")==0){
			cout << "ACTAR config file found" << endl;
			ReadingStatus=true;
		}
		else ReadingStatus=false;

		while(ReadingStatus){
			ifile >> DataBuffer;
			if(DataBuffer.compare(0,1,"#")==0){
				ifile.ignore( std::numeric_limits<std::streamsize>::max(), '\n');
			}

			else if(DataBuffer.compare(0,14,"NumberOfPadsX=")==0){
				ifile >> DataBuffer;
				fNumberOfPadsX = atoi(DataBuffer.c_str());
				check_padsX=true;
				cout << "/// Number Of Pads X= " << fNumberOfPadsX << " ///" << endl;
			}

			else if(DataBuffer.compare(0,14,"NumberOfPadsY=")==0){
				ifile >> DataBuffer;
				fNumberOfPadsY = atoi(DataBuffer.c_str());
				check_padsY=true;
				cout << "/// Number Of Pads Y= " << fNumberOfPadsY << " ///" << endl;
			}

			else if(DataBuffer.compare(0,9,"PadSizeX=")==0){
				ifile >> DataBuffer;
				fPadSizeX = atof(DataBuffer.c_str());
				check_sizeX=true;
				cout << "/// Pad Size X= " << fPadSizeX << " ///" << endl;
			}

			else if(DataBuffer.compare(0,9,"PadSizeY=")==0){
				ifile >> DataBuffer;
				fPadSizeY = atof(DataBuffer.c_str());
				check_sizeY=true;
				cout << "/// Pad Size Y= " << fPadSizeY << " ///" << endl;
			}

			else if(DataBuffer.compare(0,9,"Pressure=")==0){
				ifile >> DataBuffer;
				fPressure = atof(DataBuffer.c_str());
				check_pressure=true;
				cout << "/// Pressure= " << fPressure << " ///" << endl;
			}

			else if(DataBuffer.compare(0,14,"DriftVelocity=")==0){
				ifile >> DataBuffer;
				fDriftVelocity = atof(DataBuffer.c_str());
				check_driftvelocity=true;
				cout << "/// Drift Velocity= " << fDriftVelocity << " ///" << endl;
			}

			else if(DataBuffer.compare(0,4,"Gas=")==0){
				ifile >> DataBuffer;
				fGas = DataBuffer.c_str();
				check_gas=true;
				cout << "/// Gas Type= " << fGas << " ///" << endl;
			}

			if(check_gas && check_driftvelocity && check_pressure && check_sizeX && check_sizeY && check_padsX && check_padsY){
				ReadingStatus=false;
			}
		}

		//cout << line << endl;
	}

}
