#ifndef MDETECTORCONFIG_H
#define MDETECTORCONFIG_H

using namespace std;

class MDetectorConfig
{
	public:
	MDetectorConfig();
	~MDetectorConfig();


private:
	int fNumberOfPadsX;
	int fNumberOfPadsY;
	double fPadSizeX;
	double fPadSizeY;
	double fDriftVelocity;
	double fPressure;
	string fGas;

public:
	///////////////
	/// SETTERS ///
	///////////////
	void SetNumberOfPadsX(int NbrOfPadsX) {fNumberOfPadsX=NbrOfPadsX;}
	void SetNumberOfPadsY(int NbrOfPadsY) {fNumberOfPadsY=NbrOfPadsY;}
	void SetPadSizeX(double PadSizeX) {fPadSizeX=PadSizeX;}
	void SetPadSizeY(double PadSizeY) {fPadSizeY=PadSizeY;}
	void SetDriftVelocity(double DriftVelocity) {fDriftVelocity=DriftVelocity;}
	void SetPressure(double Pressure) {fPressure=Pressure;}
	void SetGas(string Gas) {fGas=Gas;}

	///////////////
	/// GETTERS ///
	///////////////
	int GetNumberOfPadsX() const {return fNumberOfPadsX;}
	int GetNumberOfPadsY() const {return fNumberOfPadsY;}
	double GetPadSizeX() const {return fPadSizeX;}
	double GetPadSizeY() const {return fPadSizeY;}
	double GetDriftVelocity() const {return fDriftVelocity;}
	string GetGas() const {return fGas;}

	void ReadConfigurationFile(string filename);

	//ClassDef(MDetectorConfig,1);

};

#endif
