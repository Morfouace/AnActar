#ifndef MSIMPLERANSAC_H
#define MSIMPLERANSAC_H

// ROOT Headers
#include <TLine.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraph2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>

#include <stdio.h>
#include <vector>

#include "MTrack.h"

using namespace std;

class MSimpleRansac
{
	public:
	MSimpleRansac(int NumberOfPadsX, int NumberOfPadsY, int verbose);
	~MSimpleRansac();

	void Reset();
	void Init(vector<int> v1, vector<int> v2, vector<double> v3, vector<double> v4);
	//void InitTrack(int xs, int xe, int ys, int ye);
	vector<MTrack> SimpleRansac();
	vector<double> GetChargeOfTracks();
	vector<double> GetTrackLength(double PadSizeX, double PadSizeY, double DriftVelocity);
	double Fit3D(vector<int> X, vector<int> Y, vector<double> Z, vector<double> Charge, vector<int> inliners, TVector3& V1, TVector3& V2);

private:
	TCanvas* c1;
	vector<TLine*> vline;
	vector<int> vX, vY;
	vector<double> vZ, vQ;
	TLine* L;
	//MTrack myTrack;
	TH3F* hXYZ;
	TGraph2D* pl;
	vector<MTrack> vTrack;
	vector<double> vTrackCharge;
	vector<double> vTrackXmax;
	vector<double> vTrackXmin;
	vector<double> vTrackYmax;
	vector<double> vTrackYmin;
	vector<double> vTrackZmax;
	vector<double> vTrackZmin;

private:
	float fRANSACThreshold;
	float fRANSACPointThreshold;
	float fRANSACChargeThreshold;
	float fRANSACDistance;
	float fRANSACMaxIteration;
	int fNumberOfTracksMax;
	int fOriginalCloudSize;
	double fTotalCharge;
	int fNumberOfPadsX;
	int fNumberOfPadsY;
	int fVerbose;

	ClassDef(MSimpleRansac, 1);

};

#endif
