#ifndef MHOUGH_H
#define MHOUGH_H

// ROOT Headers
#include <TLine.h>
#include <TH2F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TVector3.h>
//#include <Parameters.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <map>

#include "MTrack.h"

using namespace std;

class MHough
{
	public:
	MHough();
	~MHough();
	void Reset();
	void Init(TH2F* XY, TH2F *YZ, map<int,vector<int>> mXYZ, map<int,vector<int>> mQ);
	void FindMaxima(int panel, double dmax, bool visu);
	double FitTrack3D(vector<TVector3> vV, MTrack* T);
	void InitTrack(int xs, int xe, int ys, int ye);

private:
	vector<TLine*> vline;
	TLine* aLine;
	MTrack myTrack;
	TH2F* hHough;
	TH2F* hXY;
	TH2F* hYZ;
	map<int,vector<int>> mPads;
	map<int,vector<int>> mCharge;

	TCanvas* Ctest;
	vector<int> Hfill[128][128];
	bool visu;

	ClassDef(MHough, 1);
};

#endif
