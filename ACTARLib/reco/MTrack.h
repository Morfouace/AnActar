#ifndef MTRACK_H
#define MTRACK_H

#include <TLine.h>
#include <TVector3.h>
#include "../utils/Parameters.h"
#include <stdio.h>
#include <iostream>

class MTrack
{
	public:
	MTrack();
	~MTrack();
	void ResetLines();
	TVector3 GetDirectionVector();
	TVector3 GetVertexPostion(TVector3 vBeamDir, TVector3 vBeamPoint);
	double GetXm() {return Xm;}
	double GetXh() {return Xh;}
	double GetYm() {return Ym;}
	double GetYh() {return Yh;}
	double GetZm() {return Zm;}
	double GetZh() {return Zh;}

	int zx_s;
	int zx_e;
	int zy_s;
	int zy_e;

	double Xm;
	double Ym;
	double Zm;
	double Xh;
	double Yh;
	double Zh;

	TLine* L2DXY;
	TLine* L2DXZ;
	TLine* L2DYZ;
	TLine* L3D;

	//float ElossTable[NPADX][2];

	ClassDef(MTrack, 1);
};

#endif
