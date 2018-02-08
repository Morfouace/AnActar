///////////////////////////////////////
//                                   //
// T. Roger - GANIL 2015             //
// email: roger@ganil.fr             //
//                                   //
// MTrack class file:                //
//     - Track  graphical properties //
//                                   //
///////////////////////////////////////


#include "MTrack.h"

using namespace std;

//ClassImp(MTrack)

//////////////////////////////////////////////////////
MTrack::MTrack()
{
	L2DXY=new TLine();
	L2DXZ=new TLine();
	L2DYZ=new TLine();
	L3D=new TLine();
}

//////////////////////////////////////////////////////
MTrack::~MTrack()
{
	/*delete L2DXY;
	delete L2DXZ;
	delete L2DYZ;
	delete L3D;*/
}

//////////////////////////////////////////////////////
TVector3 MTrack::GetDirectionVector()
{
	TVector3 vTrack = TVector3(Xh-Xm, Yh-Ym, Zh-Zm);
	return vTrack;
}

//////////////////////////////////////////////////////
TVector3 MTrack::GetVertexPostion(TVector3 vBeamDir, TVector3 vBeamPoint)
{
	//y_beam = Ay_beam*x_beam + By_beam
	double Ay_beam = vBeamDir.Y()/vBeamDir.X();
	double By_beam = vBeamPoint.Y() - (vBeamDir.Y()/vBeamDir.X())*vBeamPoint.X();
	//z_beam = Az_beam*x_beam + Bz_beam
	double Az_beam = vBeamDir.Z()/vBeamDir.X();
	double Bz_beam = vBeamPoint.Z() - (vBeamDir.Z()/vBeamDir.X())*vBeamPoint.X();

	//y_track = Ay_track*x_track + By_track
	double Ay_track = GetDirectionVector().Y()/GetDirectionVector().X();
	double By_track = Ym - (GetDirectionVector().Y()/GetDirectionVector().X())*Xm;
	//z_track = Az_track*x_track + Bz_track
	double Az_track = GetDirectionVector().Z()/GetDirectionVector().X();
	double Bz_track = Zm - (GetDirectionVector().Z()/GetDirectionVector().X())*Xm;

	double DBy 		= By_beam - By_track;
	double DBz 		= Bz_beam - Bz_track;
	double alpha 	= (-Ay_beam*DBy - Az_beam*DBz)/(Ay_beam*Ay_beam + Az_beam*Az_beam + 1);
	double beta 	= (Ay_beam*Ay_track + Az_beam*Az_track + 1)/(Ay_beam*Ay_beam + Az_beam*Az_beam + 1);

	double A			= beta*(Ay_beam*Ay_beam + Az_beam*Az_beam + 1) - (Ay_beam*Ay_track + Az_beam*Az_track + 1);
	double B			= (Ay_track*Ay_track + Az_track*Az_track + 1) - beta*(Ay_beam*Ay_track + Az_beam*Az_track + 1);
	double C 			= beta*(Ay_beam*DBy + Az_beam*DBz) - Ay_track*DBy - Az_track*DBz;

	double xt			= -(A*alpha+C)/(A*beta+B);
	double xb			= beta*xt+alpha;

	double yb			= Ay_beam*xb + By_beam;
	double zb			= Az_beam*xb + Bz_beam;

	double yt			= Ay_track*xt + By_track;
	double zt			= Az_track*xt + Bz_track;

	double xvertex	= (xb+xt)/2;
	double yvertex	= (yb+yt)/2;
	double zvertex	= (zb+zt)/2;

	TVector3 vVertex = TVector3(xvertex,yvertex,zvertex);

	return vVertex;
}

//////////////////////////////////////////////////////
void MTrack::ResetLines()
{
	L2DXY->SetX1(-1);
	L2DXY->SetY1(-1);
	L2DXY->SetX2(-1);
	L2DXY->SetY2(-1);

	L2DXZ->SetX1(-1);
	L2DXZ->SetY1(-1);
	L2DXZ->SetX2(-1);
	L2DXZ->SetY2(-1);

	L2DYZ->SetX1(-1);
	L2DYZ->SetY1(-1);
	L2DYZ->SetX2(-1);
	L2DYZ->SetY2(-1);

	L3D->SetX1(-1);
	L3D->SetY1(-1);
	L3D->SetX2(-1);
	L3D->SetY2(-1);
}
