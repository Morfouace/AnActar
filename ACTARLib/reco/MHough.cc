///////////////////////////////////////
//                                   //
// P. Morfouace - GANIL 2017         //
// email: morfouace@ganil.fr         //
//                                   //
// MHough class file:                //
//     - Fit Track                   //
//                                   //
///////////////////////////////////////


#include "MHough.h"

using namespace std;
//ClassImp(MHough)
/////////////////////////////////////////////////
MHough::MHough()
{
	hHough 	= new TH2F("hHough","hHough",900,-90,90,400,-400,400);
	hXY = new TH2F("hXY","hXY",128,0,128,128,0,128);
	hYZ = new TH2F("hYZ","hYZ",128,0,128,128,0,512);
}


//////////////////////////////////////////////////
MHough::~MHough()
{
	for(unsigned int i=0; i<vline.size(); i++){
			delete vline[i];
	}
	delete hHough;
	delete hXY;
	delete hYZ;

}

//////////////////////////////////////////////////
void MHough::Init(TH2F* XY, TH2F *YZ, map<int,vector<int>> mXYZ, map<int,vector<int>> mQ)
{
	Reset();
	hXY = XY;
	hYZ = YZ;
	mPads = mXYZ;
	mCharge = mQ;

	visu=false;
	if(visu){
		Ctest=new TCanvas("Ctest","Ctest",1200,800);
		Ctest->Divide(3,1);
	}

	for(int panel=0; panel<2; panel++){
		for(int iY=0+panel*64; iY<64+panel*64; iY++){
			for(int iZ=0; iZ<128; iZ++){
				for(float angle=-90; angle<=90; angle+=0.2){
					Hfill[iY][iZ].push_back(hHough->FindBin(angle,iY*cos(angle*TMath::Pi()/180)+iZ*sin(angle*TMath::Pi()/180)));
					hHough->Fill(angle,iY*cos(angle*TMath::Pi()/180)+iZ*sin(angle*TMath::Pi()/180),hYZ->GetBinContent(iY,iZ));
				}
			}
		}
		FindMaxima(panel,6,visu);
	}
	cout << "Number of Tracks found: " << vline.size() << endl;
}
//////////////////////////////////////////////////
void MHough::FindMaxima(int panel, double dmax, bool visu)
{
	vector<TVector3> vV3;

	while(hHough->GetMaximum()>24000){
		int thm,rhm,hom;
		hHough->GetMaximumBin(thm,rhm,hom);
		double thM=hHough->GetXaxis()->GetBinCenter(thm);
		double rhM=hHough->GetYaxis()->GetBinCenter(rhm);
		//cout << thM << " " << rhM << endl;
		double x1=64*panel;
		double x2=64+64*panel;
		double y1=4*(rhM-x1*cos(thM*TMath::Pi()/180))/(sin(thM*TMath::Pi()/180));
		double y2=4*(rhM-x2*cos(thM*TMath::Pi()/180))/(sin(thM*TMath::Pi()/180));
		//aLine = new TLine(x1,y1,x2,y2);
		//vline.push_back(aLine);

		if(visu){
			Ctest->cd(1);
			hYZ->Draw("colz");
			//vline[vline.size()-1]->SetLineWidth(2);
			//vline[vline.size()-1]->Draw("same");
			Ctest->cd(2);
			hHough->Draw("colz");
			Ctest->cd(3);
			hXY->Draw("colz");
			Ctest->Update();
			Ctest->WaitPrimitive();
		}

		for(int iY=0+panel*64; iY<64+panel*64; iY++){
			for(int iZ=0; iZ<128; iZ++){
				if(hYZ->GetBinContent(iY,iZ)>0){
					double yline = (rhM-iZ*sin(thM*TMath::Pi()/180))/cos(thM*TMath::Pi()/180);
					double zline = (rhM-iY*cos(thM*TMath::Pi()/180))/sin(thM*TMath::Pi()/180);

					if(sqrt(pow(iY-yline,2)+pow(iZ-zline,2))<dmax){
						for(unsigned int i=0;i<Hfill[iY][iZ].size();i++){
								hHough->SetBinContent(Hfill[iY][iZ].at(i),hHough->GetBinContent(Hfill[iY][iZ].at(i))-hYZ->GetBinContent(iY,iZ));
							}
							int iYZ = iY + 128*iZ;
							map<int, vector<int>>::iterator it;
							map<int, vector<int>>::iterator itQ;
							it = mPads.find(iYZ);
							itQ = mCharge.find(iYZ);
							//cout << *(it->second.begin()) << "/" << *(it->second.end()-1) << endl;
							//cout << "it.size= " << it->second.size() << endl;
							for(unsigned int i=0; i<it->second.size();i++){
								int iX = *(it->second.begin()+i);
								int iQ = *(itQ->second.begin()+i);
								hXY->SetBinContent(iX,iY,hXY->GetBinContent(iX,iY)-iQ);
								vV3.push_back(TVector3(iX,iY,iZ));
							}
							hYZ->SetBinContent(iY,iZ,0);
							Hfill[iY][iZ].clear();
						}
					}
				}
			}
			//InitTrack(1,127,1,127);
			//FitTrack3D(vV3);
			//vV3.clear();
		}

}

//////////////////////////////////////////////////
void MHough::InitTrack(int xs, int xe, int ys, int ye)
{
	myTrack.zx_s=xs;
	myTrack.zx_e=xe;
	myTrack.zy_e=ye;
	myTrack.zy_s=ys;
}
//////////////////////////////////////////////////
void MHough::Reset()
{
	hHough->Reset();
	hXY->Reset();
	hYZ->Reset();
	mPads.clear();
	mCharge.clear();
	vline.clear();
}

//////////////////////////////////////////////////
/*double MHough::FitTrack3D(vector<TVector3> vV, MTrack* T)
{
	// adapted from: http://fr.scribd.com/doc/31477970/Regressions-et-trajectoires-3D

	int Rmin=T->zx_s;
	int Rmax=T->zx_e;
	int Cmi root cernn=T->zy_s;
	int Cmax=T->zy_e;
	int R, C;
	double Q,X,Y,Z;
	double Xm,Ym,Zm;
	double Sxx,Sxy,Syy,Sxz,Szz,Syz;
	double theta;
	double K11,K22,K12,K10,K01,K00;
	double c0,c1,c2;
	double p,q,r,dm2;
	double rho,phi;
	double a,b;

	Q=Xm=Ym=Zm=0.;
	Sxx=Syy=Szz=Sxy=Sxz=Syz=0.;

	for (R=Rmin;R<=Rmax;R++)
		for (C=Cmin;C<=Cmax;C++)
			if(PAD[R][C]>threshold && TIME[R][C])
			{
				X= C+0.5;// *2.+1.;
				Y= R+0.5;// *2.+1.;
				Z= TIME[R][C];// *VDRIFT;
				Q+=PAD[R][C]/10.;
				Xm+=X*PAD[R][C]/10.;
				Ym+=Y*PAD[R][C]/10.;
				Zm+=Z*PAD[R][C]/10.;
				Sxx+=X*X*PAD[R][C]/10.;
				Syy+=Y*Y*PAD[R][C]/10.;
				Szz+=Z*Z*PAD[R][C]/10.;
				Sxy+=X*Y*PAD[R][C]/10.;
				Sxz+=X*Z*PAD[R][C]/10.;
				Syz+=Y*Z*PAD[R][C]/10.;
			}
	Xm/=Q;
	Ym/=Q;
	Zm/=Q;
	Sxx/=Q;
	Syy/=Q;
	Szz/=Q;
	Sxy/=Q;
	Sxz/=Q;
	Syz/=Q;
	Sxx-=(Xm*Xm);
	Syy-=(Ym*Ym);
	Szz-=(Zm*Zm);
	Sxy-=(Xm*Ym);
	Sxz-=(Xm*Zm);
	Syz-=(Ym*Zm);

	theta=0.5*atan((2.*Sxy)/(Sxx-Syy));

	K11=(Syy+Szz)*pow(cos(theta),2)+(Sxx+Szz)*pow(sin(theta),2)-2.*Sxy*cos(theta)*sin(theta);
	K22=(Syy+Szz)*pow(sin(theta),2)+(Sxx+Szz)*pow(cos(theta),2)+2.*Sxy*cos(theta)*sin(theta);
	K12=-Sxy*(pow(cos(theta),2)-pow(sin(theta),2))+(Sxx-Syy)*cos(theta)*sin(theta);
	K10=Sxz*cos(theta)+Syz*sin(theta);
	K01=-Sxz*sin(theta)+Syz*cos(theta);
	K00=Sxx+Syy;

	c2=-K00-K11-K22;
	c1=K00*K11+K00*K22+K11*K22-K01*K01-K10*K10;
	c0=K01*K01*K11+K10*K10*K22-K00*K11*K22;


	p=c1-pow(c2,2)/3.;
	q=2.*pow(c2,3)/27.-c1*c2/3.+c0;
	r=pow(q/2.,2)+pow(p,3)/27.;


	if(r>0) dm2=-c2/3.+pow(-q/2.+sqrt(r),1./3.)+pow(-q/2.-sqrt(r),1./3.);
	if(r<0)
	{
		rho=sqrt(-pow(p,3)/27.);
		phi=acos(-q/(2.*rho));
		dm2=min(-c2/3.+2.*pow(rho,1./3.)*cos(phi/3.),min(-c2/3.+2.*pow(rho,1./3.)*cos((phi+2.*M_PI)/3.),-c2/3.+2.*pow(rho,1./3.)*cos((phi+4.*M_PI)/3.)));
	}
	a=-K10*cos(theta)/(K11-dm2)+K01*sin(theta)/(K22-dm2);
	b=-K10*sin(theta)/(K11-dm2)-K01*cos(theta)/(K22-dm2);

	T->Xm=Xm;
	T->Ym=Ym;
	T->Zm=Zm;
	T->Xh=((1.+b*b)*Xm-a*b*Ym+a*Zm)/(1.+a*a+b*b);
	T->Yh=((1.+a*a)*Ym-a*b*Xm+b*Zm)/(1.+a*a+b*b);
	T->Zh=((a*a+b*b)*Zm+a*Xm+b*Ym)/(1.+a*a+b*b);

	T->L2DXY->SetX1(T->Xm);
	T->L2DXY->SetY1(T->Ym);
	T->L2DXY->SetX2(T->Xh);
	T->L2DXY->SetY2(T->Yh);

	T->L2DXZ->SetX1(T->Xm);
	T->L2DXZ->SetY1(T->Zm);
	T->L2DXZ->SetX2(T->Xh);
	T->L2DXZ->SetY2(T->Zh);

	T->L2DYZ->SetX1(T->Ym);
	T->L2DYZ->SetY1(T->Zm);
	T->L2DYZ->SetX2(T->Yh);
	T->L2DYZ->SetY2(T->Zh);

	return(dm2/Q);
}
*/
