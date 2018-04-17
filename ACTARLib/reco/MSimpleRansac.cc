///////////////////////////////////////
//                                   //
// P. Morfouace - GANIL 2017         //
// email: morfouace@ganil.fr         //
//                                   //
// MSimpleRansac class file:                //
//     - Fit Track                   //
//                                   //
///////////////////////////////////////


#include "MSimpleRansac.h"

using namespace std;
//ClassImp(MSimpleRansac)

/////////////////////////////////////////////////
MSimpleRansac::MSimpleRansac(int NumberOfPadsX, int NumberOfPadsY, int verbose)
{
	fRANSACMaxIteration = 1000;
	fRANSACThreshold = 60;//100;
	fRANSACPointThreshold = 0.04;//0.07;
	fRANSACChargeThreshold = 0.04;//0.07;
	fRANSACDistance = 11;//7;

	fNumberOfPadsX = NumberOfPadsX;
	fNumberOfPadsY = NumberOfPadsY;
	fVerbose = verbose;
	//fNumberOfTracksMax = 10;

	if(fVerbose==1){
		hXYZ 	= new TH3F("hXYZ","hXYZ",NumberOfPadsX,0,NumberOfPadsX,NumberOfPadsY,0,NumberOfPadsY,512,0,512);
		pl		= new TGraph2D();
	}
}


//////////////////////////////////////////////////
MSimpleRansac::~MSimpleRansac()
{
	for(unsigned int i=0; i<vline.size(); i++){
			delete vline[i];
	}
	if(fVerbose==1){
		delete hXYZ;
	}
}

//////////////////////////////////////////////////
void MSimpleRansac::Init(vector<int> v1, vector<int> v2, vector<double> v3, vector<double> v4)
{
	Reset();

	vX = v1;
	vY = v2;
	vZ = v3;
	vQ = v4;

	fOriginalCloudSize = vX.size();
	double TotalCharge=0;
	for(unsigned int i=0; i< vQ.size(); i++){
		TotalCharge += vQ[i];
	}
	fTotalCharge = TotalCharge;

}

//////////////////////////////////////////////////
vector<MTrack> MSimpleRansac::SimpleRansac()
{
	if(fVerbose==1){
		c1 = new TCanvas("c1","c1",600,600);
  	c1->cd();
	}

	TRandom* Rand=new TRandom();
	double RemainingCharge = fTotalCharge;
	//vector<double> trackX, trackY, trackZ;

	int aa=0;
	//while(vX.size() > fOriginalCloudSize*fRANSACPointThreshold){
	//while(RemainingCharge > fTotalCharge*fRANSACChargeThreshold){
	while(RemainingCharge > fTotalCharge*fRANSACChargeThreshold && vX.size() > fOriginalCloudSize*fRANSACPointThreshold){
		aa++;

		double x1,x2,y1,y2,z1,z2;
	  double minErr=1E20;
		double track_charge=0;
		/*trackX.clear();
		trackY.clear();
		trackZ.clear();*/

	  TVector3 V1, V2;
	  std::vector< int> inliners;
		inliners.clear();

		if(fVerbose==1){
			hXYZ->Reset();
    	for(int i=0;i<vX.size();i++){
		   	hXYZ->Fill(vX[i],vY[i],vZ[i],vQ[i]);
			}
		}

	  for(int i=0;i<fRANSACMaxIteration;i++){
	  	int p1=(int)(Rand->Uniform(0,vX.size()));
			int p2;
	    do p2=(int)(Rand->Uniform(0,vX.size()));
	    while(p2==p1);

	    TVector3 Vu;
	    Vu.SetXYZ(vX[p2]-vX[p1],vY[p2]-vY[p1],vZ[p2]-vZ[p1]);
			x1=vX[p1];
	    x2=vX[p2];
	    y1=vY[p1];
	    y2=vY[p2];
	    z1=vZ[p1];
	    z2=vZ[p2];

	    std::vector<int> alsoin;
	   	alsoin.clear();
	    for(int p=0;p<vX.size();p++){
	    	if(p!=p1 && p!=p2 && vQ[p]){
	      	TVector3 Va;
	        Va.SetXYZ(vX[p2]-vX[p],vY[p2]-vY[p],vZ[p2]-vZ[p]);
	        TVector3 VaCVu=Va.Cross(Vu);
	        if(VaCVu.Mag()/Vu.Mag()<fRANSACDistance){
	        	alsoin.push_back(p);
						//cout << vY[p] << endl;
						if((vY[p]>0 && vY[p]<5) || (vY[p]<31 && vY[p]>26)){
							track_charge += vQ[p];
							/*trackX.push_back(vX[p]);
							trackY.push_back(vY[p]);
							trackZ.push_back(vZ[p]);*/
						}
	        }
	      }
			}

	    if(alsoin.size()>fRANSACThreshold){
	    	TVector3 v1, v2;
	      double chi2=Fit3D(vX,vY,vZ,vQ,alsoin,v1,v2);
	      if(chi2<minErr){
	      	minErr=chi2;
	        V1=v1;
	        V2=v2;
	        inliners=alsoin;
	     	}
			}
		}

		if(inliners.size()==0) break;
		//cout << "iteration " << aa << ":  with also inliners: " << inliners.size() << endl;

		TVector3 Vdir = TVector3(V2.x()-V1.x(),V2.y()-V1.y(),V2.z()-V1.z());
		Vdir = Vdir.Unit();

		MTrack myTrack = MTrack();

		myTrack.Xm=V1.x();
		myTrack.Ym=V1.y();
		myTrack.Zm=V1.z();

		myTrack.Xh=V2.x();
		myTrack.Yh=V2.y();
		myTrack.Zh=V2.z();

		vTrack.push_back(myTrack);
		//vTrackCharge.push_back(track_charge);

		/*if(trackX.size()>0){
			vTrackXmax.push_back(*max_element(trackX.begin(), trackX.end()));
			vTrackXmin.push_back(*min_element(trackX.begin(), trackX.end()));
			vTrackYmax.push_back(*max_element(trackY.begin(), trackY.end()));
			vTrackYmin.push_back(*min_element(trackY.begin(), trackY.end()));
			vTrackZmax.push_back(*max_element(trackZ.begin(), trackZ.end()));
			vTrackZmin.push_back(*min_element(trackZ.begin(), trackZ.end()));
		}*/

		/*cout << "Xh= " << myTrack.Xh << endl;
		cout << "Yh= " << myTrack.Yh << endl;
		cout << "Zh= " << myTrack.Zh << endl;*/

		if(fVerbose==1){
			hXYZ->Draw();
			pl->Clear();
			pl->SetPoint(0,V1.x(), V1.y(), V1.z());
			pl->SetPoint(1,V2.x(), V2.y(), V2.z());
			pl->SetPoint(2,V1.x() + 10*Vdir.x(), V1.y() + 10*Vdir.y(), (V1.z() + 10*Vdir.z()));
			pl->SetPoint(3,V1.x() + 100*Vdir.x(), V1.y() + 100*Vdir.y(), (V1.z() + 100*Vdir.z()));
			pl->SetPoint(4,V1.x() - 100*Vdir.x(), V1.y() - 100*Vdir.y(), (V1.z() - 100*Vdir.z()));
			pl->SetLineColor(2);
			pl->SetMarkerColor(2);
			pl->SetMarkerStyle(8);
			pl->SetMarkerSize(1);
			pl->SetLineWidth(2);
			pl->Draw("tri1 same");
			c1->WaitPrimitive();
			c1->Update();
		}

	  for(int i=inliners.size()-1 ; i>=0 ; i--){
			vQ.erase(vQ.begin()+inliners[i]);
	    vX.erase(vX.begin()+inliners[i]);
	    vY.erase(vY.begin()+inliners[i]);
			vZ.erase(vZ.begin()+inliners[i]);
	  }

		RemainingCharge = 0;
		for(unsigned int i =0; i<vQ.size(); i++){
			RemainingCharge += vQ[i];
		}
		//cout << "/// RemainingCharge= " << RemainingCharge << endl;
		//cout << "////****** " << vX.size() << " / " <<  fOriginalCloudSize*fRANSACPointThreshold << endl;
		//cout << "////****** " << RemainingCharge << " / " <<  fTotalCharge*fRANSACChargeThreshold << endl;
	}
	//cout << " //// Number of Tracks found " << vTrack.size() << " ////" << endl;

	return vTrack;
}

//////////////////////////////////////////////////
vector<double> MSimpleRansac::GetChargeOfTracks()
{
	return vTrackCharge;
}

//////////////////////////////////////////////////
vector<double> MSimpleRansac::GetTrackLength(double PadSizeX, double PadSizeY, double DriftVelocity)
{
	vector<double> length;
	double Ymin, Ymax;
	double Zmin, Zmax;
	double Xmin = 240;
	double Xmax = 256;

	for(unsigned int i=0; i<vTrack.size(); i++){
		double xh = vTrack[i].GetXh()*PadSizeX;
		double xm = vTrack[i].GetXm()*PadSizeX;
		double yh = vTrack[i].GetYh()*PadSizeY;
		double ym = vTrack[i].GetYm()*PadSizeY;
		double zh = vTrack[i].GetZh()*DriftVelocity;
		double zm = vTrack[i].GetZm()*DriftVelocity;

		Ymin = yh + (ym - yh)*(Xmin-xh)/(xm-xh);
		Ymax = yh + (ym - yh)*(Xmax-xh)/(xm-xh);

		Zmin = zh + (zm - zh)*(Xmin-xh)/(xm-xh);
		Zmax = zh + (zm - zh)*(Xmax-xh)/(xm-xh);

		length.push_back(sqrt( pow(Xmax-Xmin,2) + pow(Ymax-Ymin,2) + pow(Zmax-Zmin,2) ));
	}

	return length;
}

//////////////////////////////////////////////////
/*void MSimpleRansac::InitTrack(int xs, int xe, int ys, int ye)
{
	myTrack.zx_s=xs;
	myTrack.zx_e=xe;
	myTrack.zy_e=ye;
	myTrack.zy_s=ys;
}*/

//////////////////////////////////////////////////
void MSimpleRansac::Reset()
{
	if(fVerbose==1) hXYZ->Reset();
	vline.clear();
	vX.clear();
	vY.clear();
	vZ.clear();
	vQ.clear();
	vTrack.clear();
	vTrackCharge.clear();
	/*vTrackXmax.clear();
	vTrackXmin.clear();
	vTrackYmax.clear();
	vTrackYmin.clear();
	vTrackZmax.clear();
	vTrackZmin.clear();*/
}

//////////////////////////////////////////////////
double MSimpleRansac::Fit3D(vector<int> X, vector<int> Y, vector<double> Z, vector<double> Charge, vector<int> inliners, TVector3& V1, TVector3& V2)
{
    // adapted from: http://fr.scribd.com/doc/31477970/Regressions-et-trajectoires-3D

    int R, C;
    double Q;
    double Xm,Ym,Zm;
    double Xh,Yh,Zh;
    double a,b;
    double Sxx,Sxy,Syy,Sxz,Szz,Syz;
    double theta;
    double K11,K22,K12,K10,K01,K00;
    double c0,c1,c2;
    double p,q,r,dm2;
    double rho,phi;

    Q=Xm=Ym=Zm=0.;
		double total_charge=0;
    Sxx=Syy=Szz=Sxy=Sxz=Syz=0.;

    for (auto i : inliners)
    {
	if(X[i]>119)total_charge+=Charge[i];
        Q+=Charge[i]/10.;
        Xm+=X[i]*Charge[i]/10.;
        Ym+=Y[i]*Charge[i]/10.;
        Zm+=Z[i]*Charge[i]/10.;
        Sxx+=X[i]*X[i]*Charge[i]/10.;
        Syy+=Y[i]*Y[i]*Charge[i]/10.;
        Szz+=Z[i]*Z[i]*Charge[i]/10.;
        Sxy+=X[i]*Y[i]*Charge[i]/10.;
        Sxz+=X[i]*Z[i]*Charge[i]/10.;
        Syz+=Y[i]*Z[i]*Charge[i]/10.;
    }
    vTrackCharge.push_back(total_charge);

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
        dm2=min(-c2/3.+2.*pow(rho,1./3.)*cos(phi/3.),min(-c2/3.+2.*pow(rho,1./3.)*cos((phi+2.*TMath::Pi())/3.),-c2/3.+2.*pow(rho,1./3.)*cos((phi+4.*TMath::Pi())/3.)));
    }

    a=-K10*cos(theta)/(K11-dm2)+K01*sin(theta)/(K22-dm2);
    b=-K10*sin(theta)/(K11-dm2)-K01*cos(theta)/(K22-dm2);

    Xh=((1.+b*b)*Xm-a*b*Ym+a*Zm)/(1.+a*a+b*b);
    Yh=((1.+a*a)*Ym-a*b*Xm+b*Zm)/(1.+a*a+b*b);
    Zh=((a*a+b*b)*Zm+a*Xm+b*Ym)/(1.+a*a+b*b);

    V1.SetXYZ(Xm,Ym,Zm);
    V2.SetXYZ(Xh,Yh,Zh);

    return(fabs(dm2/Q));
}
