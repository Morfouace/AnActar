#include <iostream>
#include <fstream>
#include "TApplication.h"
#include <TROOT.h>
#include <time.h>

using namespace std;

#include"Analysis.h"

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  TApplication myApp("actaranalysis",&argc,argv);
  //gROOT->SetStyle("pierre_style");

  ReadRunToTreat("RunToTreat.txt");
  Init();


  time_t t ;
  time(&t);
  long double temps_initial = t;

  int nentries = chain->GetEntries();
  for(int i=0;i<nentries;i++){
    Clear();
    FillHistoAndMap(i);
    //cout << "/// " << hXYZ->GetEntries() << endl;
     if(vX.size()>10){
    //if(vX.size()>10 && SiMayaRawEnergy.size()>0){
      //TreatEventWithHough();
	TreatEventWithSimpleRansac();
	TVector3 vX = TVector3(1,0,0);
	TVector3 aTrack;

  	TrackMult = vTrack.size();
  
	if(TrackMult>0){
		double scalarproduct=0;
		int WhichTrack=0;
	  	for(unsigned int i=0; i<TrackMult; i++){
			TVector3 vtest = TVector3(vTrack[i].GetDirectionVector().X(),vTrack[i].GetDirectionVector().Y(),vTrack[i].GetDirectionVector().Z());
			TVector3 vunit = vtest.Unit();
			double scalar = abs(vunit.Dot(vX));
			vScalar.push_back(scalar);
			//cout << scalar << endl;
			//cout << scalarproduct << endl;
			if(scalar>scalarproduct){
				WhichTrack=i; 
				scalarproduct=scalar;
			}
  		}
		//cout << vTrack.size() << " / " << WhichTrack << " / " << scalarproduct << endl;

		double XBeam = vTrack[WhichTrack].GetDirectionVector().X()*Configurator.GetPadSizeX();
		double YBeam = vTrack[WhichTrack].GetDirectionVector().Y()*Configurator.GetPadSizeY();
		double ZBeam = vTrack[WhichTrack].GetDirectionVector().Z()*Configurator.GetDriftVelocity()*20e-3;
		TVector3 vBeam = TVector3(XBeam,YBeam,ZBeam);
	
		double XBeamPoint = vTrack[WhichTrack].GetXh()*Configurator.GetPadSizeX();
		double YBeamPoint = vTrack[WhichTrack].GetYh()*Configurator.GetPadSizeY();
		double ZBeamPoint = vTrack[WhichTrack].GetZh()*Configurator.GetDriftVelocity()*20e-3;
		TVector3 vBeamPos = TVector3(XBeamPoint,YBeamPoint,ZBeamPoint);

		for(unsigned int i=0; i<TrackMult; i++){
    			if(i!=WhichTrack){
			    double Xdir = vTrack[i].GetDirectionVector().X();
			    double Ydir = vTrack[i].GetDirectionVector().Y();
    			    double Zdir = vTrack[i].GetDirectionVector().Z();
			    aTrack = TVector3(Xdir*Configurator.GetPadSizeX(), Ydir*Configurator.GetPadSizeY(), Zdir*20e-3*Configurator.GetDriftVelocity());

		    	    XVertex.push_back(vTrack[i].GetVertexPostion(vBeam,vBeamPos).X());
    		    	    YVertex.push_back(vTrack[i].GetVertexPostion(vBeam,vBeamPos).Y());
			    ZVertex.push_back(vTrack[i].GetVertexPostion(vBeam,vBeamPos).Z());
	
			    double angle = vBeam.Angle(aTrack)*180/TMath::Pi();
			    if(angle>90) angle = 180-angle;
			    vTrackAngle.push_back(angle);

			    if(vTrackLength[i]>0){
			      dedx.push_back(vTrackCharge[i]/(16./cos(angle*TMath::Pi()/180)));
			    }
			    else dedx.push_back(-100);

			    double x1 = vTrack[i].GetXm()*Configurator.GetPadSizeX();
        		    double x2 = vTrack[i].GetXh()*Configurator.GetPadSizeX();
			    double y1 = vTrack[i].GetYm()*Configurator.GetPadSizeY()-0.5*Configurator.GetNumberOfPadsY()*Configurator.GetPadSizeY();
			    double y2 = vTrack[i].GetYh()*Configurator.GetPadSizeY()-0.5*Configurator.GetNumberOfPadsY()*Configurator.GetPadSizeY();
        		    double z1 = -(vTrack[i].GetZm()-256)*Configurator.GetDriftVelocity()*20e-3;
        		    double z2 = -(vTrack[i].GetZh()-256)*Configurator.GetDriftVelocity()*20e-3;
        		    if(vScalar[i]<0.99)GetMayaSiHitPosition(x1,x2,y1,y2,z1,z2);
			}
  		}	
  	

	
		/*for(unsigned int i=0; i<TrackMult; i++){
			if(i!=WhichTrack){
		        	double x1 = vTrack[i].GetXm()*Configurator.GetPadSizeX();
        			double x2 = vTrack[i].GetXh()*Configurator.GetPadSizeX();
			        double y1 = vTrack[i].GetYm()*Configurator.GetPadSizeY()-0.5*Configurator.GetNumberOfPadsY()*Configurator.GetPadSizeY();
			        double y2 = vTrack[i].GetYh()*Configurator.GetPadSizeY()-0.5*Configurator.GetNumberOfPadsY()*Configurator.GetPadSizeY();
        		 	double z1 = -(vTrack[i].GetZm()-256)*Configurator.GetDriftVelocity()*20e-3;
        		 	double z2 = -(vTrack[i].GetZh()-256)*Configurator.GetDriftVelocity()*20e-3;
        		 	GetMayaSiHitPosition(x1,x2,y1,y2,z1,z2);
			}
		}*/
      	}
    }

    OutputTree->Fill();
  }

  OutputFile->Write();
  OutputFile->Close();
  End();

  time(&t) ;
  long double temps_final = t;
  cout << "temps de calcul: " << (int)((temps_final-temps_initial)/60) << " min " << ((int)(temps_final-temps_initial))%60 << " sec" << endl;
  myApp.Run(kTRUE);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
void ReadRunToTreat(string PathToRun){
  ifstream ifile;
  ifile.open(PathToRun.c_str());
  chain = new TChain("ACTAR_TTree");

  /*TString root_path = "/space/morfouace/ACTAR/e750/Data/";
  //TString root_path = "../../";
  //TString root_path = "/home/morfouace/Physics/Actar/Analysis/actar_analysis/root/";
  TString total_path;
  string srun;
  getline(ifile,srun);
  TString tree_name = srun;
  total_path = root_path+tree_name;

  cout << "/// File name: " << total_path << " will be read ///"<< endl;

  TFile* f=new TFile(total_path,"read");
  //Tree=(TTree*)f->Get("Tree_Run0045");
  Tree=(TTree*)f->Get("ACTAR_TTree");*/

  TString root_path = "/space/morfouace/ACTAR/e750/Data/";

  TString total_path;
  string srun;
  while(!ifile.eof()){
    getline(ifile,srun);
    TString tree_name = srun;
    if(srun.compare(0,9,"Tree_Run_")==0){
      total_path = root_path+tree_name;

      cout << "/// File name: " << total_path << " will be read ///"<< endl;

      //TFile* f=new TFile(total_path,"read");
      //Tree=(TTree*)f->Get("ACTAR_TTree");
      chain->Add(total_path);
    }
  }
  cout << "/// Number of events to treat: " << chain->GetEntries() << " ///" << endl;

}

////////////////////////////////////////////////////////////////////////////////
void Init(){
  Configurator.ReadConfigurationFile("config/actar_config.txt");
  PadNumberX = Configurator.GetNumberOfPadsX();
  PadNumberY = Configurator.GetNumberOfPadsY();
  PadSizeX = Configurator.GetPadSizeX();
  PadSizeY = Configurator.GetPadSizeY();
  DriftVelocity = Configurator.GetDriftVelocity();

  Ransac = new MSimpleRansac(PadNumberX,PadNumberY,0);

  /*hXYZ = new TH3F("XYZ","XYZ",PadNumberX,0,PadNumberX,PadNumberY,0,PadNumberY,128,0,512);
  hXY = new TH2F("XY","XY",PadNumberX,0,PadNumberX,PadNumberY,0,PadNumberY);
	hXZ = new TH2F("XZ","XZ",PadNumberX,0,PadNumberX,128,0,512);
	hYZ = new TH2F("YZ","YZ",PadNumberY,0,PadNumberY,128,0,512);*/

  string filename = "./dat/LT.dat";
  ifstream ifile;
  ifile.open(filename.c_str());
 	cout << "Using LookupTable from: " << filename << endl;
 	for(int i=0;i<NumberOfCobo*NumberOfASAD*NumberOfAGET*NumberOfChannel;i++){
    ifile >> TABLE[0][i] >> TABLE[1][i] >> TABLE[2][i] >> TABLE[3][i] >> TABLE[4][i] >> TABLE[5][i];
	//cout << TABLE[5][16476] << endl;
  }
  ifile.close();

  string filename_si = "./config/ACTION_Si_config.dat";
  ifstream ifile_si;
  ifile_si.open(filename_si.c_str());
  string token;
  int vxi_param, si_nbr;
  for(int i=0; i<20; i++){
	ifile_si >> token >> vxi_param >> si_nbr;
	Si_map[vxi_param] = si_nbr+1;
  }

  chain->SetBranchAddress("data",&EvtRed);

  InitOutputTree();
}

////////////////////////////////////////////////////////////////////////////////
void InitOutputTree(){
  OutputFile = new TFile("FittedTracks.root","RECREATE");
  OutputTree = new TTree("PhysicsTree","PhysicsTree");

  OutputTree->Branch("Event",&Event,"Event/I");
  OutputTree->Branch("TrackMult",&TrackMult,"TrackMult/I");
  OutputTree->Branch("X",&vX);
  OutputTree->Branch("Y",&vY);
  OutputTree->Branch("Z",&vZ);
  OutputTree->Branch("vTrackLength",&vTrackLength);
  OutputTree->Branch("vTrackCharge",&vTrackCharge);
  OutputTree->Branch("vTrackAngle",&vTrackAngle);
  OutputTree->Branch("vScalar",&vScalar);
  OutputTree->Branch("dedx",&dedx);
  OutputTree->Branch("ThetaLab",&ThetaLab);
  OutputTree->Branch("PhiLab",&PhiLab);
  OutputTree->Branch("XVertex",&XVertex);
  OutputTree->Branch("YVertex",&YVertex);
  OutputTree->Branch("ZVertex",&ZVertex);

  OutputTree->Branch("SiMayaNumber",&SiMayaNumber);
  OutputTree->Branch("SiMayaRawEnergy",&SiMayaRawEnergy);
  OutputTree->Branch("SiMayaPosY",&SiMayaPosY);
  OutputTree->Branch("SiMayaPosZ",&SiMayaPosZ);
}

////////////////////////////////////////////////////////////////////////////////
void TreatEventWithSimpleRansac(){
  Ransac->Init(vX,vY,vZ,vQ);
  vTrack = Ransac->SimpleRansac();
  vTrackCharge = Ransac->GetChargeOfTracks();
  vTrackLength = Ransac->GetTrackLength(Configurator.GetPadSizeX(), Configurator.GetPadSizeY(), 20e-3*Configurator.GetDriftVelocity());
}

////////////////////////////////////////////////////////////////////////////////
void TreatEventWithHough(){
  Hough.Reset();

  cout << "/// map.size() = " << mXYZ.size() << " ///" << endl;
  //Hough.Init(hXY,hYZ,mXYZ,mQ);
}

////////////////////////////////////////////////////////////////////////////////
void FillHistoAndMap(int i){
  /*hXYZ->Reset();
  hXY->Reset();
  hYZ->Reset();
  hYZ->Reset();*/

  chain->GetEntry(i);
  Event=i;
  if(i%1000==0) cout << "//// Entry= " << i << " //// \r" << flush;
  //cout << "//// Entry= " << i << " //// \r" << endl;

  CleanPad();

  for(unsigned int it=0; it<EvtRed->CoboAsad.size(); it++){
    int co=EvtRed->CoboAsad[it].globalchannelid>>11;
    int as=(EvtRed->CoboAsad[it].globalchannelid - (co<<11))>>9;
    int ag=(EvtRed->CoboAsad[it].globalchannelid - (co<<11)-(as<<9))>>7;
    int ch=EvtRed->CoboAsad[it].globalchannelid - (co<<11)-(as<<9)-(ag<<7);
    int where=co*NB_ASAD*NB_AGET*NB_CHANNEL + as*NB_AGET*NB_CHANNEL + ag*NB_CHANNEL + ch;

    if(co!=31){
      for(unsigned int hit=0;hit<EvtRed->CoboAsad[it].peakheight.size();hit++){
        if(EvtRed->CoboAsad[it].peakheight[hit]>10 && EvtRed->CoboAsad[it].peaktime[hit]>10){
		//cout << where << " " << TABLE[4][where] << " " << TABLE[5][where] << endl; 
		if(Hit[TABLE[4][where]][TABLE[5][where]]<2){
            		if(Hit[TABLE[4][where]+1][TABLE[5][where]]>0 || Hit[TABLE[4][where]-1][TABLE[5][where]]>0 || Hit[TABLE[4][where]][TABLE[5][where]+1]>0 || Hit[TABLE[4][where]][TABLE[5][where]-1]>0){
          			vX.push_back(TABLE[4][where]);
         			vY.push_back(TABLE[5][where]);
          			vZ.push_back(EvtRed->CoboAsad[it].peaktime[hit]);
          			vQ.push_back(EvtRed->CoboAsad[it].peakheight[hit]);
          			//hXYZ->Fill(TABLE[4][where],TABLE[5][where],EvtRed->CoboAsad[it].peaktime[hit],EvtRed->CoboAsad[it].peakheight[hit]);
          			//hXY->Fill(TABLE[4][where],TABLE[5][where],EvtRed->CoboAsad[it].peakheight[hit]);
          			//hXZ->Fill(TABLE[4][where],EvtRed->CoboAsad[it].peaktime[hit],EvtRed->CoboAsad[it].peakheight[hit]);
          			//hYZ->Fill(TABLE[5][where],EvtRed->CoboAsad[it].peaktime[hit],EvtRed->CoboAsad[it].peakheight[hit]);
			}
		}
        }
      }
    }

    else if(co==31){
      for(unsigned int hit=0;hit<EvtRed->CoboAsad[it].peakheight.size();hit++){
        int vxi_parameter = EvtRed->CoboAsad[it].peaktime[hit];
	if(Si_map[vxi_parameter]<21 && Si_map[vxi_parameter]>0){
           SiMayaRawEnergy.push_back(EvtRed->CoboAsad[it].peakheight[hit]);
       	   SiMayaNumber.push_back(Si_map[vxi_parameter]);
	}
      } 
    }  
  }
}

////////////////////////////////////////////////////////////////////////////////
void CleanPad()
{
	for(unsigned int it=0; it<EvtRed->CoboAsad.size(); it++){
	    int co=EvtRed->CoboAsad[it].globalchannelid>>11;
	    int as=(EvtRed->CoboAsad[it].globalchannelid - (co<<11))>>9;
	    int ag=(EvtRed->CoboAsad[it].globalchannelid - (co<<11)-(as<<9))>>7;
	    int ch=EvtRed->CoboAsad[it].globalchannelid - (co<<11)-(as<<9)-(ag<<7);
	    int where=co*NumberOfASAD*NumberOfAGET*NumberOfChannel + as*NumberOfAGET*NumberOfChannel + ag*NumberOfChannel + ch;
	

	    if(co!=31){
	      for(unsigned int hit=0;hit<EvtRed->CoboAsad[it].peakheight.size();hit++){
	        if(EvtRed->CoboAsad[it].peakheight[hit]>10 && EvtRed->CoboAsad[it].peaktime[hit]>1){
		        Hit[TABLE[4][where]][TABLE[5][where]] += 1;
	        }
	      }
	    }
  	}
}

////////////////////////////////////////////////////////////////////////////////
void GetMayaSiHitPosition(double xm, double xh, double ym, double yh, double zm, double zh)
{
  double l = xm-xh;
  double L = SiMayaDistanceX-xm;

  double t = (l+L)/l;
  //double t = L/l;

  double zf = zh + (zm-zh)*t;
  double yf = yh + (ym-yh)*t;

  SiMayaPosY.push_back(yf);
  SiMayaPosZ.push_back(zf);

  TVector3 BeamDirection = TVector3(1,0,0);
  TVector3 HitPos = TVector3(SiMayaDistanceX,yf,zf);
  ThetaLab.push_back(BeamDirection.Angle(HitPos)*180./TMath::Pi());
  PhiLab.push_back(HitPos.Phi()*180./TMath::Pi());
}

////////////////////////////////////////////////////////////////////////////////
void Clear(){
  vTrack.clear();
  vTrackLength.clear();
  vTrackCharge.clear();
  vTrackAngle.clear();
  vScalar.clear();
  mXYZ.clear();
  mQ.clear();
  vQ.clear();
  vX.clear();
  vY.clear();
  vZ.clear();
  TrackMult = -1;
  dedx.clear();
  ThetaLab.clear();
  PhiLab.clear();

  XVertex.clear();
  YVertex.clear();
  ZVertex.clear();

  SiMayaRawEnergy.clear();
  SiMayaNumber.clear();
  SiMayaPosY.clear();
  SiMayaPosZ.clear();

  for(int i=0; i<Configurator.GetNumberOfPadsX(); i++){
    for(int j=0; j<Configurator.GetNumberOfPadsY(); j++){
      Hit[i][j]=0;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void End(){
  cout << "Analysis finished" << endl;
}
