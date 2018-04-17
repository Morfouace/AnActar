#include <iostream>
#include <fstream>
#include "TApplication.h"
#include <TROOT.h>

using namespace std;

#include"Analysis.h"

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  TApplication myApp("actaranalysis",&argc,argv);

  ReadRunToTreat("RunToTreat.txt");
  Init();

  time_t t ;
  time(&t);
  long double temps_initial = t;

  nentries = chain->GetEntries();
  for(int i=0;i<nentries;i++)
	{
    Clear();
    FillHistoAndMap(i);
    if(vX.size()>10 && SiMayaEnergy.size()>0){
      //TreatEventWithHough();ExcitationEnergy
      TreatEventWithSimpleRansac();

      for(unsigned int i=0; i<vTrack.size(); i++){
        //cout << vTrack.size() << " / "  << SiMayaEnergy.size() << endl;
        if(ThetaY[i]<70 || ThetaY[i]>110){
          for(unsigned int k=0; k<SiMayaEnergy.size(); k++){
            if(SiMayaEnergy[k]>0){
              int side = SiMayaSide[k];
              //cout << "yh = " << vTrack[k].GetYh() << endl;
              //cout << "ym = " << vTrack[k].GetYm() << endl;
              double x1 = vTrack[k].GetXm()*Configurator.GetPadSizeX();
              double x2 = vTrack[k].GetXh()*Configurator.GetPadSizeX();
              double y1 = vTrack[k].GetYm()*Configurator.GetPadSizeY();
              double y2 = vTrack[k].GetYh()*Configurator.GetPadSizeY();
              double z1 = vTrack[k].GetZm()*Configurator.GetDriftVelocity()*80e-3;
              double z2 = vTrack[k].GetZh()*Configurator.GetDriftVelocity()*80e-3;
              GetMayaSiHitPosition(side,x1,x2,y1,y2,z1,z2);
              //cout << "Maya distance= " << SiMayaDistanceY[side] << endl;
              double etot = EnergyLoss_He.EvaluateInitialEnergy(SiMayaEnergy[k]*MeV,SiMayaDistanceY[side],ThetaY[i]*deg);
              ELab.push_back(etot);
              SiMayaTrackLength.push_back(SiMayaDistanceY[side]/cos(ThetaY[i]/deg));

              double EBeamAtVertex = EnergyLoss_C.Slow(74*MeV,(55+XVertex[i])*mm,0);
              BeamEnergy.push_back(EBeamAtVertex);

              IneReaction->SetBeamEnergy(EBeamAtVertex);
              IneReaction->SetNuclei3(etot,TrackAngle[i]*deg);
              ExcitationEnergy.push_back(IneReaction->GetExcitation4());
            }
          }
        }
      }
    }
    OutputTree->Fill();
  }

  OutputFile->Write();
  OutputFile->Close();


  time(&t) ;
  long double temps_final = t;
  cout << "//****************//" << endl;
  cout << "Time: " << (int)((temps_final-temps_initial)/60) << " min " << ((int)(temps_final-temps_initial))%60 << " sec" << endl;

  End();

  myApp.Run(kTRUE);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
void ReadRunToTreat(string PathToRun){
  ifstream ifile;
  ifile.open(PathToRun.c_str());
  chain = new TChain("ACTAR_TTree");

  //TString root_path = "/space/morfouace/ACTAR/NSI-37/Data/";
  TString root_path = "/run/media/morfouace/ACTARTPC3/12C/NewFormat/";
  //TString root_path = "/home/morfouace/Physics/Actar/Analysis/actar_analysis/root/";
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
  // NPtool init //
  EnergyLoss_He = NPL::EnergyLoss("EnergyLossTable/Eloss_He_in_He-iC4H10.txt","SRIM",100);
  EnergyLoss_C = NPL::EnergyLoss("EnergyLossTable/Eloss_He_in_He-iC4H10.txt","SRIM",100);
  IneReaction = new NPL::Reaction("12C(4He,4He)12C@74");

  SiMayaDistanceY[0] = 37+63;
  //SiMayaDistanceY[1] = 37+39;
  SiMayaDistanceY[1] = 37+63;
  //SiMayaDistanceY[0] = 5+63;
  //SiMayaDistanceY[1] = 64+5+39;


  Configurator.ReadConfigurationFile("config/actar_config.txt");

  Ransac = new MSimpleRansac(Configurator.GetNumberOfPadsX(),Configurator.GetNumberOfPadsY(),0);

  //hXYZ = new TH3F("XYZ","XYZ",PadNumberX,0,PadNumberX,PadNumberY,0,PadNumberY,128,0,512);
  //hXY = new TH2F("XY","XY",PadNumberX,0,PadNumberX,PadNumberY,0,PadNumberY);
	//hXZ = new TH2F("XZ","XZ",PadNumberX,0,PadNumberX,128,0,512);
	//hYZ = new TH2F("YZ","YZ",PadNumberY,0,PadNumberY,128,0,512);

  string filename = "./dat/LT_nsi37.dat";
  ifstream ifile;
  ifile.open(filename.c_str());
 	cout << "Using LookupTable from: " << filename << endl;
  int cobo, asad, aget, channel;
  int padX, padY;
  while(!ifile.eof()){
    ifile >> cobo >> asad >> aget >> channel >> padY >> padX;
    int i=cobo*NumberOfASAD*NumberOfAGET*NumberOfChannel+asad*NumberOfAGET*NumberOfChannel+aget*NumberOfChannel+channel;
    TABLE[0][i] = cobo;
    TABLE[1][i] = asad;
    TABLE[2][i] = aget;
    TABLE[3][i] = channel;
    TABLE[5][i] = padY;
    TABLE[4][i] = padX;
  }
 	ifile.close();

  string calib_filename = "./calib/calib_total_bl_corr.txt";
  ifstream calibfile;
  calibfile.open(calib_filename.c_str());
  cout << "Using calib file from: " << calib_filename << endl;
  double slope, offset;
  while(!calibfile.eof()){
    calibfile >> cobo >> asad >> aget >> channel >> slope >> offset;
    int i=cobo*NumberOfASAD*NumberOfAGET*NumberOfChannel+asad*NumberOfAGET*NumberOfChannel+aget*NumberOfChannel+channel;
    p1[i] = slope;
    p0[i] = offset;
  }

  //Tree->SetBranchAddress("data",&EvtRed);
  chain->SetBranchAddress("data",&EvtRed);

  InitOutputTree();
}

////////////////////////////////////////////////////////////////////////////////
void InitOutputTree(){
  OutputFile = new TFile("/space/morfouace/ACTAR/NSI-37/root-files/FittedTracks.root","RECREATE");
  OutputTree = new TTree("PhysicsTree","PhysicsTree");

  OutputTree->Branch("Event",&Event,"Event/I");
  OutputTree->Branch("TrackMult",&TrackMult,"TrackMult/I");
  OutputTree->Branch("TrackAngle",&TrackAngle);
  OutputTree->Branch("ThetaY",&ThetaY);
  OutputTree->Branch("TrackCharge",&vTrackCharge);
  OutputTree->Branch("TrackLength",&vTrackLength);
  OutputTree->Branch("dedx",&dedx);
  OutputTree->Branch("XVertex",&XVertex);
  OutputTree->Branch("YVertex",&YVertex);
  OutputTree->Branch("ZVertex",&ZVertex);
  OutputTree->Branch("X",&vX);
  OutputTree->Branch("Y",&vY);
  OutputTree->Branch("Z",&vZ);
  OutputTree->Branch("DSSDNumber",&DSSDNumber);
  OutputTree->Branch("DSSDStripNumber",&DSSDStripNumber);
  OutputTree->Branch("DSSDStripFront",&DSSDStripFront);
  OutputTree->Branch("DSSDStripBack",&DSSDStripBack);
  OutputTree->Branch("DSSDRawEnergy",&DSSDRawEnergy);
  OutputTree->Branch("DSSDEnergyFront",&DSSDEnergyFront);
  OutputTree->Branch("DSSDEnergyBack",&DSSDEnergyBack);
  OutputTree->Branch("Y_DSSD",&Y_DSSD);
  OutputTree->Branch("Z_DSSD",&Z_DSSD);
  OutputTree->Branch("SiMayaSide",&SiMayaSide);
  OutputTree->Branch("SiMayaNumber",&SiMayaNumber);
  OutputTree->Branch("SiMayaRawEnergy",&SiMayaRawEnergy);
  OutputTree->Branch("SiMayaEnergy",&SiMayaEnergy);
  OutputTree->Branch("SiMayaPosX",&SiMayaPosX);
  OutputTree->Branch("SiMayaPosZ",&SiMayaPosZ);
  OutputTree->Branch("SiMayaTrackLength",&SiMayaTrackLength);
  OutputTree->Branch("BeamEnergy",&BeamEnergy);
  OutputTree->Branch("ELab",&ELab);
  OutputTree->Branch("ExcitationEnergy",&ExcitationEnergy);
}

////////////////////////////////////////////////////////////////////////////////
void TreatEventWithSimpleRansac(){
  Ransac->Init(vX,vY,vZ,vQ);

  vTrack = Ransac->SimpleRansac();
  vTrackCharge = Ransac->GetChargeOfTracks();
  vTrackLength = Ransac->GetTrackLength(Configurator.GetPadSizeX(), Configurator.GetPadSizeY(), 80e-3*Configurator.GetDriftVelocity());
  //cout << "//********* " << vTrack.size() << " / " << vTrackCharge.size() << " / " << vTrackLength.size() << endl;

  TVector3 vBeam = TVector3(1,0,0);
  TVector3 vBeamPos = TVector3(0,16,240);
  TVector3 aTrack;
  TrackMult = vTrack.size();
  for(unsigned int i=0; i<vTrack.size(); i++){
    //aTrack = TVector3(vTrack[i].X()*2, vTrack[i].Y()*2, vTrack[i].Z()*80e-3*Configurator.GetDriftVelocity());
    double Xdir = vTrack[i].GetDirectionVector().X();
    double Ydir = vTrack[i].GetDirectionVector().Y();
    double Zdir = vTrack[i].GetDirectionVector().Z();
    aTrack = TVector3(Xdir*2, Ydir*2, Zdir*80e-3*Configurator.GetDriftVelocity());

    XVertex.push_back(vTrack[i].GetVertexPostion(vBeam,vBeamPos).X());
    YVertex.push_back(vTrack[i].GetVertexPostion(vBeam,vBeamPos).Y());
    ZVertex.push_back(vTrack[i].GetVertexPostion(vBeam,vBeamPos).Z());

    double angle = vBeam.Angle(aTrack)*180/TMath::Pi();
    if(angle>90) angle = 180-angle;
    TrackAngle.push_back(angle);
    ThetaY.push_back(TVector3(0,1,0).Angle(aTrack)*180/TMath::Pi());

    if(vTrackLength[i]>0){
      dedx.push_back(vTrackCharge[i]/vTrackLength[i]);
    }
    else dedx.push_back(-100);
  }
}

////////////////////////////////////////////////////////////////////////////////
void TreatEventWithHough(){
  Hough.Reset();

  cout << "/// map.size() = " << mXYZ.size() << " ///" << endl;
  //Hough.Init(hXY,hYZ,mXYZ,mQ);
}

////////////////////////////////////////////////////////////////////////////////
void Clear(){
  mXYZ.clear();
  mQ.clear();
  vTrack.clear();
  vTrackCharge.clear();
  vTrackLength.clear();
  dedx.clear();
  TrackAngle.clear();
  ThetaY.clear();
  XVertex.clear();
  YVertex.clear();
  ZVertex.clear();
  vQ.clear();
  vX.clear();
  vY.clear();
  vZ.clear();

  DSSDNumber.clear();
  DSSDStripNumber.clear();
  DSSDStripBack.clear();
  DSSDStripFront.clear();
  DSSDEnergyFront.clear();
  DSSDEnergyBack.clear();
  DSSDRawEnergy.clear();
  DSSDSignal.clear();
  Y_DSSD.clear();
  Z_DSSD.clear();

  TrackMult = -1;

  // Intermediate variable
  DetNumber.clear();
  LTNumber.clear();
  StripFront.clear();
  StripBack.clear();

  SiMayaSignal.clear();
  SiMayaSide.clear();
  SiMayaNumber.clear();
  SiMayaRawEnergy.clear();
  SiMayaEnergy.clear();
  SiMayaPosX.clear();
  SiMayaPosZ.clear();
  SiMayaTrackLength.clear();
  LTMayaNumber.clear();

  BeamEnergy.clear();
  ELab.clear();
  ExcitationEnergy.clear();

  for(int i=0; i<Configurator.GetNumberOfPadsX(); i++){
    for(int j=0; j<Configurator.GetNumberOfPadsY(); j++){
      Hit[i][j]=0;
    }
  }

  /*hXYZ->Reset();
  hXY->Reset();
  hYZ->Reset();
  hYZ->Reset();*/
}
////////////////////////////////////////////////////////////////////////////////
void FillHistoAndMap(int i){

  chain->GetEntry(i);
  Event = i;
  if(i%1000==0){
    cout << "//// Processing... : " << 100*((double)i)/nentries << "% done //// \r" << flush;
  }

  for(unsigned int it=0; it<EvtRed->CoboAsad.size(); it++){
    int co=EvtRed->CoboAsad[it].globalchannelid>>11;
    int as=(EvtRed->CoboAsad[it].globalchannelid - (co<<11))>>9;
    int ag=(EvtRed->CoboAsad[it].globalchannelid - (co<<11)-(as<<9))>>7;
    int ch=EvtRed->CoboAsad[it].globalchannelid - (co<<11)-(as<<9)-(ag<<7);
    int where=co*NumberOfASAD*NumberOfAGET*NumberOfChannel + as*NumberOfAGET*NumberOfChannel + ag*NumberOfChannel + ch;


    if(co<2){
      for(unsigned int hit=0;hit<EvtRed->CoboAsad[it].peakheight.size();hit++){
        if(EvtRed->CoboAsad[it].peakheight[hit]>200 && EvtRed->CoboAsad[it].peaktime[hit]<450 && EvtRed->CoboAsad[it].peaktime[hit]>50){
          Hit[TABLE[4][where]][TABLE[5][where]] += 1;
        }
      }
    }
  }

  for(unsigned int it=0; it<EvtRed->CoboAsad.size(); it++){
    int co=EvtRed->CoboAsad[it].globalchannelid>>11;
    int as=(EvtRed->CoboAsad[it].globalchannelid - (co<<11))>>9;
    int ag=(EvtRed->CoboAsad[it].globalchannelid - (co<<11)-(as<<9))>>7;
    int ch=EvtRed->CoboAsad[it].globalchannelid - (co<<11)-(as<<9)-(ag<<7);
    int where=co*NumberOfASAD*NumberOfAGET*NumberOfChannel + as*NumberOfAGET*NumberOfChannel + ag*NumberOfChannel + ch;

    if(co<2){
      for(unsigned int hit=0;hit<EvtRed->CoboAsad[it].peakheight.size();hit++){
        if(EvtRed->CoboAsad[it].peakheight[hit]>200 && EvtRed->CoboAsad[it].peaktime[hit]<450 && EvtRed->CoboAsad[it].peaktime[hit]>50){
          //Hit[TABLE[4][where]][TABLE[5][where]] += 1;
          if(Hit[TABLE[4][where]][TABLE[5][where]]<2){
            if(Hit[TABLE[4][where]+1][TABLE[5][where]]>0 || Hit[TABLE[4][where]-1][TABLE[5][where]] || Hit[TABLE[4][where]][TABLE[5][where]+1] || Hit[TABLE[4][where]][TABLE[5][where]-1]){
              vX.push_back(TABLE[4][where]);
              vY.push_back(TABLE[5][where]);
              vZ.push_back(EvtRed->CoboAsad[it].peaktime[hit]);
              double Q = p0[where] + p1[where]*EvtRed->CoboAsad[it].peakheight[hit];
              vQ.push_back(Q);
              //vQ.push_back(EvtRed->CoboAsad[it].peakheight[hit]);
              //hXYZ->Fill(TABLE[4][where],TABLE[5][where],EvtRed->CoboAsad[it].peaktime[hit],EvtRed->CoboAsad[it].peakheight[hit]);
              //hXY->Fill(TABLE[4][where],TABLE[5][where],EvtRed->CoboAsad[it].peakheight[hit]);
              //hXZ->Fill(TABLE[4][where],EvtRed->CoboAsad[it].peaktime[hit],EvtRed->CoboAsad[it].peakheight[hit]);
              //hYZ->Fill(TABLE[5][where],EvtRed->CoboAsad[it].peaktime[hit],EvtRed->CoboAsad[it].peakheight[hit]);
            }
          }
        }
      }
    }

    else if(co==2){
      for(unsigned int hit=0;hit<EvtRed->CoboAsad[it].peakheight.size();hit++){
        DSSDSignal.push_back(EvtRed->CoboAsad[it].peakheight[hit]);
      }
      LTNumber.push_back(where);
      DetNumber.push_back(TABLE[5][where]);
      DSSDStripNumber.push_back(TABLE[4][where]);
      if(TABLE[4][where]<32) StripFront.push_back(TABLE[4][where]);
      else StripBack.push_back(TABLE[4][where]-32);
    }

    else if(co==3 && as==0){
      for(unsigned int hit=0;hit<EvtRed->CoboAsad[it].peakheight.size();hit++){
        SiMayaSignal.push_back(EvtRed->CoboAsad[it].peakheight[hit]);
      }
      LTMayaNumber.push_back(where);
      SiMayaSide.push_back(TABLE[5][where]);
      SiMayaNumber.push_back(TABLE[4][where]);
    }
  }
  if(DSSDSignal.size()>0){
    TreatDSSDSignal(DSSDSignal);
    ApplyCalibration();
  }
  GetDSSDImpact();

  if(SiMayaSignal.size()>0){
    TreatSiMayaSignal(SiMayaSignal);
    ApplySiMayaCalibration();
  }

  /*int p=0;
  for(int iY=0; iY<PadNumberY; iY++){
    for(int iZ=0; iZ<128; iZ++){
      if(hYZ->GetBinContent(iY,iZ)>0){
        p++;
        int iYZ = iY + 128*iZ;
        for(int iX=0; iX<PadNumberX; iX++){
          if(hXY->GetBinContent(iX,iY)>0 && hXZ->GetBinContent(iX,iZ)>0){
            int iQ = hXYZ->GetBinContent(iX,iY,iZ);
            vX.push_back(iX);
            vQ.push_back(iQ);
          }
        }
        mXYZ.insert(std::pair<int, vector<int>>(iYZ,vX));
        mQ.insert(std::pair<int, vector<int>>(iYZ,vQ));
        vX.clear();
        vQ.clear();
      }
    }
  }*/
}

////////////////////////////////////////////////////////////////////////////////
void TreatDSSDSignal(vector<double> v1)
{
  int NumberOfHit = v1.size()/512;
  double Baseline=0;
  for(int i=0; i<100; i++){
    Baseline += v1[i];
  }
  Baseline = Baseline/100;

  vector<double> vHit;
  for(int k=0; k<NumberOfHit; k++){
    vHit.clear();
    for(int i=0;i<512;i++){
      double adc_value = v1[k*512+i]-Baseline;
      vHit.push_back(adc_value);
    }

    double max = *max_element(vHit.begin(), vHit.end());
    DSSDRawEnergy.push_back(max);
  }
}

////////////////////////////////////////////////////////////////////////////////
void TreatSiMayaSignal(vector<double> v1)
{
  int NumberOfHit = v1.size()/512;
  double Baseline=0;
  for(int i=0; i<100; i++){
    Baseline += v1[i];
  }
  Baseline = Baseline/100;

  vector<double> vHit;
  for(int k=0; k<NumberOfHit; k++){
    vHit.clear();
    for(int i=0;i<512;i++){
      double adc_value = v1[k*512+i]-Baseline;
      vHit.push_back(adc_value);
    }
    double max = *max_element(vHit.begin(), vHit.end());
    SiMayaRawEnergy.push_back(max);
  }
}

////////////////////////////////////////////////////////////////////////////////
void ApplySiMayaCalibration(){
  for(int i=0; i<SiMayaRawEnergy.size(); i++){
    double Value = p0[LTMayaNumber[i]] + SiMayaRawEnergy[i]*p1[LTMayaNumber[i]];
    SiMayaEnergy.push_back(Value);
  }
}
////////////////////////////////////////////////////////////////////////////////
void ApplyCalibration(){
  vector<double> EFront;
  vector<double> EBack;
  EFront.clear();
  EBack.clear();
  for(int i=0; i<DSSDRawEnergy.size(); i++){
      double Value = p0[LTNumber[i]] + DSSDRawEnergy[i]*p1[LTNumber[i]];
      if(DSSDStripNumber[i]<32)EFront.push_back(Value);
      else EBack.push_back(Value);
  }

  if(EBack.size()==EFront.size()){
    for(int i=0; i<EFront.size(); i++){
      DSSDNumber.push_back(DetNumber[i]);

      DSSDEnergyBack.push_back(EBack[i]);
      DSSDStripBack.push_back(StripBack[i]);

      DSSDEnergyFront.push_back(EFront[i]);
      DSSDStripFront.push_back(StripFront[i]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void GetMayaSiHitPosition(int side, double xm, double xh, double ym, double yh, double zm, double zh)
{
  double x1, x2;
  double y1, y2;
  double z1, z2;


  if(ym>yh && side==1){
    x1 = xh;
    y1 = yh;
    z1 = zh;
    x2 = xm;
    y2 = ym;
    z2 = zm;
  }
  else if(ym<yh && side==1){
    x1 = xm;
    y1 = ym;
    z1 = zm;
    x2 = xh;
    y2 = yh;
    z2 = zh;
  }
  if(ym>yh && side==0){
    x1 = xm;
    y1 = ym;
    z1 = zm;
    x2 = xh;
    y2 = yh;
    z2 = zh;
  }
  else if(ym<yh && side==0){
    x1 = xh;
    y1 = yh;
    z1 = zh;
    x2 = xm;
    y2 = ym;
    z2 = zm;
  }

  double l = abs(ym-yh);
  double L;
  if(side==1) L = abs(-SiMayaDistanceY[side]-yh);
  else L = abs(SiMayaDistanceY[side]-yh);

  //double t = (l+L)/l;
  double t = L/l;

  double zf = zh + (zm-zh)*t;
  double xf = xh + (xm-xh)*t;

  SiMayaPosX.push_back(xf);
  SiMayaPosZ.push_back(zf);
}

////////////////////////////////////////////////////////////////////////////////
void GetDSSDImpact(){
  double Zval;
  double Yval;
  double StripWidth = 2;
  for(unsigned int i=0; i<DSSDStripFront.size(); i++){
    if(DSSDNumber[i]==0){
      Zval = 4 + 0.5*StripWidth + (31-DSSDStripBack[i])*StripWidth;
      Yval = -4 - 0.5*StripWidth - DSSDStripFront[i]*StripWidth;
      Y_DSSD.push_back(Yval);
      Z_DSSD.push_back(Zval);
    }
    else if(DSSDNumber[i]==1){
      Zval = -4 - 0.5*StripWidth - DSSDStripBack[i]*StripWidth;
      Yval = -4 - 0.5*StripWidth - DSSDStripFront[i]*StripWidth;
      Y_DSSD.push_back(Yval);
      Z_DSSD.push_back(Zval);
    }
    else if(DSSDNumber[i]==2){
      Zval = 4 + 0.5*StripWidth + (31-DSSDStripBack[i])*StripWidth;
      Yval = 4 + 0.5*StripWidth + (31-DSSDStripFront[i])*StripWidth;
      Y_DSSD.push_back(Yval);
      Z_DSSD.push_back(Zval);
    }
    else if(DSSDNumber[i]==3){
      Zval = -4 - 0.5*StripWidth - DSSDStripBack[i]*StripWidth;
      Yval = 4 + 0.5*StripWidth + (31-DSSDStripFront[i])*StripWidth;
      Y_DSSD.push_back(Yval);
      Z_DSSD.push_back(Zval);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void End(){
  cout << "Analysis finished" << endl;
}
