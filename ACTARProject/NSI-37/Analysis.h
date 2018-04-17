// ACTAR Header
#include "MSimpleRansac.h"
#include "MTrack.h"
#include "MHough.h"
#include "MDetectorConfig.h"
#include "MEventReduced.h"

// ROOT
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TCanvas.h"

#include "NPEnergyLoss.h"
#include "NPReaction.h"

#define NumberOfCobo 4
#define NumberOfASAD 4
#define NumberOfAGET 4
#define NumberOfChannel 64

void ReadRunToTreat(string);
void Init();
void Clear();
void InitOutputTree();
void FillHistoAndMap(int i);
void TreatEventWithSimpleRansac();
void TreatEventWithHough();
void TreatDSSDSignal(vector<double> v1);
void TreatSiMayaSignal(vector<double> v1);
void ApplyCalibration();
void ApplySiMayaCalibration();
void GetDSSDImpact();
void GetMayaSiHitPosition(int side, double xm, double xh, double ym, double yh, double zm, double zh);
void End();

//TTree* Tree;
TChain* chain;
TTree* OutputTree;
TFile* OutputFile;
//TH3F* hXYZ;
//TH2F* hXY;
//TH2F* hXZ;
//TH2F* hYZ;

int TABLE[6][NumberOfCobo*NumberOfASAD*NumberOfAGET*NumberOfChannel];
int Hit[32][64];
double p1[NumberOfCobo*NumberOfASAD*NumberOfAGET*NumberOfChannel];
double p0[NumberOfCobo*NumberOfASAD*NumberOfAGET*NumberOfChannel];

MEventReduced* EvtRed;
MSimpleRansac* Ransac;
MHough Hough;
MDetectorConfig Configurator;

map<int, vector<int>> mXYZ;
map<int, vector<int>> mQ;

// General Parameters //
int Event;
int nentries;

// TPC Paaramters //
vector<int> vX, vY;
vector<double> vZ, vQ;
int TrackMult;
vector<MTrack> vTrack;
vector<double> vTrackCharge;
vector<double> vTrackLength;
vector<double> TrackAngle;
vector<double> ThetaY;
vector<double> XVertex;
vector<double> YVertex;
vector<double> ZVertex;
vector<double> dedx;

// DSSD Parameters //
vector<int> LTNumber;
vector<int> StripFront;
vector<int> StripBack;
vector<int> DetNumber;

vector<int> DSSDStripNumber;
vector<int> DSSDStripFront;
vector<int> DSSDStripBack;
vector<int> DSSDNumber;
vector<double> DSSDRawEnergy;
vector<double> DSSDEnergyFront;
vector<double> DSSDEnergyBack;
vector<double> DSSDSignal;
vector<double> Z_DSSD;
vector<double> Y_DSSD;

// Si Maya parameters //
double SiMayaDistanceY[2];
vector<double> SiMayaSignal;
vector<int> SiMayaSide;
vector<int> SiMayaNumber;
vector<double> SiMayaRawEnergy;
vector<double> SiMayaEnergy;
vector<double> SiMayaPosX;
vector<double> SiMayaPosZ;
vector<double> SiMayaTrackLength;
vector<int> LTMayaNumber;

vector<double> BeamEnergy;
vector<double> ELab;
vector<double> ExcitationEnergy;

// NPTool object //
NPL::EnergyLoss EnergyLoss_He;
NPL::EnergyLoss EnergyLoss_C;
NPL::Reaction* IneReaction;
