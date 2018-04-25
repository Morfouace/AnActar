#include <map>

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
#include "TVector3.h"
#include "TCanvas.h"

#include "NPEnergyLoss.h"
#include "NPReaction.h"

#define NumberOfCobo 16
#define NumberOfASAD 4
#define NumberOfAGET 4
#define NumberOfChannel 68

void ReadRunToTreat(string);
void Init();
void InitOutputTree();
void FillHistoAndMap(int i);
void TreatEventWithSimpleRansac();
void TreatEventWithHough();
void CleanPad();
void Clear();
void End();
void GetMayaSiHitPosition(double xm, double xh, double ym, double yh, double zm, double zh);
void ApplySiMayaCalibration();
bool IsGoingToSilicon(int TrackID);
bool IsGoodSilicon(int ID);


TChain* chain;
TTree* OutputTree;
TFile* OutputFile;
/*TH3F* hXYZ;
TH2F* hXY;
TH2F* hXZ;
TH2F* hYZ;*/

MEventReduced* EvtRed;
MSimpleRansac* Ransac;
MHough Hough;
MDetectorConfig Configurator;

TString Base_Path="/space/morfouace/ACTAR/e750/root-files/FittedTree.root";
bool bRansacVisu=false;

map<int, vector<int>> mXYZ;
map<int, vector<int>> mQ;

double BeamAngle;


int PadNumberX;
int PadNumberY;
double PadSizeX;
double PadSizeY;
double DriftVelocity;

int TABLE[6][NumberOfCobo*NumberOfASAD*NumberOfAGET*NumberOfChannel];
int Hit[128][128];
double p0[20];
double p1[20];

vector<int> vX;
vector<int> vY;
vector<double> vZ;
vector<double> vQ;
vector<double> vTrackLength;
vector<double> vTrackCharge;
vector<double> vTrackAngle;
vector<double> vScalar;
vector<MTrack> vTrack;
int TrackMult;
int Event;
vector<double> dedx;
vector<double> ThetaLab;
vector<double> PhiLab;

vector<double> XVertex;
vector<double> YVertex;
vector<double> ZVertex;

vector<double> PositionOfInteractionX;
vector<double> PositionOfInteractionY;
vector<double> PositionOfInteractionZ;

// Si Maya parameters //;
double SiMayaDistanceX=256+58;
vector<double> SiMayaSignal;
vector<int> SiMayaNumber;
vector<double> SiMayaRawEnergy;
vector<double> SiMayaEnergy;
vector<double> SiMayaPosZ;
vector<double> SiMayaPosY;
vector<double> YMaya;
vector<double> ZMaya;
vector<double> SiMayaTrackLength;

map<int, int> Si_map;

// Physical variables //
vector<double> BeamEnergy;
vector<double> ESi;
vector<double> ELab;
vector<double> ExcitationEnergy;
vector<double> Ecm;

// NPTool //
NPL::EnergyLoss EnergyLoss_Ne;
NPL::EnergyLoss EnergyLoss_p;
NPL::Reaction* TheReaction;
