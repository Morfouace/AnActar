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

//TTree* Tree;
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

map<int, vector<int>> mXYZ;
map<int, vector<int>> mQ;


int PadNumberX;
int PadNumberY;
double PadSizeX;
double PadSizeY;
double DriftVelocity;

int TABLE[6][NumberOfCobo*NumberOfASAD*NumberOfAGET*NumberOfChannel];
int Hit[128][128];

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

// Si Maya parameters //;
double SiMayaDistanceX=256+58;
vector<double> SiMayaSignal;
vector<int> SiMayaNumber;
vector<double> SiMayaRawEnergy;
vector<double> SiMayaEnergy;
vector<double> SiMayaPosZ;
vector<double> SiMayaPosY;
vector<double> SiMayaTrackLength;

map<int, int> Si_map;



