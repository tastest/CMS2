#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TFile.h" 
#include "../../../Smurf/Core/SmurfTree.h"
#include "LeptonTreeMaker.h"
#include "SmurfDataTypes.h"
 
void processLeptonTree(TString outfileid, enum SmurfTree::DataType sample, TString file, bool realData, TString goodrunlist);

