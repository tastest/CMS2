#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "wwtypes.cc"
#include <assert.h>
#include <math.h>
#include "TCut.h"
#include "TDirectory.h"
#include "TROOT.h"
#include <sstream>
#include <stdarg.h>

// modes
// 0 - text
// 1 - twiki
// 2 - latex
int mode;
int precision;

// ===============================================================
// WARNING! Don't use plain Form from TString unless you have to.
// It's buggy.
//
const char* pm = 0;
const char* leftDivider = 0;
const char* centerDivider = 0;
const char* rightDivider = 0;
void init(){
  switch (mode){
  case 0:
    pm = "+/-";
    leftDivider = "|";
    centerDivider = "|";
    rightDivider = "|";
    break;
  case 1:
    pm = "&plusmn;";
    leftDivider = "|";
    centerDivider = "|";
    rightDivider = "|";
    break;
  case 2:
    pm = "$\\pm$";
    leftDivider = "";
    centerDivider = "&";
    rightDivider = "\\\\";
    break;
  }
}

string MForm(const char* pattern, ...){
  char buffer[1024];
  va_list args;
  va_start(args,pattern);
  vsprintf(buffer,pattern,args);
  va_end(args);
  return string(buffer);
}

string mround1(double x){
  char buffer [50];
  sprintf(buffer,"%%%d.%df",3+precision,precision);
  string patter(buffer);
  sprintf(buffer,patter.c_str(),x);
  return string(buffer);
}
string mround2(double x){
  char buffer [50];
  sprintf(buffer,"%%%d.%df",2+precision,precision);
  string patter(buffer);
  sprintf(buffer,patter.c_str(),x);
  return string(buffer);
}

TH1F* getYields(const char* file, const char* cut = 0, const char* tree = "tree")
{
  TH1F* output(0);
  TFile *f = dynamic_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(file));
  if (!f ) f = TFile::Open(file);
  if (!f) return output;
  TTree* tt = dynamic_cast<TTree*>(f->Get(tree));
  if (!tt) return output;
  if (cut)
    tt->Draw("type>>h(4,0,4)",MForm("scale1fb*(%s)",cut).c_str(),"e goff");
  else
    tt->Draw("type>>h(4,0,4)","scale1fb","e goff");
  output = dynamic_cast<TH1F*>(gDirectory->Get("h"));
  assert(output);
  output->SetDirectory(0);
  return output;
}

double getValue(const TH1F* hist, unsigned int type){
  switch (type){
  case 0:
    return hist->GetBinContent(1);
    break;
  case 1:
    return hist->GetBinContent(2);
    break;
  case 2:
    return hist->GetBinContent(3);
    break;
  case 3:
    return hist->GetBinContent(4);
    break;
  default:
    return hist->Integral(1,4);
  }
}
 
double getError(const TH1F* hist, unsigned int type){
  switch (type){
  case 0:
    return hist->GetBinError(1);
    break;
  case 1:
    return hist->GetBinError(2);
    break;
  case 2:
    return hist->GetBinError(3);
    break;
  case 3:
    return hist->GetBinError(4);
    break;
  default:
    return sqrt(pow(hist->GetBinError(1),2) + pow(hist->GetBinError(2),2)+
		pow(hist->GetBinError(3),2) + pow(hist->GetBinError(4),2));
  }
}

void showResults2(const char* cut = 0, int imode=0, int iprecision=1, bool showTotalBkg=false)
{
  mode = imode;
  precision = iprecision;
  using namespace std;
  init();
  std::vector<std::pair<TH1F*,std::string> > bkgs;
  if ( TH1F *hist  = getYields("smurf/dyee.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*DY ee*"));
  if ( TH1F *hist  = getYields("smurf/dymm.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*DY mumu*"));
  if ( TH1F *hist  = getYields("smurf/dytt.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*DY tautau*"));
  if ( TH1F *hist  = getYields("smurf/ttbar.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*ttbar*"));
  if ( TH1F *hist  = getYields("smurf/tw.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*TW*"));
  if ( TH1F *hist  = getYields("smurf/wjets.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*Wjets*"));
  if ( TH1F *hist  = getYields("smurf/wz.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*WZ*"));
  if ( TH1F *hist  = getYields("smurf/zz.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*ZZ*"));
  if ( TH1F *hist  = getYields("smurf/ggww.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*ggWW*"));
  if ( TH1F *hist  = getYields("smurf/ww.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*qqWW*"));

  const TH1F *hww120 = getYields("smurf/hww120.root",cut);
  const TH1F *hww130 = getYields("smurf/hww130.root",cut);
  const TH1F *hww140 = getYields("smurf/hww140.root",cut);
  const TH1F *hww150 = getYields("smurf/hww150.root",cut);
  const TH1F *hww160 = getYields("smurf/hww160.root",cut);
  const TH1F *hww170 = getYields("smurf/hww170.root",cut);
  const TH1F *hww180 = getYields("smurf/hww180.root",cut);
  const TH1F *hww190 = getYields("smurf/hww190.root",cut);
  const TH1F *hww200 = getYields("smurf/hww200.root",cut);
  const TH1F *hww250 = getYields("smurf/hww250.root",cut);
  const TH1F *hww300 = getYields("smurf/hww250.root",cut);
  const TH1F *vbfhww120 = getYields("smurf/vbfhww120.root",cut);
  const TH1F *vbfhww130 = getYields("smurf/vbfhww130.root",cut);
  const TH1F *vbfhww140 = getYields("smurf/vbfhww140.root",cut);
  const TH1F *vbfhww150 = getYields("smurf/vbfhww150.root",cut);
  const TH1F *vbfhww160 = getYields("smurf/vbfhww160.root",cut);
  const TH1F *vbfhww170 = getYields("smurf/vbfhww170.root",cut);
  const TH1F *vbfhww180 = getYields("smurf/vbfhww180.root",cut);
  const TH1F *vbfhww190 = getYields("smurf/vbfhww190.root",cut);
  const TH1F *vbfhww200 = getYields("smurf/vbfhww200.root",cut);
  const TH1F *vbfhww250 = getYields("smurf/vbfhww250.root",cut);
  const TH1F *vbfhww300 = getYields("smurf/vbfhww250.root",cut);
  const TH1F *data   = getYields("smurf/data.root",cut);

  const char* patternTitle = Form(" %%%ds %s",8+precision*2,centerDivider);
  const char* patternData  = Form(" %%%d.%df%s%%%d.%df %s",3+precision,precision,pm,2+precision,precision,centerDivider);

  cout << "\n" << MForm("%s %3s %s",leftDivider,"",centerDivider);
  for (unsigned int i=0; i<bkgs.size(); ++i) cout << MForm(patternTitle,bkgs.at(i).second.c_str());
  if ( showTotalBkg && bkgs.size()>0 ) cout << MForm(patternTitle,"*Total BKG*");
  cout << endl;

  double bkg[5] = {0, 0, 0, 0, 0};
  double bkgerr2[5] = {0, 0, 0, 0, 0};
  const char* names[5] = {"mm","me","em","ee","all"};
  for (int i=0; i<5; i++){
    cout << leftDivider << MForm(" %3s ",names[i]) << centerDivider;
    for (unsigned int j=0; j<bkgs.size(); ++j){
      TH1F* hist = bkgs.at(j).first;
      cout << MForm(patternData,getValue(hist,i),getError(hist,i));
      bkg[i]     += getValue(hist,i);
      bkgerr2[i] += pow(getError(hist,i),2);
    }
    if ( showTotalBkg && bkgs.size()>0 ) cout << MForm(patternData,bkg[i],sqrt(bkgerr2[i]));
    cout <<endl;
  }
  cout <<endl;

  cout << "\n" << MForm("%s %3s %s",leftDivider,"",centerDivider);
  if (hww120) cout << MForm(patternTitle,"*HWW120*");
  if (hww130) cout << MForm(patternTitle,"*HWW130*");
  if (hww140) cout << MForm(patternTitle,"*HWW140*");
  if (hww150) cout << MForm(patternTitle,"*HWW150*");
  if (hww160) cout << MForm(patternTitle,"*HWW160*");
  if (hww170) cout << MForm(patternTitle,"*HWW170*");
  if (hww180) cout << MForm(patternTitle,"*HWW180*");
  if (hww190) cout << MForm(patternTitle,"*HWW190*");
  if (hww200) cout << MForm(patternTitle,"*HWW200*");
  if (hww250) cout << MForm(patternTitle,"*HWW250*");
  cout << endl;
  flush(cout);

  for (int i=0; i<5; i++){
    cout << leftDivider << MForm(" %3s ",names[i]) << centerDivider;
    if (hww120) cout << MForm(patternData,getValue(hww120,i),getError(hww120,i));
    if (hww130) cout << MForm(patternData,getValue(hww130,i),getError(hww130,i));
    if (hww140) cout << MForm(patternData,getValue(hww140,i),getError(hww140,i));
    if (hww150) cout << MForm(patternData,getValue(hww150,i),getError(hww150,i));
    if (hww160) cout << MForm(patternData,getValue(hww160,i),getError(hww160,i));
    if (hww170) cout << MForm(patternData,getValue(hww170,i),getError(hww170,i));
    if (hww180) cout << MForm(patternData,getValue(hww180,i),getError(hww180,i));
    if (hww190) cout << MForm(patternData,getValue(hww190,i),getError(hww190,i));
    if (hww200) cout << MForm(patternData,getValue(hww200,i),getError(hww200,i));
    if (hww250) cout << MForm(patternData,getValue(hww250,i),getError(hww250,i));
    cout <<endl;
  }
  cout <<endl;

  cout << "\n" << MForm("%s %3s %s",leftDivider,"",centerDivider);
  if (vbfhww120) cout << MForm(patternTitle,"*VBFH120*");
  if (vbfhww130) cout << MForm(patternTitle,"*VBFH130*");
  if (vbfhww140) cout << MForm(patternTitle,"*VBFH140*");
  if (vbfhww150) cout << MForm(patternTitle,"*VBFH150*");
  if (vbfhww160) cout << MForm(patternTitle,"*VBFH160*");
  if (vbfhww170) cout << MForm(patternTitle,"*VBFH170*");
  if (vbfhww180) cout << MForm(patternTitle,"*VBFH180*");
  if (vbfhww190) cout << MForm(patternTitle,"*VBFH190*");
  if (vbfhww200) cout << MForm(patternTitle,"*VBFH200*");
  if (vbfhww250) cout << MForm(patternTitle,"*VBFH250*");
  cout << endl;

  for (int i=0; i<5; i++){
    cout.precision(precision);
    cout << leftDivider << MForm(" %3s ",names[i]) << centerDivider;
    if (vbfhww120) cout << MForm(patternData,getValue(vbfhww120,i),getError(vbfhww120,i));
    if (vbfhww130) cout << MForm(patternData,getValue(vbfhww130,i),getError(vbfhww130,i));
    if (vbfhww140) cout << MForm(patternData,getValue(vbfhww140,i),getError(vbfhww140,i));
    if (vbfhww150) cout << MForm(patternData,getValue(vbfhww150,i),getError(vbfhww150,i));
    if (vbfhww160) cout << MForm(patternData,getValue(vbfhww160,i),getError(vbfhww160,i));
    if (vbfhww170) cout << MForm(patternData,getValue(vbfhww170,i),getError(vbfhww170,i));
    if (vbfhww180) cout << MForm(patternData,getValue(vbfhww180,i),getError(vbfhww180,i));
    if (vbfhww190) cout << MForm(patternData,getValue(vbfhww190,i),getError(vbfhww190,i));
    if (vbfhww200) cout << MForm(patternData,getValue(vbfhww200,i),getError(vbfhww200,i));
    if (vbfhww250) cout << MForm(patternData,getValue(vbfhww250,i),getError(vbfhww250,i));
//     if (vbfhww120) cout << centerDivider << " " << mround1(getValue(vbfhww120,i)) << " " << pm << " " << mround2(getError(vbfhww120,i)) << " " << centerDivider;
//     if (vbfhww130) cout << centerDivider << " " << mround1(getValue(vbfhww130,i)) << " " << pm << " " << mround2(getError(vbfhww130,i)) << " " << centerDivider;
//     if (vbfhww140) cout << centerDivider << " " << mround1(getValue(vbfhww140,i)) << " " << pm << " " << mround2(getError(vbfhww140,i)) << " " << centerDivider;
//     if (vbfhww150) cout << centerDivider << " " << mround1(getValue(vbfhww150,i)) << " " << pm << " " << mround2(getError(vbfhww150,i)) << " " << centerDivider;
//     if (vbfhww160) cout << centerDivider << " " << mround1(getValue(vbfhww160,i)) << " " << pm << " " << mround2(getError(vbfhww160,i)) << " " << centerDivider;
//     if (vbfhww170) cout << centerDivider << " " << mround1(getValue(vbfhww170,i)) << " " << pm << " " << mround2(getError(vbfhww170,i)) << " " << centerDivider;
//     if (vbfhww180) cout << centerDivider << " " << mround1(getValue(vbfhww180,i)) << " " << pm << " " << mround2(getError(vbfhww180,i)) << " " << centerDivider;
//     if (vbfhww190) cout << centerDivider << " " << mround1(getValue(vbfhww190,i)) << " " << pm << " " << mround2(getError(vbfhww190,i)) << " " << centerDivider;
//     if (vbfhww200) cout << centerDivider << " " << mround1(getValue(vbfhww200,i)) << " " << pm << " " << mround2(getError(vbfhww200,i)) << " " << centerDivider;
//     if (vbfhww250) cout << centerDivider << " " << mround1(getValue(vbfhww250,i)) << " " << pm << " " << mround2(getError(vbfhww250,i)) << " " << centerDivider;
    cout << endl;
  }
  cout <<endl;

  cout << "\n" << MForm("%s %3s %s",leftDivider,"",centerDivider);
  if (data) cout << MForm(patternTitle,"*Data*");
  cout <<endl;
  for (int i=0; i<5; i++){
    cout << leftDivider << MForm(" %3s ",names[i]) << centerDivider;
    if (data)   cout << MForm(patternData,getValue(data,i),getError(data,i));
    cout <<endl;
  }
  cout <<endl;


}
