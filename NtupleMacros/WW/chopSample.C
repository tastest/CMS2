//
// * How to use the tool *
//
// For our current setup we have samples with sidebands, so we need to
// split samples not by expected final yield, but by yield on
// preselection level. So we need to find the number of events at
// preselection level that corresponds to the final desired yield.
//
#include "TFile.h"
#include <iostream>
#include "TRandom3.h"
#include "TTree.h"

void chopSample(const char* input, double mean, const char* out_prefix)
{
  cout << "Processing " << input << endl;
  cout << "Mean: " << mean << endl;
  TRandom3* gen = new TRandom3();
  TFile* f = TFile::Open(input);
  assert(f);
  TTree* origtree = (TTree*)f->Get("tree");
  assert(origtree);
  unsigned int nTotal = origtree->GetEntries();
  unsigned int i=0;
  unsigned int index=1;
  while (i<nTotal){
    unsigned int n = gen->Poisson(mean);
    if (i+n>=nTotal) break;
    cout << "number of entries to be stored: " << n << endl;
    TFile *output = TFile::Open(Form("%s_%u.root",out_prefix,index),"RECREATE");
    assert(output);
    TTree* t_clone=origtree->CloneTree(0);
    for(unsigned int iev=i; iev <i+n ; iev++ ){
      origtree->GetEntry( iev );
      t_clone->Fill();
    }
    output->Write();
    output->Close();
    i+=n;
    index++;
  }
}

void chopAllSamples(){
  //   chopSample("/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets_3500ipb/qqww.root",1.43862277221679688e+03*3.6*1.2,"/smurf/dmytro/ichep2012/pseudo-data/qqww");
  chopSample("/smurf/dmytro/ichep2012/pseudo-data/input/gghww125.root",8.75360107421875000e+01*3.6,"/smurf/dmytro/ichep2012/pseudo-data/gghww125");
  chopSample("/smurf/dmytro/ichep2012/pseudo-data/input/gghww130.root",1.23496063232421875e+02*3.6,"/smurf/dmytro/ichep2012/pseudo-data/gghww130");

//   chopSample("/smurf/data/EPS/tas/qqww.root",   407/0.46896,     "/smurf/data/EPS/pseudo-data/qqww");
//   chopSample("/smurf/data/EPS/tas/dytt.root",   1.6/0.000243996, "/smurf/data/EPS/pseudo-data/dytt");
//   chopSample("/smurf/data/EPS/tas/ggww.root",   17.2/0.472701,   "/smurf/data/EPS/pseudo-data/ggww");
//   chopSample("/smurf/data/EPS/tas/wgamma.root", 8.7/0.0949153,   "/smurf/data/EPS/pseudo-data/wgamma");
//   chopSample("/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets/wjets-monster.root", 106.9/0.0522066, "/smurf/data/EPS/pseudo-data/wjets");
//   chopSample("/smurf/data/EPS/tas/dyee.root", 7.2/3/1.02176e-05, "/smurf/data/EPS/pseudo-data/dyee");
//   chopSample("/smurf/data/EPS/tas/dymm.root", 7.2/3*2/2.53514e-05, "/smurf/data/EPS/pseudo-data/dymm");
//   chopSample("/smurf/data/EPS/tas/ttbar.root", 42.38/0.00524255, "/smurf/data/EPS/pseudo-data/ttbar");
//   chopSample("/smurf/data/EPS/tas/tw.root", 21.38/0.0379785, "/smurf/data/EPS/pseudo-data/tw");
//   chopSample("/smurf/data/EPS/tas/wz.root", 9.22/0.050117, "/smurf/data/EPS/pseudo-data/wz");
//   chopSample("/smurf/data/EPS/tas/zz.root", 4.27/0.0547364, "/smurf/data/EPS/pseudo-data/zz");

}
