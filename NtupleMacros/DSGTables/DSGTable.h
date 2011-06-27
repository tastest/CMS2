// -*- C++ -*-

// $Id: DSGTable.h,v 1.9 2010/08/14 20:11:03 jmuelmen Exp $

#ifndef DSGTABLE_H
#define DSGTABLE_H

#include <map>
#include <vector>
#include "TNamed.h"
#include "Tools/Sample.h"
#include "DSGSearchWindow.h"

class TH1F;
class DSGTable : public TNamed {
public:
     static const int	nZcat 		= 2;
     static const int	nMETcat 	= 5;
     static const int	nSumJetcat	= 3;
     static const int	nJetcat 	= 3;
     static const int	nBuckets 	= 10;

public: 
#ifndef __CINT__
     DSGTable (Sample s) : TNamed(s.name.c_str(), s.name.c_str()) 
	  {
	       memset(events_, 0, sizeof(events_));
	       memset(w2s_, 0, sizeof(w2s_));
	       for (int i = 0; i <= nZcat; ++i) {
		    for (int j = 0; j < nMETcat; ++j) {
		    for (int jj = 0; jj < nSumJetcat; ++jj) {
			 for (int k = 0; k < nJetcat; ++k) {
			      for (int l = 0; l < nBuckets; ++l) {
				   hmet_[i][j][jj][k][l] = new TH1F(Form("hmet%s%d%d%d%d", s.name.c_str(), i, j, k, l), "MET;MET", 10, 0, 500);
				   hmet_[i][j][jj][k][l]->SetFillColor(s.histo_color);
				   hmet_[i][j][jj][k][l]->SetFillStyle(1001);

				   hmll_[i][j][jj][k][l] = new TH1F(Form("hmll%s%d%d%d%d", s.name.c_str(), i, j, k, l), "Mll;Mll", 10, 0, 500);
				   hmll_[i][j][jj][k][l]->SetFillColor(s.histo_color);
				   hmll_[i][j][jj][k][l]->SetFillStyle(1001);

				   hht_[i][j][jj][k][l] = new TH1F(Form("hht%s%d%d%d%d", s.name.c_str(), i, j, k, l), "HT;HT", 20, 0, 1000);
				   hht_[i][j][jj][k][l]->SetFillColor(s.histo_color);
				   hht_[i][j][jj][k][l]->SetFillStyle(1001);

				   hjsumet_[i][j][jj][k][l] = new TH1F(Form("hjsumet%s%d%d%d%d", s.name.c_str(), i, j, k, l), "jetSumEt;jetSumEt", 20, 0, 1000);
				   hjsumet_[i][j][jj][k][l]->SetFillColor(s.histo_color);
				   hjsumet_[i][j][jj][k][l]->SetFillStyle(1001);

				   hmaxjetpt_[i][j][jj][k][l] = new TH1F(Form("hmaxjetpt%s%d%d%d%d", s.name.c_str(), i, j, k, l), "maxJetPt;maxJetPt", 10, 0, 500);
				   hmaxjetpt_[i][j][jj][k][l]->SetFillColor(s.histo_color);
				   hmaxjetpt_[i][j][jj][k][l]->SetFillStyle(1001);

				   hmaxleppt_[i][j][jj][k][l] = new TH1F(Form("hmaxleppt%s%d%d%d%d", s.name.c_str(), i, j, k, l), "maxLepPt;maxLepPt", 10, 0, 250);
				   hmaxleppt_[i][j][jj][k][l]->SetFillColor(s.histo_color);
				   hmaxleppt_[i][j][jj][k][l]->SetFillStyle(1001);

				   hlepdphi_[i][j][jj][k][l] = new TH1F(Form("hlepdphi%s%d%d%d%d", s.name.c_str(), i, j, k, l), "lepDphi;lepDphi", 10, -3.2, 3.2);
				   hlepdphi_[i][j][jj][k][l]->SetFillColor(s.histo_color);
				   hlepdphi_[i][j][jj][k][l]->SetFillStyle(1001);
			      }
			 }
		    }
		    }
	       }
	  }
     DSGTable (const DSGTable &other) : TNamed(other)
	  {
	       memcpy(events_, other.events_, sizeof(events_));
	       memset(w2s_, 0, sizeof(w2s_));
	       for (int i = 0; i <= nZcat; ++i) {
		    for (int j = 0; j < nMETcat; ++j) {
		    for (int jj = 0; jj < nSumJetcat; ++jj) {
			 for (int k = 0; k < nJetcat; ++k) {
			      for (int l = 0; l < nBuckets; ++l) {
				   hmet_[i][j][jj][k][l] = new TH1F(*other.hmet_[i][j][jj][k][l]);
				   hmll_[i][j][jj][k][l] = new TH1F(*other.hmll_[i][j][jj][k][l]);
				   hht_[i][j][jj][k][l] = new TH1F(*other.hht_[i][j][jj][k][l]);
				   hjsumet_[i][j][jj][k][l] = new TH1F(*other.hjsumet_[i][j][jj][k][l]);
				   hmaxjetpt_[i][j][jj][k][l] = new TH1F(*other.hmaxjetpt_[i][j][jj][k][l]);
				   hmaxleppt_[i][j][jj][k][l] = new TH1F(*other.hmaxleppt_[i][j][jj][k][l]);
				   hlepdphi_[i][j][jj][k][l] = new TH1F(*other.hlepdphi_[i][j][jj][k][l]);
			      }
			 }
		    }
		    }
	       }
	  }
     DSGTable () : TNamed()
	  {
	       memset(events_, 0, sizeof(events_));
	       memset(w2s_, 0, sizeof(w2s_));
	       for (int i = 0; i <= nZcat; ++i) {
		    for (int j = 0; j < nMETcat; ++j) {
		    for (int jj = 0; jj < nSumJetcat; ++jj) {
			 for (int k = 0; k < nJetcat; ++k) {
			      for (int l = 0; l < nBuckets; ++l) {
				   hmet_     [i][j][jj][k][l] = new TH1F;
				   hmll_     [i][j][jj][k][l] = new TH1F;
				   hht_      [i][j][jj][k][l] = new TH1F;
				   hjsumet_  [i][j][jj][k][l] = new TH1F;
				   hmaxjetpt_[i][j][jj][k][l] = new TH1F;
				   hmaxleppt_[i][j][jj][k][l] = new TH1F;
				   hlepdphi_ [i][j][jj][k][l] = new TH1F;


			      }
			 }
		    }
		    }
	       }
	  }
#endif
     double 	Increment (int zcat, int metcat, int sumjetcat,  int jetcat, int bucket, 
			   double weight)
	  {
	       w2s_[zcat][metcat][sumjetcat][jetcat][bucket] += weight * weight;
	       return events_[zcat][metcat][sumjetcat][jetcat][bucket] += weight;
	  }
     DSGTable &operator *= (double scale)
	  {
	       for (int i = 0; i <= nZcat; ++i) {
		    for (int j = 0; j < nMETcat; ++j) {
		    for (int jj = 0; jj < nSumJetcat; ++jj) {
			 for (int k = 0; k < nJetcat; ++k) {
			      for (int l = 0; l < nBuckets; ++l) {
				   events_[i][j][jj][k][l] 	*= scale;
				   w2s_[i][j][jj][k][l] 	*= scale * scale;
				   hmet_[i][j][jj][k][l]      	->Scale(scale);
				   hmll_[i][j][jj][k][l]      	->Scale(scale);
				   hht_[i][j][jj][k][l]       	->Scale(scale);
				   hjsumet_[i][j][jj][k][l]   	->Scale(scale);
				   hmaxjetpt_[i][j][jj][k][l] 	->Scale(scale);
				   hmaxleppt_[i][j][jj][k][l] 	->Scale(scale);
				   hlepdphi_[i][j][jj][k][l]  	->Scale(scale);
			      }
			 }
		    }
		    }
	       }
	       for (unsigned int i = 0; i < search_windows_.size(); ++i) {
		    *search_windows_[i] *= scale;
	       }
	       return *this;
	  }
//      void	FillMET (int zcat, int metcat, int jetcat, int bucket, 
// 			 double met, double weight) 
// 	  {
	       
// 	  }
//      void	FillMll (int zcat, int metcat, int jetcat, int bucket, 
// 			 double mll, double weight) 
// 	  {
	       
// 	  }

public:
     typedef TH1F* table_t     [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     typedef std::vector<DSGSearchWindow *>	sw_t;
     double		events_[nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     double		w2s_   [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     TH1F		*hmet_ [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     TH1F		*hmll_ [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     TH1F		*hht_       [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     TH1F		*hjsumet_   [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     TH1F		*hmaxjetpt_ [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     TH1F		*hmaxleppt_ [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     TH1F		*hlepdphi_  [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     sw_t		search_windows_;

public:
//      struct chunk_of_runs {
	  int		min_run_;
	  int		max_run_;
	  double 	lumi_;
//      };
//      std::vector<chunk_of_runs> chunks_of_runs_;
     void	SetMinMaxLumi (int min, int max, double lumi) 
	  {
// 	       chunk_of_runs chunk = { min, max, lumi };
// 	       chunks_of_runs_.clear();
// 	       chunks_of_runs_.push_back(chunk);
	       min_run_ = min;
	       max_run_ = max;
	       lumi_ = lumi;
	  }

public:
     ClassDef(DSGTable, 1)
};

class DSGTables : public TNamed {
public:
     DSGTables () : TNamed() { }
     DSGTables (const char *name) : TNamed(name, name) { }

public:
     std::vector<DSGTable>	tables_;
     
public: 
     ClassDef(DSGTables, 1)
};

#endif
