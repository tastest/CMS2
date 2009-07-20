// -*- C++ -*-

// $Id: DSGTable.h,v 1.5 2009/07/20 13:53:21 avi Exp $

#ifndef DSGTABLE_H
#define DSGTABLE_H

#include <vector>
#include "TNamed.h"
#include "Tools/Sample.h"

class TH1F;
class DSGTable : public TNamed {
public:
     static const int	nZcat 		= 2;
     static const int	nMETcat 	= 3;
     static const int	nSumJetcat	= 3;
     static const int	nJetcat 	= 3;
     static const int	nBuckets 	= 10;

public: 
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
				   hmet_[i][j][jj][k][l] = new TH1F(Form("hmet%d%d%d%d", i, j, k, l), "MET;MET", 10, 0, 500);
				   hmll_[i][j][jj][k][l] = new TH1F(Form("hmll%d%d%d%d", i, j, k, l), "Mll;Mll", 10, 0, 500);
			      }
			 }
		    }
		    }
	       }
	  }
     double 	Increment (int zcat, int metcat, int sumjetcat,  int jetcat, int bucket, 
			   double weight)
	  {
	       w2s_[zcat][metcat][sumjetcat][jetcat][bucket] += weight * weight;
	       return events_[zcat][metcat][sumjetcat][jetcat][bucket] += weight;
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
     double		events_[nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     double		w2s_   [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     TH1F		*hmet_ [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];
     TH1F		*hmll_ [nZcat + 1   ][nMETcat ][nSumJetcat][nJetcat ][nBuckets];

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
