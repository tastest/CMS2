// -*- C++ -*-

// $Id: DSGTable.h,v 1.1 2009/07/02 15:37:49 jmuelmen Exp $

#ifndef DSGTABLE_H
#define DSGTABLE_H

#include "TNamed.h"
#include "Tools/Sample.h"

class TH1F;
class DSGTable : public TNamed {
public:
     static const int	nZcat 		= 2;
     static const int	nMETcat 	= 3;
     static const int	nJetcat 	= 3;
     static const int	nBuckets 	= 10;

public: 
     DSGTable () : TNamed() { }
     DSGTable (Sample s) : TNamed(s.name.c_str(), s.name.c_str()) 
	  {
	       memset(events_, 0, sizeof(events_));
	       for (int i = 0; i < nZcat; ++i) {
		    for (int j = 0; j < nMETcat; ++j) {
			 for (int k = 0; k < nJetcat; ++k) {
			      for (int l = 0; l < nBuckets; ++l) {
				   hmet_[i][j][k][l] = new TH1F(Form("hmet%d%d%d%d", i, j, k, l), "MET;MET", 100, 0, 500);
				   hmet_[i][j][k][l]->SetFillColor(s.histo_color);
				   hmet_[i][j][k][l]->SetFillStyle(1001);
				   hmll_[i][j][k][l] = new TH1F(Form("hmll%d%d%d%d", i, j, k, l), "Mll;Mll", 100, 0, 500);
				   hmll_[i][j][k][l]->SetFillColor(s.histo_color);
				   hmll_[i][j][k][l]->SetFillStyle(1001);
			      }
			 }
		    }
	       }
	  }
     double 	Increment (int zcat, int metcat, int jetcat, int bucket, 
			   double weight)
	  {
	       return events_[zcat][metcat][jetcat][bucket] += weight;
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
     double		events_[nZcat   ][nMETcat ][nJetcat ][nBuckets];
     TH1F		*hmet_ [nZcat   ][nMETcat ][nJetcat ][nBuckets];
     TH1F		*hmll_ [nZcat   ][nMETcat ][nJetcat ][nBuckets];

public: 
     ClassDef(DSGTable, 1)
};

#endif
