// -*- C++ -*-

// $Id: jetSelections.h,v 1.1 2010/02/19 20:23:33 jmuelmen Exp $

#ifndef JETSELECTIONS_H
#define JETSELECTIONS_H

#include <vector>
#include "CMS2.h"

enum JetType {
     JETS_TYPE_JPT,
     JETS_TYPE_PF_CORR,
     JETS_TYPE_PF_UNCORR,
     JETS_TYPE_CALO_CORR,
     JETS_TYPE_CALO_UNCORR
};

enum CleaningType {
     JETS_CLEAN_NONE,			// dirty, dirty jets
     JETS_CLEAN_HYP_E_MU,		// e or mu from hypothesis
     JETS_CLEAN_HYP_E,			// e only from hypothesis
     JETS_CLEAN_SINGLE_E,		// e in single-lepton final state (QCD)
};

#define JET_DEFAULT_TYPE 	JETS_TYPE_JPT
#define JET_DEFAULT_CLEANING 	JETS_CLEAN_HYP_E_MU
#define JET_DEFAULT_DR		0.4
#define JET_DEFAULT_PT		30
#define JET_DEFAULT_ETA		2.4

std::vector<LorentzVector> getJets (unsigned int i_hyp,  // hyp or single-e to use for cleaning
				    bool sort = false,
				    enum JetType = JET_DEFAULT_TYPE,
				    enum CleaningType = JET_DEFAULT_CLEANING,
				    double deltaR = JET_DEFAULT_DR,
				    double min_pt = JET_DEFAULT_PT,
				    double max_eta = JET_DEFAULT_ETA);

int nJets (unsigned int i_hyp,  // hyp or single-e to use for cleaning
	   bool sort = false,
	   enum JetType = JET_DEFAULT_TYPE,
	   enum CleaningType = JET_DEFAULT_CLEANING,
	   double deltaR = JET_DEFAULT_DR,
	   double min_pt = JET_DEFAULT_PT,
	   double max_eta = JET_DEFAULT_ETA);

double sumPt (unsigned int i_hyp,  // hyp or single-e to use for cleaning
	      bool sort = false,
	      enum JetType = JET_DEFAULT_TYPE,
	      enum CleaningType = JET_DEFAULT_CLEANING,
	      double deltaR = JET_DEFAULT_DR,
	      double min_pt = JET_DEFAULT_PT,
	      double max_eta = JET_DEFAULT_ETA);

#endif // SEL_JETS_H
