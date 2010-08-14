// -*- C++ -*-

// $Id: DSGSearchWindow.h,v 1.2 2010/08/14 20:10:32 jmuelmen Exp $

#ifndef DSGSEARCHWINDOW_H
#define DSGSEARCHWINDOW_H

#include <string>
#include <vector>
#include "TNamed.h"
#include "TH1.h"
#include "Tools/Sample.h"

class TH1F;
class DSGTable;
class DSGSearchWindow : public TNamed {
public:
     double events_;
     double w2s_;
     std::string name_;

public: 
     TH1F *h_;
     DSGTable *dsgTable_;

public:
     Sample s_; //!

public: 
     DSGSearchWindow (Sample s, const char *name) : TNamed(s.name.c_str(), s.name.c_str()),
						    events_(0),
						    w2s_(0),
						    name_(name),
						    h_(0),
						    dsgTable_(0),
						    s_(s)
	  {

	  }
     DSGSearchWindow () : TNamed(),
			  events_(0),
			  w2s_(0),
			  h_(0),
			  dsgTable_(0)
	  {

	  }
     DSGSearchWindow &operator *= (double scale);
     void 	SetH (TH1F *h) 
	  {
	       h_ = h;
	       h_->SetFillColor(s_.histo_color);
	       h_->SetFillStyle(1001);
	  }
     double 	Increment (double value, double weight)
	  {
	       h_->Fill(value, weight);
	       w2s_ += weight * weight;
	       return events_ += weight;
	  }

public:
     ClassDef(DSGSearchWindow, 1)
};

#endif
