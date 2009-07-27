
#ifndef DATASOURCE_H
#define DATASOURCE_H

#include "TString.h"

typedef UInt_t sources_t;

enum {
        H_WW      ,
        H_TTBAR   ,
        H_DYMM    ,
        H_DYEE    ,
        H_DYTT    ,
        H_WJETS   ,
        H_TW      ,
        H_ZZ      ,
        H_WZ,
	H_WENU,
	H_EM30_80,
	H_BC30_80,
};

const static sources_t sources_all =
 (1ll << H_WW)          |
 (1ll << H_TTBAR)       |
 (1ll << H_DYMM)        |
 (1ll << H_DYEE)        |
 (1ll << H_DYTT)        |
 (1ll << H_WJETS)       |
 (1ll << H_TW)          |
 (1ll << H_ZZ)          |
 (1ll << H_WZ);

const static sources_t sources_dy =
 (1ll << H_DYMM)        |
 (1ll << H_DYEE);

const static sources_t sources_peaking =
 (1ll << H_DYMM)        |
 (1ll << H_DYEE)        |
 (1ll << H_ZZ)          |
 (1ll << H_WZ);

const static sources_t sources_nonpeaking =
 (1ll << H_WW)          |
 (1ll << H_TTBAR)       |
 (1ll << H_DYTT)        |
 (1ll << H_WJETS)       |
 (1ll << H_TW);

class DataSource {

        public:
		DataSource() {}
                DataSource(TString name, sources_t source) {
                        sourceName_ = name;
                        source_ = source;
                }
                ~DataSource() {}

                TString         getName()       { return sourceName_; }
                sources_t       getSource()     { return source_; }
		sources_t 	getBit()	{ return 1ll << source_; }

        private:
                TString         sourceName_;
                sources_t       source_;

};

DataSource fH_WW();
DataSource fH_TTBAR();
DataSource fH_DYMM();
DataSource fH_DYEE();
DataSource fH_DYTT();
DataSource fH_WJETS();
DataSource fH_TW();
DataSource fH_ZZ();
DataSource fH_WZ();

DataSource fH_WENU();
DataSource fH_EM30_80();
DataSource fH_BC30_80();


#endif

