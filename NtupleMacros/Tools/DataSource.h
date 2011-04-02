
#ifndef DATASOURCE_H
#define DATASOURCE_H

#include "TString.h"
#include "TColor.h"
#include <iostream>

typedef ULong64_t   uint64;
typedef uint64      sources_t;

class DataSource {

	public:
		DataSource();
		DataSource(TString name, sources_t source, Color_t color = 0, TString legendName = "");
		~DataSource();

		TString         getName()       { return sourceName_; }
        TString         getLegendName();
		sources_t       getSource()     { return source_; }
		sources_t 	getBit()	{ return 1ll << source_; }
		Color_t		getColor()	{ return color_; }
	private:
		TString         sourceName_;
        TString         legendName_;
		sources_t       source_;
		Color_t		color_;

};

#endif

