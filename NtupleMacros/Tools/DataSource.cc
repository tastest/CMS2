
#include "DataSource.h"

DataSource::DataSource()
{
}

DataSource::DataSource(TString name, sources_t source, Color_t color, Int_t fillStyle, TString legendName) 
{
	sourceName_ = name;
    legendName_ = legendName;
	source_ = source;
	color_ = color;
    fillStyle_ = fillStyle;
}

DataSource::~DataSource()
{
}

TString DataSource::getLegendName() 
{ 
    if (legendName_ != "") return legendName_;
    return sourceName_;
}

