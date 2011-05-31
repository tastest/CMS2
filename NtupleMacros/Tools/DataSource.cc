
#include "DataSource.h"

DataSource::DataSource()
{
}

DataSource::DataSource(TString name, sources_t source, Color_t color, TString legendName) 
{
	sourceName_ = name;
    legendName_ = legendName;
	source_ = source;
	color_ = color;
}

DataSource::~DataSource()
{
}

TString DataSource::getLegendName() 
{ 
    if (legendName_ != "") return legendName_;
    return sourceName_;
}

