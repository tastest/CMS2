
#include "DataSource.h"

DataSource::DataSource()
{
}

DataSource::DataSource(TString name, sources_t source, Color_t color) 
{
	sourceName_ = name;
	source_ = source;
	color_ = color;
}

DataSource::~DataSource()
{
}

