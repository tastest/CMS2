
//
// electron selections
//

// CMS2 includes
#include "electronSelections.h"
#include "CMS2.h"

//
// candidate electron id function
//
bool electronId_cand01()
{
	return false;
}

//
// class based electron id that we have used before
//
bool electronId_classBased()
{
	return false;
}

//
// candidate electron isolation function
//
bool electronIsolation_cand01()
{
	return false;
}

//
// functions for individual subdet electron isolations
//
bool electronTkIsolation(const unsigned int index, const float threshold)
{
	return false;
}

bool electronEcalIsolation(const unsigned int index, const float threshold)
{
        return false;
}

bool electronHcalIsolation(const unsigned int index, const float threshold)
{
        return false;
}

