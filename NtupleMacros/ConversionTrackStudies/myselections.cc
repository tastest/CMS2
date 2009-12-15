#include "myselections.h"

bool isTrackQuality( int index, int cuts ) {
     return ( ( cms2.trks_qualityMask().at(index) & cuts ) == cuts );
}


