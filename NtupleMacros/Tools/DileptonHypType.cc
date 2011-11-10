#include "DileptonHypType.h"
#include <assert.h>

enum DileptonHypType hyp_typeToHypType (int hyp_type)
{
    switch (hyp_type) {
    case 0:
        return DILEPTON_MUMU;
    case 1: case 2:
        return DILEPTON_EMU;
    case 3:
        return DILEPTON_EE;
    default:
        assert(hyp_type < 4);
    }
    return DILEPTON_ALL;
}
