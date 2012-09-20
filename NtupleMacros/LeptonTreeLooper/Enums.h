#ifndef LEPTONTREELOOPER_ENUMS_H
#define LEPTONTREELOOPER_ENUMS_H

namespace plotbin {
enum PlotPtBin {
    EB      = (1<<0),
    EE      = (1<<1),
    CENTRAL     = (1<<2),
    TRANSITION  = (1<<3),
    EELO       = (1<<4),
    EEHI        = (1<<5),
    ALLETA = (1<<6),

    PT1015  = (1<<7),
    PT1020  = (1<<8),
    PT20UP  = (1<<9),
    PT10UP  = (1<<10)

};
}

#endif

