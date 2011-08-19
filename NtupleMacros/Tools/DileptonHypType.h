// -*- C++ -*-
#ifndef DILEPTONHYPTYPE_H
#define DILEPTONHYPTYPE_H

enum DileptonHypType {
     DILEPTON_ALL,
     DILEPTON_MUMU,
     DILEPTON_EMU,
     DILEPTON_EE
};

static const char dilepton_hypo_names[][128] = { "all", "mm", "em", "ee" };

enum DileptonHypType hyp_typeToHypType (int hyp_type);

#endif
