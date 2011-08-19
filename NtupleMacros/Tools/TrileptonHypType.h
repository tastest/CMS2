// -*- C++ -*-
#ifndef TRILEPTONHYPTYPE_H
#define TRILEPTONHYPTYPE_H

enum TrileptonHypType {
     TRILEPTON_MPMPMP,
     TRILEPTON_MPMPMM,
     TRILEPTON_MPMMMM,
     TRILEPTON_MMMMMM,
     TRILEPTON_MPMPEP,
     TRILEPTON_MPMPEM,
     TRILEPTON_MPMMEP,
     TRILEPTON_MPMMEM,
     TRILEPTON_MMMMEP,
     TRILEPTON_MMMMEM,
     TRILEPTON_MPEPEP,
     TRILEPTON_MPEPEM,
     TRILEPTON_MPEMEM,
     TRILEPTON_MMEPEP,
     TRILEPTON_MMEPEM,
     TRILEPTON_MMEMEM,
     TRILEPTON_EPEPEP,
     TRILEPTON_EPEPEM,
     TRILEPTON_EPEMEM,
     TRILEPTON_EMEMEM,
     TRILEPTON_ALL,
};

static const char trilepton_hypo_names[][128] = {"mpmpmp","mpmpmm","mpmmmm","mmmmmm","mpmpep","mpmpem","mpmmep","mpmmem","mmmmep","mmmmem","mpepep","mpepem","mpemem","mmepep","mmepem","mmemem","epepep","epepem","epemem","ememem","all"};

#endif
