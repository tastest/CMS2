// -*- C++ -*-

#ifndef TRK_SELECTION_H
#define TRK_SELECTION_H

#include <utility>

enum TrackAlgorithm { undefAlgorithm=0, ctf=1, rs=2, cosmics=3, iter0=4, 
                      iter1=5, iter2=6, iter3=7, iter4=8, iter5=9, iter6=10, iter7=11, iter8=12, iter9=13,iter10=14,
                      outInEcalSeededConv=15, inOutEcalSeededConv=16, 
                      nuclInter=17,
                      standAloneMuon=18,globalMuon=19,cosmicStandAloneMuon=20,cosmicGlobalMuon=21,
                      iter1LargeD0=22,iter2LargeD0=23,iter3LargeD0=24,iter4LargeD0=25,iter5LargeD0=26,
                      bTagGhostTracks=27,
                      beamhalo=28, 
                      gsf=29,
                      algoSize=30 };

enum TrackQuality { undefQuality=-1, loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4, qualitySize=5};

std::pair<double, double> trks_d0_pv (int itrk, int ipv);
bool isTrackQuality( int index, int cuts );

#endif
