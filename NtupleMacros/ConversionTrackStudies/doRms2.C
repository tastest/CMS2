//
{
  gSystem->Load("/code/osgcode/cmssoft/cms/slc5_ia32_gcc434/cms/cmssw/CMSSW_3_4_0/lib/slc5_ia32_gcc434/libValidationRecoParticleFlow.so");
  gROOT->ProcessLine(".L warrenRms2.C++");

  //data first arg, mc 2nd, third is is2tev
  //doRms("FlatTree_jan23_12.root", "FlatTree_mc_3.root");

  doRms("FlatTree_2tev_5.root", "FlatTree_mc_2tev_1.root", true);

}
