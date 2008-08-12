// -*- C++ -*-
//Add includes for your classes here
#include <vector>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace {
   namespace {
      std::vector<std::vector<int> > vi2d;
      edm::Wrapper<std::vector<std::vector<int> > > wvi2d;
	
      std::vector<std::vector<std::vector<int> > > vi3d;
      edm::Wrapper<std::vector<std::vector<std::vector<int> > > > wvi3d;
      
      std::vector<std::vector<float> > vf2d;
      edm::Wrapper<std::vector<std::vector<float> > > wvf2d;
      
      std::vector<std::vector<std::vector<float> > > vf3d;
      edm::Wrapper<std::vector<std::vector<std::vector<float> > > > wvf3d;
      
      std::vector<std::vector<math::XYZTLorentzVector> > vp2d;
      edm::Wrapper<std::vector<std::vector<math::XYZTLorentzVector> > > wvp2d;
      
      std::vector<std::vector<std::vector<math::XYZTLorentzVector> > > vp3d;
      edm::Wrapper<std::vector<std::vector<std::vector<math::XYZTLorentzVector> > > > wvp3d;
   }
}
