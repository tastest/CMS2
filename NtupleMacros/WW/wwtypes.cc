#include "wwtypes.h"
#include <iostream>
const char* SampleName(Sample sample){
  switch (sample){
  case WW:
    return "ww";
    break;
  case WZ:
    return "wz";
    break;
  case ZZ:
    return "zz";
    break;
  case Wjets:
    return "wjets";
    break;
  case DYee:
    return "dyee";
    break;
  case DYmm:
    return "dymm";
    break;
  case DYtt:
    return "dytt";
    break;
  case ttbar:
    return "ttbar";
    break;
  case tW:
    return "tw";
    break;
  case qcd:
    return "qcd";
    break;
  case Data:
    return "data";
    break;
  }
  std::cout << "ERROR: unknown sample name is requested: " << sample << std::endl;
  return "";
}

const char* HypothesisTypeName(HypothesisType type){
  switch (type){
  case MM:
    return "mm";
    break;
  case EM:
    return "em";
    break;
  case EE:
    return "ee";
    break;
  case ALL:
    return "all";
    break;
  }
  std::cout << "ERROR: unknown hypothesis name is requested: " << type << std::endl;
  return "";
}

const char* HypothesisTypeName(unsigned int type){
  return HypothesisTypeName( HypothesisType(type) );
}

HypothesisType getHypothesisType( HypTypeInNtuples type ){
  switch (type){
  case ElEl:
    return EE;
    break;
  case MuMu:
    return MM;
    break;
  case ElMu: 
  case MuEl:
    return EM;
    break;
  }
  return EM; //to stop warnings
}

HypothesisType getHypothesisType( unsigned int type ){
  return getHypothesisType( HypTypeInNtuples(type) );
}
