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
  case ggWW:
    return "ggww";
    break;
  case hWW120:
    return "hww120";
    break;
  case hWW130:
    return "hww130";
    break;
  case hWW140:
    return "hww140";
    break;
  case hWW150:
    return "hww150";
    break;
  case hWW160:
    return "hww160";
    break;
  case hWW170:
    return "hww170";
    break;
  case hWW180:
    return "hww180";
    break;
  case hWW190:
    return "hww190";
    break;
  case hWW200:
    return "hww200";
    break;
  case hWW210:
    return "hww210";
    break;
  case hWW220:
    return "hww220";
    break;
  case hWW230:
    return "hww230";
    break;
  case hWW240:
    return "hww240";
    break;
  case hWW250:
    return "hww250";
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
