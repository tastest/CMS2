#include "wwtypes.h"
#include <iostream>

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
