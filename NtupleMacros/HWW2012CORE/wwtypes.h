#ifndef WW_Types_h
#define WW_Types_h
//
// hypothesis type as it's stored in the ntuples
//
enum HypTypeInNtuples {MuMu, MuEl , ElMu, ElEl}; 

//
// hypothesis types as they are used for analysis
//
enum HypothesisType {MM, ME, EM, EE, ALL}; 
const char* HypothesisTypeName(HypothesisType type); 
const char* HypothesisTypeName(unsigned int type); 
HypothesisType getHypothesisType( HypTypeInNtuples type );
HypothesisType getHypothesisType( unsigned int );
HypothesisType getHypothesisTypeNew( unsigned int i_hyp );

#endif
