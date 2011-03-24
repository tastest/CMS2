#ifndef WW_Types_h
#define WW_Types_h
//
// signal samples
//
enum Sample {WW=0, WZ=1, ZZ=2, Wjets=3, DYee=4, DYmm=5, DYtt=6, ttbar=7, tW=8, qcd=9, Data=10, ggWW=11, 
	     hWW120=12, hWW130=13, hWW140=14, hWW150=15, hWW160=16, hWW170=17, hWW180=18, hWW190=19, hWW200=20,
	     hWW210=21, hWW220=22, hWW230=23, hWW240=24, hWW250=25
};

const char* SampleName(Sample sample);

//
// hypothesis type as it's stored in the ntuples
//
enum HypTypeInNtuples {MuMu, MuEl , ElMu, ElEl}; 

//
// hypothesis types as they are used for analysis
// (em and me counted as same)
//
enum HypothesisType {MM, EM, EE, ALL}; 
const char* HypothesisTypeName(HypothesisType type); 
const char* HypothesisTypeName(unsigned int type); 
HypothesisType getHypothesisType( HypTypeInNtuples type );
HypothesisType getHypothesisType( unsigned int );

#endif
