#include "Generator.hh"
#include "KinematicSolver.hh"
#include "PhaseSpace.hh"
#include "TEvtProb.hh"
#include "TMCFM.hh"
#include "TMatrixElement.hh"
#include "TUtil.hh"
#include "TVar.hh"

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// #pragma link C++ function ScanChain(TChain*, enum Sample, double, bool);
#pragma link C++  nestedclasses;
#pragma link C++  nestedtypedefs;


#pragma link C++  function  GetRandom;
#pragma link C++  function WWL1L2Sol_MHiggsMw1;
#pragma link C++  class  TVar;
#pragma link C++  class  TEvtProb;
#pragma link C++  class  TMatrixElement;
#pragma link C++  function breitw;
#pragma link C++  function SetTGCParameter;
#pragma link C++  function GetTGCParameter;
#pragma link C++  function SumMatrixElementPDF;



#endif
