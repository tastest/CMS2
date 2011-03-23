#include "TVar.hh"
#include "TEvtProb.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include <iostream>

double ERRORthreshold=1.0;
using namespace std;

enum process {
  proc_WW,
  proc_Wpj,
  proc_Wmj,
  proc_ZZ,
  proc_HWW160,
  kNProc
};

enum { kNDilep=3}; 


float lumi=1000.0;
float yield [kNProc][kNDilep];
float acceptance [kNProc][kNDilep];

float BR [kNProc][kNDilep] = {{ (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1784)/3, (2.0+0.1736+0.1784)/3, (1+0.1784)/3},
			      { (1.0+0.1784)/3, (2.0+0.1736+0.1784)/3, (1+0.1784)/3},
			      { 2*0.2*(0.0363+0.0370*0.1736*0.1736), 0.0, 2*0.2*(0.0363+0.0370*0.1784*0.1784) },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
};

float NLOXsec[kNProc] = {  4.5*0.919, 31314.0/2.0, 31314.0/2.0, 5.9, 0.8664429};
float MCFMXsec[kNProc] = { 28.4, 11270, 11270.0, 4.3, 3.25};




void CalculateAcceptance(){

    for (int i=0; i<kNProc; i++){
      for (int j=0; j<kNDilep; j++){
      float denominator = lumi* BR[i][j] * NLOXsec[i];
      if (denominator!=0) acceptance [i][j]= yield [i][j]/denominator; else acceptance [i][j]=1E+20;
      }
    }

}

/*
void InitializeYields(){

  float temp_yield[kNProc][kNDilep] = {
    { 36.1939, 110.983, 23.1874},
    { 0.0001, 11.35, 7.22 }, // MM 0; EM 22.7089; EE 14.4511                                                                                             
    { 0.0001, 11.35, 7.22 },
    { 1.0,  0.158, 0.6 },
    { 18.2133, 34.9456, 12.5212},
  };

  for (int i=0; i<kNProc; i++){
    for (int j=0; j<kNDilep; j++){
      yield[i][j]=temp_yield[i][j];
	}
  }
 
}

*/

void InitializeYields(TString Higgs){

  TString FileName;
  FileName = "../Util.root";
  TFile *f = TFile::Open(FileName, "READ");
  assert(f);
  gROOT->cd();

  TString histName;
  
  TString procName[kNProc] = {  "WW" , "Wpj", "Wmj", "ZZ", Higgs };
  TString dilName[kNDilep]={"mm", "em", "ee" };

 for (int i=0; i<kNProc; i++){                                                                                                                           
   for (int j=0; j<kNDilep; j++){                                                                                                                       
     histName=procName[i]+"_dilmass_"+dilName[j];
     TH1F* h = (TH1F*) f->Get(histName);
     if (h!= 0x0) yield[i][j]=h->Integral(); 
     else yield[i][j]=0;
   }       
   cout<<procName[i]<<":   "<<yield[i][0]<<"   "<<yield[i][1]<<"   "<<yield[i][2]<<"\n"; 
  }        
  f->Close();
}



void NeutrinoIntegration(int process,TString inputFileName,int seed, int SmearLevel,int ncalls,int maxevt, int Mass);




//###################
//# main function
//###################
void runME_WWHWW(TString inputFileName,int seed,int SmearLevel,int ncalls,double Error, int Mass, int nev){
  /*
    gSystem->SetIncludePath("-I$ROOTSYS/include"); 
    gSystem->Load("./libgfortran.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("./libmcfm.so");
    gSystem->Load("./libME.so");

  */

  ERRORthreshold=Error;
  int process=TVar::HWW;
  int maxevt=nev;


  cout <<"=== Neutrino Integration ==========" <<endl;  
  NeutrinoIntegration(process,inputFileName,seed, SmearLevel,ncalls,maxevt, Mass); 
 
}


void NeutrinoIntegration(int process,TString inputFileName,int seed, int SmearLevel,int ncalls,int maxevt, int Mass){

    cout << "Input File: "<< inputFileName <<" seed " << seed <<" SmearLevel " << SmearLevel <<" ncalls "<< ncalls<<endl;
    TFile* fin = new TFile(inputFileName);
    TString outFileName=inputFileName;
    outFileName.ReplaceAll(".root","_ME.root");

    TFile *newfile = new TFile(outFileName,"recreate");
   	

    TEvtProb Xcal2;  
 
    Xcal2.SetApplyFake(1);
    Xcal2.SetSmearLevel(SmearLevel);
    Xcal2.SetSeed(seed);
    Xcal2.SetMatrixElement(TVar::MCFM);
    Xcal2.SetNcalls(ncalls);
    

    cout <<"Integration Seed= "<< Xcal2._seed << " SmearLevel= "<< Xcal2._smearLevel << " Ncalls = " << Xcal2._ncalls <<  endl;  

    TTree* ch=(TTree*)fin->Get("Events"); 
    if (ch==0x0) ch=(TTree*)fin->Get("ev");
    TTree* evt_tree=(TTree*) ch->CloneTree(0);
    ch->Print();
  
    int ProcIdx=process;
    int nProc=60;
    int runNumber;
    int eventNumber;
    double dXsecList   [60];
    double dXsecErrList[60];
    double LR[60];
    
    evt_tree->Branch("ProcIdx"    ,&ProcIdx            ,"ProcIdx/I");
    evt_tree->Branch("nProc"      ,&nProc              ,"nProc/I");
    evt_tree->Branch("dXsec"      ,dXsecList           ,"dXsec[nProc]/D");
    evt_tree->Branch("dXsecErr"   ,dXsecErrList        ,"dXsecErr[nProc]/D");
    evt_tree->Branch("LR"   ,      LR                  ,"LR[nProc]/D");
    
    // SM Models to Process
    int NProcessCalculate=0;
    TVar::Process processList[10];
    
    processList[ NProcessCalculate++]=TVar::WW;  
    //processList[ NProcessCalculate++]=TVar::ZZ;  
    //processList[ NProcessCalculate++]=TVar::Wp_1jet;
    //processList[ NProcessCalculate++]=TVar::Wm_1jet;
  

    TVar::Process ProcInt;
    TVar::Process Vproc;
    
    cdf_event_type cms_event;
    
    //==========================================
    // Loop All Events
    //==========================================
    
    int Ntot = ch->GetEntries();
    if(maxevt<Ntot) Ntot=maxevt;
    printf("Total number of events = %d\n", Ntot);
    
    for(int ievt=0;ievt<Ntot;ievt++){
      
      for(int idx=0;idx<nProc;idx++) {
	dXsecList   [idx] = 0;
	dXsecErrList[idx] = 0;
      }
      
      ch->GetEntry(ievt);            
      
      runNumber      = (int)ch->GetLeaf("runNumber"  )->GetValue();
      eventNumber    = (int)ch->GetLeaf("eventNumber")->GetValue();
      
      
      int Njets             = ch->GetLeaf("Njets"       )->GetValue();
      double MetX   	 = ch->GetLeaf("MetX"	    )->GetValue();
      double MetY	         = ch->GetLeaf("MetY"	    )->GetValue();
      
      double lep1_Px	 = ch->GetLeaf("lep1_Px"      )->GetValue();
      double lep1_Py	 = ch->GetLeaf("lep1_Py"      )->GetValue();
      double lep1_Pz	 = ch->GetLeaf("lep1_Pz"      )->GetValue();
      double lep1_E	 = ch->GetLeaf("lep1_E"       )->GetValue();
      double lep1_Charge    = (int)ch->GetLeaf("lep1_Charge"  )->GetValue();
      double lep1_Type      = (int)ch->GetLeaf("lep1_Type"    )->GetValue();
      
      double lep2_Px	 = ch->GetLeaf("lep2_Px"      )->GetValue();
      double lep2_Py	 = ch->GetLeaf("lep2_Py"      )->GetValue();
      double lep2_Pz	 = ch->GetLeaf("lep2_Pz"      )->GetValue();
      double lep2_E	 = ch->GetLeaf("lep2_E"      )->GetValue();
      double lep2_Charge    = (int)ch->GetLeaf("lep2_Charge"  )->GetValue();
      double lep2_Type      = (int)ch->GetLeaf("lep2_Type"    )->GetValue();
      
      double weight      = (double)ch->GetLeaf("weight"    )->GetValue();
      int dilflavor      = (int)ch->GetLeaf("dilflavor"    )->GetValue();
      
      cms_event.PdgCode[0]=lep1_Type;
      cms_event.PdgCode[1]=lep2_Type;    
      cms_event.MetX = MetX;
      cms_event.MetY = MetY;
      
      if (lep1_Charge>lep2_Charge){
	cms_event.p[0].SetXYZM(lep1_Px,lep1_Py,lep1_Pz,0.0);
	cms_event.p[1].SetXYZM(lep2_Px,lep2_Py,lep2_Pz,0.0);
      }
      else{
	cms_event.p[0].SetXYZM(lep2_Px,lep2_Py,lep2_Pz,0.0);
	cms_event.p[1].SetXYZM(lep1_Px,lep1_Py,lep1_Pz,0.0);
      }
      
      double Mll  =(cms_event.p[0]+cms_event.p[1]).M();
      double Phill=TVector2::Phi_0_2pi(cms_event.p[0].DeltaPhi(cms_event.p[1]));
      if(Phill>TMath::Pi()) Phill=2*TMath::Pi()-Phill;			   
      
      cout<<endl<<endl<<endl;
      cout <<"========================================================="<<endl;
      cout<<"Entry: "<<ievt<<"   Run: "<<runNumber<<"   Event: "<<eventNumber<<endl;
      //      cout<<ProcIdx<<"   "<<TVar::ProcessName(ProcIdx).Data()<<"   "<<lep1_Type<<"   "<<lep2_Type<<endl;

      cout <<"Input: =================================================="<<endl;
      for(int lep=0;lep<2;lep++){
	printf("lep%d : %d %4.4f %4.4f %4.4f %4.4f %4.4f \n",
	       lep, cms_event.PdgCode[lep], cms_event.p[lep].Px(),
	       cms_event.p[lep].Py(), cms_event.p[lep].Pz(),
	       cms_event.p[lep].Energy(), cms_event.p[lep].M());
      } 
      cout <<"========================================================="<<endl;
      
      double Xsec=0;
      double XsecErr=0;
      double Ratio=0;


      for(int iproc=0; iproc<NProcessCalculate; iproc++){ 
	
	ProcInt=processList[iproc];
	Xcal2.SetNcalls(ncalls);
	Xcal2.SetMCHist(ProcInt); 
	
	printf(" Calculate Evt %4i Run %9i Evt %8i Proc %4i %s Lep %4i %4i\n", ievt, runNumber, eventNumber, ProcIdx, TVar::ProcessName(ProcIdx).Data(),lep1_Type,lep2_Type);
	
	Xcal2.NeutrinoIntegrate(ProcInt,cms_event, &Xsec, &XsecErr);
	dXsecList[ProcInt]=Xsec;  
	dXsecErrList[ProcInt]=XsecErr;
	if (Xsec>0) Ratio = XsecErr/Xsec;
	
	if (Ratio>ERRORthreshold){
	  cout <<"IntegrateNTrials "<<" Ncalls " << Xcal2._ncalls<<" "<<Vproc<<" "<< TVar::ProcessName(Vproc)
	       <<" Xsec = " <<  Xsec << " +- " << XsecErr << " ( " << Ratio << " ) " << endl;
	}
      }
            
      //The Higgs calculations now : 
      ProcInt =  Vproc = TVar::HWW;
      
      int i;
      
      int l1 = TVar::HWW160;
      int l2 = TVar::HWW160;
      
      for(i = l1; i<=l2; i++){
	Xcal2.SetMCHist(i);     
	InitializeYields(TVar::ProcessName(i).Data());
	CalculateAcceptance();
	int HiggsMASS;
	
	if (i==TVar::HWW110) HiggsMASS = 110;
	if (i==TVar::HWW120) HiggsMASS = 120;
	if (i==TVar::HWW130) HiggsMASS = 130;
	if (i==TVar::HWW140) HiggsMASS = 140;
	if (i==TVar::HWW150) HiggsMASS = 150;
	if (i==TVar::HWW160) HiggsMASS = 160;
	if (i==TVar::HWW170) HiggsMASS = 170;
	if (i==TVar::HWW180) HiggsMASS = 180;
	if (i==TVar::HWW190) HiggsMASS = 190;
	if (i==TVar::HWW200) HiggsMASS = 200;
	
	if (i==TVar::HWW115) HiggsMASS = 115;
	if (i==TVar::HWW125) HiggsMASS = 125;
	if (i==TVar::HWW135) HiggsMASS = 135;
	if (i==TVar::HWW145) HiggsMASS = 145;
	if (i==TVar::HWW155) HiggsMASS = 155;
	if (i==TVar::HWW165) HiggsMASS = 165;
	if (i==TVar::HWW175) HiggsMASS = 175;
	if (i==TVar::HWW185) HiggsMASS = 185;
	if (i==TVar::HWW195) HiggsMASS = 195;
	
	if (i==TVar::HWW210) HiggsMASS = 210;
	if (i==TVar::HWW220) HiggsMASS = 220;
	if (i==TVar::HWW230) HiggsMASS = 230;
	if (i==TVar::HWW240) HiggsMASS = 240;
	if (i==TVar::HWW250) HiggsMASS = 250;
	if (i==TVar::HWW260) HiggsMASS = 260;
	if (i==TVar::HWW270) HiggsMASS = 270;
	if (i==TVar::HWW280) HiggsMASS = 280;
	if (i==TVar::HWW290) HiggsMASS = 290;
	if (i==TVar::HWW300) HiggsMASS = 300;
	
	
	cout<<i<<"   "<<HiggsMASS<<"   "<<"\n";
	Xcal2.SetNcalls(ncalls);
	Xcal2.SetProcess(TVar::HWW);
	Xcal2.SetHiggsMass(HiggsMASS); 
	
	if ( HiggsMASS<160.8 ) 
	  Xcal2.SetHWWPhaseSpace(TVar::MH);  
	else	
	  Xcal2.SetHWWPhaseSpace(TVar::MHMW);
	
	printf(" Calculate Evt %4i Run %9i Evt %8i Proc %4i %s Lep %4i %4i\n", ievt, runNumber, eventNumber, ProcIdx, TVar::ProcessName(ProcIdx).Data(),lep1_Type,lep2_Type);
	
	Xcal2.NeutrinoIntegrate(ProcInt,cms_event, &Xsec, &XsecErr);
	dXsecList[i]=Xsec; 
	dXsecErrList[i]=XsecErr; 
	if (Xsec>0) Ratio = XsecErr/Xsec;
      }
      
      dXsecList[0]=dXsecList[0];
            
      cout << "START summary ====================" <<endl;
      printf(" Evt %4i/%4i Run %9i Evt %8i Proc %4i %s lep %4i %4i Njets %d \n", ievt,Ntot, runNumber, eventNumber, ProcIdx, TVar::ProcessName(ProcIdx).Data(),lep1_Type,lep2_Type,Njets);
      printf(" MetX %8.8f MetY %8.8f Mll %8.8f Phill %8.8f\n",cms_event.MetX, cms_event.MetY,Mll,Phill);
      printf(" L1 %d ( %8.8f , %8.8f , %8.8f , %8.8f)\n",cms_event.PdgCode[0],cms_event.p[0].Px(),cms_event.p[0].Py(),cms_event.p[0].Pz(),cms_event.p[0].Energy());
      printf(" L2 %d ( %8.8f , %8.8f , %8.8f , %8.8f)\n",cms_event.PdgCode[1],cms_event.p[1].Px(),cms_event.p[1].Py(),cms_event.p[1].Pz(),cms_event.p[1].Energy());
      
      for(int j=0;j<NProcessCalculate;j++){
	TVar::Process proc = processList[j];
	printf("%2i %8s  :", proc,TVar::ProcessName(proc).Data());
	double ratio=0;
	if(dXsecList[processList[j]]>0) ratio=dXsecErrList[processList[j]]/dXsecList[processList[j]];
	cout <<  dXsecList[processList[j]] << " +- " << dXsecErrList[processList[j]] << " ( " << ratio <<" ) "<<endl;
      }
      
      
      for( int ii = l1; ii<=l2; ii++){
	double ratio=0;
	printf("%2i  %8s  :", ii,TVar::ProcessName(ii).Data());
	if(dXsecList[ii]>0) ratio=dXsecErrList[ii]/dXsecList[ii];
	cout <<  dXsecList[ii] << " +- " << dXsecErrList[ii]<< " ( " << ratio <<" ) "<<endl;
      }
      
      cout << "END summary =====================" <<endl;
      
      
      double denom=0.0, numer=0.0;
      double yield_bg=0.0;
      
      // Begin HWW likelihood

      yield_bg = yield[proc_WW][dilflavor];

      cout<<"bg_yield="<<yield_bg<<"\n";
      
      numer=1/(MCFMXsec[proc_HWW160]*acceptance[proc_HWW160][dilflavor]) * dXsecList[TVar::HWW160];
      denom  = numer;
      cout<<" PHHWW160= "<<numer<<"\n";
      
      denom += 1/(MCFMXsec[proc_WW ]*acceptance[proc_WW ][dilflavor]) * dXsecList[proc_WW ] * yield[proc_WW][dilflavor]/yield_bg;
      cout<<" PWW= "<< 1/(MCFMXsec[proc_WW]*acceptance[proc_WW][dilflavor]) * dXsecList[proc_WW] * yield[proc_WW][dilflavor]/yield_bg <<"\n";	
      
      if(denom!=0)
	LR[proc_HWW160]=numer/denom;
      cout<<"LR_HWW160= "<<LR[proc_HWW160]<<"\n";
      evt_tree->Fill();
      
    }//nevent
    
    
    cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
    
    newfile->cd(); 
    evt_tree->Write(); 
    newfile->Close();
    
}  

