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
  proc_Wpg,
  proc_Wmg,
  proc_ZZ,
  proc_HWW120,
  proc_HWW130,
  proc_HWW140,
  proc_HWW150,
  proc_HWW160,
  proc_HWW170,
  proc_HWW180,
  proc_HWW190,
  proc_HWW200,
  proc_HWW210,
  proc_HWW220,
  proc_HWW230,
  proc_HWW250,
  proc_HWW300,
  kNProc
};

enum { kNDilep=3 }; 


float lumi=1000.0;
float yield [kNProc][kNDilep];
float acceptance [kNProc][kNDilep];

float BR [kNProc][kNDilep] = {{ (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1784)/3, (2.0+0.1736+0.1784)/3, (1+0.1784)/3},
			      { (1.0+0.1784)/3, (2.0+0.1736+0.1784)/3, (1+0.1784)/3},
			      { (1.0+0.1784)/3, (2.0+0.1736+0.1784)/3, (1+0.1784)/3}, //Wpgamma placeholder
                              { (1.0+0.1784)/3, (2.0+0.1736+0.1784)/3, (1+0.1784)/3}, //Wmgamma placefolder 
			      { 2*0.2*(0.0363+0.0370*0.1736*0.1736), 0.0, 2*0.2*(0.0363+0.0370*0.1784*0.1784) },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 },
			      { (1.0+0.1736)*(1.0+0.1736)/9, (1.0+0.1736)*(1.0+0.1784)/9*2, (1.0+0.1784)*(1.0+0.1784)/9 }
};

// The NLO xsec includes the W->l BR, where l includes e/mu/tau
// The values are obtained through http://ceballos.web.cern.ch/ceballos/hwwlnln/sigma_hww_2l2n_ecm7tev.txt
float NLOXsec[kNProc] = {  4.5*0.919, 31314.0/2.0, 31314.0/2.0, 31314.0/2.0, 31314.0/2.0, 5.9, 0.249642, 0.452090, 0.641773, 0.770471, 0.866443, 0.782962, 0.659328, 0.486486, 0.408305, 0.358465, 0.321398, 0.290454, 0.243724, 0.175652};
// Note that MCFMXsec is obtained from a standalone MCFM calculations with no generator cuts applied, excluding the W->l B.R.
float MCFMXsec[kNProc] = { 28.4, 11270, 11270.0,  11270, 11270.0, 4.3, 0.6619, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25, 3.25};



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

void InitializeYields(){

  TString FileName;
  FileName = "../Util.root";
  TFile *f = TFile::Open(FileName, "READ");
  assert(f);
  gROOT->cd();

  TString histName;
  
  TString procName[kNProc] = {  "WW" , "Wpj", "Wmj", "Wpg", "Wmg", "ZZ", "HWW120", "HWW130", "HWW140", "HWW150", "HWW160",
				"HWW170", "HWW180", "HWW190", "HWW200", "HWW210", "HWW220", "HWW230", "HWW250", "HWW300"};
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

  ERRORthreshold=Error;
  int process=TVar::HWW;
  int maxevt=nev;

  InitializeYields();
  CalculateAcceptance();

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
      double MetX   	    = ch->GetLeaf("MetX"	    )->GetValue();
      double MetY	    = ch->GetLeaf("MetY"	    )->GetValue();
      
      double lep1_Px	 = ch->GetLeaf("lep1_Px"      )->GetValue();
      double lep1_Py	 = ch->GetLeaf("lep1_Py"      )->GetValue();
      double lep1_Pz	 = ch->GetLeaf("lep1_Pz"      )->GetValue();
      double lep1_E	 = ch->GetLeaf("lep1_E"       )->GetValue();
      double lep1_Charge    = (int)ch->GetLeaf("lep1_Charge"  )->GetValue();
      int lep1_Type      = (int)ch->GetLeaf("lep1_Type"    )->GetValue();
      
      double lep2_Px	 = ch->GetLeaf("lep2_Px"      )->GetValue();
      double lep2_Py	 = ch->GetLeaf("lep2_Py"      )->GetValue();
      double lep2_Pz	 = ch->GetLeaf("lep2_Pz"      )->GetValue();
      double lep2_E	 = ch->GetLeaf("lep2_E"      )->GetValue();
      double lep2_Charge    = (int)ch->GetLeaf("lep2_Charge"  )->GetValue();
      int lep2_Type      = (int)ch->GetLeaf("lep2_Type"    )->GetValue();
      
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

      cout <<"Input: =================================================="<<endl;
      for(int lep=0;lep<2;lep++){
	printf("lep%d : %d %4.4f %4.4f %4.4f %4.4f %4.4f \n",
	       lep, cms_event.PdgCode[lep], cms_event.p[lep].Px(),
	       cms_event.p[lep].Py(), cms_event.p[lep].Pz(),
	       cms_event.p[lep].Energy(), cms_event.p[lep].M());
      } 
      cout <<"========================================================="<<endl;
      

      // SM Models to Process                                                                                                                            
      int NProcessCalculate=0;
      TVar::Process processList[10];
      processList[ NProcessCalculate++]=TVar::WW;
      processList[ NProcessCalculate++]=TVar::HWW160;

      //processList[ NProcessCalculate++]=TVar::ZZ;                                                                                                  
      //processList[ NProcessCalculate++]=TVar::Wp_1jet;                                                                                             
      //processList[ NProcessCalculate++]=TVar::Wm_1jet;                                                                                              

      TVar::Process ProcInt;
      TVar::Process Vproc;

      double Xsec=0;
      double XsecErr=0;
      double Ratio=0;

      int HiggsMASS[20]={0,0,0,0,0,0,120,130,140,150,160,170,180,190,200,210,220,230,250,300};

      for(int iproc=0; iproc<NProcessCalculate; iproc++){ 
	
	ProcInt=processList[iproc];
	Xcal2.SetNcalls(ncalls);
	Xcal2.SetMCHist(ProcInt); 
	
	printf(" Calculate Evt %4i Run %9i Evt %8i Proc %4i %s Lep %4i %4i\n", ievt, runNumber, eventNumber, ProcIdx, TVar::ProcessName(ProcIdx).Data(),lep1_Type,lep2_Type);

	if ((processList[iproc]>=TVar::HWW120) && (processList[iproc]<=TVar::HWW300)){
	  Xcal2.SetProcess(TVar::HWW);
	  Xcal2.SetHiggsMass(HiggsMASS[processList[iproc]]);
	  cout<<HiggsMASS[processList[iproc]]<<"\n";
	  
	  if ( HiggsMASS[processList[iproc]]<160.8 )
	    Xcal2.SetHWWPhaseSpace(TVar::MH);
	  else
	    Xcal2.SetHWWPhaseSpace(TVar::MHMW);

	  Xcal2.NeutrinoIntegrate(TVar::HWW,cms_event, &Xsec, &XsecErr);
	} 
	else Xcal2.NeutrinoIntegrate(ProcInt,cms_event, &Xsec, &XsecErr);
	 
	dXsecList[ProcInt]=Xsec;  
	dXsecErrList[ProcInt]=XsecErr;
	if (Xsec>0) Ratio = XsecErr/Xsec;
	
	if (Ratio>ERRORthreshold){
	  cout <<"IntegrateNTrials "<<" Ncalls " << Xcal2._ncalls<<" "<<Vproc<<" "<< TVar::ProcessName(Vproc)
	       <<" Xsec = " <<  Xsec << " +- " << XsecErr << " ( " << Ratio << " ) " << endl;
	}
      }
            

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
      cout << "END summary =====================" <<endl;
      
      
      double denom=0.0, numer=0.0;
      double yield_bg=0.0;
      
      // Begin HWW likelihood

      for(int k=proc_HWW120; k<proc_HWW300; k++){

        yield_bg = yield[proc_WW][dilflavor];

	cout<<"bg_yield="<<yield_bg<<"\n";
      
	numer=1/(MCFMXsec[k]*acceptance[k][dilflavor]) * dXsecList[k];
	denom  = numer;
	cout<<" PHWW= "<<numer<<"\n";
      
	denom += 1/(MCFMXsec[proc_WW ]*acceptance[proc_WW ][dilflavor]) * dXsecList[proc_WW ] * yield[proc_WW][dilflavor]/yield_bg;
	cout<<" PWW= "<< 1/(MCFMXsec[proc_WW]*acceptance[proc_WW][dilflavor]) * dXsecList[proc_WW] * yield[proc_WW][dilflavor]/yield_bg <<"\n";	
	
	if(denom!=0)
	  LR[k]=numer/denom;
	//cout<<"LR_HWW["<<k<<"]= "<<LR[k]<<"\n";
      }      


	evt_tree->Fill();
      
    }//nevent
    
    
    cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
    
    newfile->cd(); 
    evt_tree->Write(); 
    newfile->Close();
    
}  

