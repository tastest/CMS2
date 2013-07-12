#include "Analysis/user.h"
#include "TLatex.h"

void user::Initialize()
{
  //Initializing Physics services 
  PHYSICS->mcConfig().Reset();
  // Initializing the histogram 
  double xbins[7];
  xbins[0]=acos(1.0);  xbins[1]=acos(0.6);  xbins[2]=acos(0.3); xbins[3]=acos(0.0); 
  xbins[4]=acos(-0.3); xbins[5]=acos(-0.6); xbins[6]=acos(-1.0);
  myHisto = new TH1F("myHisto","d#phi_{ll}", 6, xbins);
}

void user::Execute(const SampleFormat& sample, const EventFormat& event)
{

  const MCParticleFormat* top	= 0; 
  const MCParticleFormat* w	= 0; 
  const MCParticleFormat* lepton1 = 0;
  const MCParticleFormat* lepton2 = 0;
  
  for (unsigned int i=0;i<event.mc()->particles().size();i++) {

    const MCParticleFormat* prt=&(event.mc()->particles()[i]);
    // The lepton is a final state particle 
    if (!PHYSICS->IsFinalState(prt)) continue;

    const MCParticleFormat* mother = prt->mother1();   
    const MCParticleFormat* grandmother = mother->mother1();

    if (mother==0) continue; 
    if (std::abs(mother->pdgid())!=24) continue;

    //if (( prt->pdgid()==11) ||  ( prt->pdgid()==13) ||  ( prt->pdgid()==15))   lepton1 = prt; 
    //if (( prt->pdgid()==-11) ||  ( prt->pdgid()==-13) ||  ( prt->pdgid()==-15))   lepton2 = prt; 

    if (( prt->pdgid()==11) ||  ( prt->pdgid()==13) )   lepton1 = prt; 
    if (( prt->pdgid()==-11) ||  ( prt->pdgid()==-13) )   lepton2 = prt; 

    /*
    // Lepton selection based on the PDG-id 
    if ( std::abs(prt->pdgid())!=11 && std::abs(prt->pdgid())!=13 ) continue;
    const MCParticleFormat* mother = prt->mother1();    
    if (mother==0) continue; 
    if (std::abs(mother->pdgid())!=24) continue;
    // Getting the grand-mother of the lepton and 
    // checking if it is a top quark 
     const MCParticleFormat* grandmother = mother->mother1(); 
     if (grandmother==0) continue; 
     if (std::abs(grandmother->pdgid())!=6) continue;
     // Saving the selected particles 
     lepton = prt; 
     w	= mother; 
     top	= grandmother;
     // Particles are found: breaking the loop 
     break;
    */
  }

  // Rejection of the event if the decay chain is not found 
  if ((lepton1==0)||(lepton2==0)) {
    //WARNING << "dilepton decay chain not found!" << std::endl;
    return;
  }

// Computing the observable; filling the histogram
  //myHisto -> Fill(fabs(TVector2::Phi_mpi_pi(lepton1->phi()-lepton2->phi())));
  
  //std::cout<<event.mc()->weight()<<std::endl;
  myHisto -> Fill(fabs(TVector2::Phi_mpi_pi(lepton1->phi()-lepton2->phi())), event.mc()->weight() );

}


void user::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  // Color of the canvas background 
  gStyle->SetCanvasColor(0);
// Turning off the border lines of the canvas 
  gStyle->SetCanvasBorderMode(0);
// Configuring statistics printing 
  gStyle->SetOptStat(111110);
// Creating the output root file 
  TCanvas* myCanvas = new TCanvas("myCanvas","");
// Setting background color 
  myHisto->SetFillColor(kBlue);
// Normalization of the histogram: L = 10 fb-1 
//  double nrm = summary.mc()->xsection() * 10000. / static_cast<float>(summary.nevents()); 
  for (int ibin=1; ibin<=myHisto->GetNbinsX(); ibin++)
    myHisto->SetBinContent(ibin, myHisto->GetBinContent(ibin)/myHisto->GetBinWidth(ibin));

  //double Afb=(myHisto->Integral(4,6)-myHisto->Integral(1,3))/(myHisto->Integral(4,6)+myHisto->Integral(1,3));
  double Afb=(myHisto->Integral(4,6,"width")-myHisto->Integral(1,3,"width"))/(myHisto->Integral(4,6,"width")+myHisto->Integral(1,3,"width"));
  std::cout<< "Asymmetry = " << Afb<<"\n";

// Setting axis title 
  myHisto->GetXaxis()->SetTitle("d#phi_{ll}");
// Drawing histogram 
  myHisto->SetMinimum(0.0);
  myHisto->SetMaximum(1.5*myHisto->GetMaximum());
  myHisto->Draw();
  
  TLatex* lat=new TLatex(0.5,5000, Form("Asym = %f",Afb));
  lat->SetTextFont(42);
  lat->Draw();

// Saving plot 
  myCanvas->SaveAs((outputName_+".pdf").c_str());
  myCanvas->SaveAs((outputName_+".C").c_str());

}
