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
  myHisto = new TH1F("myHisto","d#phi_{t#bar{t}}", 600, 0.,3.1415926536);
  myHisto_neg = new TH1F("myHisto_neg","d#phi_{t#bar{t}}", 600, 0.,3.1415926536);
  myHisto2 = new TH1F("myHisto2","t#bar{t} p_{T}", 600, 0.,1200);
  myHisto2_neg = new TH1F("myHisto2_neg","t#bar{t} p_{T}", 600, 0.,1200);
  myHisto3 = new TH1F("myHisto3","#theta_{t#bar{t}} ", 600, 0.,3.1415926536);
  myHisto4 = new TH1F("myHisto4","y_{t#bar{t}} ", 600, -3.,3.);

}

void user::Execute(const SampleFormat& sample, const EventFormat& event)
{

  const MCParticleFormat* top	= 0; 
  const MCParticleFormat* w	= 0; 
  const MCParticleFormat* top1 = 0;
  const MCParticleFormat* top2 = 0;
  
  for (unsigned int i=0;i<event.mc()->particles().size();i++) {

    const MCParticleFormat* prt=&(event.mc()->particles()[i]);
    // The top is not a final state particle 
    if (PHYSICS->IsFinalState(prt)) continue;

    //const MCParticleFormat* mother = prt->mother1();   
    //const MCParticleFormat* grandmother = mother->mother1();

    if (( prt->pdgid()==6) )   top1 = prt; 
    if (( prt->pdgid()==-6) )   top2 = prt; 
    
  }

  // Rejection of the event if top(s) not found 
  if ((top1==0)||(top2==0)) {
    WARNING << "top(s) not found!" << std::endl;
    return;
  }

// Computing the observable; filling the histogram

  TLorentzVector topplus_genp_p4(0,0,0,0), topminus_genp_p4(0,0,0,0);
  topplus_genp_p4.SetXYZT( top1->px(), top1->py(), top1->pz(), top1->e() );
  topminus_genp_p4.SetXYZT( top2->px(), top2->py(), top2->pz(), top2->e() );

  //std::cout<<event.mc()->weight()<<std::endl;
  myHisto -> Fill(fabs(TVector2::Phi_mpi_pi(top1->phi()-top2->phi())), event.mc()->weight()/fabs(event.mc()->weight()));
  //myHisto2 -> Fill(  sqrt( (top1->px() + top2->px())*(top1->px() + top2->px()) + (top1->py() + top2->py())*(top1->py() + top2->py()) ) , event.mc()->weight()/fabs(event.mc()->weight())  );
  myHisto2 -> Fill( (topplus_genp_p4+topminus_genp_p4).Pt() , event.mc()->weight()/fabs(event.mc()->weight())  );
  myHisto3 -> Fill( (top1->angle(top2)) , event.mc()->weight()/fabs(event.mc()->weight()) );
  myHisto4 -> Fill( (topplus_genp_p4+topminus_genp_p4).Rapidity() , event.mc()->weight()/fabs(event.mc()->weight())  );

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
  TCanvas* myCanvas = new TCanvas("myCanvas","myCanvas");
  myCanvas->cd(1);
  myCanvas->SetLogy();
// Setting background color 
  myHisto->SetFillColor(kBlue);
// Normalization of the histogram: L = 10 fb-1 
//  double nrm = summary.mc()->xsection() * 10000. / static_cast<float>(summary.nevents()); 
  //for (int ibin=1; ibin<=myHisto->GetNbinsX(); ibin++) myHisto->SetBinContent(ibin, fabs(myHisto->GetBinContent(ibin)));
  for (int ibin=1; ibin<=myHisto->GetNbinsX(); ibin++) if(myHisto->GetBinContent(ibin)<0) myHisto_neg->SetBinContent(ibin, fabs(myHisto->GetBinContent(ibin)));
  myHisto_neg->SetFillColor(kBlue+3);

  //double Afb=(myHisto->Integral(4,6)-myHisto->Integral(1,3))/(myHisto->Integral(4,6)+myHisto->Integral(1,3));
  //double Afb=(myHisto->Integral(31,60,"width")-myHisto->Integral(1,30,"width"))/(myHisto->Integral(31,60,"width")+myHisto->Integral(1,30,"width"));
  //std::cout<< "Asymmetry = " << Afb<<"\n";

// Setting axis title 
  myHisto->GetXaxis()->SetTitle("d#phi_{t#bar{t}}");
// Drawing histogram 
  myHisto->SetMinimum(0.9);
  myHisto->SetMaximum(30.*myHisto->GetMaximum());
  myHisto->Draw();
  myHisto_neg->Draw("same");
  
  //TLatex* lat=new TLatex(0.5,5000, Form("Asym = %f",Afb));
  //lat->SetTextFont(42);
  //lat->Draw();

// Saving plot 
  myCanvas->SaveAs((outputName_+"_dPhi.pdf").c_str());
  myCanvas->SaveAs((outputName_+"_dPhi.C").c_str());

  
  TCanvas* myCanvas2 = new TCanvas("myCanvas2","myCanvas2");
  myCanvas2->cd(1);
  myCanvas2->SetLogy();
  myHisto2->SetFillColor(kGreen);
  //for (int ibin=1; ibin<=myHisto2->GetNbinsX(); ibin++) myHisto2->SetBinContent(ibin, fabs(myHisto2->GetBinContent(ibin)));
  for (int ibin=1; ibin<=myHisto2->GetNbinsX(); ibin++) if(myHisto2->GetBinContent(ibin)<0) myHisto2_neg->SetBinContent(ibin, fabs(myHisto2->GetBinContent(ibin)));
  myHisto2_neg->SetFillColor(kGreen+3);
  myHisto2->GetXaxis()->SetTitle("t#bar{t} p_{T}");
  myHisto2->SetMinimum(0.9);
  myHisto2->SetMaximum(30.*myHisto2->GetMaximum());  
  myHisto2->Draw();
  myHisto2_neg->Draw("same");
  
  myCanvas2->SaveAs((outputName_+"_pTtt.pdf").c_str());
  myCanvas2->SaveAs((outputName_+"_pTtt.C").c_str());
  

  TCanvas* myCanvas3 = new TCanvas("myCanvas3","myCanvas3");
  myCanvas3->cd(1);
  myHisto3->SetFillColor(kRed);
  //for (int ibin=1; ibin<=myHisto3->GetNbinsX(); ibin++) myHisto3->SetBinContent(ibin, fabs(myHisto3->GetBinContent(ibin)));
  myHisto3->GetXaxis()->SetTitle("#theta_{t#bar{t}}");
  myHisto3->SetMinimum(0.);
  myHisto3->SetMaximum(1.5*myHisto3->GetMaximum());  
  myHisto3->Draw();
  
  myCanvas3->SaveAs((outputName_+"_thetatt.pdf").c_str());
  myCanvas3->SaveAs((outputName_+"_thetatt.C").c_str());

  
  TCanvas* myCanvas4 = new TCanvas("myCanvas4","myCanvas4");
  myCanvas4->cd(1);
  myHisto4->SetFillColor(kYellow);
  //for (int ibin=1; ibin<=myHisto4->GetNbinsX(); ibin++) myHisto4->SetBinContent(ibin, fabs(myHisto4->GetBinContent(ibin)));
  myHisto4->GetXaxis()->SetTitle("y_{t#bar{t}}");
  //std::cout<<"min "<<myHisto4->GetMinimum()<<std::endl;
  myHisto4->SetMinimum(0.);
  myHisto4->SetMaximum(1.5*myHisto4->GetMaximum());  
  myHisto4->Draw();
  
  myCanvas4->SaveAs((outputName_+"_ytt.pdf").c_str());
  myCanvas4->SaveAs((outputName_+"_ytt.C").c_str());

}
