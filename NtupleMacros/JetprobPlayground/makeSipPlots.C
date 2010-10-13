#include <assert.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"

void makeSipPlots ()
{
     TFile *f_trkqual = TFile::Open("Results.root");
     assert(f_trkqual != 0);
     TH1 *h_ww_trkd0_good 		= (TH1 *)f_trkqual->Get("ww_trkd00_all");
     TH1 *h_ww_trkDeltaz0_good 		= (TH1 *)f_trkqual->Get("ww_trd0errorGood_all");
     TH1 *h_ww_trknchi2_good 		= (TH1 *)f_trkqual->Get("ww_trknchi20_all");
     TH1 *h_ww_trkvalidhits_good 	= (TH1 *)f_trkqual->Get("ww_trkvalidhits0_all");
     TH1 *h_ww_trkd0_bad 		= (TH1 *)f_trkqual->Get("ww_trkd02_all");
     TH1 *h_ww_trkDeltaz0_bad 		= (TH1 *)f_trkqual->Get("ww_trd0errorBad_all");
     TH1 *h_ww_trknchi2_bad 		= (TH1 *)f_trkqual->Get("ww_trknchi22_all");
     TH1 *h_ww_trkvalidhits_bad 	= (TH1 *)f_trkqual->Get("ww_trkvalidhits2_all");
     h_ww_trkd0_good               ->SetLineColor(kBlack);
     h_ww_trkDeltaz0_good          ->SetLineColor(kBlack);
     h_ww_trknchi2_good            ->SetLineColor(kBlack);
     h_ww_trkvalidhits_good        ->SetLineColor(kBlack);
     h_ww_trkd0_bad                ->SetLineColor(kBlack);
     h_ww_trkDeltaz0_bad           ->SetLineColor(kBlack);
     h_ww_trknchi2_bad             ->SetLineColor(kBlack);
     h_ww_trkvalidhits_bad         ->SetLineColor(kBlack);
     h_ww_trkd0_good               ->SetMarkerStyle(26);
     h_ww_trkDeltaz0_good          ->SetMarkerStyle(26);
     h_ww_trknchi2_good            ->SetMarkerStyle(26);
     h_ww_trkvalidhits_good        ->SetMarkerStyle(26);
     h_ww_trkd0_bad                ->SetMarkerStyle(29);
     h_ww_trkDeltaz0_bad           ->SetMarkerStyle(29);
     h_ww_trknchi2_bad             ->SetMarkerStyle(29);
     h_ww_trkvalidhits_bad         ->SetMarkerStyle(29);
     h_ww_trkvalidhits_good->SetMinimum(0.1);
     TCanvas *c_trkqual = new TCanvas();
     int i = 0;
     c_trkqual->Divide(2,2);
     c_trkqual->cd(++i); 	gPad->SetLogy(); h_ww_trkd0_good	->Draw();	h_ww_trkd0_bad         ->Draw("same");       
     c_trkqual->cd(++i); 	gPad->SetLogy(); h_ww_trkDeltaz0_good	->Draw();	h_ww_trkDeltaz0_bad    ->Draw("same");       
     c_trkqual->cd(++i); 	gPad->SetLogy(); h_ww_trknchi2_good	->Draw();	h_ww_trknchi2_bad      ->Draw("same");       
     c_trkqual->cd(++i); 	gPad->SetLogy(); h_ww_trkvalidhits_good	->Draw();	h_ww_trkvalidhits_bad  ->Draw("same");       
     c_trkqual->Print("trkqual.eps");
     c_trkqual->Print("trkqual.root");
     
}
