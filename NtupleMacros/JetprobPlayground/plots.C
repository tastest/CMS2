#include "TCanvas.h"
#include "TDirectory.h"

void plots_tracks ()
{
     TCanvas *c1 = new TCanvas;
     c1->Divide(4, 1);
     for (int i = 0; i < 4; ++i) {
	  c1->cd(i + 1);
	  gPad->SetLogy();
// 	  printf("%s\n", Form("ttbar_trd0sigByPt%d_all", i));
// 	  TObject *obj = gDirectory->FindObject(Form("ttbar_trd0sigByPt%d_all", i));
	  gDirectory->Get(Form("ttbar_trd0sigByPt%d_all", i))->Draw();
	  gDirectory->Get(Form("ww_trd0sigByPt%d_all", i))->Draw("same");
     }
     TCanvas *c2 = new TCanvas;
     c2->Divide(4, 1);
     for (int i = 0; i < 4; ++i) {
	  c2->cd(i + 1);
	  gPad->SetLogy();
	  gDirectory->Get(Form("ttbar_trd0sigByNtrks%d_all", i))->Draw();
	  gDirectory->Get(Form("ww_trd0sigByNtrks%d_all", i))->Draw("same");
     }
}

void plots_jets ()
{
     TCanvas *c1 = new TCanvas;
     c1->Divide(4, 1);
     for (int i = 0; i < 4; ++i) {
	  c1->cd(i + 1);
	  gPad->SetLogy();
// 	  printf("%s\n", Form("ttbar_jetd0sigByPt%d_all", i));
// 	  TObject *obj = gDirectory->FindObject(Form("ttbar_jetd0sigByPt%d_all", i));
	  gDirectory->Get(Form("ttbar_jetd0sigByPt%d_all", i))->Draw();
	  gDirectory->Get(Form("ww_jetd0sigByPt%d_all", i))->Draw("same");
     }
     TCanvas *c2 = new TCanvas;
     c2->Divide(4, 1);
     for (int i = 0; i < 4; ++i) {
	  c2->cd(i + 1);
	  gPad->SetLogy();
	  gDirectory->Get(Form("ttbar_jetd0sigByNtrks%d_all", i))->Draw();
	  gDirectory->Get(Form("ww_jetd0sigByNtrks%d_all", i))->Draw("same");
     }
}
