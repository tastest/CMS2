
#include "histtools.h"

typedef TH1F H;

H cumulate (const H &in, bool increasing) 
{
     H h_out(in.GetName() + TString("tmp"), in.GetTitle(), in.GetNbinsX(), 
	     in.GetBinLowEdge(1), in.GetBinLowEdge(in.GetNbinsX() + 1));
     h_out.Sumw2();
     h_out.SetFillColor(in.GetFillColor());
     h_out.SetFillStyle(in.GetFillStyle());
     h_out.SetLineStyle(in.GetLineStyle());
     h_out.SetLineColor(in.GetLineColor());
     double sum = 0;
     double err2 = 0;
     if (increasing) {
	  for (int j = 0; j <= in.GetNbinsX() + 1; ++j) {
	       sum += in.GetBinContent(j);
           err2 += in.GetBinError(j)*in.GetBinError(j);
	       h_out.SetBinContent(j, sum);
           h_out.SetBinError(j, sqrt(err2));
	  }
     } else {
	  for (int j = in.GetNbinsX() + 1; j >= 0; --j) {
	       sum += in.GetBinContent(j);
           err2 += in.GetBinError(j)*in.GetBinError(j);
           h_out.SetBinContent(j, sum);
           h_out.SetBinError(j, sqrt(err2));
	  }
     }
     return h_out;
}

TGraph eff_rej (const H &signal, H &background, bool normalize, bool increasing)
{
     H sig = *(TH1F*)signal.Clone("h_tmp_s");
     if (normalize)
	  sig.Scale(1 / sig.Integral(0, sig.GetNbinsX() + 1));
     H bg = *(TH1F*)background.Clone("h_tmp_bg");
     if (normalize)
	  bg.Scale(1 / bg.Integral(0, bg.GetNbinsX() + 1));
     H sig_cum = cumulate(sig, increasing);
     H bg_cum = cumulate(bg, increasing);
     TGraph ret(signal.GetNbinsX());
     for (int i = 1; i <= signal.GetNbinsX(); ++i) {
	  const double x = sig_cum.GetBinCenter(i);
 	  const double sig = sig_cum.GetBinContent(i);
	  const double bg = bg_cum.GetBinContent(bg_cum.FindBin(x));
	  ret.SetPoint(i - 1, sig, bg); // gotta love offsets
//  	  printf("point %d: %f sig, %f bg\n", i, sig, bg);
     }
     return ret;
}
