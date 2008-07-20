// Make the stacks and then browse them interactively
// if (makePictures==true), saves all plots to out/stacks.ps
// Also: make out_stacks_xx.png where xx is the page number 
// on the stack.ps file
//
// keep2D=false to skip the annoying 2D histograms
//
void browseStacks( bool makePictures=false, bool wait=true ) {

  bool keep2D=false;

// Find out what the names of the existing histograms are
// The histogram names are XX_YY_ZZ, where XX is the sample,
// eg, "tt", YY is the actual name, ZZ is the final state, eg, "ee"
  TObjArray* myNames = getMyHistosNames("wz","all",keep2D);


// Now loop over histograms, and make stacks
   TCanvas *c = new TCanvas();
   c->Divide(3,2);
  // define old categories for WZ
  const int num_suffix = 5;
  char *suffix[num_suffix];
  suffix[0] = "zee_e";
  suffix[1] = "zee_m";
  suffix[2] = "zmm_e";
  suffix[3] = "zmm_m";
  suffix[4] = "all";
   if (makePictures) c->Print("out/stacks.ps[");
   for (int i=0; i<myNames->GetEntries(); i++) {
     for (int sample=0; sample<num_suffix; sample++) {
       hist::stack(Form("st_%s_%s",myNames->At(i)->GetName(),suffix[sample]),
                    Form("%s_%s$",myNames->At(i)->GetName(), suffix[sample]));
        THStack* thisStack = (THStack*) gROOT->FindObjectAny(Form("st_%s_%s", myNames->At(i)->GetName(), suffix[sample]));
	TLegend* thisLeg = hist::legend(thisStack, "lpf", 0, 0);
        c->cd(sample+1);
        thisStack->Draw("hist");
        thisLeg->Draw();
        c->Update();
     }
     if (makePictures) {
       c->Print("out/stacks.ps");
       //       c->Print(Form("out/stacks_%d.png",i+1));
       c->Print(Form("out/stacks_%s.png",myNames->At(i)->GetName()));
     }
     if (wait) {
       cout << "Enter carriage return for the next set of plots....q to quit" << endl;
       char in = getchar();
       if (in == 'q') break;
     }
   }
   if (makePictures) c->Print("out/stacks.ps]");
}
