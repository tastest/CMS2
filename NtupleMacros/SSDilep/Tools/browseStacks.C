// Make the stacks and then browse them interactively
// if (makePictures==true), saves all plots to out/stacks.ps
// Also: make out_stacks_xx.png where xx is the page number 
// on the stack.ps file
//
// keep2D=false to skip the annoying 2D histograms
//
void browseStacks( bool makePictures=false, bool wait=true, const char* dir = "out" ) {

  // check if dir exists, if not, create it
  gSystem->mkdir(dir);

  bool keep2D=false;

// Find out what the names of the existing histograms are
// The histogram names are XX_YY_ZZ, where XX is the sample,
// eg, "tt", YY is the actual name, ZZ is the final state, eg, "ee"
  TObjArray* myNames = getMyHistosNames("tt","ee",keep2D);


// Now loop over histograms, and make stacks
   TCanvas *c = new TCanvas();
   c->Divide(2,2);
   char* suffix[4];
   suffix[0] = "ee";
   suffix[1] = "mm";
   suffix[2] = "em";
   suffix[3] = "all";
   if (makePictures) c->Print(Form("%s/stacks.ps[",dir));
   for (int i=0; i<myNames->GetEntries(); i++) {
     //  cout << myNames->At(i)->GetName() << endl;
     for (int sample=0; sample<4; sample++) {
       hist::stack(Form("st_%s_%s",myNames->At(i)->GetName(),suffix[sample]),
                    Form("%s_%s$",myNames->At(i)->GetName(), suffix[sample]));
        THStack* thisStack = (THStack*) gROOT->FindObjectAny(
                              Form("st_%s_%s", myNames->At(i)->GetName(), suffix[sample]));
        c->cd(sample+1);
	if ( thisStack ) {
	  TLegend* thisLeg = hist::legend(thisStack, "lpf", 0, 0);
	  thisStack->Draw("hist");
	}
        thisLeg->Draw();
        c->Update();
     }
     if (makePictures) {
       c->Print(Form("%s/stacks.ps",dir));
       //       c->Print(Form("out/stacks_%d.png",i+1));
       c->Print(Form("%s/stacks_%s.png",dir,myNames->At(i)->GetName()));
       c->Print(Form("%s/stacks_%s.root",dir,myNames->At(i)->GetName()));
     }
     if (wait) {
       cout << "Enter carriage return for the next set of plots....q to quit" << endl;
       char in = getchar();
       if (in == 'q') break;
     }
   }
   if (makePictures) c->Print(Form("%s/stacks.ps]",dir));
}

void browseStacksTrilep( bool makePictures=false, bool wait=true, const char* dir = "out" ) {

  bool keep2D=false;

// Find out what the names of the existing histograms are
// The histogram names are XX_YY_ZZ, where XX is the sample,
// eg, "tt", YY is the actual name, ZZ is the final state, eg, "ee"
  TObjArray* myNames = getMyHistosNames("tt","mpmpmp",keep2D);


// Now loop over histograms, and make stacks
   TCanvas *c = new TCanvas();
   c->Divide(6,4);
  const unsigned int allBuckets = 20;
  char *suffix[allBuckets+1];
  suffix[0]  = "mpmpmp";
  suffix[1]  = "mpmpmm";
  suffix[2]  = "mpmmmm";
  suffix[3]  = "mmmmmm";
  suffix[4]  = "mpmpep";
  suffix[5]  = "mpmpem";
  suffix[6]  = "mpmmep";
  suffix[7]  = "mpmmem";
  suffix[8]  = "mmmmep";
  suffix[9]  = "mmmmem";
  suffix[10] = "mpepep";
  suffix[11] = "mpepem";
  suffix[12] = "mpemem";
  suffix[13] = "mmepep";
  suffix[14] = "mmepem";
  suffix[15] = "mmemem";
  suffix[16] = "epepep";
  suffix[17] = "epepem";
  suffix[18] = "epemem";
  suffix[19] = "ememem";
  suffix[20] = "all";
  if (makePictures) c->Print(Form("%s/stacks.ps[",dir));
   for (int i=0; i<myNames->GetEntries(); i++) {
     //  cout << myNames->At(i)->GetName() << endl;
     for (int sample=0; sample<=allBuckets; sample++) {
       hist::stack(Form("st_%s_%s",myNames->At(i)->GetName(),suffix[sample]),
                    Form("%s_%s$",myNames->At(i)->GetName(), suffix[sample]));
        THStack* thisStack = (THStack*) gROOT->FindObjectAny(
                              Form("st_%s_%s", myNames->At(i)->GetName(), suffix[sample]));
	TLegend* thisLeg = hist::legend(thisStack, "lpf", 0, 0);
        c->cd(sample+1);
        thisStack->Draw("hist");
        thisLeg->Draw();
        c->Update();
     }
     if (makePictures) {
       c->Print(Form("%s/stacks.ps",dir));
       //       c->Print(Form("out/stacks_%d.png",i+1));
       c->Print(Form("%s/stacks_%s.png",dir,myNames->At(i)->GetName()));
     }
     if (wait) {
       cout << "Enter carriage return for the next set of plots....q to quit" << endl;
       char in = getchar();
       if (in == 'q') break;
     }
   }
   if (makePictures) c->Print(Form("%s/stacks.ps]",dir));
}
