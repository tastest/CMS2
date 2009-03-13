void printNJets( bool latex=false){
  char* suffix[4];
  suffix[0] = "ee";
  suffix[1] = "mm";
  suffix[2] = "em";
  suffix[3] = "all";


  if (latex) {
    std::cout << "\\begin{table}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|l|c|c|c|} \\hline" << std::endl;
  }

     
  for (int sample=0; sample<4; sample++) {
    
        hist::stack(Form("st_hnJet_%s", suffix[sample]), Form("hnJet_%s$", suffix[sample]));
    THStack* thisStack = (THStack*) 
      gROOT->FindObjectAny(Form("st_hnJet_%s", suffix[sample]));
    //    std::cout<<"Found "<<thisStack->GetName()<<std::endl;

    int nHists = thisStack->GetHists()->GetSize();
    if (latex) {
      std::cout << "\\hline \\hline" << std::endl;
      std::cout << "\\multicolumn{4}{|l|}{" << suffix[sample] << " final state} \\"<<"\\  " << std::endl;
      std::cout << "   & $N_{jets}=0$       & $N_{jets}=1$      & $N_{jets} \\geq 2$ \\"<<"\\ \\hline" <<std::endl;
    } else {
      std::cout<<"===================================================="<<std::endl;
      std::cout<<suffix[sample]<<std::endl;
      std::cout<<" |   sample  |        nJet = 0        |       nJet = 1         |       nJet >= 2       |"<<std::endl;
    }
    for(int iH=0; iH< nHists; ++iH){
      TH1F* h1F = (TH1F*)(thisStack->GetHists()->At(iH));
      TObjArray* sampleNs = TString(h1F->GetName()).Tokenize("_");
      double n0 = h1F->GetBinContent(1);
      double n0E = h1F->GetBinError(1);
      double n1 = h1F->GetBinContent(2);
      double n1E = h1F->GetBinError(2);
      int nBins = h1F->GetNbinsX();
      double n2 = 0; 
      for (int i=3; i<= nBins+1; ++i) n2+= h1F->GetBinContent(i);
      double n2E = 0; 
      for (int i=3; i<= nBins+1; ++i) n2E += h1F->GetBinError(i)*h1F->GetBinError(i);
      n2E = sqrt(n2E);
      if (latex) {
	std::cout< <Form("%9s",sampleNs->At(0)->GetName()) << " & $"
		 <<Form("%6.1f",n0) <<" \\pm "<< Form("%6.1f",n0E) <<"$ & $"
	         <<Form("%6.1f",n1) <<" \\pm "<< Form("%6.1f",n1E) <<"$ & $"
		 <<Form("%6.1f",n2) <<" \\pm "<< Form("%6.1f",n2E) <<"$ \\"<<"\\"
	         <<std::endl;
      } else {
	std::cout<<" | "<<Form("%9s",sampleNs->At(0)->GetName())
	         <<" | "<<Form("%6.1f",n0) <<" &plusmn; "<<Form("%6.1f",n0E)
	         <<" | "<<Form("%6.1f",n1) <<" &plusmn; "<<Form("%6.1f",n1E)
	         <<" | "<<Form("%6.1f",n2) <<" &plusmn; "<<Form("%6.1f",n2E)
	         <<" | "<<std::endl;
      }
    }
   if(!latex)
	std::cout << " " << std::endl;
  }
 
  if (latex) {
    std::cout << "\\hline " << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\caption{ Put your caption here}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\label{tab:dummyLabel}" <<std::endl;
    std::cout << "\\end{table}" << std::endl;
  }
}
