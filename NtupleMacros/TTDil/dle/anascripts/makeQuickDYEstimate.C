
void printR(TString fileName, TString sourceName, TString nJets, TString hyp_type)
{

    TFile f(fileName, "READ");

    TH1F *h1_met_in = f.Get(sourceName + "_dyest_met_in_" + nJets + "_" + hyp_type);
    TH1F *h1_met_out = f.Get(sourceName + "_dyest_met_out_" + nJets + "_" + hyp_type);

    Double_t R = -1.0;
    Double_t err = -1.0;
    Double_t err_nIn = 0.0;
    Double_t err_nOut = 0.0;
    Double_t nOut = h1_met_out->IntegralAndError(0, h1_met_out->GetNbinsX() + 1, err_nOut);
    Double_t nIn = h1_met_in->IntegralAndError(0, h1_met_in->GetNbinsX() + 1, err_nIn);

    if (nIn != 0 && nOut != 0) {
        R = nOut/nIn;
        err = R*sqrt( ((err_nOut*err_nOut)/(nOut*nOut)) + ((err_nIn*err_nIn)/(nIn*nIn)) );
    }

    printf("out: %.2f, in: %.2f (R = %.2f $\\pm$ %.4f) \t", nOut, nIn, R, err);

    delete h1_met_in;
    delete h1_met_out;

}

void makeTable(TString fileName)
{

    std::cout << "doing ratio for cuts in file " << fileName << std::endl << std::endl;

    //
    // ee
    //

    std::cout << "\\begin{table}[ht]" << std::endl;
    std::cout << "\\caption{In-Out ratio for a veriety of cuts (ee)}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|l*{3}{|c}|r|}\\hline" << std::endl;
    std::cout << "Cuts   & Data & MC  \\\\ \\hline" << std::endl;

    std::cout << fileName << " & ";
    printR(fileName, "data", "allj", "ee");
    std::cout << "& ";
    printR(fileName, "DYee", "allj", "ee");
    std::cout << "\\\\ \\hline" << std::endl;


    std::cout <<"\\end{tabular}" << std::endl;
    std::cout <<"\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;

    std::cout << std::endl << std::endl;

    //
    // mm
    //

    std::cout << "\\begin{table}[ht]" << std::endl;
    std::cout << "\\caption{In-Out ratio for a veriety of cuts (mm)}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|l*{3}{|c}|r|}\\hline" << std::endl;
    std::cout << "Cuts   & Data & MC  \\\\ \\hline" << std::endl;

    std::cout << fileName << " & ";
    printR(fileName, "data", "allj", "mm");
    std::cout << "& ";
    printR(fileName, "DYmm", "allj", "mm");
    std::cout << "\\\\ \\hline" << std::endl;


    std::cout <<"\\end{tabular}" << std::endl;
    std::cout <<"\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;

}

void makeQuickDYEstimate()
{
    makeTable("hist_usePtGt1010_applyFOv2Cuts_applyTriggers_usetcMET_usejptJets_useOS.root");
}

