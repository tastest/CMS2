
Double_t getError2(const TH1F *hist, const Int_t binMin, const Int_t binMax)
{

	//Double_t integral = 0.0;
	Double_t error2 = 0.0;
	for (Int_t i = binMin; i <= binMax; ++i)
	{
		//integral += hist->GetBinContent(i);
		error2 += pow(hist->GetBinError(i), 2);
	}
        //std::cout << hist->Integral(binMin, binMax) << std::endl;
	//std::cout << integral << std::endl;
	return error2;
}

void saveCanvas(const TCanvas *c1, const TString name)
{

	c1->SaveAs(name + ".eps");
        c1->SaveAs(name + ".png");
        c1->SaveAs(name + ".root");
        c1->SaveAs(name + ".C");

}

