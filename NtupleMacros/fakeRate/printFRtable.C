#include <iostream>
#include <string>

#include "TROOT.h"
#include "TH2.h"

void printFRtable(TH2F* hist) {

	 // some commonly used charactors for printing table
	 std::string colsep = " & ";
	 std::string pmSign = " $\\pm$ ";
	 std::string endL   = " \\\\ \\hline";

	 // print out header stuff for latex table
	 std::cout << "\\begin{table}[htb]" << std::endl;
	 std::cout << "\\begin{center}" << std::endl;
	 std::cout << "\\caption{}" << std::endl;
	 std::cout << "\\begin{tabular}{c|ccccc}" << std::endl;
	 std::cout << "\\hline" << std::endl;
	 std::cout << "\\backslashbox{$|\\eta|$}{$p_T$}";

	 int nbinsx = hist->GetNbinsX(); // number of bins along x axis
	 int nbinsy = hist->GetNbinsY(); // number of bins along y axis

	 // first, print pt ranges
	 for (int ny = 1; ny < nbinsy+1; ny++) {
		  float lowedge  = hist->GetYaxis()->GetBinLowEdge(ny);
		  float width    = hist->GetYaxis()->GetBinWidth(ny);
		  
		  printf("%s%.3f%s%.3f", colsep.c_str(), lowedge, " -- ", lowedge+width);
	 }

	 std::cout << endL << "\\hline" << std::endl;
	 
	 // loop over bins
	 for (int nx = 1; nx < nbinsx+1; nx++) {
		  for (int ny = 1; ny < nbinsy+1; ny++) {

			   int nbin = hist->GetBin(nx, ny);

			   // if first bin in column, print out eta range
			   if (ny == 1) {
					float lowedge  = hist->GetXaxis()->GetBinLowEdge(nx);
					float width    = hist->GetXaxis()->GetBinWidth(nx);		 

					printf("%.3f%s%.3f", lowedge, " -- ", lowedge+width);
			   }

			   float fr    = hist->GetBinContent(nbin);
			   float frerr = hist->GetBinError(nbin);

			   printf("%s%.4f%s%.4f", colsep.c_str(), fr, pmSign.c_str(), frerr);
		  }

		  std::cout << endL << std::endl;
	 }

	 // print out final header stuff for table
	 std::cout << "\\end{tabular}" << std::endl;
	 std::cout << "\\end{center}" << std::endl;
	 std::cout << "\\end{table}" << std::endl;
}
