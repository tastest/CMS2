#include "TH1D.h"
#include "TH2D.h"
#include <stdio.h>  
#include <iostream>
#include <vector>
using namespace std;

void checkZeroBin2D( TH2D* hist2d ) {
  Double_t bin_buffer = -999;
  bool emptyBins = false;
 
  std::vector<TH1D*> slices;
  slices.reserve( hist2d->GetNbinsX() +1 );
  for(int binx = 1; binx <= hist2d->GetNbinsX(); ++binx ) {
    slices[binx] = hist2d->ProjectionY("",binx,binx);
    
    //    if( slices[binx]->Integral() != 0.0 ) slices[binx]->Scale( 1./slices[binx]->Integral() ); // temp: disabled norm!

    for(int binx1d = 1; binx1d <= slices[binx]->GetNbinsX(); ++binx1d) {
      bin_buffer = slices[binx]->GetBinContent(binx1d);
      if( bin_buffer == 0 ) {
        emptyBins = true;
      }; 
      //      hist2d->SetBinContent( binx, binx1d, bin_buffer );
      //      hist2d->SetBinError(   binx, binx1d, slices[binx]->GetBinError(binx1d)   );
    }

  }
  if(emptyBins)  std::cout<<"ALARM!! 2D histo "<<hist2d->GetName()<<" has bin with 0 entries!!"<<std::endl;

}
