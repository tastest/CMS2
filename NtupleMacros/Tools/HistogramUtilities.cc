
#include "HistogramUtilities.h"
#include "AllDataSources.h"
#include "histtools.h"

using namespace std;

HistogramUtilities::HistogramUtilities(TString fileName, std::vector<DataSource> potentialSources, Double_t lumiNorm) 
{
    setOrder(potentialSources);

    // open root file
    file_ = new TFile(fileName, "READ");
    //std::cout << "opening file " << fileName << "  " << file_->GetName() << std::endl;
    //will i want to initialize file2_ so that i can check if it's null???

    // leave the luminosity norm as in the root file
    lumiNorm_ = lumiNorm;
    verbose_ = false;
}

HistogramUtilities::HistogramUtilities(TString fileName, TString fileName2, std::vector<DataSource> potentialSources, Double_t lumiNorm) 
{
    setOrder(potentialSources);

    // open root file
    file_  = new TFile(fileName , "READ");
    file2_ = new TFile(fileName2, "READ");

    // leave the luminosity norm as in the root file
    lumiNorm_ = lumiNorm;
    verbose_ = false;
}


void HistogramUtilities::setOrder(std::vector<DataSource> potentialSources)
{
    sources_ = potentialSources;
}


TH1F* HistogramUtilities::getHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin, TString nameprefix) 
{
    TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
    TH1F *h1_data = 0;
    for (unsigned int i = 0; i < sources_.size(); ++i)
    {
        if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
            //std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
            if (!h1_data) h1_data = (TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix)->Clone());
            else h1_data->Add((TH1F*)file_->Get(sources_[i].getName() + histNameSuffix));
        }
    }       

    // set the last bin to be the overflow
    Int_t lastBin = h1_data->GetNbinsX();
    h1_data->SetBinContent(lastBin, h1_data->GetBinContent(lastBin) + h1_data->GetBinContent(lastBin+1));

    h1_data->Scale(lumiNorm_);
    if( rebin != 1 )
        h1_data->Rebin(rebin);
    h1_data->SetName( nameprefix + var + "_" + hyp_type );
    h1_data->SetFillColor( 0 );
    return h1_data;
}

//hyp2 is a required argument. Ex. use case: get sum of all samples ee + mm
TH1F* HistogramUtilities::getHistogramSum(sources_t theSources, TString var, TString nJets, TString hyp1, TString hyp2, Int_t rebin) 
{
    TString histNameSuffix1 = "_" + var + nJets + "_" + hyp1;
    TString histNameSuffix2 = "_" + var + nJets + "_" + hyp2;

    TH1F* h_result = 0;
    for (int i = sources_.size() - 1; i >= 0; --i) {
        if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
            //std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
            if( !h_result )	  
                h_result = (TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix1)->Clone());
            else
                h_result->Add( (TH1F*)file_->Get(sources_[i].getName() + histNameSuffix1) );
            TH1F *h2_temp = (TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix2)->Clone());
            h_result->Add(h2_temp);
        }
    }
    if( lumiNorm_ != 1. )
        h_result->Scale(lumiNorm_);
    if( rebin != 1 )
        h_result->Rebin(rebin);

    return h_result;
}


TH2F* HistogramUtilities::get2dHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin) 
{
    //TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
    TString histNameSuffix = "_" + var + nJets + "_" + hyp_type;

    TH2F *h_data = 0;
    for (unsigned int i = 0; i < sources_.size(); ++i)
    {
        if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
            //std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
            if (!h_data) h_data = (TH2F*)(file_->Get(sources_[i].getName() + histNameSuffix)->Clone());
            else h_data->Add((TH2F*)file_->Get(sources_[i].getName() + histNameSuffix));
        }
    }       
    if( lumiNorm_ != 1. )
        h_data->Scale(lumiNorm_);
    if( rebin != 1 )
        h_data->Rebin(rebin);
    return h_data;
}

//get a stack of single var, hyp for given sources
THStack* HistogramUtilities::getStack(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin, bool plotCumulated, bool cumulateAscending) 
{
    // create a new stack object
    //THStack *st_temp = new THStack("st_temp", ""); //instead, use a sensible name:
    TString name = var + nJets + "_" + hyp_type;
    THStack *st_temp = new THStack(name, name);  

    // get each constituent in turn and add to the stack
    TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
    for (int i = sources_.size() - 1; i >= 0; --i)
    {
        if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
            std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
            //std::cout << "reading file " << file_->GetName() << std::endl;
            //file_->cd(); //we used this at one point to debug why we couldn't find a hist
            //gDirectory->ls();
            TH1F *h1_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix)->Clone()));
            //std::cout << h1_temp->GetBinContent(1) << std::endl;
            Int_t lastBin = h1_temp->GetNbinsX();
            h1_temp->SetBinContent(lastBin, h1_temp->GetBinContent(lastBin) + h1_temp->GetBinContent(lastBin+1));
            if( lumiNorm_ != 1. )
                h1_temp->Scale(lumiNorm_);
            if( rebin != 1 )
                h1_temp->Rebin(rebin);
            if (plotCumulated) {
                h1_temp = (TH1F*)cumulate(*h1_temp, cumulateAscending).Clone();
            }
            if (sources_[i].getColor() != 0) {
                h1_temp->SetFillColor(sources_[i].getColor());
                h1_temp->SetLineColor(sources_[i].getColor());
            }
            st_temp->Add(h1_temp);
        }
    }
    return st_temp;
}

//for combining two hyps or two vars, and scaling one--scale by -1 to subtract
THStack* HistogramUtilities::getSumStack(sources_t theSources, TString var1, TString nJets, TString hyp1, TString hyp2, Int_t rebin, TString var2, double scale) {

    if( var2 == "" ) var2 = var1; //default is same var

    // create a new stack object
    TString name1 = var1 + nJets + "_" + hyp1 + hyp2; //name includes both
    THStack *st_temp = new THStack(name1, name1);

    // get each constituent in turn and add to the stack
    TString histNameSuffix1 = "_" + var1 + nJets + "_" + hyp1;
    TString histNameSuffix2 = "_" + var2 + nJets + "_" + hyp2;

    for (int i = sources_.size() - 1; i >= 0; --i) {
        if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
            //std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
            TH1F *h1_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix1)->Clone()));
            TH1F *h2_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix2)->Clone()));
            if (sources_[i].getColor() != 0) {
                h1_temp->SetFillColor(sources_[i].getColor());
                h1_temp->SetLineColor(sources_[i].getColor());
                h2_temp->SetFillColor(sources_[i].getColor());
                h2_temp->SetLineColor(sources_[i].getColor());
            }
            if( lumiNorm_ != 1. ) {
                h1_temp->Scale(lumiNorm_);
                h2_temp->Scale(lumiNorm_);
            }
            if( rebin != 1 ) {
                h1_temp->Rebin(rebin);
                h2_temp->Rebin(rebin);
            }
            h1_temp->Add(h2_temp, scale); //now h1 is h1+h2
            st_temp->Add(h1_temp); //put the sum in the stack
        }
    }
    return st_temp;
}

//for combining same hists from two files in one stack--note had to add second tfile to class
THStack* HistogramUtilities::get2fileStack(sources_t theSources, TString var, TString nJets, TString hyp1, Int_t rebin) {

    // create a new stack object
    TString name = var + nJets + "_" + hyp1;
    THStack *st_temp = new THStack(name, name);

    // get each constituent in turn and add to the stack
    TString histNameSuffix1 = "_" + var + nJets + "_" + hyp1;

    for (int i = sources_.size() - 1; i >= 0; --i) {
        if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
            //std::cout << "getting " << sources_[i].getName() + histNameSuffix1 << std::endl;
            TH1F *h1_temp = ((TH1F*)(file_ ->Get(sources_[i].getName() + histNameSuffix1)->Clone()));
            TH1F *h2_temp = ((TH1F*)(file2_->Get(sources_[i].getName() + histNameSuffix1)->Clone())); //suffix1
            if (sources_[i].getColor() != 0) {
                h1_temp->SetFillColor(sources_[i].getColor());
                h1_temp->SetLineColor(sources_[i].getColor());
                h2_temp->SetFillColor(sources_[i].getColor());
                h2_temp->SetLineColor(sources_[i].getColor());
            }
            if( lumiNorm_ != 1. ) {
                h1_temp->Scale(lumiNorm_);
                h2_temp->Scale(lumiNorm_);
            }
            if( rebin != 1 ) {
                h1_temp->Rebin(rebin);
                h2_temp->Rebin(rebin);
            }
            h1_temp->Add(h2_temp); //now h1 is h1+h2
            st_temp->Add(h1_temp); //put the sum in the stack
        }
    }
    return st_temp;
}


//for combining three hyps, sum of the first two minus the third
//note: last argument (Int_t posneg) is my hack for fixing the negative histogram stacking bug in root:
// root doesn't stack properly if there is a mixture of positive and negative bin content for teh same bin,
// so i'm making one of two hists for each hist--one with only positive bins, other only negative.
// Obviously, this is after subtraction. Return a stack of either positive (1), negative(-1), or normal(0)
// based on the value of posneg
THStack* HistogramUtilities::getSumDifStack(sources_t theSources, TString var, TString nJets, TString hyp1, TString hyp2, TString hyp3, Int_t rebin, Int_t posneg) {

    // create a new stack object
    TString name = var + nJets + "_" + hyp1 + hyp2 + "-" + hyp3; //name includes all
    THStack *st_temp = new THStack(name, name);

    // get each constituent in turn and add to the stack
    //TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
    TString histNameSuffix1 = "_" + var + nJets + "_" + hyp1;
    TString histNameSuffix2 = "_" + var + nJets + "_" + hyp2;
    TString histNameSuffix3 = "_" + var + nJets + "_" + hyp3;

    for (int i = sources_.size() - 1; i >= 0; --i) {
        if( (theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
            //std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
            TH1F *h1_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix1)->Clone()));
            TH1F *h2_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix2)->Clone()));
            TH1F *h3_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix3)->Clone()));
            if (sources_[i].getColor() != 0) h1_temp->SetFillColor(sources_[i].getColor());
            if (sources_[i].getColor() != 0) h2_temp->SetFillColor(sources_[i].getColor());
            if (sources_[i].getColor() != 0) h3_temp->SetFillColor(sources_[i].getColor());
            if( lumiNorm_ != 1. ) {
                h1_temp->Scale(lumiNorm_);
                h2_temp->Scale(lumiNorm_);
                h3_temp->Scale(lumiNorm_);
            }
            if( rebin != 1 ) {
                h1_temp->Rebin(rebin);
                h2_temp->Rebin(rebin);
                h3_temp->Rebin(rebin);
            }

            if( verbose_ && sources_[i].getName() == "dytt" ) {
                int bin = 40; //80 bins, the 0 bin should be 40, this should be -5 to 0
                //cout << h1_temp->GetName() << "  " << sources_[i].getName() << "  bin " << bin << "  bincenter " << h1_temp->GetXaxis()->GetBinCenter(bin) << "  content " << h1_temp->GetBinContent( bin ) << endl;
                cout << sources_[i].getName() << "  bin " << bin << "  bincenter " << h1_temp->GetXaxis()->GetBinCenter(bin) << endl;
                cout << h1_temp->GetName() << "  content " << h1_temp->GetBinContent( bin ) << endl;
                cout << h2_temp->GetName() << "  content " << h2_temp->GetBinContent( bin ) << endl;
                cout << h3_temp->GetName() << "  content " << h3_temp->GetBinContent( bin ) << endl;
                cout << "  1+2-3 " << h1_temp->GetBinContent( bin ) + h2_temp->GetBinContent( bin ) - h3_temp->GetBinContent( bin ) << endl;
            }

            h1_temp->Add(h2_temp); //now h1 is h1+h2
            h1_temp->Add(h3_temp, -1.0); //this adds -1*h3, ie, subtracts h3

            h1_temp = GetPosNeg( h1_temp, posneg );

            st_temp->Add(h1_temp); //put the sum in the stack
        }
    }
    return st_temp;
}


TLegend* HistogramUtilities::getLegend(sources_t theSources, TString var, TString nJets, TString hyp_type) 
{

    // create a new legend object and make it look nice
    TLegend *lg = new TLegend(0.7, 0.5, 0.9, 0.9);
    lg->SetFillColor(kWhite);
    lg->SetLineColor(kWhite);
    lg->SetShadowColor(kWhite);

    // get each constituent in turn and add to the legend
    TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
    for (unsigned int i = 0; i < sources_.size(); ++i)
    {
        if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
            TH1F *h1_temp = (TH1F*)file_->Get(sources_[i].getName() + histNameSuffix)->Clone();
            if (sources_[i].getColor() != 0) {
                h1_temp->SetFillColor(sources_[i].getColor());
                h1_temp->SetLineColor(sources_[i].getColor());
            }
            lg->AddEntry(h1_temp, sources_[i].getLegendName(), "f");
        }
    }
    return lg;
}

//if posneg is -1, set all positive bins to zero, keep all negative
//if posneg is +1, set all negative bins to zero, keep all positive
TH1F* GetPosNeg( TH1F* h1_temp, Int_t posneg ) {

    if( posneg != 0 && posneg != 1 && posneg != -1 )
        cout << "HistogramUtilities Error: bad posneg value: give either -1, 0, or 1." << endl;
    else if( posneg == 0 )
        return h1_temp; //no change if arg is 0

    TH1F* cpy = new TH1F( *h1_temp ); //copy may not be necessary...
    for( int i=0; i < h1_temp->GetNbinsX(); i++ ) {
        if( posneg == -1 && h1_temp->GetBinContent( i ) > 0 )
            cpy->SetBinContent( i, 0 ); //kill bin
        else if( posneg == 1 && h1_temp->GetBinContent( i ) < 0 )
            cpy->SetBinContent( i, 0 ); //kill bin
    }

    return cpy;
}
