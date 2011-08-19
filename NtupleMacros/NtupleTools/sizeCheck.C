/*
  
Code that gives you a breakdown of branch sizes. Works for anything that has a tree structure 
similar to an EDM file (so CMS2 ntuples and EDM files). The groupByMaker option is mostly for 
CMS2 ntuples and lets you combine branches made from the same maker. If you groupByMaker, you'll 
get three beautiful pie charts, showing the sizes of the top 6 branches and then the size of 
everything else. 

The sizes are determined in three ways:

TTree::GetTotalSize -> Gets the uncompressed size
TTree::GetZipBytes  -> Get the compressed size (this is what you want for large files)
TTree::GetTotBytes  -> Similar to GetTotalSize but does something funny with the buffers 
                       Not sure if this is really useful, but its there

*/


#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPie.h"
#include "TPieSlice.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <map>
//#include <pair>
#include <algorithm> 
using namespace std;


double sum_zipBytes;
double sum_totBytes;
double sum_totSize;


void printBranches(const TObjArray *objArray, double& totSize, double& totBytes, double& zipBytes) {
  
  if(objArray == NULL)
    return;

  for(int i = 0; i < objArray->GetSize(); i++) {
    TBranch *b = (TBranch*)objArray->At(i);
    if(b == NULL)
      continue;

    totSize += b->GetTotalSize("*"); 
    totBytes+= b->GetTotBytes("*");
    zipBytes+= b->GetZipBytes("*");
    printBranches(((TBranch*)objArray->At(i))->GetListOfBranches(),
		  totSize, totBytes, zipBytes);
  }
}


TString Append(const TString s, const TString toAppend, const int numchars) {
  
  TString temp = s;
  for(unsigned int i = 0; i < numchars; i ++) 
    temp = temp + toAppend;

  return temp;
}

TPie *makeTPie(string name, string title, const map<double, TString> m, float norm) {
  
  
  Int_t colors[7] = {2,1,3,4,5,6,7};
  unsigned int i = 0;
  Float_t sum = 0;
  int msize = m.size();
  TPie *p = new TPie(name.c_str(), title.c_str(), 7);

  for(map<double, TString>::const_reverse_iterator rit = m.rbegin();
      rit != m.rend(); rit++, i++) {
    if(i<6) {
      p->GetSlice(i)->SetValue(rit->first/norm);
      p->GetSlice(i)->SetTitle(rit->second);
      p->SetEntryFillColor(i, colors[i]);
    } else 
      sum += rit->first;
  }

  p->GetSlice(6)->SetValue(sum/norm);
  p->GetSlice(6)->SetTitle("Everything Else");
  p->SetEntryFillColor(6, colors[6]);
  p->SetRadius(0.33);
  p->SetAngle3D(50.);
  p->SetTextSize(0.037);
  //p->SetTextFont(62);
  p->SetX(0.45);

  return p;
}

void sizeCheck(string filename, bool groupByMaker = true) {

  TFile *f = TFile::Open(filename.c_str(), "READ");
  if(f==NULL) 
    cout << "Cannot find " << filename << ". Exiting." << endl;
  TTree *t = (TTree*)f->Get("Events");
  if(t==NULL)
    cout << "Cannot find tree named \"Events\". Exiting." << endl;
  TObjArray *fullarray =  t->GetListOfBranches();

  
  map<TString, double> m_bs_zipBytes; //branch, size
  map<TString, double> m_bs_totBytes; //branch, size
  map<TString, double> m_bs_totSize; //branch, size

  
  //for sorting by size
  map<double, TString> m_sb_zipBytes; //size, branch
  map<double, TString> m_sb_totBytes; //size, branch
  map<double, TString> m_sb_totSize; //size, branch


  int maxStringSize = 0;



   sum_zipBytes = 0.; 
   sum_totBytes = 0.;
   sum_totSize  = 0.;
  
  
  for(Int_t i = 0; i < fullarray->GetSize(); ++i) {
    TString bname(fullarray->At(i)->GetName());
    TBranch *b =  (TBranch*)t->GetBranch(bname);
    TObjArray *objArray = b->GetListOfBranches();
    double size_zipBytes = 0.;
    double size_totBytes = 0.;
    double size_totSize =  0.;
    /*
      double size_zipBytes = b->GetZipBytes("*");
      double size_totBytes = b->GetTotBytes("*");
      double size_totSize  = b->GetTotalSize("*");
    */
    //if(TString(b->GetName()).Contains("trkstrkp4"))
      printBranches(objArray, size_totSize,size_totBytes, size_zipBytes);

    //continue;
      
      
    
    sum_zipBytes += size_zipBytes;
    sum_totBytes += size_totBytes;
    sum_totSize  += size_totSize;

    
    if(!groupByMaker) {

      map<double, TString>::iterator it_zipBytes = m_sb_zipBytes.find(size_zipBytes);
      map<double, TString>::iterator it_totBytes = m_sb_totBytes.find(size_totBytes);
      map<double, TString>::iterator it_totSize  = m_sb_totSize.find(size_totSize);

      if(bname.Sizeof() > maxStringSize)
	maxStringSize = bname.Sizeof();

      if(it_zipBytes == m_sb_zipBytes.end())
	m_sb_zipBytes.insert(make_pair(size_zipBytes, bname));
      else 
	m_sb_zipBytes.insert(make_pair(size_zipBytes + 0.0000001, bname));


      if(it_totBytes == m_sb_totBytes.end())
	m_sb_totBytes.insert(make_pair(size_totBytes, bname));
      else 
	m_sb_totBytes.insert(make_pair(size_totBytes + 0.0000001, bname));


      if(it_totSize == m_sb_totSize.end())
	m_sb_totSize.insert(make_pair(size_totSize, bname));
      else 
	m_sb_totSize.insert(make_pair(size_totSize + 0.0000001, bname));

    } else {//group by Maker

      TObject *obj = bname.Tokenize("_")->At(1);
      TString makerName;
      if(obj == NULL )
	makerName = "other";
      else 
	makerName = obj->GetName();
      
      if(makerName.Sizeof() > maxStringSize)
	maxStringSize = makerName.Sizeof();

      map<TString, double>::iterator it_zipBytes = m_bs_zipBytes.find(makerName);
      map<TString, double>::iterator it_totBytes = m_bs_totBytes.find(makerName);	
      map<TString, double>::iterator it_totSize  = m_bs_totSize.find(makerName);
      
      if(it_zipBytes == m_bs_zipBytes.end())
	m_bs_zipBytes.insert(make_pair(makerName, size_zipBytes));
      else {
	double temp = it_zipBytes->second;
	m_bs_zipBytes[makerName] = size_zipBytes + temp;
      }

      if(it_totBytes == m_bs_totBytes.end())
	m_bs_totBytes.insert(make_pair(makerName, size_totBytes));
      else {
	double temp = it_totBytes->second;
	m_bs_totBytes[makerName] = size_totBytes + temp;
      }

      if(it_totSize == m_bs_totSize.end())
	m_bs_totSize.insert(make_pair(makerName, size_totSize));
      else {
	double temp = it_totSize->second;
	m_bs_totSize[makerName] = size_totSize + temp;
      }      
    }//group by maker

  }//loop over branch names


  //open the log files
  ofstream zipBytes, totBytes, totSize;
  if(groupByMaker) {
    zipBytes.open("zipBytes_groupByMaker.txt");
    totBytes.open("totBytes_groupByMaker.txt");
    totSize.open("totSize_groupByMaker.txt");
  } else { 
    zipBytes.open("zipBytes.txt");
    totBytes.open("totBytes.txt");
    totSize.open("totSize.txt");
  }



  if(groupByMaker) {


    //now create a map sorted by the size from m_bs
    map<double, TString> m_bs_zipBytes_sizeSort;
    map<double, TString> m_bs_totBytes_sizeSort;
    map<double, TString> m_bs_totSize_sizeSort;


    

    //sort
    for(map<TString, double>::iterator it = m_bs_zipBytes.begin(); 
	it!= m_bs_zipBytes.end(); it++) {   
    
      if(m_bs_zipBytes_sizeSort.find(it->second) != m_bs_zipBytes_sizeSort.end())
	m_bs_zipBytes_sizeSort.insert(make_pair(it->second + 0.00001, it->first));
      else
	m_bs_zipBytes_sizeSort.insert(make_pair(it->second, it->first));
    }
  

    for(map<TString, double>::iterator it = m_bs_totBytes.begin(); 
	it!= m_bs_totBytes.end(); it++) {   
    
      if(m_bs_totBytes_sizeSort.find(it->second) != m_bs_totBytes_sizeSort.end())
	m_bs_totBytes_sizeSort.insert(make_pair(it->second + 0.00001, it->first));
      else
	m_bs_totBytes_sizeSort.insert(make_pair(it->second, it->first));
    }
  
  
    for(map<TString, double>::iterator it = m_bs_totSize.begin(); 
	it!= m_bs_totSize.end(); it++) {   
    
      if(m_bs_totSize_sizeSort.find(it->second) != m_bs_totSize_sizeSort.end())
	m_bs_totSize_sizeSort.insert(make_pair(it->second + 0.00001, it->first));
      else
	m_bs_totSize_sizeSort.insert(make_pair(it->second, it->first));
    }
  
    
    //make some TPies ->show only the top 5, and lump in the rest with the rest
    TPie *p_zipBytes = makeTPie("p_zipBytes", "Breakdown according to zipBytes",
				m_bs_zipBytes_sizeSort, sum_zipBytes);
    TPie *p_totBytes = makeTPie("p_totBytes", "Breakdown according to totBytes",
				m_bs_totBytes_sizeSort, sum_totBytes);
    TPie *p_totSize = makeTPie("p_totSize",  "Breakdown according to totSize",
			       m_bs_totSize_sizeSort, sum_totSize);
    TCanvas *c1 = new TCanvas();
    p_zipBytes->Draw("3d>");
    c1->SaveAs("zipBytes.png");
    TCanvas *c2 = new TCanvas();
    p_totBytes->Draw("3d>");
    c2->SaveAs("totBytes.png");
    TCanvas *c3 = new TCanvas();
    p_totSize->Draw("3d>");
    c3->SaveAs("totSize.png");


    //cout to logs
    for(map<double, TString>::reverse_iterator rit = m_bs_zipBytes_sizeSort.rbegin(); 
	rit!= m_bs_zipBytes_sizeSort.rend(); rit++) {   
      
      TString tempString = rit->second;
      zipBytes << Append(tempString, " ", maxStringSize - tempString.Sizeof() + 2)  
	       << (rit->first)/sum_zipBytes  << "    " << rit->first << endl;
    }

    for(map<double, TString>::reverse_iterator rit = m_bs_totBytes_sizeSort.rbegin(); 
	rit!= m_bs_totBytes_sizeSort.rend(); rit++) {   
      
      TString tempString = rit->second;
      totBytes << Append(tempString, " ", maxStringSize - tempString.Sizeof() + 2)  
	       << rit->first/sum_totBytes  << "    " << rit->first << endl;
    }
    
    for(map<double, TString>::reverse_iterator rit = m_bs_totSize_sizeSort.rbegin(); 
	rit!= m_bs_totSize_sizeSort.rend(); rit++) {   
      
      TString tempString = rit->second;
      
      totSize <<  Append(tempString, " ", maxStringSize - tempString.Sizeof() + 2)   
	      << rit->first/sum_totSize  << "    " << rit->first << endl;
    }


    //summarize
    zipBytes << "Sum of Branches:    " << sum_zipBytes << " Bytes" << endl;
    totBytes << "Sum of Branches:    " << sum_totBytes << " Bytes" << endl;
    totSize << "Sum of Branches:    " << sum_totSize << " Bytes" << endl;



  } else {
  
    for(map<double, TString>::reverse_iterator rit = m_sb_zipBytes.rbegin(); 
	rit != m_sb_zipBytes.rend(); rit++) {

      TString tempString = rit->second;
      zipBytes << Append(tempString, " ", maxStringSize - tempString.Sizeof() + 2)  
	       << rit->first/sum_zipBytes << "    " << rit->first << endl;
    }


    for(map<double, TString>::reverse_iterator rit = m_sb_totBytes.rbegin(); 
	rit != m_sb_totBytes.rend(); rit++) {

      TString tempString = rit->second;
      totBytes << Append(tempString, " ", maxStringSize - tempString.Sizeof() + 2)  
	       << rit->first/sum_totBytes  << "    " << rit->first << endl;
    }

    for(map<double, TString>::reverse_iterator rit = m_sb_totSize.rbegin(); 
	rit != m_sb_totSize.rend(); rit++)  {

      TString tempString = rit->second;      
      totSize << Append(tempString, " ", maxStringSize - tempString.Sizeof() + 2) 
	      << rit->first/sum_totSize  << "    " << rit->first << endl;

    }
  
    

    zipBytes << "Sum of Branches:    " << sum_zipBytes << " Bytes" << endl;
    totBytes << "Sum of Branches:    " << sum_totBytes << " Bytes" << endl;
    totSize << "Sum of Branches:    " << sum_totSize << " Bytes" << endl;


  }


  zipBytes.close();
  totBytes.close();
  totSize.close();

  f->Close();
  
  

}

