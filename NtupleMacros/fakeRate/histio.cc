histio()
{
}

saveHist(const char* filename, const char* pat="*") 
{
  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;

  TFile outf(filename,"RECREATE") ;
  while(obj=iter->Next()) {    
    if (TString(obj->GetName()).Index(re)>=0) {
      obj->Write() ;
      cout << "." ;
      cout.flush() ;
    }
  }
  cout << endl ;
  outf.Close() ;

  delete iter ;
}


loadHist(const char* filename, const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE) 
{
  TFile inf(filename) ;
  //inf.ReadAll() ;
  TList* list = inf.GetListOfKeys() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;
  cout << "pat = " << pat << endl ;

  gDirectory->cd("Rint:") ;

  TObject* obj ;
  TKey* key ;
  cout << "doAdd = " << (doAdd?"T":"F") << endl ;
  cout << "loadHist: reading." ;
  while(key=(TKey*)iter->Next()) {
   
    Int_t ridx = TString(key->GetName()).Index(re) ;    
    if (ridx==-1) {
      continue ;
    }

    obj = inf->Get(key->GetName()) ;
    TObject* clone ;
    if (pfx) {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(Form("%s_%s",pfx,obj->GetName())) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }
      if (oldObj) {
	clone = oldObj ;
	((TH1*)clone)->Add((TH1*)obj) ;
      } else {
	clone = obj->Clone(Form("%s_%s",pfx,obj->GetName())) ;
      }


    } else {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(key->GetName()) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }

      if (oldObj) {
	clone = oldObj ;
	((TH1*)clone)->Add((TH1*)obj) ;
      } else {
	clone = obj->Clone() ;
      }
    }
    if (!gDirectory->GetList()->FindObject(clone)) {
      gDirectory->Append(clone) ;
    }
    cout << "." ;
    cout.flush() ;
  }
  cout << endl;
  inf.Close() ;
  delete iter ;
}
