
{

    //
    TChain *chain_whunt_skim = new TChain("Events");
    chain_whunt_skim->Add("/tas03/disk01/whunt/skim/emuskim_*.root");
    chain_whunt_skim->Merge("/tmp/emuskim_merged.root", "fast");

    //
    delete chain_whunt_skim;

}

