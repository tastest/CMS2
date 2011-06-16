{
  gSystem->CompileMacro("../../../Smurf/Core/SmurfTree.h","k");
  gSystem->CompileMacro("smurfAnalysis.C","k");

  SmurfAnalysis a(0.187*1.1,"/smurf/data/Run2011_Spring11_SmurfV6/tas-TightLooseFullMET-alljets");
  // SmurfAnalysis a(0.187*1.1,"/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/");

  // SmurfAnalysis a(0.191,"/smurf/dmytro/samples/smurf", "files/ww_el_fr_smurfV4.root", "el_fr_v4", "files/ww_mu_fr_smurfV4.root", "mu_fr_m2");
  // a.setGoodRunList("Cert_160404-163869_7TeV_PromptReco_Collisions11_JSON.txt");

  cout << "\n****************** WW 0-jet Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::qqww;
  a.measurement_.type = SmurfAnalysis::WW0Jet;
  a.processSamples();
  a.showYields();
  // a.estimateFakeBackground();
  // a.estimateDYBackground();
  a.estimateTopBackground();
  // a.showYields(1.0);

  return;

  cout << "\n****************** WW 1-jet Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::qqww;
  a.measurement_.type = SmurfAnalysis::WW1Jet;
  a.processSamples();
  a.showYields();
  // a.estimateFakeBackground();
  a.estimateTopBackground();

  cout << "\n****************** WW 2-jet Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::qqww;
  a.measurement_.type = SmurfAnalysis::WW2Jets;
  a.processSamples();
  a.showYields();
  // a.estimateFakeBackground();
  a.estimateTopBackground();
  return;
//   cout << "\n****************** HWW115 Selection *******************\n" << endl;
//   a.measurement_.sig_type = SmurfTree::hww115;
//   a.measurement_.type = SmurfAnalysis::HWWCutBased0Jet;
//   a.processSamples();
//   a.showYields(1.0);
//   cout << "\n****************** HWW120 Selection *******************\n" << endl;
//   a.measurement_.sig_type = SmurfTree::hww120;
//   a.measurement_.type = SmurfAnalysis::HWWCutBased0Jet;
//   a.processSamples();
//   a.showYields(1.0);

  // a.estimateFakeBackground();

  a.estimateTopBackground();
  return;
  a.estimateDYBackground();

  cout << "\n****************** WW 1-jet Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::qqww;
  a.measurement_.type = SmurfAnalysis::WW1Jet;
  a.estimateFakeBackground();

  cout << "\n****************** WW 2-jet Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::qqww;
  a.measurement_.type = SmurfAnalysis::WW2Jets;
  a.estimateFakeBackground();

  cout << "\n****************** HWW120 Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::hww120;
  a.measurement_.type = SmurfAnalysis::HWWCutBased0Jet;
  a.estimateFakeBackground();

  cout << "\n****************** HWW130 Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::hww130;
  a.measurement_.type = SmurfAnalysis::HWWCutBased0Jet;
  a.estimateFakeBackground();

  cout << "\n****************** HWW140 Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::hww140;
  a.measurement_.type = SmurfAnalysis::HWWCutBased0Jet;
  a.estimateFakeBackground();

  cout << "\n****************** HWW150 Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::hww150;
  a.measurement_.type = SmurfAnalysis::HWWCutBased0Jet;
  a.estimateFakeBackground();

  cout << "\n****************** HWW160 Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::hww160;
  a.measurement_.type = SmurfAnalysis::HWWCutBased0Jet;
  a.estimateFakeBackground();

  cout << "\n****************** HWW200 Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::hww200;
  a.measurement_.type = SmurfAnalysis::HWWCutBased0Jet;
  a.estimateFakeBackground();

  cout << "\n****************** HWW250 Selection *******************\n" << endl;
  a.measurement_.sig_type = SmurfTree::hww250;
  a.measurement_.type = SmurfAnalysis::HWWCutBased0Jet;
  a.estimateFakeBackground();
}
