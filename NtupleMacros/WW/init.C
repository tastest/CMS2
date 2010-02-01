{
TH1::AddDirectory(true);
gSystem->AddIncludePath("-I../");
gSystem->Load("libPhysics.so");  
gSystem->Load("libEG.so");
}
