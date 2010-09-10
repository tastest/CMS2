// ROOT includes
#include "TSystem.h"
#include <iostream>

// TAS includes
#include "fakerates.h"

// muons
bool isFakeableMuon (int index, enum fakeRateVersion version){
	std::cout << std::endl  << "fakerates.cc is depricated. Fakeable object definitions are now in muonSelections.cc" <<  std::endl << std::endl;
	gSystem->Exit(1);
	return false;
} 
double muFakeProb (int i_mu, enum fakeRateVersion version){
	std::cout << std::endl  << "fakerates.cc is depricated. Fakeable object definitions are now in muonSelections.cc" <<  std::endl << std::endl;
	gSystem->Exit(1);
	return 0;
}
double muFakeProbErr (int i_mu, enum fakeRateVersion version){
	std::cout << std::endl  << "fakerates.cc is depricated. Fakeable object definitions are now in muonSelections.cc" <<  std::endl << std::endl;
	gSystem->Exit(1);
	return 0;
}
void SetMuonFile(char * filename){
	std::cout << std::endl  << "fakerates.cc is depricated. Fakeable object definitions are now in muonSelections.cc" <<  std::endl << std::endl;
	gSystem->Exit(1);
	return;
}

// electrons
bool isFakeableElectron (int index, enum fakeRateVersion version){
	std::cout << std::endl  << "fakerates.cc is depricated. Fakeable object definitions are now in electronSelections.cc" <<  std::endl << std::endl;
	gSystem->Exit(1);
	return false;
}
double elFakeProb (int i_el, enum fakeRateVersion version, char* sample ){
	std::cout << std::endl  << "fakerates.cc is depricated. Fakeable object definitions are now in electronSelections.cc" <<  std::endl << std::endl;
	gSystem->Exit(1);
	return 0;
}
double elFakeProbErr (int i_el, enum fakeRateVersion version, char* sample ){
	std::cout << std::endl  << "fakerates.cc is depricated. Fakeable object definitions are now in electronSelections.cc" <<  std::endl << std::endl;
	gSystem->Exit(1);
	return 0;
}
void SetElectronFile(char * filename){
	std::cout << std::endl  << "fakerates.cc is depricated. Fakeable object definitions are now in electronSelections.cc" <<  std::endl << std::endl;
	gSystem->Exit(1);
	return;
}

