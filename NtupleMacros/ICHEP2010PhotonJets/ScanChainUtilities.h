//#include "CMS2.cc"
#include "TH1F.h"
#include "TH2F.h"


//global var for tty magic
static int i_permille_old = 0;


bool isGoodTrk( unsigned int i );
bool passesTrigger(bool isGEN);
bool passesTrackCuts();
std::pair<float, float> getConversionInfo(int idx1, int idx2, float bfiled); 
float getTwrHFSwiss( int seedidx, bool em );
float deltaPhi(float phi1,float phi2);
int CaloTwr_ieta( int detid );
int CaloTwr_iphi( int detid );
void progressbar( int nEventsTotal, int nEventsChain );

void NewHist(TH1F*& h, const char* name, const char* title, int bins, double min, double max, bool isGEN=false);
void NewHist(TH2F*& h, const char* name, const char* title, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax, bool isGEN=false);
TH1F* MakeHist(const char* name, const char* title, int bins, double min, double max, bool isGEN=false);
TH2F* MakeHist(const char* name, const char* title, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax, bool isGEN=false);
