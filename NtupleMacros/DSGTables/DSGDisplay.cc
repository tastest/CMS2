#include <assert.h>
#include <math.h>
#include <ncurses.h>
#include <stdlib.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "DSGTable.h"

/*  cheat sheet ...
i   Z/notZ     major row
j   MET        major column
jj  sumJetEt   alt major row
k   njets      minor column
l   bucket     minor row

*/

enum layout {ZMET, JETMET};

static DSGTable **dsgs_;
static int n_dsgs_;

int nBucketsGroups[4] = {10, 6, 4, 3};
int nBucketsPerGroup[4][10] = {
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
  {2, 2, 2, 2, 1, 1},
  {4, 2, 2, 2},
  {6, 2, 2}
};
int buckets[4][10][10] = {
  { {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9} },
  { {0,1},    {2,3},    {4,5},    {6,7},    {8}, {9} },   
  { {0,1,2,3},          {4,5},    {6,7},    {8,9} },    
  { {0,1,2,3,4,5},    {6,7},    {8,9} },   
};
const static char bucketStrs[4][10][1280] = {
  {"e+e+", "e-e-", "m+m+", "m-m-", "e+m+", "e-m-", "e+m-", "m+e-", "e+e-", "m+m-"},
  {"e+e+ e-e-",    "m+m+ m-m-",    "e+m+ e-m-",    "e+m- m+e-",    "e+e-", "m+m-"},
  {"SS SF",    "SS OF",    "OS OF",    "OS SF"},
  {"SS *F",    "OS OF",    "OS SF"},
};

static void printNumbers (int i, int j, int k, int l, enum layout ZJ, int whichBucketGrouping)
{
  double sm = 0;
  double smw2 = 0;
  double data = 0; 
  double dataw2 = 0;
  for (int n = 0 ; n < nBucketsPerGroup[whichBucketGrouping][l]; n++) {
    for (int m = 0; m < n_dsgs_ - 1; ++m) {
      if (ZJ == ZMET) {
	sm += dsgs_[m]->events_[i][j][0][k][buckets[whichBucketGrouping][l][n]];
	smw2 += dsgs_[m]->w2s_[i][j][0][k][buckets[whichBucketGrouping][l][n]];
      } else {
	sm += dsgs_[m]->events_[2][j][i][k][buckets[whichBucketGrouping][l][n]];
	smw2 += dsgs_[m]->w2s_[2][j][i][k][buckets[whichBucketGrouping][l][n]];
      }
    }
    if (ZJ == ZMET) {
      data += dsgs_[n_dsgs_ - 1]->events_[i][j][0][k][buckets[whichBucketGrouping][l][n]];
      dataw2 += dsgs_[n_dsgs_ - 1]->w2s_[i][j][0][k][buckets[whichBucketGrouping][l][n]];
    } else {
      data += dsgs_[n_dsgs_ - 1]->events_[2][j][i][k][buckets[whichBucketGrouping][l][n]];
      dataw2 += dsgs_[n_dsgs_ - 1]->w2s_[2][j][i][k][buckets[whichBucketGrouping][l][n]];
    }
  }
  double sig = (data - sm) / sqrt(smw2 + dataw2);
  if (sig > 3)
    attron(COLOR_PAIR(1));
  mvprintw(3 + i * 25 + (l + 1) * 2,
	   10 + j * 50 + (k + 1) * 10,
	   "%6.2f ", sm);
  mvprintw(4 + i * 25 + (l + 1) * 2,
	   10 + j * 50 + (k + 1) * 10,
	   "%6.2f", data);
  if (sig > 3)
    attroff(COLOR_PAIR(1));
}

// wow, pointer to member!
// we use this so we can specify which of the histograms we want plotted
static void plotTheDistribution (DSGTable::table_t DSGTable::* h, 
				 int i, int j, int k, int l, enum layout ZJ, int whichBucketGrouping)
{

  THStack *stack = new THStack((dsgs_[0]->*h)[0][0][0][0][0]->GetName(), (dsgs_[0]->*h)[0][0][0][0][0]->GetTitle() );
  
  for (int m = 0; m < n_dsgs_ - 1; ++m) {
    TH1F *hist;
    if (ZJ == ZMET) {
      hist = dynamic_cast<TH1F*> ((dsgs_[m]->*h)[i][j][0][k][buckets[whichBucketGrouping][l][0]]->Clone());
    } else {
      hist = dynamic_cast<TH1F*> ((dsgs_[m]->*h)[2][j][i][k][buckets[whichBucketGrouping][l][0]]->Clone());
    }
    for (int n = 1 ; n < nBucketsPerGroup[whichBucketGrouping][l]; n++) {
      if (ZJ == ZMET) {
	hist->Add((dsgs_[m]->*h)[i][j][0][k][buckets[whichBucketGrouping][l][n]]);
      } else {
	hist->Add((dsgs_[m]->*h)[2][j][i][k][buckets[whichBucketGrouping][l][n]]);
      }
    }
    stack->Add(hist);
  }
  TH1F *data;
  if (ZJ == ZMET) {
    data = dynamic_cast<TH1F*> ((dsgs_[n_dsgs_ - 1]->*h)[i][j][0][k][buckets[whichBucketGrouping][l][0]]->Clone());
  } else {
    data = dynamic_cast<TH1F*> ((dsgs_[n_dsgs_ - 1]->*h)[2][j][i][k][buckets[whichBucketGrouping][l][0]]->Clone());
  }
  for (int n = 1 ; n < nBucketsPerGroup[whichBucketGrouping][l]; n++) {
    if (ZJ == ZMET) {
      data->Add((dsgs_[n_dsgs_ - 1]->*h)[i][j][0][k][buckets[whichBucketGrouping][l][n]]);
    } else {
      data->Add((dsgs_[n_dsgs_ - 1]->*h)[2][j][i][k][buckets[whichBucketGrouping][l][n]]);
    }
  }
  data->SetFillStyle(0);
  data->SetMarkerStyle(30);
  data->Draw("pee");
  stack->Draw("same");
  data->Draw("same pee");
}

static void displayHistos (int i, int j, int k, int l, enum layout ZJ, int whichBucketGrouping)

{
     static TCanvas *c = 0;
     if (c == 0) {
	  c = new TCanvas("c", "c", 1200, 800);
// 	  c->Update();
// 	  c->SetWindowSize(c->GetWindowWidth() * 2, c->GetWindowHeight() * 2);
// 	  c->Update();
     }
     c->Clear();
     c->Divide(4, 2);
     int i_pad = 0;
     c->cd(++i_pad);
     plotTheDistribution (&DSGTable::hmet_, i, j, k, l, ZJ, whichBucketGrouping) ;
     c->cd(++i_pad);
     plotTheDistribution (&DSGTable::hmll_, i, j, k, l, ZJ, whichBucketGrouping) ;
     c->cd(++i_pad);
     plotTheDistribution (&DSGTable::hht_, i, j, k, l, ZJ, whichBucketGrouping) ;
     c->cd(++i_pad);
     plotTheDistribution (&DSGTable::hjsumet_, i, j, k, l, ZJ, whichBucketGrouping) ;
     c->cd(++i_pad);
     plotTheDistribution (&DSGTable::hmaxjetpt_, i, j, k, l, ZJ, whichBucketGrouping) ;
     c->cd(++i_pad);
     plotTheDistribution (&DSGTable::hmaxleppt_, i, j, k, l, ZJ, whichBucketGrouping) ;
     c->cd(++i_pad);
     plotTheDistribution (&DSGTable::hlepdphi_, i, j, k, l, ZJ, whichBucketGrouping) ;
     c->Update();
}


void printTable(enum layout ZJ, int whichBucketGrouping)
{
     
  clear();
  mvprintw(0, 20, "  MET > 0");
  mvprintw(0, 70, "  MET > 20");
  mvprintw(0, 120, "  MET > 30");
  mvprintw(0, 170, "  MET > 100");
  mvprintw(0, 220, "  MET > 175");
  if (ZJ == ZMET) {  
    mvprintw(15, 0, "no Z");
    mvprintw(40, 0, "Z");
    for (int i = 0; i < DSGTable::nZcat; ++i) {
      for (int j = 0; j < DSGTable::nMETcat; ++j) {
	mvprintw(2 + i * 25, 30 + j * 50, "Njet");
	mvprintw(3 + i * 25, 22 + j * 50, "<2        2-4       >4");
	for (int k = 0; k < DSGTable::nJetcat; ++k) { 
	  for (int l = 0; l < nBucketsGroups[whichBucketGrouping]; ++l) { // # of elements in current grp: 10, 5, 4 elements
	    printNumbers(i, j, k, l, ZJ, whichBucketGrouping);
	    if (k == 0) {
	      mvprintw(3 + i * 25 + (l + 1) * 2,
		       10 + j * 50, bucketStrs[whichBucketGrouping][l]);
	    }
	  }
	}
      }
    }
  } else {
    mvprintw(15, 0, "SJPT>0");
    mvprintw(40, 0, "SJPT>300");
    mvprintw(65, 0, "SJPT>500");
    for (int i = 0; i < DSGTable::nSumJetcat; ++i) {
      for (int j = 0; j < DSGTable::nMETcat; ++j) {
	mvprintw(2 + i * 25, 30 + j * 50, "Njet");
	mvprintw(3 + i * 25, 22 + j * 50, "<2        2-4       >4");
	for (int k = 0; k < DSGTable::nJetcat; ++k) { 
	  for (int l = 0; l < nBucketsGroups[whichBucketGrouping]; ++l) { // # of elements in current grp: 10, 5, 4 elements
	    printNumbers(i, j, k, l, ZJ, whichBucketGrouping);
	    if (k == 0) {
	      mvprintw(3 + i * 25 + (l + 1) * 2,
		       10 + j * 50, bucketStrs[whichBucketGrouping][l]);
	    }
	  }
	}
      }
    }
  }
  refresh();
}

static void printSWs (int i_sw_)
{
     mvprintw(1, 200, "SEARCH WINDOWS");
     for (int i_sw = 0; i_sw < dsgs_[0]->search_windows_.size(); ++i_sw) {
	  double sm = 0;
	  double smw2 = 0;
	  double data = 0; 
	  double dataw2 = 0;
	  for (int m = 0; m < n_dsgs_ - 1; ++m) {
	       sm += dsgs_[m]->search_windows_[i_sw]->events_;
	       smw2 += dsgs_[m]->search_windows_[i_sw]->w2s_;
	  } 
	  data += dsgs_[n_dsgs_ - 1]->search_windows_[i_sw]->events_;
	  dataw2 += dsgs_[n_dsgs_ - 1]->search_windows_[i_sw]->w2s_;
	  double sig = (data - sm) / sqrt(dataw2 + smw2);
	  mvprintw(8 + i_sw * 10, 200, "%s", dsgs_[0]->search_windows_[i_sw]->name_.c_str());
	  if (i_sw == i_sw_)
	       attron(A_REVERSE);
	  if (sig > 3)
	       attron(COLOR_PAIR(1));
	  mvprintw(10 + i_sw * 10, 200, "%6.2f", sm);
	  mvprintw(11 + i_sw * 10, 200, "%6.2f", data);
	  if (sig > 3)
	       attroff(COLOR_PAIR(1));
	  if (i_sw == i_sw_)
	       attroff(A_REVERSE);
     }
}

static void displaySW (int i_sw)
{
     static TCanvas *c = 0;
     if (c == 0) {
	  c = new TCanvas("csw", "csw", 300, 200);
     }
     THStack *stack = new THStack(dsgs_[0]->search_windows_[i_sw]->h_->GetName(), 
				  dsgs_[0]->search_windows_[i_sw]->h_->GetTitle());
     for (int m = 0; m < n_dsgs_ - 1; ++m) {
	  TH1F *hist = dynamic_cast<TH1F*>(dsgs_[m]->search_windows_[i_sw]->h_->Clone());
	  stack->Add(hist);
     }
     TH1F *data = dynamic_cast<TH1F*>(dsgs_[n_dsgs_ - 1]->search_windows_[i_sw]->h_->Clone());
     data->SetFillStyle(0);
     data->SetMarkerStyle(30);
     data->Draw("pee");
     stack->Draw("same");
     data->Draw("same pee");
     c->Update();
}

void DSGDisplay ()
{       
     initscr();                       /* Start curses mode            */
     raw();                           /* Line buffering disabled      */
     keypad(stdscr, TRUE);            /* We get F1, F2 etc..          */
     noecho();                        /* Don't echo() while we do getch */
     // get the file with the tables
     TFile *f = TFile::Open("Results.root", "read");
     assert(f != 0);
     DSGTable *dsg_vv 		= dynamic_cast<DSGTable *>(f->Get("vv"));
//      DSGTable *dsg_wz   	= dynamic_cast<DSGTable *>(f->Get("wz"));
//      DSGTable *dsg_zz   	= dynamic_cast<DSGTable *>(f->Get("zz"));
     DSGTable *dsg_we  		= dynamic_cast<DSGTable *>(f->Get("we"));
     DSGTable *dsg_wmu 		= dynamic_cast<DSGTable *>(f->Get("wmu"));
     DSGTable *dsg_wtau		= dynamic_cast<DSGTable *>(f->Get("wtau"));
     DSGTable *dsg_zee 	= dynamic_cast<DSGTable *>(f->Get("zeejets"));
     DSGTable *dsg_zmm 	= dynamic_cast<DSGTable *>(f->Get("zmmjets"));
     DSGTable *dsg_ztt 	= dynamic_cast<DSGTable *>(f->Get("zttjets"));
     DSGTable *dsg_ttbar	= dynamic_cast<DSGTable *>(f->Get("ttbar"));
     DSGTable *dsg_tw	= dynamic_cast<DSGTable *>(f->Get("tw"));
//      DSGTable *dsg_tw   	= dynamic_cast<DSGTable *>(f->Get("tw"));
     DSGTable *dsg_data 	= dynamic_cast<DSGTable *>(f->Get("data"));
     DSGTable *dsgs[] = { 
// 	  dsg_ww          ,
//  	  dsg_wz          ,
//  	  dsg_zz          ,
  	  dsg_vv          ,
 	  dsg_we         ,
 	  dsg_wmu        ,
 	  dsg_wtau       ,
 	  dsg_zee        ,
 	  dsg_zmm        ,
 	  dsg_ztt        ,
 	  dsg_ttbar       ,
 	  dsg_tw          ,
  	  dsg_data        ,
     };
     const int n_dsgs = sizeof(dsgs) / sizeof(DSGTable *);
     dsgs_ = dsgs;
     n_dsgs_ = n_dsgs;
     for (int i = 0; i < n_dsgs - 1; ++i) {
//	  *dsgs[i] *= 254e-6;
     }

     // make a placeholder array to be used for search window tables
     DSGTable *dsgs_sw[n_dsgs];

     // set up subtractions
     DSGTable *dsg_ww_emu 		= dynamic_cast<DSGTable *>(f->Get("ww_emu"));
     DSGTable *dsg_wz_emu   	= dynamic_cast<DSGTable *>(f->Get("wz_emu"));
     DSGTable *dsg_zz_emu   	= dynamic_cast<DSGTable *>(f->Get("zz_emu"));
     DSGTable *dsg_we_emu  		= dynamic_cast<DSGTable *>(f->Get("we_emu"));
     DSGTable *dsg_wmu_emu 		= dynamic_cast<DSGTable *>(f->Get("wmu_emu"));
     DSGTable *dsg_wtau_emu		= dynamic_cast<DSGTable *>(f->Get("wtau_emu"));
     DSGTable *dsg_zee_emu 	= dynamic_cast<DSGTable *>(f->Get("zee_emu"));
     DSGTable *dsg_zmm_emu 	= dynamic_cast<DSGTable *>(f->Get("zmm_emu"));
     DSGTable *dsg_ztt_emu 	= dynamic_cast<DSGTable *>(f->Get("ztt_emu"));
     DSGTable *dsg_ttbar_emu	= dynamic_cast<DSGTable *>(f->Get("ttbar_emu"));
//      DSGTable *dsg_tw_emu   	= dynamic_cast<DSGTable *>(f->Get("tw_emu"));
     DSGTable *dsg_data_emu 	= dynamic_cast<DSGTable *>(f->Get("data_emu"));
     DSGTable *dsgs_emu[] = { 
	  dsg_ww_emu          ,
 	  dsg_wz_emu          ,
 	  dsg_zz_emu          ,
 	  dsg_we_emu         ,
 	  dsg_wmu_emu        ,
 	  dsg_wtau_emu       ,
 	  dsg_zee_emu        ,
 	  dsg_zmm_emu        ,
 	  dsg_ztt_emu        ,
 	  dsg_ttbar_emu       ,
// 	  dsg_tw_emu          ,
  	  dsg_data_emu        ,
     };

     // set up search windows
     static const unsigned int n_sws = dsgs_[0]->search_windows_.size();
//      std::vector<std::string> sw_names;
//      sw_names.reserve(n_sws);
//      for (DSGTable::sw_t::const_iterator i = dsgs[0]->search_windows_.begin(),
// 	       i_end = dsgs[0]->search_windows_.end();
// 	  i != i_end; ++i) {
// 	  sw_names.push_back(i->first);
//      }

     // print table
     start_color();
     init_pair(1, COLOR_RED, COLOR_BLACK);
     printTable(JETMET, 0);

     refresh();
     int i = 0;
     int j = 0;
     int k = 0;
     int l = 0;
     enum layout ZJ = JETMET;
     int iBucketGrouping = 0;
     attron(A_REVERSE);
     printNumbers(i, j, k, l, ZJ, iBucketGrouping);
     while (1) {
	  // read input
	  int ch = getch();
	  if (tolower(ch) == 'q')
	       break;
	  // unhighlight current cell
	  attroff(A_REVERSE);
	  printNumbers(i, j, k, l, ZJ, iBucketGrouping );
	  switch (ch) {
	  case KEY_UP: case 'k': case 'K':
	       if (l == 0) {
		    if (i != 0) {
			 i--;
			 l = nBucketsGroups[iBucketGrouping] - 1;
		    }
	       } else {
		    l--;
	       }
	       break;
	  case KEY_DOWN: case 'j': case 'J':
	       if (l == nBucketsGroups[iBucketGrouping] - 1) {
		    if (ZJ == ZMET && i != DSGTable::nZcat - 1 ||
			ZJ == JETMET && i != DSGTable::nSumJetcat - 1) {
			 i++;
			 l = 0;
		    }
	       } else {
		    l++;
	       }
	       break;
	  case KEY_LEFT: case 'h': case 'H':
	       if (k == 0) {
		    if (j != 0) {
			 j--;
			 k = DSGTable::nJetcat - 1;
		    }
	       } else {
		    k--;
	       }
	       break;
	  case KEY_RIGHT: case 'l': case 'L':
	       if (k == DSGTable::nJetcat - 1) {
		    if (j != DSGTable::nMETcat - 1) {
			 j++;
			 k = 0;
		    }
	       } else {
		    k++;
	       }
	       break;
	  case 'b': case 'B':
	    iBucketGrouping++;
	    iBucketGrouping %= 4;
	    // for now, if we've gone beyond the end of the new
	    // grouping, we get put at the new end; ideally, we would
	    // be put into the bucket group that swallowed the old
	    // bucket group
	    if (l > nBucketsGroups[iBucketGrouping] - 1)
		 l = nBucketsGroups[iBucketGrouping] - 1;
	    printTable(ZJ, iBucketGrouping);
	    break;
	  case 'z': case 'Z':
	    if (ZJ == ZMET) {
	      ZJ = JETMET;
	    } else {
	      ZJ = ZMET;
	    }
	    // if we are in the 3rd major row, switch to 2 major rows.
	    if (ZJ == ZMET) {
	      if (i == DSGTable::nSumJetcat - 1) {
		i =  DSGTable::nZcat - 1;
	      }
	    }
	    printTable(ZJ, iBucketGrouping);
	    break;
	  case '\n':
	    displayHistos(i, j, k, l, ZJ, iBucketGrouping);
	    break;
	  case 'w': case 'W':
	  {
	       // search window menu
	       mvprintw(80, 0, "Search windows are active, press 'x' to exit");
	       unsigned int i = 0;
	       while (1) {
		    printSWs(i);
		    refresh();
		    // read input
		    int ch = getch();
		    if (tolower(ch) == 'x')
			 break;
		    switch (ch) {
		    case KEY_UP: case 'k': case 'K':
			 if (i > 0)
			      i--;
			 break;
		    case KEY_DOWN: case 'j': case 'J':
			 if (i < n_sws - 1)
			      i++;
			 break; 
		    case '\t':
			 if (dsgs_ == dsgs) {
			      for (int j = 0; j < n_dsgs; ++j) {
				   dsgs_sw[j] = dsgs[j]->search_windows_[i]->dsgTable_;
			      }
			      dsgs_ = dsgs_sw;
			 } else dsgs_ = dsgs;
			 goto end_windows;
		    case '\n':
			 displaySW(i);
			 break;
		    }
	       }
	       mvprintw(80, 0, "Table mode is active                        ");
	       break;
	  end_windows: 
	       printTable(JETMET, iBucketGrouping);
	       mvprintw(80, 0, "Table mode is active for search window ");
	       attron(A_REVERSE);
	       printw("%s", dsgs[0]->search_windows_[i]->name_.c_str());
	       attroff(A_REVERSE);
	       printw(", <TAB> to exit");
	       break;
	  }
	  case 's': case 'S':
	  {
	       // subtractions menu
	       mvprintw(80, 0, "Subtractions menu: ");
	       attron(A_REVERSE); printw("e");
	       attroff(A_REVERSE); printw("/mu ");
	       attron(A_REVERSE); printw("f");
	       attroff(A_REVERSE); printw("ake rate ");
	       printw("e");
	       attron(A_REVERSE); printw("x");
	       attroff(A_REVERSE); printw("it ");
	       while (1) {
		    refresh();
		    // read input
		    int ch = getch();
		    if (tolower(ch) == 'x')
			 break;
		    switch (ch) {
		    case 'e': case 'E':
			 if (dsgs_ == dsgs)
			      dsgs_ = dsgs_emu;
			 else dsgs_ = dsgs;
			 goto end;
		    case '\n':
			 break;
		    }
	       }
	  end:
	       printTable(JETMET, iBucketGrouping);
	       mvprintw(80, 0, "Table mode is active                        ");
	       break;

	  }
	  case '\t':
	       // tab gets us out of whatever table we're in and puts
	       // us into the default one
	       dsgs_ = dsgs;
	       printTable(JETMET, iBucketGrouping);
	       mvprintw(80, 0, "Table mode is active                        ");
	       break;
	  default:
	       break;
	  }
	  // highlight new cell
	  attron(A_REVERSE);
	  printNumbers(i, j, k, l, ZJ, iBucketGrouping);
	  attroff(A_REVERSE);
	  refresh();
     } 
     endwin(); 
     exit(0);
     return;
}

