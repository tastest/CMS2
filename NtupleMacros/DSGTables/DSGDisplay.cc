#include <assert.h>
#include <ncurses.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "DSGTable.h"

void DSGDisplay ()
{       
     initscr();                       /* Start curses mode            */
     raw();                           /* Line buffering disabled      */
     keypad(stdscr, TRUE);            /* We get F1, F2 etc..          */
     noecho();                        /* Don't echo() while we do getch */
     // get the file with the tables
     TFile *f = TFile::Open("Results.root", "read");
     assert(f != 0);
     DSGTable *dsg_ww 		= dynamic_cast<DSGTable *>(f->Get("ww"));
#if 0
     DSGTable *dsg_wz   	= dynamic_cast<DSGTable *>(f->Get("wz"));
     DSGTable *dsg_zz   	= dynamic_cast<DSGTable *>(f->Get("zz"));
     DSGTable *dsg_wjets	= dynamic_cast<DSGTable *>(f->Get("wjets"));
     DSGTable *dsg_dyee 	= dynamic_cast<DSGTable *>(f->Get("dyee"));
     DSGTable *dsg_dymm 	= dynamic_cast<DSGTable *>(f->Get("dymm"));
     DSGTable *dsg_dytt 	= dynamic_cast<DSGTable *>(f->Get("dytt"));
     DSGTable *dsg_ttbar	= dynamic_cast<DSGTable *>(f->Get("ttbar"));
     DSGTable *dsg_tw   	= dynamic_cast<DSGTable *>(f->Get("tw"));
     DSGTable *dsg_data 	= dynamic_cast<DSGTable *>(f->Get("LM2"));
#endif
     const DSGTable *dsgs[] = { 
	  dsg_ww          ,
// 	  dsg_wz          ,
// 	  dsg_zz          ,
// 	  dsg_wjets       ,
// 	  dsg_dyee        ,
// 	  dsg_dymm        ,
// 	  dsg_dytt        ,
// 	  dsg_ttbar       ,
// 	  dsg_tw          ,
// 	  dsg_data        ,
	  0
     };
     const int n_dsgs = sizeof(dsgs) / sizeof(DSGTable *);

     // print table
     mvprintw(0, 20, "MET > 0");
     mvprintw(0, 70, "MET > 45");
     mvprintw(0, 120, "MET > 175");
     mvprintw(15, 0, "no Z");
     mvprintw(40, 0, "Z");
     for (int i = 0; i < DSGTable::nZcat; ++i) {
	  for (int j = 0; j < DSGTable::nMETcat; ++j) {
	       mvprintw(2 + i * 25, 30 + j * 50, "Njet");
	       mvprintw(3 + i * 25, 22 + j * 50, "<2        2-4       >4");
	       for (int k = 0; k < DSGTable::nJetcat; ++k) {
		    for (int l = 0; l < DSGTable::nBuckets; ++l) {
 			 double sm = 0;
 			 for (int m = 0; m < n_dsgs - 1; ++m) {
 			      sm += dsgs[m]->events_[i][j][k][l];
 			 }
// 			 double data = dsgs[n_dsgs - 1]->events_[i][j][k][l];
			 mvprintw(3 + i * 25 + (l + 1) * 2,
				  10 + j * 50 + (k + 1) * 10,
				  "%6.2f", sm);
			 if (k == 0) {
			      const static char bucketStr[10][1280] = {
				   "e+e+", "m+m+", "e+m+",
				   "e-e-", "m-m-", "e-m-",
				   "e+m-", "m+e-",
				   "e+e-", "m+m-"
			      };
			      mvprintw(3 + i * 25 + (l + 1) * 2,
				       10 + j * 50, bucketStr[l]);
			 }
		    }
	       }
	  }
     }
     refresh();
     int i = 0;
     int j = 0;
     int k = 0;
     int l = 0;
     attron(A_REVERSE);
     double sm = 0;
     for (int m = 0; m < n_dsgs - 1; ++m) {
	  sm += dsgs[m]->events_[i][j][k][l];
     }
     mvprintw(3 + i * 25 + (l + 1) * 2,
	      10 + j * 50 + (k + 1) * 10,
	      "%6.2f", sm);
     while (1) {
	  // read input
	  int ch = getch();
	  if (tolower(ch) == 'q')
	       break;
	  // unhighlight current cell
	  attroff(A_REVERSE);
	  double sm = 0;
	  for (int m = 0; m < n_dsgs - 1; ++m) {
	       sm += dsgs[m]->events_[i][j][k][l];
	  }
	  mvprintw(3 + i * 25 + (l + 1) * 2,
		   10 + j * 50 + (k + 1) * 10,
		   "%6.2f", sm);
	  switch (ch) {
	  case KEY_UP:
	       if (l == 0) {
		    if (i != 0) {
			 i--;
			 l = DSGTable::nBuckets - 1;
		    }
	       } else {
		    l--;
	       }
	       break;
	  case KEY_DOWN:
	       if (l == DSGTable::nBuckets - 1) {
		    if (i != DSGTable::nZcat - 1) {
			 i++;
			 l = 0;
		    }
	       } else {
		    l++;
	       }
	       break;
	  case KEY_LEFT:
	       if (k == 0) {
		    if (j != 0) {
			 j--;
			 k = DSGTable::nJetcat - 1;
		    }
	       } else {
		    k--;
	       }
	       break;
	  case KEY_RIGHT:
	       if (k == DSGTable::nJetcat - 1) {
		    if (j != DSGTable::nMETcat - 1) {
			 j++;
			 k = 0;
		    }
	       } else {
		    k++;
	       }
	       break;
	  case 'h': case 'H':
	  {
	       static TCanvas *c = new TCanvas;
	       c->Clear();
	       c->Divide(2, 1);
	       c->cd(1);
	       dsgs[0]->hmet_[i][j][k][l]->Draw();
	       c->cd(2);
	       dsgs[0]->hmll_[i][j][k][l]->Draw();
	       c->Update();
	       break;
	  }
	  default:
	       break;
	  }
	  // highlight new cell
	  attron(A_REVERSE);
	  sm = 0;
	  for (int m = 0; m < n_dsgs - 1; ++m) {
	       sm += dsgs[m]->events_[i][j][k][l];
	  }
	  mvprintw(3 + i * 25 + (l + 1) * 2,
		   10 + j * 50 + (k + 1) * 10,
		   "%6.2f", sm);
     } 
     endwin(); 
     return;
}
