#include <stdio.h>
#include <string>
#include "TChain.h"

TChain *makeChain (const char *glob, const char *treename = "Events")
{
     std::string cmd = "ls ";
     cmd += glob;
     FILE *f = popen(cmd.c_str(), "r");
     if (!f) {
	  perror("Opening pipe");
	  return 0;
     }
     TChain *c = new TChain(treename);
     int s;
     do {
	  char fname[1024];
	  s = fscanf(f, " %1024s\n", fname);
	  if (s != 1) {
	       if (s != EOF)
		    perror("scanning file list");
	  } else {
	       printf("Adding %s\n", fname);
	       c->Add(fname);
	  }
     } while (s == 1);
     if (pclose(f) == -1) 
	  perror("Closing pipe");
     return c;
}

