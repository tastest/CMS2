// $Id: bx.cc,v 1.1 2010/01/29 19:39:34 jmuelmen Exp $

// CINT is allowed to see this, but nothing else:
bool bx (unsigned int evt_run, int evt_bunchCrossing);

#ifndef __CINT__

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <set>

struct run_and_bx {
     unsigned int run;
     std::set<int> bxs;
     run_and_bx (unsigned int run_) : run(run_) { }
};

bool operator < (const struct run_and_bx &r1, const struct run_and_bx &r2)
{
     return r1.run < r2.run;
}

typedef std::multiset<struct run_and_bx> set_t; 
static set_t bxs_;
static bool bxs_loaded_ = false;

static void load_bxs (const char *fname)
{
     FILE *file = fopen(fname, "r");
     if (file == 0) {
	  perror("opening bx list");
	  return;
     }
     int s;
     int line = 0;
     do {
	  int n;
	  char buf[1024] = "";
	  // read a line from the file, not including the newline (if
	  // there is a newline)
	  s = fscanf(file, "%1024[^\n]%n", buf, &n);
	  assert(n < 1023);
	  if (s != 1) {
	       if (s != EOF) {
		    perror("reading bx list");
		    return;
	       } else {
		    if (ferror(file)) {
			 perror("reading bx list");
			 return;
		    }
	       }
	  } else {
	       line++;
	       // printf("Read a line from the bx list: %s\n", buf);
	       unsigned int run;
	       char *pbuf = buf;
	       int s = sscanf(pbuf, " %u%n", &run, &n);
	       if (s != 1) {
		    fprintf(stderr, "Expected a run number (unsigned int)"
			    " in the first position of line: %s\n", buf);
		    assert(s == 1);
	       }
	       pbuf += n;
	       // keep reading ints from the line, interpret them as a list of bx numbers 
	       struct run_and_bx new_entry(run);
	       do { 
		    int bx;
		    s = sscanf(pbuf, " %d%n", &bx, &n);
		    if (s == 1) {
			 pbuf += n;
			 new_entry.bxs.insert(bx);
		    }
	       } while (s == 1);
	       char trail[1024] = "";
	       s = sscanf(pbuf, " %s", trail);
	       if (strlen(trail) != 0) {
		    fprintf(stderr, "Unexpected trailing junk (%s) on line: %s\n", trail, buf);
		    assert(s == 0);
	       }
	       bxs_.insert(new_entry);
	       // advance past the newline 
	       char newlines[1024] = "";
	       s = fscanf(file, "%[ \f\n\r\t\v]", newlines); 
	       if (strlen(newlines) != 1) {
		    fprintf(stderr, "Warning: unexpected white space following line %d\n", line);
		    // but that's just a warning
	       }
	  }
     } while (s == 1);
     fclose(file);
}

bool bx (unsigned int run, int bx)
{
     if (not bxs_loaded_) {
	  load_bxs("bxs.txt");
	  bxs_loaded_ = true;
     }
     // we assume that an empty list means accept anything
     if (bxs_.size() == 0)
	  return true;
     // find all blocks with this run number
     struct run_and_bx r_a_b(run);
     std::pair<set_t::const_iterator, set_t::const_iterator> blocks = 
	  bxs_.equal_range(r_a_b);
     for (set_t::const_iterator i = blocks.first; 
	  i != blocks.second; ++i) {
	  if (i->bxs.find(bx) != i->bxs.end())
	       return true;
     }
     return false;
}

#endif // __CUNT__
