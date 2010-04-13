// $Id: goodrun.cc,v 1.4 2010/04/13 21:39:48 warren Exp $

// CINT is allowed to see this, but nothing else:
bool goodrun (unsigned int run, unsigned int lumi_block);

#ifndef __CINT__

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <set>

struct run_and_lumi {
     unsigned int run;
     long long int lumi_min;
     long long int lumi_max;
};

bool operator < (const struct run_and_lumi &r1, const struct run_and_lumi &r2)
{
     return r1.run < r2.run;
}

typedef std::multiset<struct run_and_lumi> set_t; 
static set_t good_runs_;
static bool good_runs_loaded_ = false;

static int load_runs (const char *fname)
{
     good_runs_.clear();
     FILE *file = fopen(fname, "r");
     if (file == 0) {
	  perror("opening good run list");
	  return 0;
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
		    perror("reading good run list");
		    return 0;
	       } else {
		    if (ferror(file)) {
			 perror("reading good run list");
			 return 0;
		    }
	       }
	  } else if (strlen(buf) != 0 && buf[0] == '#') {
	       line++;
	       // printf("Read a comment line (line %d) from the good run list: %s\n", line, buf);
	  } else {
	       line++;
	       // printf("Read a line from the good run list: %s\n", buf);
	       unsigned int run;
	       char *pbuf = buf;
	       int s = sscanf(pbuf, " %u%n", &run, &n);
	       if (s != 1) {
		    fprintf(stderr, "Expected a run number (unsigned int)"
			    " in the first position of line %d: %s\n", line, buf);
		    assert(s == 1);
	       }
	       pbuf += n;
	       long long int lumi_min = -1;
	       long long int lumi_max = -1;
	       s = sscanf(pbuf, " %lld%n", &lumi_min, &n);
	       // if there is no lumi_min specified, that means the
	       // entire run is good
	       if (s == 1) {
		    pbuf += n;
		    s = sscanf(pbuf, " %lld%n", &lumi_max, &n);
		    if (s != 1) {
			 fprintf(stderr, "Expected a max lumi section in a lumi section range"
				 " (int) in the third position of line %d: %s\n", line, buf);
			 assert(s == 1);
		    }
		    pbuf += n;
	       }
	       char trail[1024] = "";
	       s = sscanf(pbuf, " %s", trail);
	       if (strlen(trail) != 0) {
		    fprintf(stderr, "Unexpected trailing junk (%s) on line %d: %s\n", trail, line, buf);
		    assert(s == 0);
	       }
	       // printf("Read line: run %u, min lumi %lld, max lumi %lld\n", run, lumi_min, lumi_max);
	       struct run_and_lumi new_entry = { run, lumi_min, lumi_max };
	       good_runs_.insert(new_entry);
	  }
	  // advance past the newline 
	  char newlines[1024] = "";
	  s = fscanf(file, "%[ \f\n\r\t\v]", newlines); 
	  if (s != -1 && strlen(newlines) != 1) {
		fprintf(stderr, "Warning: unexpected white space following line %d\n", line);
		// but that's just a warning
	  } 
	 } while (s == 1);
	 fclose(file);
	 return line;
}	       

bool goodrun (unsigned int run, unsigned int lumi_block)
{
     if (not good_runs_loaded_) {
	  load_runs("goodruns.txt");
	  good_runs_loaded_ = true;
     }
     // we assume that an empty list means accept anything
     if (good_runs_.size() == 0)
	  return true;
     // find all blocks with this run number
     struct run_and_lumi r_a_l = { run, 0, 0 };
     std::pair<set_t::const_iterator, set_t::const_iterator> good_blocks = 
	  good_runs_.equal_range(r_a_l);
     for (set_t::const_iterator i = good_blocks.first; 
	  i != good_blocks.second; ++i) {
// 	  printf("considering run %u, min %lld, max %lld\n", 
// 		 i->run, i->lumi_min, i->lumi_max);
	  if (i->lumi_min <= lumi_block) {
	       if (i->lumi_max == -1 || i->lumi_max >= lumi_block)
		    return true;
	  }
     }
     return false;
}

void set_goodrun_file (const char* filename)
{
  int ret = load_runs(filename);
  assert(ret != 0);
  good_runs_loaded_ = true;
}

#endif // __CUNT__
