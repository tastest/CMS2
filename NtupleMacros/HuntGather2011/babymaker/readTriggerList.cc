
// CINT is allowed to see this, but nothing else:
#include "readTriggerList.h"

#ifndef __CINT__

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>

struct run_and_trig {
    unsigned int run_min;
    unsigned int run_max;
    std::string trig;
    std::string type;
};

typedef std::vector<struct run_and_trig> set_t;
static set_t trigger_list_;
static bool trigger_list_loaded_ = false;

//bool operator < (const struct run_and_trig &r1, const struct run_and_trig &r2)
//{
//    return r1.run_min < r2.run_min;
//}

static int load_triggers (const char *fname)
{
    trigger_list_.clear();

    FILE *file = 0;
    file = fopen(fname, "r");
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
                perror("reading trigger list");
                return 0;
            } else {
                if (ferror(file)) {
                    perror("reading trigger list");
                    return 0;
                }
            }
        } else if (strlen(buf) != 0 && buf[0] == '#') {
            line++;
            // printf("Read a comment line (line %d) from the trigger list: %s\n", line, buf);
        } else {
            line++;
            //printf("Read a line from the trigger list: %s\n", buf);
            unsigned int run_min = -1;
            unsigned int run_max = -1;
            char trig[1024] = "";
            char type[1024] = "";

            // get run min
            char *pbuf = buf;
            int s = sscanf(pbuf, "%u%n", &run_min, &n);
            if (s != 1) {
                fprintf(stderr, "Expected a run number (unsigned int)"
                        " in the first position of line %d: %s\n", line, buf);
                assert(s == 1);
            }
            pbuf += n;

            // get run max
            s = sscanf(pbuf, " %u%n", &run_max, &n);
            if (s != 1) {
                fprintf(stderr, "Expected a run number (unsigned int)"
                        " in the second position of line %d: %s\n", line, buf);
                assert(s == 1);
            }
            pbuf += n;

            // get trigger name
            s = sscanf(pbuf, " %s%n", trig, &n);
            if (s != 1) {
                fprintf(stderr, "Expected a trigger name (char)"
                        " in the third position of line %d: %s\n", line, buf);
                assert(s == 1);
            }
            pbuf += n;

            // get trigger type
            s = sscanf(pbuf, " %s%n", type, &n);
            if (s != 1) {
                fprintf(stderr, "Expected a trigger type (char)"
                        " in the third position of line %d: %s\n", line, buf);
                assert(s == 1);
            }
            pbuf += n;

            // check for trailing characters
            char trail[1024] = "";
            s = sscanf(pbuf, " %s", trail);
            if (strlen(trail) != 0) {
                fprintf(stderr, "Unexpected trailing junk (%s) on line %d: %s\n", trail, line, buf);
                assert(s == 0);
            }

            // if we got here all is ok
            struct run_and_trig new_entry = { run_min, run_max, trig, type };
            //printf("read %u %u %s %s\n", run_min, run_max, trig, type);
            trigger_list_.push_back(new_entry);
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

void set_trigger_file (const char* filename)
{
    int ret = load_triggers(filename);
    assert(ret != 0);
    trigger_list_loaded_ = true;
}

std::vector<std::string> get_trigger_names(unsigned int run, const char *type)
{
    std::vector<std::string> triggers;
     if (trigger_list_loaded_) {
        for (unsigned int i = 0; i < trigger_list_.size(); ++i) 
        {
            //printf("%u\n", trigger_list_[i].run_min);
            if (run >= trigger_list_[i].run_min && run <= trigger_list_[i].run_max) {
                if (trigger_list_[i].type == type) triggers.push_back(trigger_list_[i].trig);
            }
        }
    }

    return triggers;

} 

#endif // __CUNT__

