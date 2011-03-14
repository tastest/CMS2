
#ifndef READTRIGGERLIST_H
#define READTRIGGERLIST_H

#include <vector>
#include <set>
#include <string>

void set_trigger_file (const char* filename);
std::vector<std::string> get_trigger_names(unsigned int run, const char *type);

#endif

