#!/bin/bash
#need to make this into makefile

gcc `root-config --cflags` -c -shared CORE/CMS2.cc -o libCMS2.so -O3 -I./ -I./CORE
gcc `root-config --cflags` -c -shared loader.C -o libloader.so -O3 -I./ -I./CORE
rm libgetMyHistosNames.so
gcc `root-config --cflags` -c -shared getMyHistosNames.C -o libgetMyHistosNames.so -O3 -I./ -I./CORE
rm libhisttools.so
gcc `root-config --cflags` -c -shared histtools.C -o libhisttools.so -O3 -I./ -I./CORE
rm libbrowseStacks.so
gcc `root-config --cflags` -c -shared browseStacks.C -o libbrowseStacks.so -O3 -I./ -I./CORE
rm libttDilCounts_looper.so
gcc `root-config --cflags`  -c -shared ttDilCounts_looper.C  -o libttDilCounts_looper.so  -O3 -I./ -I./CORE
gcc `root-config --cflags` doAll_oneM.C -o doAll_oneM `root-config --libs`  -L./ -lloader -lbrowseStacks -lgetMyHistosNames -lhisttools  -lttDilCounts_looper  `root-config --libs` -lPhysics -lEG -lGenVector -lCMS2 `root-config --ldflags` -z muldefs
gcc `root-config --cflags` doAllOctX.C -o doAllOctX `root-config --libs`  -L./ -lloader -lbrowseStacks -lgetMyHistosNames -lhisttools  -lttDilCounts_looper  `root-config --libs` -lPhysics -lEG -lGenVector -lCMS2 `root-config --ldflags` -z muldefs
gcc `root-config --cflags` doAllOctX7TeV.C -o doAllOctX7TeV `root-config --libs`  -L./ -lloader -lbrowseStacks -lgetMyHistosNames -lhisttools  -lttDilCounts_looper  `root-config --libs` -lPhysics -lEG -lGenVector -lCMS2 `root-config --ldflags` -z muldefs
