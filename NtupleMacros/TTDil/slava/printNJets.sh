f=$1
[ ! -f "${f}" ] && echo exit 33
fmt=%6.3f
[ "x$2" != "x" ] && fmt=$2
echo -e ".x setup.C \n .L histscripts/printHists.C \n hist::loadHist(\"${f}\",0,\"*hnJet*\", \"^LM*\");\n printNJets(false, \"${fmt}\");\n gSystem->Exit(0);" | root -l -b

