#!/bin/bash

HADOOPDIR=$1

echo "
void mergeFiles() {
  TChain *chain = new TChain(\"tree\");
  chain->SetMaxTreeSize(5000000000LL); //default is 100000000000LL
" > mergeFiles.C

PREFIX=""
for FILE in `ls -lh smurf_19*/data.root | awk '{print $9}'`; do
    echo "  chain->Add(\"${FILE}\");" >> mergeFiles.C
done

echo "
  chain->Merge(\"merged_data_nojson.root\", \"fast\");
}" >> mergeFiles.C

#root -b -q mergeFiles.C

