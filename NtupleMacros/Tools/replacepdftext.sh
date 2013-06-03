#!/bin/bash
#!-----------------------------------
#! Replace "oldtext" with "newtext"
#! inside all pdf files in the directory.
#! The old files are saved as *.pdf.bak
#!
#! Needs the pdftk package
#! http://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/
#!-------------------------------------
#
oldtext=Preliminary
newtext=

for pdffile in *.pdf
do
   echo "Processing $pdffile..."
   cp $pdffile $pdffile.bak
   pdftk $pdffile output $pdffile.tmp uncompress
   rm $pdffile
   sed -i -e "s/$oldtext/$newtext/g" $pdffile.tmp
   pdftk $pdffile.tmp output $pdffile compress
   rm $pdffile.tmp
   rm $pdffile.tmp-e
done