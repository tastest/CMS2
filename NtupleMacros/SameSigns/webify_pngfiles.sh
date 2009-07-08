#!/bin/bash
if [ ! $1 ] 
then
 echo "usage: webify_rootfile.sh <folder with pngs>"
else

x=$1
y=${x%.root}
#echo "Creating target directory Histos_"${y##*/}

$RAWfilename=${y##*/}
directoryname=$1
#mkdir $directoryname
#cp $1 $directoryname
cd    $directoryname
current_dir=$PWD
#root -b -q -l '/uscms/home/ibloch/afshome/CMS/dumpAllPlots.C("'$1'","'$2'")'

htmlname="index.html"
files=`ls *.png`

echo "If not run the first time first delete *SM.png"
echo "If not run the first time first delete *SM.png"
echo "If not run the first time first delete *SM.png"

echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">" > $htmlname
echo "<html>"  >> $htmlname
echo "  <head>"  >> $htmlname
echo "    <meta http-equiv=\"Content-Script-Type\" content=\"type\">"  >> $htmlname
echo "    <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">"  >> $htmlname
echo "    <meta name=\"Author\" content=\"Ingo Bloch\">"  >> $htmlname
echo "    <meta name=\"description\" content=\"\">"  >> $htmlname
echo "    <meta name=\"keywords\" content=\"quick html page for image display\">"  >> $htmlname
echo ""  >> $htmlname
echo "    <title>Quick picture display</title>"  >> $htmlname
echo "  </head>"  >> $htmlname
echo "  "  >> $htmlname
echo "  <body>"  >> $htmlname

echo "    <center style=\"font-family: helvetica,arial,sans-serif;\">" >> $htmlname
echo "      <table border=\"0\" cellspacing=\"1\" cellpadding=\"2\" width=\"950\" >" >> $htmlname
echo "	<tbody>" >> $htmlname
echo "	  <tr>" >> $htmlname
echo "" >> $htmlname
echo "	    <td valign=\"top\"><b><font size=\"+3\">Histograms from file:</font></b>" >> $htmlname
echo "	      <br><br>" >> $htmlname
echo "	    <font size=\"+0\">"$current_dir"/"$1"</font>" >> $htmlname
echo "	      <br><br>" >> $htmlname
echo "	    <font size=\"+0\">Only pictures with names containing: \"<b>"$2"</b>\" are shown.</font>" >> $htmlname
echo "	      <br><br>" >> $htmlname
echo "	      " >> $htmlname
echo "	      <br>&nbsp;" >> $htmlname
echo "	      <br><br><br><br>" >> $htmlname
echo "" >> $htmlname
echo "" >> $htmlname
echo "	      <!-- room for the page begin: -->" >> $htmlname
## dump all the gif files that have been produced by dumpAllPlots.C into the html
for file in $files
do
  fil=${file%%"_sm.gif"}
#  bigfil=$fil".gif"
  bigfil=$fil""
  echo $bigfil
  smfilbuf=${file%%".png"}
  smfil=$smfilbuf"_sm.PNG"
  echo $smfil
  convert -scale 233 $bigfil $smfil
#echo "<a href=\""$fil"\"><img src=\""$file"\" width=\"400\"></a>"  >> $htmlname
echo "<a href=\""$bigfil"\"><img src=\""$smfil"\" width=\"233\" alt=\"\"></a>"  >> $htmlname

done  
## continue with html body
echo "	      <!-- room for the page end: -->" >> $htmlname
echo "" >> $htmlname
echo "" >> $htmlname
echo "" >> $htmlname
echo "	    </td>" >> $htmlname
echo "	  </tr>" >> $htmlname
echo "" >> $htmlname
echo "" >> $htmlname
echo "	</tbody>" >> $htmlname
echo "      </table>" >> $htmlname
echo "    </center>" >> $htmlname
echo "  </body>"  >> $htmlname
echo "</html>"  >> $htmlname

cd -
fi
