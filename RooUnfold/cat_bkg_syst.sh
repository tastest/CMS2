#!/bin/bash
for file in summary_1Dunfolding*formated.txt 
do
        cat $file | awk '{print $4"\t"$6}' > ${file}.temp
done
paste *_*_10_10_10_10_formated.txt.temp *_10_*_10_10_10_formated.txt.temp *_10_10_*_10_10_formated.txt.temp *_10_10_10_*_10_formated.txt.temp *_10_10_10_10_*_formated.txt.temp > output_1D
rm *.temp

for file in summary_2Dunfolding_formated.txt summary_2Dunfolding_1*formated.txt summary_2Dunfolding_0*formated.txt summary_2Dunfolding_2*formated.txt
do
        cat $file | awk '{print $4"\t"$6}' > ${file}.temp
done
paste *_*_10_10_10_10_formated.txt.temp *_10_*_10_10_10_formated.txt.temp *_10_10_*_10_10_formated.txt.temp *_10_10_10_*_10_formated.txt.temp *_10_10_10_10_*_formated.txt.temp > output_2D
rm *.temp

for file in summary_2Dunfolding_ttpt*formated.txt
do
        cat $file | awk '{print $4"\t"$6}' > ${file}.temp
done
paste *_*_10_10_10_10_formated.txt.temp *_10_*_10_10_10_formated.txt.temp *_10_10_*_10_10_formated.txt.temp *_10_10_10_*_10_formated.txt.temp *_10_10_10_10_*_formated.txt.temp > output_2D_ttpt
rm *.temp

for file in summary_2Dunfolding_ttrapidity2*formated.txt
do
        cat $file | awk '{print $4"\t"$6}' > ${file}.temp
done
paste *_*_10_10_10_10_formated.txt.temp *_10_*_10_10_10_formated.txt.temp *_10_10_*_10_10_formated.txt.temp *_10_10_10_*_10_formated.txt.temp *_10_10_10_10_*_formated.txt.temp > output_2D_ttrapidity2
rm *.temp

cat output_1D
echo ""
echo ""
echo ""
cat output_2D
echo ""
echo ""
echo ""
cat output_2D_ttpt
echo ""
echo ""
echo ""
cat output_2D_ttrapidity2
echo ""
echo ""
echo ""
