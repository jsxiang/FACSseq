df='../Data/Intensities/BaseCalls/MFSI_peared.assembled.fastq'

# run the following first
#while read p; do
#  echo $p
#  grep $p $df >> MFSI_alllib.seq
#done <libraryregex2.txt

#libregex=$(<libraryregex.txt)
#libname=$(<libnames.txt)
fext=".fasta"
#for a in 0 1 2 3 4 5 6 7 8
#do
#    echo $a
#    echo ${libregex[a]}
#    echo ${libnames[$a]}
#done

while IFS= read -r lineA && IFS= read -r lineB <&3; do
  echo "$lineA"; echo "$lineB"
  grep $lineB $df > $lineA.seq
done <t3libnames.txt 3<t3libregex.txt 