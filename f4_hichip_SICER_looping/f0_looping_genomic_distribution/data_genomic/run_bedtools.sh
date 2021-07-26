

for ii in *.bed
do
    echo $ii
    bedtools sort -i $ii > $ii.sorted
    bedtools merge -c 4 -o distinct -i $ii.sorted > $ii.sorted.merged
done