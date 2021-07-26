
mkdir merged_islands

for factor in UTX UTXFEB H3K4me1 H3K4me2 H3K4me3 H3K27ac MLL4 H3K27me3
do 
  echo $factor
  cat sicer_out/*_${factor}_*.scoreisland  > merged_islands/${factor}.bed
  bedtools sort -i merged_islands/${factor}.bed > merged_islands/${factor}.sort.bed
  bedtools merge -i merged_islands/${factor}.sort.bed > merged_islands/${factor}.merged.bed
done


## version 2
for factor in H3K4me2APR
do 
  echo $factor
  cat sicer_out/*_${factor}_*.scoreisland  > merged_islands/${factor}.bed
  bedtools sort -i merged_islands/${factor}.bed > merged_islands/${factor}.sort.bed
  bedtools merge -i merged_islands/${factor}.sort.bed > merged_islands/${factor}.merged.bed
done