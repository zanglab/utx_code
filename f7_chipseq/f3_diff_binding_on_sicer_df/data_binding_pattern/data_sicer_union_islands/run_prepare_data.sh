cp ../../data_sicer_df/sicer_df_results/*/*union.island .

for ii in *.island
do
  echo $ii
  cat $ii |awk '{OFS="\t";print$0,NR,"0","+"}' >${ii}.bed
done
