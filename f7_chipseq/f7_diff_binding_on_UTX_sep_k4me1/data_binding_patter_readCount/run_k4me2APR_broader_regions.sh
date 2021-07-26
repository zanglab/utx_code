

outdir=readCount_csv_k4me2APR_broader_regions

regions=10000
regions_name=es10kb

regions=20000
regions_name=es20kb

regions=50000
regions_name=es50kb

regions=100000
regions_name=es100kb

python get_pattern_near_site_readcount_write_out_revised.py \
-i /Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/UTX_peaks.bed \
-t /Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202104_UTX_H3K4me2/process_qc_out/DEL_H3K4me2APR/DEL_H3K4me2APR_treat.bam \
-s hg38 -f bam -w ${regions} -b 200 -m -o ${outdir}/DEL_H3K4me2APR_on_UTX_peaks_${regions_name}_bin200.csv

python get_pattern_near_site_readcount_write_out_revised.py \
-i /Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/UTX_peaks.bed \
-t /Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202104_UTX_H3K4me2/process_qc_out/EIF_H3K4me2APR/EIF_H3K4me2APR_treat.bam \
-s hg38 -f bam -w ${regions} -b 200 -m -o ${outdir}/EIF_H3K4me2APR_on_UTX_peaks_${regions_name}_bin200.csv

python get_pattern_near_site_readcount_write_out_revised.py \
-i /Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/UTX_peaks.bed \
-t /Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202104_UTX_H3K4me2/process_qc_out/Vector_H3K4me2APR/Vector_H3K4me2APR_treat.bam \
-s hg38 -f bam -w ${regions} -b 200 -m -o ${outdir}/Vector_H3K4me2APR_on_UTX_peaks_${regions_name}_bin200.csv

python get_pattern_near_site_readcount_write_out_revised.py \
-i /Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f7_chipseq/f7_diff_binding_on_UTX_sep_k4me1/data_utx_binding/UTX_peaks.bed \
-t /Volumes/zanglab/zw5j/since2019_projects/UTX_HaoJiang/f0_data_process/chip_seq/final_chipseq/202104_UTX_H3K4me2/process_qc_out/WT_H3K4me2APR/WT_H3K4me2APR_treat.bam \
-s hg38 -f bam -w ${regions} -b 200 -m -o ${outdir}/WT_H3K4me2APR_on_UTX_peaks_${regions_name}_bin200.csv



