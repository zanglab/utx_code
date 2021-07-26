
python find_overlap_keep_info_NOT_sep_strand.py -a ../data/WT_peaks.bed  -b active_mark.bed -s hg38 -p WT_on_active_mark.bed 
python find_overlap_keep_info_NOT_sep_strand.py -a ../data/WT_peaks.bed  -b active_mark.bed -s hg38 -p WT_on_active_mark_extend1kb.bed -e1 1000 
python find_overlap_keep_info_NOT_sep_strand.py -a ../data/WT_peaks.bed  -b active_mark.bed -s hg38 -p WT_on_active_mark_extend5kb.bed -e1 5000 
python find_overlap_keep_info_NOT_sep_strand.py -a ../data/WT_peaks.bed  -b active_mark.bed -s hg38 -p WT_on_active_mark_extend10kb.bed -e1 10000 


python find_overlap_keep_info_NOT_sep_strand.py -a ../data/DEL_NOT_overlap_with_WT.bed  -b active_mark.bed -s hg38 -p DEL_on_active_mark.bed 
python find_overlap_keep_info_NOT_sep_strand.py -a ../data/DEL_NOT_overlap_with_WT.bed  -b active_mark.bed -s hg38 -p DEL_on_active_mark_extend1kb.bed -e1 1000 
python find_overlap_keep_info_NOT_sep_strand.py -a ../data/DEL_NOT_overlap_with_WT.bed  -b active_mark.bed -s hg38 -p DEL_on_active_mark_extend5kb.bed -e1 5000 
python find_overlap_keep_info_NOT_sep_strand.py -a ../data/DEL_NOT_overlap_with_WT.bed  -b active_mark.bed -s hg38 -p DEL_on_active_mark_extend10kb.bed -e1 10000 

    
