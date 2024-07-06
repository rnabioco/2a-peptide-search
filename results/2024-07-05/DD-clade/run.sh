grep "#=GR" imgvr.class-2.dd-clade.sto | cut -f2 -d' ' | cut -f1 -d"|" > DD-ids.tab

for id in `cat DD-ids.tab`; do
    zgrep -m1 $id ../IMGVR_all_Sequence_information-high_confidence.tsv.gz
done > imgvr-sequence-info.dd-clade.tab
