
# Original IMGVR file throws on hmmsearch:
#
# Parse failed (sequence file IMGVR_all_proteins-high_confidence.faa.gz):
# Line 207760398: illegal character -

pigz -dc IMGVR_all_proteins-high_confidence.faa.gz | sed 's/-/A/g' | pigz -c > IMGVR_all_proteins-high_confidence.fixed.faa.gz
