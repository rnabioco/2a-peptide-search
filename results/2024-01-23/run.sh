zcat lit-annotations.tsv.gz  | grep -v '^#' | awk 'BEGIN {FS="\t"} {print ">"$1" " $2"\n"$3}' | gzip -c > lit-annotations.faa.gz

hmmsearch --tblout lit-annotations.class-1.hmmsearch.tab ../../curated-models/2A-class-1.hmm lit-annotations.faa.gz | gzip -c > lit-annotations.class-1.hmmsearch.gz
hmmsearch --tblout lit-annotations.class-2.hmmsearch.tab ../../curated-models/2A-class-2.hmm lit-annotations.faa.gz | gzip -c > lit-annotations.class-2.hmmsearch.gz
