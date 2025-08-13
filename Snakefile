# Snakemake pipeline for 2A peptide search
# Replaces ad hoc shell scripts with a structured workflow

import os
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Global variables
PROJECT_DIR = Path(workflow.basedir)
MODELS = ["class-1", "class-2"]
DATABASES = ["uniprot", "mgnify", "ref-proteomes", "imgvr"]

# Helper function to get database files
def get_db_files():
    db_files = {}
    for db in DATABASES:
        if db == "uniprot":
            db_files[db] = config["databases"]["uniprot"]
        elif db == "mgnify":
            db_files[db] = expand(config["databases"]["mgnify_pattern"], chunk=range(1, 31))
        elif db == "ref-proteomes":
            db_files[db] = config["databases"]["ref_proteomes"]
        elif db == "imgvr":
            db_files[db] = config["databases"]["imgvr"]
    return db_files

DB_FILES = get_db_files()

# Rule all - defines final outputs
rule all:
    input:
        # HMM models
        expand("curated-models/2A-{model}.hmm", model=MODELS),

        # Database search results
        expand("db-searches/{db}.{model}.sto.gz", db=DATABASES, model=MODELS),

        # Combined alignments
        expand("db-searches/combined-{model}.sto.gz", model=MODELS),

        # InterProScan download
        "tools/interproscan-5.62-94.0-64-bit.tar.gz"

# Rule to download InterProScan
rule download_interproscan:
    output:
        "tools/interproscan-5.62-94.0-64-bit.tar.gz"
    params:
        url="http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.62-94.0/interproscan-5.62-94.0-64-bit.tar.gz"
    shell:
        """
        mkdir -p tools
        wget -O {output} {params.url}
        """

# Rule to build HMM models from Stockholm alignments
rule build_hmm_models:
    input:
        "curated-models/2A-{model}.sto.gz"
    output:
        "curated-models/2A-{model}.hmm"
    params:
        name="2A-{model}"
    threads: 1
    shell:
        "hmmbuild -n {params.name} {output} {input}"

# Rule for hmmsearch against single database files (uniprot, ref-proteomes, imgvr)
rule hmmsearch_single:
    input:
        hmm="curated-models/2A-{model}.hmm",
        db=lambda wildcards: DB_FILES[wildcards.db]
    output:
        aln="db-searches/{db}.{model}.sto.gz",
        out="results/{db}.{model}.out"
    params:
        cpu=config["hmmsearch"]["cpu"],
        evalue=config["hmmsearch"]["evalue"]
    threads: config["hmmsearch"]["cpu"]
    wildcard_constraints:
        db="uniprot|ref-proteomes|imgvr"
    shell:
        """
        hmmsearch --cpu {params.cpu} -E {params.evalue} \
            -o {output.out} \
            -A {output.aln}.tmp \
            {input.hmm} {input.db}

        gzip {output.aln}.tmp
        mv {output.aln}.tmp.gz {output.aln}
        """

# Rule for hmmsearch against MGnify (multiple files)
rule hmmsearch_mgnify:
    input:
        hmm="curated-models/2A-{model}.hmm",
        db=DB_FILES["mgnify"]
    output:
        aln="db-searches/mgnify.{model}.sto.gz",
        out="results/mgnify.{model}.out"
    params:
        cpu=config["hmmsearch"]["cpu"],
        evalue=config["hmmsearch"]["evalue"]
    threads: config["hmmsearch"]["cpu"]
    shell:
        """
        # Create temporary directory for individual searches
        mkdir -p tmp/mgnify_{wildcards.model}

        # Search each MGnify file
        for db_file in {input.db}; do
            base=$(basename ${{db_file}} .fa.gz)
            hmmsearch --cpu {params.cpu} -E {params.evalue} \
                -A tmp/mgnify_{wildcards.model}/${{base}}.sto \
                {input.hmm} ${{db_file}}
        done

        # Merge all alignments
        find tmp/mgnify_{wildcards.model} -name "*.sto" -size +0c > tmp/mgnify_{wildcards.model}_files.txt

        if [ -s tmp/mgnify_{wildcards.model}_files.txt ]; then
            esl-alimerge --list tmp/mgnify_{wildcards.model}_files.txt > {output.aln}.tmp
            gzip {output.aln}.tmp
            mv {output.aln}.tmp.gz {output.aln}

            # Create combined output file
            cat tmp/mgnify_{wildcards.model}/*.out > {output.out}
        else
            # Create empty files if no hits found
            touch {output.aln}.tmp
            gzip {output.aln}.tmp
            mv {output.aln}.tmp.gz {output.aln}
            touch {output.out}
        fi

        # Clean up
        rm -rf tmp/mgnify_{wildcards.model}
        rm -f tmp/mgnify_{wildcards.model}_files.txt
        """

# Rule to merge alignments from all databases for each model
rule merge_alignments:
    input:
        lambda wildcards: expand("db-searches/{db}.{model}.sto.gz",
                                db=DATABASES, model=wildcards.model)
    output:
        "db-searches/combined-{model}.sto.gz"
    shell:
        """
        # Create list of input files
        ls {input} > tmp/{wildcards.model}_files.txt

        # Merge alignments
        esl-alimerge --list tmp/{wildcards.model}_files.txt > {output}.tmp

        # Compress
        gzip {output}.tmp
        mv {output}.tmp.gz {output}

        # Clean up
        rm -f tmp/{wildcards.model}_files.txt
        """

# Rule to expand hits - extract sequences with flanking regions
rule expand_hits:
    input:
        aln="db-searches/{db}.{model}.sto.gz",
        db_file=lambda wildcards: DB_FILES[wildcards.db] if wildcards.db != "mgnify" else DB_FILES["mgnify"][0]
    output:
        "expanded-hits/{db}.{model}.expanded.fa"
    shell:
        "python3 src/expand-hits.py {input.aln} {input.db_file} > {output}"

# Rule for maximum sensitivity search (equivalent to hmmsearch.max.sh)
rule hmmsearch_max_sensitivity:
    input:
        hmm="curated-models/2A-{model}.hmm",
        db=lambda wildcards: DB_FILES[wildcards.db]
    output:
        aln="results/max/{db}.{model}.sto.gz",
        out="results/max/{db}.{model}.out"
    params:
        cpu=config["hmmsearch"]["cpu"]
    threads: config["hmmsearch"]["cpu"]
    wildcard_constraints:
        db="uniprot|ref-proteomes|imgvr"
    shell:
        """
        mkdir -p results/max

        hmmsearch --cpu {params.cpu} --max \
            -o {output.out} \
            -A {output.aln}.tmp \
            {input.hmm} {input.db}

        gzip {output.aln}.tmp
        mv {output.aln}.tmp.gz {output.aln}
        """

# Rule to create sequence logos (requires external tools)
rule create_logos:
    input:
        "db-searches/combined-{model}.sto.gz"
    output:
        "img/{model}.logo.png"
    shell:
        """
        # This would require integration with logo generation tools
        # Placeholder for now - you'll need to implement based on your preferred method
        echo "Logo generation for {wildcards.model} - implement with your preferred tool"
        touch {output}
        """

# Utility rule to clean intermediate files
rule clean:
    shell:
        """
        rm -rf tmp/
        rm -rf results/
        rm -f db-searches/*.out
        """

# Utility rule to clean everything
rule clean_all:
    shell:
        """
        rm -rf tmp/
        rm -rf results/
        rm -rf tools/
        rm -rf expanded-hits/
        rm -f db-searches/*.sto.gz
        rm -f db-searches/*.out
        rm -f curated-models/*.hmm
        rm -f img/*.png
        """
