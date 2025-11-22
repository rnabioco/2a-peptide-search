# Prokaryotic 2A-like Peptide Discovery

## Overview

This pipeline discovers prokaryotic ribosomal stalling peptides using principles similar to eukaryotic 2A peptide discovery. While eukaryotic 2A peptides cause ribosomal skipping at H/NPGP motifs, prokaryotic stalling peptides have related but distinct sequences (e.g., RAGP in SecM).

## Strategy

### Two Parallel Approaches

#### 1. Domain-Guided Discovery (Primary)
- Extract all GP-containing sequences from prokaryotic proteomes
- Annotate protein domains (Pfam/InterPro)
- **Focus on inter-domain GP motifs** (hypothesis: stalling peptides occur at domain boundaries)
- Cluster by sequence context
- Build consensus patterns and HMMs

#### 2. Comprehensive GP Discovery (Control)
- Extract ALL GP motifs regardless of domain position
- Should recover known stalling peptides (SecM, TnaC, MifM)
- Provides validation of domain-guided approach

### Key Principles

1. **Conservation with Variation**: Like eukaryotic 2A, prokaryotic stalling peptides are conserved but variable
2. **Inter-Domain Location**: Hypothesis that functional stalling occurs between protein domains
3. **GP Core Motif**: GP dipeptide is common feature (RAGP, NPGP, etc.)
4. **Upstream Context**: N-terminal sequence interacts with ribosome exit tunnel

## Workflow

```
┌─────────────────────────────────────────────────────────────┐
│ Data Collection                                             │
├─────────────────────────────────────────────────────────────┤
│ • Bacterial/Archaeal proteomes (UniProt)                    │
│ • Phage genomes (INPHARED, IMG/VR)                          │
│ • Predict ORFs from phage genomes (Prodigal)                │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ GP Motif Extraction                                         │
├─────────────────────────────────────────────────────────────┤
│ • Find all GP-containing sequences                          │
│ • Extract context (30 upstream, 15 downstream)              │
│ • Record position in protein                                │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ Domain Annotation                                           │
├─────────────────────────────────────────────────────────────┤
│ • Run Pfam/InterProScan on GP-containing proteins           │
│ • Classify GP positions:                                    │
│   - Inter-domain (between domains)                          │
│   - Intra-domain (within domains)                           │
│   - Unstructured (no domains)                               │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ Pattern Discovery                                           │
├─────────────────────────────────────────────────────────────┤
│ • Cluster inter-domain GP motifs by sequence context        │
│ • Analyze conservation within clusters                      │
│ • Identify consensus patterns                               │
│ • Build sequence logos                                      │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ HMM Building                                                │
├─────────────────────────────────────────────────────────────┤
│ • Build HMMs for major clusters                             │
│ • Validate against known stalling peptides:                 │
│   - SecM (FSTPVWISQAQGIRAGP)                                │
│   - TnaC (WDPXXXIGP)                                        │
│   - MifM (XXXGPXXGIAGP)                                     │
│   - CydA (RAGP)                                             │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ Genome-wide Search                                          │
├─────────────────────────────────────────────────────────────┤
│ • Search all prokaryotic proteomes                          │
│ • Search predicted phage ORFs                               │
│ • Compare phage vs bacterial patterns                       │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ Analysis & Validation                                       │
├─────────────────────────────────────────────────────────────┤
│ • Compare domain-guided vs all-GP approaches                │
│ • Analyze phage vs bacterial differences                    │
│ • Generate comprehensive report                             │
└─────────────────────────────────────────────────────────────┘
```

## Known Prokaryotic Stalling Peptides

### SecM (Secretion Monitor, *E. coli*)
- **Sequence**: FSTPVWISQAQGIRAGP
- **Stalling motif**: RAGP
- **Function**: Monitors SecA-dependent protein export
- **Mechanism**: Tryptophan (W) interacts with ribosome exit tunnel

### TnaC (Tryptophanase Leader, *E. coli*)
- **Sequence**: WDP...IGP
- **Function**: Regulates tryptophanase expression
- **Mechanism**: Tryptophan-dependent stalling

### MifM (YidC2 Monitor, *B. subtilis*)
- **Sequence**: ...GPXXGIAGP
- **Function**: Monitors membrane protein insertion
- **Mechanism**: Stalls during YidC2 insertion

### Key Differences from Eukaryotic 2A
| Feature | Eukaryotic 2A | Prokaryotic Stalling |
|---------|---------------|----------------------|
| Core motif | H/NPGP | RAGP, GP variants |
| Mechanism | Ribosome skipping | Ribosome stalling |
| N-terminal | Leucine-rich | Aromatic-rich (W, F) |
| Central motif | GDVE | Variable |
| Function | Polyprotein processing | Gene regulation |

## Configuration

Edit `workflow/config-prokaryotic.yaml`:

```yaml
# Select phage databases
phage_databases_to_use:
  - inphared     # Curated phage genomes (~5k)
  # - imgvr      # Large database (~600k)

# Clustering parameters
prokaryotic:
  clustering_identity: 0.7    # 70% identity
  clustering_coverage: 0.8     # 80% coverage
  min_cluster_size: 5         # Minimum sequences per cluster
```

## Usage

### Run Prokaryotic Discovery

```bash
# Include prokaryotic rules
snakemake --configfile workflow/config-prokaryotic.yaml \
          --use-conda --cores 12 \
          --snakefile workflow/rules/prokaryotic.smk \
          prokaryotic_discovery_report
```

### With SLURM

```bash
# Create submission script for prokaryotic pipeline
sbatch submit-prokaryotic.sh
```

### Specific Targets

```bash
# Extract GP motifs from bacteria
snakemake --configfile workflow/config-prokaryotic.yaml \
          extract_gp_motifs

# Predict ORFs from phage genomes
snakemake --configfile workflow/config-prokaryotic.yaml \
          merge_predicted_orfs

# Cluster motifs
snakemake --configfile workflow/config-prokaryotic.yaml \
          cluster_gp_motifs

# Validate against known peptides
snakemake --configfile workflow/config-prokaryotic.yaml \
          validate_against_known_peptides
```

## Output Structure

```
results/prokaryotic/
├── gp_motifs/
│   ├── all_gp_motifs.tsv.gz              # All GP motifs
│   ├── interdomain_gp_motifs.tsv.gz      # Inter-domain GP
│   ├── intradomain_gp_motifs.tsv.gz      # Intra-domain GP
│   └── gp_motif_stats.tsv                # Statistics
├── clusters/
│   ├── gp_clusters.tsv.gz                # Cluster assignments
│   ├── cluster_representatives.fasta      # Representative sequences
│   ├── cluster_conservation.tsv.gz        # Conservation analysis
│   └── logos/                             # Sequence logos per cluster
├── consensus/
│   ├── consensus_patterns.tsv             # Identified patterns
│   └── cluster_*.sto                      # Alignments per cluster
├── models/
│   └── initial/
│       └── cluster_*.hmm                  # HMMs per cluster
├── searches/
│   └── cluster_*.{hmmsearch,tblout,sto}.gz # Search results
├── validation/
│   ├── known_peptide_hits.tsv             # Hits to known peptides
│   └── validation_summary.txt             # Summary
└── reports/
    └── prokaryotic_discovery.html         # Comprehensive report

results/phage/
├── {db}/
│   ├── orfs/
│   │   └── predicted_orfs.faa.gz          # Predicted ORFs
│   └── gp_motifs/
│       ├── all_gp_motifs.tsv.gz
│       └── interdomain_gp_motifs.tsv.gz
└── analysis/
    ├── gp_distribution.tsv                # Distribution analysis
    └── distribution_plots/                # Visualizations

results/combined/
├── interdomain_gp_motifs_all.tsv.gz       # Merged prok + phage
├── clusters/                              # Combined clustering
└── analysis/
    ├── phage_vs_bacterial.tsv             # Comparative analysis
    └── comparison_plots/
```

## Expected Results

### Pattern Classes

Based on known stalling peptides, expect to find:

1. **RAGP Family** (SecM-like)
   - Strong conservation of RAGP
   - Aromatic upstream (W, F)
   - Hydrophobic context

2. **[RK]AGP Family** (Broader)
   - R or K at -3 position
   - Variable central region

3. **NPGP Family** (Eukaryotic-like)
   - Present in some prokaryotes
   - May represent horizontal transfer

4. **Novel Patterns**
   - Conserved GP variants not yet characterized

### Validation Criteria

1. **Known Peptide Recovery**: Should identify all known stalling peptides
2. **Inter-domain Enrichment**: Inter-domain GP motifs should be enriched for stalling signals
3. **Phage vs Bacterial**: May show distinct patterns reflecting different regulation needs
4. **Conservation**: True stalling peptides should show strong sequence conservation

## Interpretation

### High-Confidence Candidates
- Large clusters (>50 sequences)
- High conservation (>70%)
- Inter-domain localization
- Match to known stalling patterns

### Novel Candidates
- Conserved GP variants not matching known patterns
- Enriched in specific organisms/phages
- Consistent domain junction localization

### Validation Experiments
For interesting candidates:
1. Toeprinting assays for ribosome stalling
2. Reporter constructs for functional testing
3. Mutagenesis of conserved residues
4. Structural modeling of exit tunnel interactions

## References

- Ito, K., et al. (2010). Nascentome analysis uncovers novel translational arrest sequences. *Mol. Cell*.
- Ramu, H., et al. (2011). Nascent peptide in the ribosome exit tunnel affects functional properties. *Mol. Cell*.
- Woolstenhulme, C.J., et al. (2013). Nascent peptides that block protein synthesis in bacteria. *PNAS*.
