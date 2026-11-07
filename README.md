# Bacterial WGS Workflow — Illumina + Oxford Nanopore

A hands-on **bacterial whole-genome sequencing (WGS)** training workflow combining
**Illumina short reads** and **Oxford Nanopore long reads**, from raw data download
through de novo assembly, reference-based variant calling, structural-variant
calling, and functional annotation.

This repository was created as a **completed training assignment** (see
*Attribution* below). It uses only open-source bioinformatics tools and publicly
available NCBI SRA data.

> **This is a reference / lab-notebook of the commands run during the training**,
> not a single turn-key script. Sample names, paths and accessions should be
> adapted to your own working directory. Remote server addresses, logins and
> personal file paths have been **redacted** (`<REMOTE_HOST>`, `<user>`).

---

## Workflow overview

```
                    public NCBI SRA (kingfisher / prefetch)
                                   │
            ┌──────────────────────┴──────────────────────┐
            ▼                                              ▼
   Illumina short reads                          Nanopore long reads
            │                                              │
   QC: FastQC / fastp                             QC: NanoPlot
   taxonomy: BBMap sendsketch                     taxonomy: sendsketch
            │                                              │
   assembly: SPAdes                               assembly: Flye
   graph: Bandage                                 rotate: dnaapler
   stats: assemblathon / QUAST                    stats: flye.log
            │                                              │
            └──────────────────────┬──────────────────────┘
                                   ▼
        Reference-based variant calling (Bowtie2 → SAMtools)
                                   │
          ┌────────────────┬───────┴────────┬───────────────┐
          ▼                ▼                ▼               ▼
      GATK             BCFtools         FreeBayes      SV: Manta / DELLY / CNVkit
   HaplotypeCaller    mpileup/call
          │
          ▼
   Filtering (QUAL/DP)  →  Annotation (snpEff / VEP)  →  Typing (Prokka/Bakta/AMRFinder/abricate)
```

## Datasets (public NCBI SRA)

| Strain label    | Platform  | SRA accession |
|-----------------|-----------|---------------|
| Agro_illumina   | Illumina  | SRR24759617   |
| Agro_nanopore   | Nanopore  | SRR24759618   |
| Shig1_illumina  | Illumina  | SRR10506616   |
| Shig1_nanopore  | Nanopore  | SRR32486128   |
| Shig2_illumina  | Illumina  | SRR28839724   |
| Shig2_nanopore  | Nanopore  | SRR32486129   |
| E. coli (ref-based variant calling) | Illumina | SRR5384356 |

Reference genome: *E. coli* K-12 MG1655 (RefSeq **GCF_000005845.2 / ASM584v2**).

All data are downloaded from the public NCBI SRA / ENA, so the workflow is
reproducible.

## Main tools

| Stage | Tools |
|-------|-------|
| Download | kingfisher, sra-tools (`prefetch`, `fastq-dump`) |
| QC / trimming | FastQC, MultiQC, fastp, Trimmomatic, BBMap (`reformat.sh`, `sendsketch.sh`) |
| Short-read assembly | SPAdes, Bandage, assemblathon_stats / QUAST |
| Long-read assembly | Flye, NanoPlot, dnaapler |
| Alignment | Bowtie2, SAMtools, Picard |
| Variant calling | GATK4 (HaplotypeCaller), BCFtools, FreeBayes |
| Structural variants | Manta, DELLY, CNVkit |
| Annotation / typing | snpEff, Ensembl-VEP, Prokka, Bakta, Kraken2, NCBI-AMRFinderPlus, abricate, antiSMASH |

See [`wgs_pipeline.sh`](wgs_pipeline.sh) for the full, section-by-section commands.

## Environment setup (conda / Bioconda)

```bash
conda config --add channels bioconda
conda config --add channels conda-forge

# assembly / QC environment
conda create -n wgs_env -y
conda install -n wgs_env -y kingfisher bbmap fastqc fastp spades bandage \
    flye nanoplot dnaapler bakta prokka kraken2 abricate ncbi-amrfinderplus

# a separate environment for GATK is recommended
conda create -n gatk_env -c bioconda -c conda-forge gatk4 openjdk=8 -y
```
Manta (v1.6.0) uses Python 2.7 and is installed in its own environment, as noted
in the script.

## Repository contents

| File | Description |
|------|-------------|
| `wgs_pipeline.sh` | Full documented workflow (download → assembly → variants → annotation) |
| `README.md` | This file |
| `LICENSE` | MIT license |
| `.gitignore` | Excludes large sequencing files (FASTQ/BAM/VCF/assemblies/indices) |

Large intermediate files (reads, alignments, assemblies, VCFs, indices) are **not**
committed — see [`.gitignore`](.gitignore). They regenerate by re-running the workflow.

## Attribution

> **Please update the bracketed fields below before sharing.**

This work was completed as a hands-on assignment during a bacterial-genomics /
WGS bioinformatics training:

- **Training / workshop:** _[training or workshop name]_
- **Conducted by / instructor:** _[instructor or organizing institution]_
- **Date:** _[month, year]_

The workflow design and tool choices follow what was taught during the training.
The commands in this repository are my own worked implementation of the assigned
exercise. All datasets used are publicly available (NCBI SRA); no proprietary
course material, server credentials, or other participants' personal information
is redistributed here.

## License

Code in this repository is released under the [MIT License](LICENSE). The
underlying bioinformatics tools and SRA datasets retain their own respective
licenses and terms of use.
