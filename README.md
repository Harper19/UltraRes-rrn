# UltraRes-rrn


## Workflow

| Step | Script | What it does |
|---|---|---|
| 1 | `03_Preprocess_basecalling_QC.sh` | dorado SUP basecalling and demultiplexing, NanoComp QC, chopper length/quality filtering |
| 2 | `04_UMI-based_trim_filter.sh` | Extract UMI pairs from both read ends, cluster (usearch), build per-barcode UMI reference, remove chimeric UMIs |
| 3 | `05_Reads_binning.sh` | Map reads to the UMI reference (bwa) and bin reads by UMI with strand, error and bin/cluster-size filters |
| 4 | `06_Consensus_Polishing_parallel.sh` | Per-bin seed selection, racon polishing, medaka polishing (SLURM array) |
| 5 | `07_minimap2.sh` | Extract 16S / 23S / ITS / rrn / V3V4 regions (ipcress) and map to the reference database (CoverM, minimap2) |
| 6 | `07_classifier_EMU.sh` | Train and run RDP classifiers on the extracted regions |

## Software

dorado, NanoComp, chopper, cutadapt, seqtk, gawk, GNU parallel, usearch v11.0.667, bwa, samtools,
minimap2, racon, medaka, ipcress, CoverM, seqkit, RDP classifier, conda/anaconda.
Reference database: GROND (GTDB-derived rrn database).
Record the exact versions you use; the dorado model used here was `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`.

## Attribution

The UMI binning (`05_*`) and consensus polishing (`06_*`) steps are adapted from
**longread_umi** (https://github.com/SorenKarst/longread_umi, GPL-3.0):

> Karst SM, Ziels RM, Kirkegaard RH, et al. High-accuracy long-read amplicon sequences using unique
> molecular identifiers with Nanopore or PacBio sequencing. *Nat Methods* 18, 165-169 (2021).
> https://doi.org/10.1038/s41592-020-01041-y

Modifications include adaptation to R10.4.1 data and medaka 2.x, SLURM array parallelisation,
and changed clustering and filtering parameters.

## License

GPL-3.0, see `LICENSE`. This follows the license of longread_umi, from which parts of this code are derived.

## Citation

If you use these scripts, please cite the preprint:

> Bian K, Sudarshan AS, Meng Z, Yang J, Bachmann M, Bott C, Graham KE, Pinto AJ. High-accuracy hierarchical rRNA operon profiling resolves genomovar-level taxonomy and microdiversity using Nanopore sequencing. *bioRxiv* (2026). https://doi.org/10.64898/2026.06.01.729381

and longread_umi (see Attribution).
