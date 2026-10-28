# UltraRes-rrn

## Workflow

| Step | Script | What it does |
|---|---|---|
| 1 | `Preprocess_basecalling_QC.sh` | dorado SUP basecalling and demultiplexing (GPU), chopper length/quality filtering |
| 2 | `UMI-based_trim_filter.sh` | Extract UMI pairs from both read ends, cluster (usearch), build per-barcode UMI reference, remove chimeric UMIs |
| 3 | `Reads_binning.sh` | Map reads to the UMI reference (bwa) and bin reads by UMI with strand, error and bin/cluster-size filters |
| 4 | `Consensus_Polishing_parallel.sh` | Per-bin seed selection, racon polishing, medaka polishing (GPU; run per barcode) |

## Software

dorado, chopper, cutadapt, seqtk, gawk, GNU parallel, usearch v11.0.667, bwa, samtools,
minimap2, racon, medaka, conda/anaconda.
Record the exact versions you use; the dorado model used here was `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`.

## Attribution

The UMI binning and consensus polishing steps are adapted from
**longread_umi** (https://github.com/SorenKarst/longread_umi, GPL-3.0):

> Karst SM, Ziels RM, Kirkegaard RH, et al. High-accuracy long-read amplicon sequences using unique
> molecular identifiers with Nanopore or PacBio sequencing. *Nat Methods* 18, 165-169 (2021).
> https://doi.org/10.1038/s41592-020-01041-y

Modifications include adaptation to R10.4.1 data and medaka 2.x, per-barcode parallelisation,
and changed clustering and filtering parameters.

## License

GPL-3.0, see `LICENSE`. This follows the license of longread_umi, from which parts of this code are derived.

## Citation

If you use these scripts, please cite the preprint:

> Bian K, Sudarshan AS, Meng Z, Yang J, Bachmann M, Bott C, Graham KE, Pinto AJ. High-accuracy hierarchical rRNA operon profiling resolves genomovar-level taxonomy and microdiversity using Nanopore sequencing. *bioRxiv* (2026). https://doi.org/10.64898/2026.06.01.729381

and longread_umi (see Attribution).
