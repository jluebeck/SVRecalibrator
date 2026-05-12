# SVRecalibrator

Tools to collect reads around structural-variant (SV) breakpoints and refine predictions using split-read evidence and optional de-novo scaffolds.

## Prerequisites

* Linux / macOS (bash)
* **git**
* **conda**

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/tshneour/SVRecalibrator.git
```

---

## 2. Create conda environment and install requirements

Use conda with strict channel priority to ensure consistent dependency resolution across conda-forge and bioconda:

```bash

conda create -n sv-analysis "python>=3.11,<=3.14.0a0"
conda activate sv-analysis
conda install -c bioconda biopython pysam samtools bwa spades=4.2.0
conda install pandas numpy natsort marisa-trie

conda activate sv-analysis
```

---

## Quickstart

1. **Setup**

```bash
conda activate sv-analysis
```

2. **Collect reads** for all AA breakpoints into an alignment table + per-breakpoint read pairs:

```
python /path/to/collect.py 350 path/to/AA_summaries/ path/to/sample.bam \
  --strict \
  -v \
  -f alignments.tsv
```

This writes:

* `alignments.tsv` (or `<bam_basename>.tsv` if `-f/--file` not provided)
* gzip’d FASTQs per breakpoint in `fastq/`:

  * `fastq/b_<chr1>_<pos1+1>_<chr2>_<pos2+1>_1.fastq.gz`
  * `fastq/b_<chr1>_<pos1+1>_<chr2>_<pos2+1>_2.fastq.gz`

3. **Run Scaffold and Split Read** and get a combined table:

```bash
python /path/to/refine.py alignments.tsv \
  --mode both \
  --fasta /path/to/genome.fa \
  --out-table final_augmented \
  -v
# Produces final_augmented.tsv
```

> Tip: **Do not** add “.tsv” to `--out-table`; the program appends “.tsv” automatically.

---

## Input format

### SV summary (`sum`) files

`collect.py` expects one or more **tab-separated (****`.tsv`****) files** describing structural-variant breakpoints with a fixed set of required columns.

Each TSV in the `sum/` directory represents a collection of SV breakpoints which are required to be **coordinate-sorted**. The filename is used only to derive a `sample` identifier or `amplicon` (taken as `filename.split("_")[1]`) if present.

#### Required columns

Each input TSV **must** contain the following columns (case-sensitive):

| Column name         | Type   | Description                                                                                                                  |
| ------------------- | ------ | -----------------------------------------------------------------------------------------------------------------------------|
| `chrom1`            | string | Chromosome of the first breakpoint end (e.g. `chr8`)                                                                         |
| `pos1`              | int    | 0-based genomic coordinate of the first breakpoint                                                                           |
| `chrom2`            | string | Chromosome of the second breakpoint end                                                                                      |
| `pos2`              | int    | 0-based genomic coordinate of the second breakpoint                                                                          |
| `sv_type`           | string | Structural variant type (i.e. `deletion`, `duplication`, `interchromosomal`, `inversion`, `foldback`)                        |
| `orientation`       | string | Breakpoint orientation as a 2-character string (`++`, `--`, `+-`, `-+`)                                                      |

Coordinates in the output tables are reported as **1-based**, but the input `pos1` / `pos2` values are treated as **0-based** internally. Note that `sv_type` is a case-sensitive field.

---

#### Optional columns

If provided, these columns (case-sensitive) will be included in the final output for comparison purposes.

| Column name         | Type   | Description                                                                               |
| ------------------- | ------ | ----------------------------------------------------------------------------------------- |
| `features`          | string | Arbitrary annotation or metadata for the breakpoint                                       |
| `read_support`      | int    | Number of reads supporting the breakpoint (used for reporting only)                       |
| `homology_length`   | int    | Length of homology reported for the breakpoint (may be 0)                                 |
| `homology_sequence` | string | Homology or inserted sequence (may be empty)                                              |

## Usage

### `collect.py`

Collect reads around SV breakpoints for refinement and write per-breakpoint paired FASTQs plus a combined TSV of read evidence.

```
usage: collect.py [-h] [-v] [--strict] [-f FILE] refine sum bam
```

**Positional arguments**

* `refine` (int): Radius (bp) around each breakpoint to fetch alignments (e.g., 500).
* `sum` (path): Directory containing SV summary `.tsv` files in the format described above.
* `bam` (path): Coordinate-sorted BAM with index (`.bai`).

**Options**

* `-v, --verbose`           Verbose output.
* `--strict`                Only keep read pairs whose alignments fully fall within the region(s) of interest.
* `-f, --file FILE`         Output TSV path (default: `<bam_basename>.tsv`).

**Outputs**

* TSV of reads (`split` and `nonsplit`) with columns including:

  * `break_chrom1`, `break_pos1`, `break_chrom2`, `break_pos2`
  * `break_sv_type`, `break_orientation`
  * `query_*` (read-level alignment details)
  * `proper_pair`, `split`, `amplicon`
  * `break_read_support`, `break_features`, `homology_len`, `homology_seq` (if present in the input SV summary TSVs)

* Per-breakpoint FASTQs in `fastq/`, gzipped:

  * `fastq/b_<chr1>_<pos1>_<chr2>_<pos2>_1.fastq.gz`
  * `fastq/b_<chr1>_<pos1>_<chr2>_<pos2>_2.fastq.gz`

### `refine.py`

Refine SV breakpoints using split-read evidence and/or local de-novo scaffold reconstruction.

```
usage: refine.py FILE [--mode {split,scaffold,both}]
                     [--out-table OUT] [--split-log PATH]
                     [--scaffold-log PATH] [--outdir DIR]
                     [-l | --list] [-b IDX [IDX ...]] [-v ...]
                     [--fasta FASTA]
```

**Required**

* `FILE` (path): TSV produced by `collect.py`.

**Modes**

* `split`: infer refined breakpoint positions and homology using split-read evidence only.
* `scaffold`: perform local assembly (SPAdes) around each breakpoint and align scaffolds to the reference.
* `both`: run scaffold refinement first, then augment with split-read results.

**Outputs**

* **Split mode** columns added:

  * `split_support`, `soft_support`
  * `sp_left_sv`, `sp_right_sv`
  * `sp_hom_len` (positive = homology, negative = insertion)
  * `hom` (homology or insertion sequence)
  * `tst_cords` (if `hom` is a templated insertion, this will show the relevant reference coordinates)

* **Scaffold mode** columns added:

  * `sc_pos1`, `sc_pos2` (refined breakpoint coordinates)
  * `sc_hom_len` (Int64; negative = insertion, 0 = abutting, positive = homology)
  * `sc_hom` (sequence)

* **Both**: combined scaffold and split-read columns in a single TSV.

---

## Batch processing

`batch_run.sh` runs `collect.py` → `refine.py` for every sample in a shared SV summary directory.

### Sample list format

A tab-separated file with at minimum two columns — sample name and BAM path (any additional columns are ignored):

```
ACHN    /data/bam/ACHN/ACHN.cs.rmdup.bam
BT474M1 /data/bam/BT474M1/BT474M1.cs.rmdup.bam
```

### SV summary file naming

Files in the SV summary directory must match either:
- `<SAMPLE>_amplicon<N>_SV_summary.tsv`  (one file per amplicon)
- `<SAMPLE>_SV_summary.tsv`              (single file per sample)

### Running

```bash
# All samples
bash batch_run.sh \
  --samples phase2_samples.txt \
  --sv-dir /data/SV_summaries \
  --outdir /data/SVRecalibrator_outputs \
  --fasta /data/ref/hg38.fa \
  --strict

# Single sample
bash batch_run.sh \
  --samples phase2_samples.txt \
  --sv-dir /data/SV_summaries \
  --outdir /data/SVRecalibrator_outputs \
  --fasta /data/ref/hg38.fa \
  --sample ACHN
```

**Options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-s / --samples` | required | Tab-separated sample list |
| `-v / --sv-dir` | required | Directory of SV summary TSV files |
| `-o / --outdir` | required | Output root directory |
| `-f / --fasta` | required (scaffold/both) | Indexed reference FASTA |
| `-m / --mode` | `both` | `split`, `scaffold`, or `both` |
| `-r / --radius` | `350` | Collect radius (bp) around each breakpoint |
| `-t / --threads` | `16` | SPAdes threads per sample |
| `--spades-timeout HOURS` | `2` | Per-breakpoint SPAdes time limit in hours; `0` = no limit |
| `--strict` | off | Only keep reads fully within the breakpoint region |
| `--sample NAME` | — | Run a single sample instead of all |

Each sample writes to `<outdir>/<SAMPLE>/`. A `done.flag` is created on success; re-runs skip completed samples automatically. If a sample directory exists without a `done.flag` (i.e., a previously interrupted run), the directory is cleared and the sample is re-run from scratch.

---

## Summarizing homology improvement

`summarize_homology.py` loads all `final_augmented.tsv` files from a batch output directory and reports how many SVs gained microhomology by split-read and/or scaffold approaches that were absent from the original AA calls.

```bash
# Print summary to stdout only
python summarize_homology.py /data/SVRecalibrator_outputs

# Save all outputs under a common prefix
python summarize_homology.py /data/SVRecalibrator_outputs -o /data/results/summary
# produces:
#   summary_summary.csv          — per-sample table
#   summary_venn.png/.pdf        — proportional Euler diagrams: Homology / Insertion / Blunt (AA vs SVRecalibrator)
#   summary_venn_combined.png/.pdf — single Euler diagram: any junction call (AA vs SVRecalibrator), with unrefined count
#   summary_hist.png/.pdf        — stacked junction-length histograms (AA / split / scaffold)

# Restrict to specific samples or raise the minimum homology length
python summarize_homology.py /data/SVRecalibrator_outputs \
  --sample ACHN BT474M1 \
  --min-hom-len 3 \
  -o /data/results/summary
```

**Options:**

| Flag | Default | Description |
|------|---------|-------------|
| `outdir` | required | Batch output directory (one subdirectory per sample) |
| `--sample NAME [...]` | all | Restrict to these sample(s) |
| `--min-hom-len N` | `1` | Minimum homology length to count as detected |
| `-o / --output-prefix PREFIX` | — | Prefix for all output files (CSV + plots × two formats) |

**Output sections (stdout):**

1. **Dataset overview** — samples loaded, total SVs, AA homology prevalence
2. **SVRecalibrator detection (all SVs)** — split / scaffold / both / either, with percentages
3. **Gain over AA** — among SVs with no AA homology: how many gained it and by which approach
4. **Concordance** — on SVs AA did call, how many SVRecalibrator confirmed
5. **Homology length distributions** — median / mean / min / max / IQR per approach
6. **Per-sample breakdown** — per-row counts for every sample

**Plots:**

- **Venn diagram** (`_venn`) — 1×3 proportional Euler diagrams showing AA vs SVRecalibrator overlap separately for Homology, Insertion, and Blunt junctions. Circle areas and overlaps are geometrically correct; labels shown as a legend. When AA data is absent, falls back to split-read vs scaffold.
- **Combined Venn diagram** (`_venn_combined`) — single Euler diagram showing overlap of any junction call (homology, insertion, or blunt) between AA and SVRecalibrator, with the number of unrefined SVs labelled prominently in the rectangle outside both circles.
- **Stacked histograms** (`_hist`) — junction lengths (negative = insertion, positive = homology) for AA, split reads, and scaffold on separate tracks with shared axes and uniform y-limits. Values outside ±50 bp are clipped into the ≤−50 and ≥50 bins. The top track includes insertion/homology region shading and legend.

Blunt junctions (`*_hom_len == 0`) are included in the histograms but excluded from all "detected" counts. Samples with empty `final_augmented.tsv` files are silently skipped.

---

## Working assumptions & tips

* **BAM**: Coordinate-sorted, indexed (`.bai` present).
* **SV summaries**: Input TSVs must conform exactly to the column specification above.
* **Mapping quality**: Reads with MAPQ > 15 (or mapped status) are retained.
* **FASTQs**: Reconstructed from reads overlapping each breakpoint pair and used by scaffold mode.
* **FASTA**: Required for scaffold mode; must be indexed first:

  ```bash
  samtools faidx /path/to/genome.fa
  ```
