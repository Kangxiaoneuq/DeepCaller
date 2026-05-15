# DeepCaller

<p align="center">
  <img src="https://img.shields.io/badge/version-1.0.0-blue" alt="Version">
  <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  <img src="https://img.shields.io/badge/python-3.9-blue" alt="Python">
  <img src="https://img.shields.io/badge/platform-linux-lightgrey" alt="Platform">
</p>

**DeepCaller** is a deep learning–based variant caller for the accurate detection of SNPs and small indels in polyploid genomes from short-read sequencing data. It provides pre-trained models for tetraploid and hexaploid crops and supports both speed-optimized and accuracy-optimized inference modes.

> **Note**: This repository accompanies a manuscript currently under review. Full use of the software is permitted upon publication. See [LICENSE](LICENSE) for details.

---

## Table of Contents

- [Features](#features)
- [Supported Species](#supported-species)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Output](#output)
- [Citation](#citation)

---

<p align="center">
  <img src="docs/flow.png" alt="DeepCaller Workflow" width="800">
</p>

---

## Features

- **Polyploid-aware genotyping** — native support for tetraploid and hexaploid genomes with ploidy-specific genotype encoding
- **Five pre-trained models** — covering potato, alfalfa, modern rose, sweetpotato, and a potato-derived synthetic hexaploid
- **Dual inference modes** — `speed` mode for rapid genome-wide screening; `performance` mode for maximum accuracy
- **Whole-genome depth normalisation** — sequencing depth is estimated from all chromosomes regardless of the target region, ensuring robust variant filtering
- **GPU-optional** — runs on CPU-only servers; GPU acceleration is supported when available
- **Standard output** — bgzip-compressed, tabix-indexed VCF 4.2

---

## Supported Species

| `--species`   | Common name                | Ploidy | Training dataset |
|---------------|----------------------------|--------|------------------|
| `potato`      | Tetraploid potato          | 4x     | C88              |
| `alfalfa`     | Alfalfa                    | 4x     | Bolivia          |
| `rose`        | Modern rose                | 4x     | Samantha         |
| `sweetpotato` | Sweetpotato                | 6x     | Tanzania         |
| `syn_potato`  | Synthetic hexaploid potato | 6x     | SyntheticPotato  |

---

## Installation

### Requirements

- Linux (x86_64)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) ≥ 4.10
- samtools, mosdepth, bgzip, tabix (all installed automatically via the conda environment)

### Steps

```bash
# 1. Clone the repository
git clone https://github.com/JiaoLab2021/DeepCaller.git
cd DeepCaller

# 2. Create and activate the conda environment
conda env create -f DeepCaller_env.yml
conda activate DeepCaller_env

# 3. Install DeepCaller
pip install -e .

# 4. Verify installation
DeepCaller --version
```

---

## Quick Start

A small demo dataset (chromosome 10, 1 Mb region; tetraploid potato C88 at ~20× coverage) is provided in the `Demo/` directory.

```bash
cd Demo

DeepCaller \
    -r DM8.1_chr10_100000_1100000.fa \
    -b C88_20x_chr10_100000_1100000.bam \
    -p 4 \
    --mode speed \
    -o demo_output.vcf
```

> `--species` is omitted here — DeepCaller automatically selects `potato` as the default model when `--ploidy 4` is set.

Expected runtime: < 2 minutes on a 24-core server.

---

## Usage

```
DeepCaller -r <REF> -b <BAM> -p <PLOIDY> [options]
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `-r`, `--ref` | Reference FASTA file (must be indexed) |
| `-b`, `--bam` | Input BAM file (must be coordinate-sorted) |
| `-p`, `--ploidy` | Ploidy level: `4` or `6` |

### Input/output configuration

| Argument | Default | Description |
|----------|---------|-------------|
| `-o`, `--output` | `output.vcf` | Output VCF file (will be bgzip-compressed) |
| `-c`, `--chroms` | all | Chromosomes to process (space-separated) |
| `-l`, `--bed` | — | BED file restricting variant calling to target regions; overrides `--chroms` |

### Processing options

| Argument | Default | Description |
|----------|---------|-------------|
| `-s`, `--species` | auto | Species model — defaults to `potato` for ploidy 4, `sweetpotato` for ploidy 6; see [Supported Species](#supported-species) |
| `-m`, `--mode` | `speed` | Inference mode: `speed` or `performance` |
| `-t`, `--cpus` | `24` | CPU threads; use `-1` for all available |
| `--min_vaf` | `0.10` | Minimum allele frequency at candidate sites |
| `--rd_floor` | `10` | Minimum read depth at candidate sites |
| `--no_gpu` | — | Disable GPU acceleration |

### Example commands

```bash
# Tetraploid potato, whole genome, performance mode
DeepCaller -r ref.fa -b sample.bam -p 4 --mode performance -o out.vcf -t 32

# Hexaploid sweetpotato, specific chromosomes
DeepCaller -r ref.fa -b sample.bam -p 6 -c chr1 chr2 chr3 -o out.vcf

# Alfalfa, target regions only (BED file)
DeepCaller -r ref.fa -b sample.bam -p 4 --species alfalfa -l targets.bed -o out.vcf
```

---

## Output

DeepCaller produces a bgzip-compressed, tabix-indexed VCF file (`<output>.gz` and `<output>.gz.tbi`).

### FORMAT fields

| Field | Description |
|-------|-------------|
| `GT`  | Polyploid genotype (e.g. `0/0/0/1` for tetraploid simplex) |
| `GQ`  | Genotype quality (Phred-scaled) |
| `DP`  | Read depth at the site |
| `AD`  | Allelic depth (ref, alt) |
| `AF`  | Allele frequency |

---

## Citation

If you use DeepCaller in your research, please cite:

> 

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.  
Full use is permitted upon official publication of the accompanying manuscript.

---

## Contact

Kang Xiao · [xiaokangneuq@163.com](mailto:xiaokangneuq@163.com)