# aDNA_Damage

**Description**  
This repository contains scripts designed to handle ancient DNA (aDNA) damage patterns, specifically targeting common transition mutations associated with postmortem degradation in ancient DNA samples. The tools here aim to improve the accuracy of downstream analyses by addressing artefactual C->T and G->A transitions commonly observed in ancient samples.

---

## Contents

- **adna_sslib_damage_removal.py**  
  This script processes single-stranded ancient DNA libraries to remove damage-specific substitutions. Taking into account the characteristics of single-stranded libraries, this script reconstructs the reference sequence using the CIGAR and MD tags and then removes C->T transitions on the forward strand and G->A transitions on the reverse strand.

- **maskDamage.py**
  This script processes both single-stranded and double-stranded libraries undistinctly to remove damage-specific substitutions. It removes C->T and G->A transitions, which are common damage patterns in ancient DNA, while retaining transitions at positions likely representing true mutations based on the transition frequency and the minimum read depth at a position. 

---

## Usage

### adna_sslib_damage_removal.py

**Arguments**
- `-b`, `--bamfile`: Path to the input BAM file.
- `-o`, `--output`: Path to the output BAM file where processed reads will be saved.

**Example Command**
```bash
python adna_sslib_damage_removal.py -b input.bam -o output.bam
```

### maskDamage.py

**Arguments**  
- `-b`, `--bamfile`: Path to the input BAM file.
- `-o`, `--output`: Path to the output BAM file where processed reads will be saved.
- `-r`, `--reference`: Path to the reference FASTA file.
- `-t`, `--threshold`: Threshold for transition frequency to consider as a true mutation. (Default: 0.9)
- `-d`, `--min_depth`: Minimum read depth required to consider a position as a true mutation. (Default: 3)

**Example Command**
```bash
python maskDamage.py -b input.bam -o output.bam -r reference.fasta -t 0.9 -d 3
```
---

## Installation

Both scripts require Python 3, and the following Python libraries must be installed:
- `pysam` (for BAM file handling)
- `tqdm` (for progress display)

You can install these dependencies via `pip`:
```bash
pip install pysam tqdm
```
---

## Authors

- **Bilal Sharif** - Developer of `adna_sslib_damage_removal.py`  
  [GitHub Profile](https://github.com/bilalbioinfo)

- **Benjamin Guinet** - Contributor to `adna_sslib_damage_removal.py`  
  [GitHub Profile](https://github.com/BenjaminGuinet)

- **Zoé Pochon** - Developer of `maskDamage.py`  
  [GitHub Profile](https://github.com/ZoePochon)
