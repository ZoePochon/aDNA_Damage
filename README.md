# aDNA_Damage

**Description**  
This repository contains scripts designed to handle ancient DNA (aDNA) damage patterns, specifically targeting common transition mutations associated with postmortem degradation in ancient DNA samples. The tools here aim to improve the accuracy of downstream analyses by addressing artefactual C->T and G->A transitions commonly observed in ancient samples.

---

## Contents

- **adna_sslib_damage_removal.py**  
  This script processes single-stranded ancient DNA libraries to remove C->T substitutions, which are often indicative of ancient DNA damage, rather than true genetic mutations. The script reconstructs the reference sequence using the CIGAR and MD tags and then removes C->T substitutions on the forward strand and G->A substitutions on the reverse strand. This script is especially useful for single-stranded library preparations typical of ancient DNA studies.

- **rmTrans_and_keep_trueMut_for_minDepth.py**
  This script removes C->T and G->A transitions that are typical of ancient DNA damage, while retaining transitions at positions likely representing true mutations based on user-defined thresholds. It allows for the identification of genuine mutations in low-quality, ancient DNA samples by applying a minimum read depth and a transition frequency threshold.

---

## Usage

### adna_sslib_damage_removal.py

**Description**
This script removes C->T substitutions from single-stranded ancient DNA libraries, which are indicative of ancient DNA damage rather than genuine mutations. It reconstructs the reference sequence based on the CIGAR and MD tags, allowing for accurate identification and removal of damage-specific transitions. This is particularly useful for researchers working with single-stranded library preparations of ancient DNA.

**Arguments**
- `-b`, `--bamfile`: Path to the input BAM file.
- `-o`, `--output`: Path to the output BAM file where processed reads will be saved.

**Example Command**
```bash
python adna_sslib_damage_removal.py -b input.bam -o output.bam
```

### rmTrans_and_keep_trueMut_for_minDepth.py

**Description**  
This script processes a BAM file to remove C->T and G->A transitions unless they appear as true mutations, based on specific thresholds. The script can be customised to set minimum read depth and transition frequency thresholds, making it flexible for different data qualities and research needs.

**Arguments**  
- `-b`, `--bamfile`: Path to the input BAM file.
- `-o`, `--output`: Path to the output BAM file where processed reads will be saved.
- `-r`, `--reference`: Path to the reference FASTA file.
- `-t`, `--threshold`: Threshold for transition percentage to consider as a true mutation. (Default: 0.9)
- `-d`, `--min_depth`: Minimum read depth required to consider a position as a true mutation. (Default: 3)

**Example Command**
```bash
python rmTrans_and_keep_trueMut_for_minDepth.py -b input.bam -o output.bam -r reference.fasta -t 0.9 -d 3
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

- **Zoé Pochon** - Developer of `rmTrans_and_keep_trueMut_for_minDepth.py`  
  [GitHub Profile](https://github.com/ZoePochon)
