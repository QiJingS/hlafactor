# 🧬 HLAfactor

<img src="Qijing_TCR.jpeg" width="600">

**HLAfactor** is a command-line tool for **HLA (Human Leukocyte Antigen) quality control (QC)** and **downstream immunogenetic analysis**.

It supports multiple analytical modules including heterozygosity, expression, supertypes, and NK-related ligand markers.

---

## 👤 Author
**Qijing Shen**

---

## ⏬ Installation

### 🔹 Download 

### Linux (x86_64)
```bash
wget https://github.com/QiJingS/hlafactor/releases/download/v0.0.0/hlafactor-0.0.0-Linux-x86_64.tar.gz
tar -xzf hlafactor-0.0.0-Linux-x86_64.tar.gz
cd hlafactor-0.0.0-Linux-x86_64
chmod +x hlafactor
```
### macOS 
```bash
wget https://github.com/QiJingS/hlafactor/releases/download/v0.0.0/hlafactor-0.0.0-Darwin-arm64.tar.gz
tar -xzf hlafactor-0.0.0-Darwin-arm64.tar.gz
cd hlafactor-0.0.0-Darwin-arm64
chmod +x hlafactor
```
---

## 🚀 Basic Usage

```bash
hlafactor [options]
```

---

## ⚙️ Required Arguments

| Option | Description |
|--------|------------|
| `-i`   | Input file (.txt, .dosage, .vcf, .vcf.gz) |
| `-o`   | Output directory or file prefix |

---

## 💻 1. HLA Quality Control (QC)

```bash
hlafactor -i <input_file> -o <output_folder> -hlarg38|-hlarg19 \
          [-af <min> [<max>]] [-impr2 <min> [<max>]]
```

### Options

| Option | Description |
|--------|------------|
| `-hlarg38` | Use hg38 (GRCh38) coordinates |
| `-hlarg19` | Use hg19 (GRCh19) coordinates |
| `-af`      | Filter variants by allele frequency |
| `-impr2`   | Filter variants by imputation R² |

---

## 🧪 2. HLA Downstream Analysis

```bash
hlafactor -i <input_file> -o <output_folder> [option]
```

### Available analyses

| Option | Description |
|--------|------------|
| `-tapasin` | Compute tapasin dependence score |
| `-het`     | Compute general HLA heterozygosity |
| `-hetf`    | Compute functional heterozygosity (HLA-A/B/C) |
| `-exp`     | Estimate HLA-A and HLA-C expression levels |
| `-sup`     | Assign HLA supertypes (Class I & II) |
| `-lignk`   | Compute NK-related HLA markers (KIR ligands) |

---

## ⚠️ Important Notes

- ❌ QC and analysis options **cannot be combined**
- 📂 Output files are automatically named based on:
  - input filename  
  - analysis type suffix  
- 📁 Output will be written to:
  - specified folder, or  
  - current directory if not provided  

## ⚠️ Error Handling and Return Values
- HLAfactor performs strict validation during parsing and downstream computation to improve robustness and reproducibility.
- Supported missing / invalid values
  - Missing values such as NAN, ., or other invalid entries will be treated as missing and propagated as NA in the final output.
  - If an input value cannot be matched to a valid HLA allele or genotype, the corresponding result will be returned as NA instead of forcing an incorrect calculation.
- Imputed HLA probabilities
  - For imputation outputs that are reported as probabilities or dosages, HLAfactor supports both hard-call and probabilistic formats.
  - For feature calculation, probabilistic genotypes can be rounded to the nearest call when required by the downstream analysis.
  - For supertype assignment, dosage/probability values would be preserved in the output to retain more original imputation information.
- Locus-specific genotype requirements
  - For loci that must follow a strict diploid format, such as HLA-A, HLA-B, and HLA-C, HLAfactor will only compute results after valid allele pairs are matched.
  - If an imputation platform returns multiple candidate alleles, ambiguous pairs, or unmatched references, the corresponding sample-locus result will be returned as NA.
- Reference matching report
  - When alleles or allele pairs cannot be matched to the reference panel used by HLAfactor, the software will generate an additional report listing:
  - sample ID
  - gene / locus
  - unmatched allele or genotype
---

## ✅ Multi-platform imputation support
- HLAfactor is designed to accept HLA imputation results from multiple platforms and output formats, with consistent parsing and standardized downstream processing.
---
## 📥 Supported Input Formats
- `.txt`
- `.dosage`
- `.vcf`
- `.vcf.gz`

---

## 📤 Output

Output files are generated automatically with informative suffixes:

| Analysis | Example Output |
|----------|--------------|
| QC       | `sample_qc.txt` |
| Heterozygosity | `sample_het.txt` |
| Expression | `sample_expression.txt` |
| Supertypes | `sample_supertypes.txt` |

---

## 🧩 Example Commands

### QC with allele frequency filtering
```bash
hlafactor -i data.vcf.gz -o results/ -hlarg38 -af 0.01 0.5 -impr2 0.9 1.0
```

### HLA general heterozygosity
```bash
hlafactor -i data.txt -o results/ -het
```
### HLA functional heterozygosity
```bash
hlafactor -i data.txt -o results/ -hetf
```

### Expression analysis
```bash
hlafactor -i data.txt -o results/ -exp
```

### NK ligand markers
```bash
hlafactor -i data.txt -o results/ -lignk
```
---

## 🛠️ Features

- ✔ Flexible input formats  
- ✔ Modular analysis pipeline  
- ✔ Built-in QC filtering  
- ✔ Automatic output naming  
- ✔ Multiple immunogenetic metrics  

---

## 📜 License
(MIT / GPL ??????)

---

## 🙌 Acknowledgements
**** ****
