# 🧬 HLAfactor

<img src="Qijing_TCR.jpeg" width="600">

**HLAfactor** is a command-line tool for **HLA (Human Leukocyte Antigen) quality control (QC)** and **downstream immunogenetic analysis**.

It supports multiple analytical modules including heterozygosity, expression, supertypes, and NK-related ligand markers.

---

## 👤 Author
**Qijing Shen**

---

## ⏬ Installation

### 🔹 Download binaries

#### Linux (x86_64)
```bash
wget http://*****.gz
tar -xzf hlafactor_linux_x86_64.gz
chmod +x hlafactor
```

#### macOS
```bash

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
hlafactor -i <input_file> -o <output_folder> [analysis_option]
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
hlafactor -i data.vcf.gz -o results/ -hlarg38 -af 0.01 0.5
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
