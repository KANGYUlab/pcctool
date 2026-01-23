# PAFconfidenceChain(PCC)

![pccgenechain](config/pccgenechain.svg)
**PAFconfidenceChain** is a high-confidence genomic alignment filtering and chaining tool constrained by **gene anchors**. Through rigorous **synteny screening** and sequence identity verification, it extracts high-confidence **anchored genes** from complex whole-genome alignments and utilizes them to construct stable **syntenic block** mappings. This provides the most reliable underlying coordinate mapping for **variant calling (VCF) mapping**, functional annotation transfer, and comparative genomics research.


# Introduction to pcctool

**pcctool** is a specialized toolkit designed for high-precision genomic coordinate transformation. It offers the following core advantages:

1. **High-Accuracy Anchor Filtering**
   The tool identifies anchor genes based on strict synteny and sequence identity. By effectively filtering out multi-copy genes and paralogs, it ensures that only the most reliable genomic anchors are used for chaining.

2. **Streamlined Coordinate Transformation**
   Converting coordinates between the **yao** assembly and the **hg38** reference is highly simplified. The tool provides a clean and concise workflow, making complex liftover tasks straightforward and efficient.

3. **Reliability Grading via chainQ**
   Each mapping result is assigned a **chainQ score**, which provides a quantitative grade for the alignment's reliability. This allows users to easily filter and categorize results based on their specific confidence requirements.


# pcctool
##### ***Warning:*** This is a ***0-based*** system.
### PCCfile Release 1.0
##### Data Path Declaration
The high-confidence mapping files are located at:
* **Maternal Haplotype:** data/yaomat2hg38v1.0.pcc
* **Paternal Haplotype:** data/yaopat2hg38v1.0.pcc

---

##### Data Summary

The pcc for both haplotypes are summarized below:

| Metric | yaomat (Maternal) | yaopat (Paternal) |
| :--- | :--- | :--- |
| **Blocks Coverage** | **42.36%** | **42.34%** |
| **Protein-coding Genes** | 14,351 (74.7%) | 14,387 (74.7%) |
| **lncRNA Genes** | 11,717 (74.4%) | 11,717 (74.0%) |




### Usage and Options
```bash
git clone https://github.com/KANGYUlab/pcctool.git
cd pcctool
python setup.py install
```
`pcctool` provides a command-line interface for coordinate liftover and mapping between **hg38** and **yao** genomes.

```bash
usage: pcctool [-h] -src {hg38,yao} -dest {hg38,yao} -p {mat,pat} -chr CHR 
                        [-pos POS] [-Seg START END] 
                        [-posfile POSFILE] [-Segfile SEGFILE] [-V]
                        (> outfile)
```

Below is a summary of the available options for `pcctool`:
| Argument | Description |
| :--- | :--- |
| `-h, --help` | Show this help message and exit. |
| **`-src {hg38,yao}`** | **Source genome** for mapping. |
| **`-dest {hg38,yao}`** | **Destination genome** for mapping. |
| **`-p, --parent {mat,pat}`** | **Parental haplotype**: `mat` (maternal) or `pat` (paternal). |
| **`-chr CHR`** | **Chromosome ID** (e.g., `chr1`, `chr2`...). |
| `-pos POS` | A single **genomic position** to map. |
| `-Seg START END` | A **genomic region** defined by two integers (e.g., `100 200`). |
| `-posfile POSFILE` | Path to a `.txt` file containing a **list of coordinates** (e.g., `chr1 110124436`). |
| `-Segfile SEGFILE` | Path to a `.txt` file containing a **list of regions** (e.g., `chr1 110124436 110124578`). |
| `-V, --version` | Show program version. |

---
##### 1. Locus Mapping

Convert a specific genomic segment from the reference genome to the target genome:

```bash
pcctool -src hg38 -dest yao -p mat -chr chr2 -Seg 110124436 110124578
```

**Output:** hg38 chr2:[110124436, 110124578] → yaomat chr2:[110610091, 110610233] (strand=-, ChainQ=4.9991)

##### 2. Position Mapping

Pinpoint the orthologous position of a single nucleotide:

```bash
pcctool -src hg38 -dest yao -p mat -chr chr2 -pos 110124436
```
**Output:** hg38 chr2:110124436 → yaomat chr2:110610233 (strand=-, ChainQ=4.9991)
##### 3. FILE Mapping
use -posfile or -locusfile to map a list of positions or loci.
```bash
pcctool -src hg38 -dest yao -p mat -chr chr2 -posfile ./data/test.txt 
```
**Output:**
hg38 chr2:110124436 → yaomat chr2:110610233 (strand=-, ChainQ=4.9991)
hg38 chr2:110124578 → yaomat chr2:110610091 (strand=-, ChainQ=4.9991)
##### 4. Reverse Mapping

Perform reverse mapping by simply swapping the `-src` and `-dest` parameters:

```bash
pcctool -src yao -dest hg38 -p mat -chr chr2_hap1 -Seg 110610091 110610233
```


# PCC Format
| Col  | Type   | Description                                      |
| ---- | ------ | ------------------------------------------------ |
| 1    | string | Query sequence name                              |
| 2    | int    | Query sequence block length                      |
| 3    | int    | Query start (0-based; BED-like; closed)          |
| 4    | int    | Query end (0-based; BED-like; open)              |
| 5    | char   | Relative strand: "+" or "-"                      |
| 6    | string | Target sequence name                             |
| 7    | int    | Target sequence length                           |
| 8    | int    | Target start on original strand (0-based)        |
| 9    | int    | Target end on original strand (0-based)          |
| 10   | int    | Number of residue matches                        |
| 11   | int    | Alignment block length                           |
| 12   | int    | Mapping quality (0-255; 255 for missing)         |
| 13   | float  | Number of nucleotides matches in alignment block |
| 14   | float  | Alignment block synteny conservation             |
| 15   | float  | Alignment block Completeness                     |
| 16   | float  | chainQ(Confidence score  0-5)                    |
| 17   | char   | Cigar list (e.g., "10M2I4M1D3M")                 |
| 18   | char   | block name                                       |



## chainQ(Confidence Score Calculation) 

The confidence score is comprehensively evaluated using the following formula, reflecting sequence quality, local structural stability, and alignment continuity:

$$
chainQ = 3 \times \left( 1 - \frac{\text{Edit Distance}}{\text{Alignment Length}} \right) + 1 \times S_{\text{conservation}} + 1 \times B_{\text{completeness}}
$$

---

#### 1. Flanking Conservation (S)

This metric scores stability by evaluating the conservation of the 4 neighboring blocks (2 upstream and 2 downstream) relative to the current block:

* **1.00**: All 4 neighboring blocks remain syntenic (Conserved)
* **0.75**: 3 neighboring blocks remain syntenic
* **0.50**: 2 neighboring blocks remain syntenic
* **0.25**: 1 neighboring block remains syntenic



#### 2. Block Alignment Completeness (B)

This metric reflects the physical continuity (degree of fragmentation) of the alignment block across the genome:

$$
B_{\text{completeness}} = 
\begin{cases} 
1.0 & \text{Entire block aligned with full integrity (Integrity)} \\
0.8 & \text{Block split into two segments} \\
0.6 & \text{Block split into multiple segments (Multi-segment)} 
\end{cases}
$$


# hg38 Whole-Genome High-chainQ Chains

![hg38linkyaomatpccConfidencemap](config/yaomat2hg38_pcc.svg)

<center>Figure 1. high chainQ Mapping of GRCh38 Genome (from yaomat)</center>

![hg38linkyaopatpccConfidencemap](config/yaopat2hg38_pcc.svg)

<center>Figure 2. high chainQ Mapping of GRCh38 Genome (from yaopat)</center>

