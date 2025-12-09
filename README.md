# Gene Sequence Analyzer 🧬

A lightweight, beginner-friendly bioinformatics tool that analyzes DNA sequences from FASTA files and generates beautiful statistics and visualizations.

**Built with:**
- Biopython → sequence parsing
- Pandas → data handling & statistics
- Matplotlib → publication-ready plots

### Features
- Calculates sequence length, A/T/C/G counts, GC%
- Stores everything in a clean Pandas DataFrame
- Generates four insightful plots:
  - Nucleotide composition (stacked bar)
  - GC% across sequences (horizontal bar)
  - Sequence length distribution (histogram)
  - Average base composition (pie chart)
- Exports results to CSV + PNG files

### Data
Uses the classic orchid chloroplast dataset (`ls_orchid.fasta`) containing 94 sequences — perfect for testing and comparison across species.

### How to Run
```bash
git clone https://github.com/GhashiaMukhtar/gene-sequence-analyzer.git
cd gene-sequence-analyzer
python analyzer.py
