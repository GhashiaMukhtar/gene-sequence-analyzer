from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd
import matplotlib.pyplot as plt
import os


os.makedirs("results", exist_ok=True)
file_path = input("Enter FASTA file path: ")

data = []
for record in SeqIO.parse(file_path, "fasta"):
    seq = record.seq
    species = " ".join(record.description.split()[1:3])
    
    data.append({
        "ID": record.id,
        "Species": species,
        "Length": len(seq),
        "A": seq.count("A"),
        "T": seq.count("T"),
        "C": seq.count("C"),
        "G": seq.count("G"),
        "GC%": gc_fraction(seq) * 100,
    })

df = pd.DataFrame(data)
df.to_csv("results/sequence_statistics.csv", index=False)

summary = f"""
GENE SEQUENCE ANALYSIS SUMMARY
==============================

Dataset: ls_orchid.fasta (94 orchid chloroplast genes)

Total sequences          : {len(df)}
Average length           : {df['Length'].mean():.1f} bp
Longest sequence         : {df['Length'].max()} bp
Shortest sequence        : {df['Length'].min()} bp
Average GC content       : {df['GC%'].mean():.2f}%
Highest GC%              : {df['GC%'].max():.2f}%  → {df.loc[df['GC%'].idxmax(), 'Species']}
Lowest GC%               : {df['GC%'].min():.2f}%  → {df.loc[df['GC%'].idxmin(), 'Species']}
Standard deviation (length): {df['Length'].std():.1f} bp

Generated on: {pd.Timestamp('today').strftime('%Y-%m-%d')}
By: Ghashia Mukhtar | QAU Bioinformatics
"""

with open("results/summary.txt", "w",encoding='utf-8') as f:
    f.write(summary)

# 1. GC% — sorted, top/bottom 30
df_sorted = df.sort_values("GC%", ascending=False)
top_bottom = pd.concat([df_sorted.head(15), df_sorted.tail(15)])

plt.figure(figsize=(10, 8))
plt.barh(top_bottom["Species"], top_bottom["GC%"], color="teal")
plt.xlabel("GC %")
plt.title("GC Content Across Orchid Species (Top & Bottom 15)")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig("results/gc_content.png", dpi=300)
plt.close()

# 2. Nucleotide composition — top 10
top10 = df.nlargest(10, "Length")
top10.plot(x="Species", y=["A", "T", "C", "G"], kind="bar", stacked=True, figsize=(10, 6))
plt.title("Nucleotide Composition (Top 10 Longest)")
plt.ylabel("Count")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("results/base_composition.png", dpi=300)
plt.close()

# 3. Length histogram
plt.figure(figsize=(8, 6))
plt.hist(df["Length"], bins=20, color="skyblue", edgecolor="black")
plt.title("Distribution of Sequence Lengths")
plt.xlabel("Length (bp)")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig("results/length_distribution.png", dpi=300)
plt.close()

# 4. Average pie
avg = [df["A"].mean(), df["T"].mean(), df["C"].mean(), df["G"].mean()]
labels = ["A", "T", "C", "G"]
plt.figure(figsize=(7, 7))
plt.pie(avg, labels=labels, autopct="%1.1f%%", colors=["#ec1212","#0f77de","#13ce13","#e27e1a"])
plt.title("Average Nucleotide Composition")
plt.savefig("results/average_composition.png", dpi=300)
plt.close()