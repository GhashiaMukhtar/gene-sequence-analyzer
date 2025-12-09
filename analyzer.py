from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd
import sys
try:
    file_name = input("Enter file path\n")
    with open(file_name,"r") as handle:
        data = []
        for record in SeqIO.parse(handle,"fasta"):
            name = record.description.split(" ")
            my_seq = Seq(record.seq)
            A = my_seq.count("A")
            G = my_seq.count("G")
            C = my_seq.count("C")
            T = my_seq.count("T")
            count = A+T+C+G
            gc_percenrage = (gc_fraction(my_seq))*100
            at_percentage = ((A+T)/count)*100
            newdata = {
                "Sequence ID": record.id,
                "Specie Name": name[1],
                "Length": len(record.seq),
                "Base count": count,
                "GC percentage": gc_percenrage,
                "AT percentage" : at_percentage
                        }
            data.append(newdata)

except FileNotFoundError:
    print("File doesn't exist")
    sys.exit(1)
except :
    print("Unexpected error")
    sys.exit(1)

df = pd.DataFrame(data)
with open ("results/sequence_statistics.csv","w") as f:
    f.write(f"""Total number of sequences: {len(df)}\n
Average sequence length: {df["Length"].mean()}\n
Longest sequence: {df["Length"].max()}\n
Shortest sequence: {df['Length'].min()}\n
Average GC percentage: {df["GC percentage"].mean()}\n
Highest GC%: {df['GC percentage'].max()}\n
Lowest GC%: {df['GC percentage'].min()}\n
Median sequence length: {df['Length'].median()}\n
Standard deviation of sequence lengths: {df['Length'].std()}\n
Top 3 species with highest GC content:\n {df.sort_values(by='Length', ascending=False).head(3)}""")