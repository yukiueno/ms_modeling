from Bio import SeqIO
from Bio.SeqUtils import GC
for seq_record in SeqIO.parse("out.fa", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    print("GC:")
    print(GC(seq_record.seq) / 200 * len(seq_record))
    print("AT:")
    print((100.0 - GC(seq_record.seq)) / 200 * len(seq_record))
