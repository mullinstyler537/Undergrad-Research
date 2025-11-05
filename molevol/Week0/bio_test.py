# bio_test.py
from Bio.Seq import Seq

print("Biopython import OK")

s = Seq("ATGCGT")
print("Sequence:", s)
print("Reverse complement:", s.reverse_complement())
print("Length:", len(s))