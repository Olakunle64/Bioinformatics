# from Bio.Seq import Seq
# dna = "TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG"
# seq = Seq(dna)
# print(seq.translate())

from Bio.Blast import NCBIWWW, NCBIXML
from io import StringIO

query_sequence = "TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG"
result_seq = NCBIWWW.qblast("blastn", "nt", query_sequence)

# Read the result and store it as a string
result_data = result_seq.read()

# Convert the string into a file-like object using StringIO
result_handle = StringIO(result_data)

# Parse the BLAST output
blast_record = NCBIXML.read(result_handle)

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:  # High-scoring segment pairs
        print(f"****Alignment****")
        print(f"Sequence: {alignment.title}")
        print(f"Length: {alignment.length}")
        print(f"e value: {hsp.expect}")
        print(f"Query: {hsp.query}")
        print(f"Match: {hsp.match}")
        print(f"Subject: {hsp.sbjct}")
    break  # Only print the top hit
