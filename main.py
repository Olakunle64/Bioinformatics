from fasta_functions import (
    longest_orf,
    get_longest_orf_of_sequence,
    record_count_in_fasta,
    repeated_substring,
    find_maximum_values,
    find_all_maximum_seq,
    find_all_minimum_seq,
    get_sequence,
)
print(record_count_in_fasta("dna2.fasta"))
seq_dict = get_sequence("dna2.fasta")
# print(find_all_maximum_seq(seq_dict))
# print(find_all_minimum_seq(seq_dict))
# print(longest_orf("dna2.fasta", frame=3))
# print(get_longest_orf_of_sequence("gi|142022655|gb|EQ086233.1|16", "dna2.fasta"))
# print(repeated_substring(6, "dna2.fasta"))
repeated_substring(7, "dna2.fasta")

