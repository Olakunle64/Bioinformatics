"""This module has functions that perform various process on a
    fasta file from next generation sequencing
"""

def record_count_in_fasta(filename):
    """return the number of records in a fasta file"""
    count = 0
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(">"):
                count += 1
        # print("record count", count)
        return count

def find_all_minimum_seq(seq_dict):
    """Find the shortest sequence(s) in a fasta file
    
        return both the identifier and the sequence as a key/value pair
    """
    minimum_values = {}
    if not seq_dict:
        return minimum_values
    new_dict = {key: len(value) for key, value in seq_dict.items()}
    keys = new_dict.keys()
    values = new_dict.values()
    min_value = min(values)
    for key, value in new_dict.items():
        if value == min_value:
            minimum_values[key] = seq_dict.get(key)
    
    return minimum_values

def find_all_maximum_seq(seq_dict):
    """Find the longest sequence(s) in a fasta file
    
        return both the identifier and the sequence as a key/value pair
    """
    maximum_values = {}
    if not seq_dict:
        return maximum_values
    new_dict = {key: len(value) for key, value in seq_dict.items()}
    keys = new_dict.keys()
    values = new_dict.values()
    max_value = max(values)
    for key, value in new_dict.items():
        if value == max_value:
            maximum_values[key] = seq_dict.get(key)
    
    return maximum_values

def get_sequence(filename):
    """Get all sequence in a fasta file
    
        use the identifer as the key and the sequence as
        the value
        
        return a dictionary
    """
    seq_dict = {}
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(">"):
                identifier = line.strip('\n')
                seq_dict[identifier] = ""
                # print(line)
                continue
            # print(line)
            seq_dict[identifier] += line.strip('\n')
        return seq_dict

def find_maximum_values(my_dict):
    """return maximum value(s) of a dict including their key"""
    maximums = {}
    keys = my_dict.keys()
    values = my_dict.values()
    max_value = max(values)
    for key, value in my_dict.items():
        if value == max_value:
            maximums[key] = my_dict.get(key)
    return maximums

def longest_orf(filename):
    """find the longest ORF in a fasta file"""
    seq_dict = get_sequence("dna.example.fasta")
    length_of_orfs = {}
    for identifier, sequence, in seq_dict.items():
        start_codon = {}
        stop_codon = {}
        index = 0
        while index < len(sequence):
            if sequence[index: index + 3] == "ATG":
                index_of_start_codon = index
                index += 3
                j = index
                while j < len(sequence):
                    if sequence[j: j + 3] in ["TAA", "TAG", "TGA"]:
                        index_of_stop_codon = j
                        length_of_orfs[identifier] = len(sequence[index_of_start_codon:index_of_stop_codon + 3])
                        index += 3
                        break
                    j += 1
            index += 1
    print("all ORF in the fasta file", length_of_orfs)
    return find_maximum_values(length_of_orfs)

def get_longest_orf_of_sequence(identifier, filename):
    """find the longest of sequence that belongs to identifier"""
    seq_dict = get_sequence(filename)
    length_of_orfs = []
    maximum_orfs = {}
    index = 0
    if not seq_dict.get(identifier):
        return 0
    sequence = seq_dict.get(identifier)
    while index < len(sequence):
        if sequence[index: index + 3] == "ATG":
            index_of_start_codon = index
            index += 3
            j = index
            while j < len(sequence):
                if sequence[j: j + 3] in ["TAA", "TAG", "TGA"]:
                    index_of_stop_codon = j
                    length_of_orfs.append(len(sequence[index_of_start_codon:index_of_stop_codon + 3]))
                    index += 3
                    break
                j += 1
        index += 1
    return max(length_of_orfs)


def repeated_substring(n, filename):
    """find all repeated substring of length n of a sequence in a fasta file"""
    seq_dict = get_sequence(filename)
    reads = seq_dict.values()
    dict_track = {}
    index = 0
    for read in reads:
        while index < len(read):
            sub_string = read[index: index + n]
            if len(sub_string) != n:
                index += 1
                continue
            # print(len(sub_string), n)
            if sub_string in dict_track:
                dict_track[sub_string] += 1
            else:
                dict_track[sub_string] = 1
            index += 1
    print("Most frequent repeat", find_maximum_values(dict_track))
    return dict_track

