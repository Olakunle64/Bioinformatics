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
            minimum_values[key] = [seq_dict.get(key), len(seq_dict.get(key))]
    
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
            maximum_values[key] = [seq_dict.get(key), len(seq_dict.get(key))]
    
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
                identifier = line.strip('\n').split()[0].lstrip(">")
                # print(identifier)
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

def longest_orf(filename, frame=0):
    """Find the longest ORF in a fasta file for the specified reading frame."""
    seq_dict = get_sequence(filename)  # Assumes this function reads the fasta file
    longest_orf_info = {"identifier": None, "start": None, "length": 0}
    
    for identifier, sequence in seq_dict.items():
        index = frame  # Start at the specified frame
        while index < len(sequence) - 2:  # Ensure we have full codons
            if sequence[index:index + 3] == "ATG":  # Start codon found
                index_of_start_codon = index
                index += 3
                while index < len(sequence) - 2:
                    codon = sequence[index:index + 3]
                    if codon in ["TAA", "TAG", "TGA"]:  # Stop codon found
                        index_of_stop_codon = index
                        orf_length = index_of_stop_codon - index_of_start_codon + 3
                        if orf_length > longest_orf_info["length"]:
                            longest_orf_info = {
                                "identifier": identifier,
                                "start": index_of_start_codon,
                                "length": orf_length
                            }
                        break  # Exit the inner loop to find the next start codon
                    index += 3  # Continue looking for stop codons
            index += 1  # Continue searching for start codons
    
    return longest_orf_info  # Return the longest ORF's info

def get_longest_orf_of_sequence(identifier, filename):
    """Find the longest ORF (Open Reading Frame) of the sequence that belongs to the identifier."""
    seq_dict = get_sequence(filename)
    length_of_orfs = []
    
    if identifier not in seq_dict:
        return 0
    
    sequence = seq_dict[identifier]
    index = 0
    
    while index < len(sequence):
        if sequence[index:index + 3] == "ATG":
            index_of_start_codon = index
            index += 3
            
            while index < len(sequence):
                if sequence[index:index + 3] in ["TAA", "TAG", "TGA"]:
                    index_of_stop_codon = index
                    length_of_orfs.append(index_of_stop_codon - index_of_start_codon + 3)
                    index = index + 3  # Move past the stop codon
                    break  # Exit the inner loop to find the next start codon
                index += 3  # Move to the next codon (3 bases)
        else:
            index += 1  # Move to the next base if not a start codon
            
    return max(length_of_orfs) if length_of_orfs else 0



def repeated_substring(n, filename):
    """find all repeated substring of length n of a sequence in a fasta file"""
    seq_dict = get_sequence(filename)
    reads = seq_dict.values()
    dict_track = {}
    
    for read in reads:
        index = 0
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
