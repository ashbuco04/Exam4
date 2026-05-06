import sys

def validate_sequence(sequence, k):
    """
    Check whether a DNA sequence is long enough for k-mer analysis

	valid sequence must contain no numeric characters and least length k

    Parameters:
        sequence (str): DNA sequence to check
        k (int): Length of the k-mer

    Returns:
        bool: True if the sequence is valid, False otherwise
    """
# Check if sequence is shorter than k
    if len(sequence) < k:
        return False
    for nucleotide in sequence:
        if nucleotide in '1234567890':
            return False
    return True

def update_kmer_count(kmer_data, kmer, next_char):
    """
    Updates the frequency information for a k-mer and the character that follows it

    Parameters:
        kmer_data (dict): Dictionary storing k-mer counts and characters
        kmer (str): The current k-mer substring
        next_char (str): Character immediately following the k-mer

    Returns:
        dict: Updated k-mer dictionary with modified counts
    """
# Create a new dictionary for a k-mer if it does not exist
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}
# Increase the total count for this k-mer
    kmer_data[kmer]['count'] += 1
    
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0
    kmer_data[kmer]['next_chars'][next_char] += 1

    return kmer_data

def count_kmers_with_context(sequence, k):
    """
    Counts all k-mers in a sequence and records the character that follows each k-mer

    Parameters:
        sequence (str): DNA sequence to analyze
        k (int): Length of the k-mer

    Returns:
        dict: Dictionary containing k-mer counts and following character frequencies
    """

    kmer_data = {}
# Loop through every possible k-mer position in the sequence
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        next_char = sequence[i+k]
        
        kmer_data = update_kmer_count(kmer_data, kmer, next_char)
    
    return kmer_data


def write_results_to_file(kmer_data, output_filename):
    """
    Writes k-mer frequency results to an output file

    The output file contains each k-mer followed by the frequency
    of characters that occur immediately after it

    Parameters:
        kmer_data (dict): Dictionary of k-mer frequency data
        output_filename (str): Name of the output file

    Returns:
        None
    """
    sorted_kmers = sorted(kmer_data.keys())
    
    with open(output_filename, 'w') as f:
        for kmer in sorted_kmers:
            next_chars = kmer_data[kmer]['next_chars']
            
            next_char_str = " ".join(
                f"{char}:{freq}" 
                for char, freq in sorted(next_chars.items())
            )
# Write formatted k-mer information to the output file            
            f.write(f"{kmer} {next_char_str}\n")


def main():
    """
    Runs the k-mer analysis program from the command line.
    """

    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]

    print(f"Reading sequences from {sequence_file}...")

 # Create one dictionary for all sequences in the file
    kmer_data = {}

    with open(sequence_file, 'r') as f:
        for sequence in f:
            sequence = sequence.strip()

            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence")
                continue

            # Count k-mers for this sequence
            sequence_kmers = count_kmers_with_context(sequence, k)

            # Add this sequence's k-mer counts to the total data
            for kmer in sequence_kmers:
                for next_char in sequence_kmers[kmer]['next_chars']:
                    freq = sequence_kmers[kmer]['next_chars'][next_char]

                    for i in range(freq):
                        kmer_data = update_kmer_count(kmer_data, kmer, next_char)

# Write all combined results after reading the whole file
    write_results_to_file(kmer_data, output_file)
   
if __name__ == '__main__':
    main()
