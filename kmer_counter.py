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

    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}
    
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
            
            f.write(f"{kmer} {next_char_str}\n")


def main():
    """
    Runs the k-mer analysis program from the command line.

    The program:
    1. Reads sequence from an input file
    2. Validates each sequence
    3. Counts k-mers and following characters
    4. Writes results to an output file

    Command line arguments:
        sys.argv[1]: Input sequence filename
        sys.argv[2]: k-mer length
        sys.argv[3]: Output filename

    Returns:
        None
    """
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]
    
    print(f"Reading sequences from {sequence_file}...")

    with open(sequence_file, 'r') as f:
        for sequence in f:
            sequence = sequence.strip()

            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence")
                continue
            
            kmer_data = count_kmers_with_context(sequence, k) 
            
            write_results_to_file(kmer_data, output_file)

if __name__ == '__main__':
    main()
