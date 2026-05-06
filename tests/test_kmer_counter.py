from kmer_counter import validate_sequence, update_kmer_count, count_kmers_with_context

def test_validate_sequence_good_sequence():
    assert validate_sequence("ATGC", 2) == True

def test_validate_sequence_too_short():
    assert validate_sequence("A", 2) == False

def test_validate_sequence_numbers_invalid():
    assert validate_sequence("AT3C", 2) == False

def test_update_kmer_count_new_kmer():
    result = update_kmer_count({}, "AT", "G")
    assert result["AT"]["count"] == 1
    assert result["AT"]["next_chars"]["G"] == 1

def test_update_kmer_count_existing_kmer_same_next_char():
    data = {"AT": {"count": 1, "next_chars": {"G": 1}}}
    result = update_kmer_count(data, "AT", "G")
    assert result["AT"]["count"] == 2
    assert result["AT"]["next_chars"]["G"] == 2

def test_count_kmers_with_context():
    result = count_kmers_with_context("ATGT", 2)
    assert result["AT"]["count"] == 1
    assert result["AT"]["next_chars"]["G"] == 1
    assert result["TG"]["count"] == 1
    assert result["TG"]["next_chars"]["T"] == 1
