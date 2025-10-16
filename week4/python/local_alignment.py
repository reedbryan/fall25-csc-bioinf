"""
Local Sequence Alignment using Smith-Waterman Algorithm
Python implementation for bioinformatics sequence analysis
"""
import time

def smith_waterman(seq1: str, seq2: str, match_score: int = 2, mismatch_score: int = -1, gap_penalty: int = -2):
    """
    Perform local sequence alignment using Smith-Waterman algorithm.
    
    Parameters
    ----------
    seq1 : str
        First sequence to align
    seq2 : str
        Second sequence to align
    match_score : int
        Score for matching characters (default: 2)
    mismatch_score : int
        Score for mismatching characters (default: -1)
    gap_penalty : int
        Penalty for gaps/indels (default: -2)
        
    Returns
    -------
    tuple
        (aligned_seq1, aligned_seq2, alignment_score)
    """
    len1 = len(seq1)
    len2 = len(seq2)
    
    # Initialize scoring matrix
    score_matrix = [[0 for _ in range(len2 + 1)] for _ in range(len1 + 1)]
    
    # Initialize traceback matrix
    traceback = [["" for _ in range(len2 + 1)] for _ in range(len1 + 1)]
    
    # Track maximum score and position
    max_score = 0
    max_i = 0
    max_j = 0
    
    # Fill the scoring matrix (Smith-Waterman modification)
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            # Calculate scores for three possible moves
            if seq1[i-1] == seq2[j-1]:
                diagonal_score = score_matrix[i-1][j-1] + match_score
            else:
                diagonal_score = score_matrix[i-1][j-1] + mismatch_score
            
            up_score = score_matrix[i-1][j] + gap_penalty
            left_score = score_matrix[i][j-1] + gap_penalty
            
            # Choose the maximum score, but never below 0 (key difference from global)
            max_val = diagonal_score
            direction = "diagonal"
            
            if up_score > max_val:
                max_val = up_score
                direction = "up"
            
            if left_score > max_val:
                max_val = left_score
                direction = "left"
            
            # Smith-Waterman: set negative scores to 0
            if max_val < 0:
                max_val = 0
                direction = "stop"
            
            score_matrix[i][j] = max_val
            traceback[i][j] = direction
            
            # Track the maximum score position
            if max_val > max_score:
                max_score = max_val
                max_i = i
                max_j = j
    
    # Traceback from maximum score position
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = max_i
    j = max_j
    
    while i > 0 and j > 0 and traceback[i][j] != "stop" and score_matrix[i][j] > 0:
        if traceback[i][j] == "diagonal":
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback[i][j] == "up":
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif traceback[i][j] == "left":
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        else:
            break
    
    return (aligned_seq1, aligned_seq2, max_score)

def read_fasta(filename: str) -> str:
    """
    Read a FASTA file and return the sequence.
    
    Parameters
    ----------
    filename : str
        Path to FASTA file
        
    Returns
    -------
    str
        DNA/RNA/protein sequence
    """
    sequence = ""
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line.upper()
    return sequence

def run_alignment_test(file1: str, file2: str) -> float:
    """
    Run alignment test on two FASTA files and return runtime.
    
    Parameters
    ----------
    file1 : str
        Path to first FASTA file
    file2 : str  
        Path to second FASTA file
        
    Returns
    -------
    float
        Runtime in milliseconds
    """
    start_time = time.time()
    
    # Read sequences
    seq1 = read_fasta(file1)
    seq2 = read_fasta(file2)
    
    # Truncate sequences for reasonable runtime (first 1000 chars)
    seq1 = seq1[:1000]
    seq2 = seq2[:1000]
    
    # Perform local alignment
    aligned1, aligned2, score = smith_waterman(seq1, seq2)
    
    end_time = time.time()
    runtime_ms = (end_time - start_time) * 1000
    
    return runtime_ms

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python local_alignment.py <file1.fa> <file2.fa>")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    runtime = run_alignment_test(file1, file2)
    print(f"Python runtime: {runtime:.2f}ms")