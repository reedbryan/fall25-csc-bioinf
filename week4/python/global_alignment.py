"""
Global Sequence Alignment using Needleman-Wunsch Algorithm
Python implementation for bioinformatics sequence analysis
"""
import time

def needleman_wunsch(seq1: str, seq2: str, match_score: int = 2, mismatch_score: int = -1, gap_penalty: int = -2):
    """
    Perform global sequence alignment using Needleman-Wunsch algorithm.
    
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
    
    # Initialize first row and column
    for i in range(len1 + 1):
        score_matrix[i][0] = i * gap_penalty
        if i > 0:
            traceback[i][0] = "up"
    
    for j in range(len2 + 1):
        score_matrix[0][j] = j * gap_penalty
        if j > 0:
            traceback[0][j] = "left"
    
    traceback[0][0] = "done"
    
    # Fill the scoring matrix
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            # Calculate scores for three possible moves
            if seq1[i-1] == seq2[j-1]:
                diagonal_score = score_matrix[i-1][j-1] + match_score
            else:
                diagonal_score = score_matrix[i-1][j-1] + mismatch_score
            
            up_score = score_matrix[i-1][j] + gap_penalty
            left_score = score_matrix[i][j-1] + gap_penalty
            
            # Choose the maximum score
            max_score = diagonal_score
            direction = "diagonal"
            
            if up_score > max_score:
                max_score = up_score
                direction = "up"
            
            if left_score > max_score:
                max_score = left_score
                direction = "left"
            
            score_matrix[i][j] = max_score
            traceback[i][j] = direction
    
    # Traceback to construct alignment
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = len1
    j = len2
    
    while traceback[i][j] != "done":
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
    
    final_score = score_matrix[len1][len2]
    return (aligned_seq1, aligned_seq2, final_score)

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
    
    # Perform alignment
    aligned1, aligned2, score = needleman_wunsch(seq1, seq2)
    
    end_time = time.time()
    runtime_ms = (end_time - start_time) * 1000
    
    return runtime_ms

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python global_alignment.py <file1.fa> <file2.fa>")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    runtime = run_alignment_test(file1, file2)
    print(f"Python runtime: {runtime:.2f}ms")