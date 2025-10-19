"""
Affine Gap Penalty Global Alignment Algorithm
Python implementation for bioinformatics sequence analysis
"""
import time

def affine_global_alignment(seq1: str, seq2: str, match_score: int = 2, mismatch_score: int = -1, gap_open: int = -3, gap_extend: int = -1):
    """
    Perform global sequence alignment with affine gap penalties.
    
    Uses three matrices to track different alignment states:
    - M: Match/mismatch state
    - I: Insertion state (gap in seq1)
    - D: Deletion state (gap in seq2)
    
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
    gap_open : int
        Penalty for opening a gap (default: -3)
    gap_extend : int
        Penalty for extending a gap (default: -1)
        
    Returns
    -------
    tuple
        (aligned_seq1, aligned_seq2, alignment_score)
    """
    len1 = len(seq1)
    len2 = len(seq2)
    
    # Initialize three scoring matrices
    # M: match/mismatch, I: insertion (gap in seq1), D: deletion (gap in seq2)
    M = [[float('-inf') for _ in range(len2 + 1)] for _ in range(len1 + 1)]
    I = [[float('-inf') for _ in range(len2 + 1)] for _ in range(len1 + 1)]
    D = [[float('-inf') for _ in range(len2 + 1)] for _ in range(len1 + 1)]
    
    # Initialize traceback matrices
    M_trace = [["" for _ in range(len2 + 1)] for _ in range(len1 + 1)]
    I_trace = [["" for _ in range(len2 + 1)] for _ in range(len1 + 1)]
    D_trace = [["" for _ in range(len2 + 1)] for _ in range(len1 + 1)]
    
    # Base case
    M[0][0] = 0
    M_trace[0][0] = "done"
    
    # Initialize first row and column
    for i in range(1, len1 + 1):
        D[i][0] = gap_open + (i - 1) * gap_extend
        D_trace[i][0] = "D"
        
    for j in range(1, len2 + 1):
        I[0][j] = gap_open + (j - 1) * gap_extend
        I_trace[0][j] = "I"
    
    # Fill the matrices
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            # Calculate match/mismatch score
            if seq1[i-1] == seq2[j-1]:
                match_val = match_score
            else:
                match_val = mismatch_score
            
            # M matrix: can come from M, I, or D
            m_from_m = M[i-1][j-1] + match_val
            m_from_i = I[i-1][j-1] + match_val
            m_from_d = D[i-1][j-1] + match_val
            
            M[i][j] = m_from_m
            M_trace[i][j] = "M"
            
            if m_from_i > M[i][j]:
                M[i][j] = m_from_i
                M_trace[i][j] = "I"
            
            if m_from_d > M[i][j]:
                M[i][j] = m_from_d
                M_trace[i][j] = "D"
            
            # I matrix: gap in seq1 (insertion in seq2)
            i_from_m = M[i][j-1] + gap_open
            i_from_i = I[i][j-1] + gap_extend
            
            I[i][j] = i_from_m
            I_trace[i][j] = "M"
            
            if i_from_i > I[i][j]:
                I[i][j] = i_from_i
                I_trace[i][j] = "I"
            
            # D matrix: gap in seq2 (deletion from seq2)
            d_from_m = M[i-1][j] + gap_open
            d_from_d = D[i-1][j] + gap_extend
            
            D[i][j] = d_from_m
            D_trace[i][j] = "M"
            
            if d_from_d > D[i][j]:
                D[i][j] = d_from_d
                D_trace[i][j] = "D"
    
    # Find the best final score
    final_score = M[len1][len2]
    final_state = "M"
    
    if I[len1][len2] > final_score:
        final_score = I[len1][len2]
        final_state = "I"
    
    if D[len1][len2] > final_score:
        final_score = D[len1][len2]
        final_state = "D"
    
    # Traceback
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = len1
    j = len2
    current_state = final_state
    
    while i > 0 or j > 0:
        if current_state == "M":
            if M_trace[i][j] == "done":
                break
            
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            
            next_state = M_trace[i][j]
            i -= 1
            j -= 1
            current_state = next_state
            
        elif current_state == "I":
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            
            next_state = I_trace[i][j]
            j -= 1
            current_state = next_state
            
        elif current_state == "D":
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            
            next_state = D_trace[i][j]
            i -= 1
            current_state = next_state
    
    return (aligned_seq1, aligned_seq2, int(final_score))

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
    
    # Perform affine gap penalty global alignment
    aligned1, aligned2, score = affine_global_alignment(seq1, seq2)
    
    end_time = time.time()
    runtime_ms = (end_time - start_time) * 1000
    
    return runtime_ms

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python affine_global_alignment.py <file1.fa> <file2.fa>")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    runtime = run_alignment_test(file1, file2)
    print(f"Python runtime: {runtime:.2f}ms")