# 427 Deliverable 4
### Reed Bryan

## Step 1: Understanding the Assignment Requirements
- Reviewed the assignment specification requiring implementation of four sequence alignment algorithms
- Identified the need for implementing:
  - Global alignment (Needleman-Wunsch)
  - Local alignment (Smith-Waterman) 
  - Semi-global alignment
  - Affine gap penalty global alignment
- Set up performance comparison between Python and Codon implementations

## Step 2: Implementing Global Alignment Algorithm
- Created global alignment using dynamic programming
- Aligns sequences from end-to-end
- Used scoring matrices and traceback to find optimal alignment
- Built both Python and Codon versions with FASTA file reading

## Step 3: Implementing Local Alignment Algorithm  
- Created local alignment to find similar regions within sequences
- Key difference: negative scores reset to zero
- Starts traceback from highest score instead of bottom corner
- Good for finding conserved domains or motifs

## Step 4: Implementing Semi-Global Alignment Algorithm
- Created algorithm for fitting one sequence inside another
- No penalty for gaps at the ends of the longer sequence
- Useful for finding genes within larger DNA sequences
- Traceback starts from best score in the last column

## Step 5: Implementing Affine Gap Penalty Global Alignment
- Improved gap penalty model with different costs for opening vs extending gaps
- Uses three matrices instead of one (more complex but more realistic)
- Gap opening costs -3, gap extension costs -1
- Takes longer to run but gives better biological results

## Step 6: Creating Evaluation Framework
- Built script to test all four algorithms automatically
- Used real DNA data (human vs orangutan mitochondrial sequences)
- Measured timing for both Python and Codon versions
- Generated comparison table showing runtime differences

## Step 7: Testing and Results
- All algorithms work in both Python and Codon
- Created separate files for each algorithm to avoid conflicts
- Automated testing shows performance differences between languages
- Framework makes it easy to add new algorithms later

## Key Lessons Learned
- Different alignment algorithms solve different biological problems
- More complex algorithms (like affine gaps) are slower but more accurate
- Codon generally runs faster than Python for these computational tasks
- Good project organization helps when comparing multiple implementations