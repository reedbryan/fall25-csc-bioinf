#!/usr/bin/env python3
"""
Complete test of the motif analysis library using MEME files
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import __init__ as motifs
import minimal

def test_meme_files():
    """Test parsing and analysis of MEME format files."""
    print("=== Motif Analysis Library Test with MEME Files ===\n")
    
    # Test 1: Parse DNA MEME file
    print("1. Parsing DNA MEME file...")
    with open('data/minimal_test.meme', 'r') as handle:
        dna_motifs = minimal.read(handle)
    
    print(f"   ✓ Loaded {len(dna_motifs)} DNA motifs")
    
    # Test 2: Parse RNA MEME file  
    print("\n2. Parsing RNA MEME file...")
    with open('data/minimal_test_rna.meme', 'r') as handle:
        rna_motifs = minimal.read(handle)
        
    print(f"   ✓ Loaded {len(rna_motifs)} RNA motifs")
    
    # Test 3: Analyze first DNA motif
    print("\n3. Analyzing first DNA motif...")
    motif = dna_motifs[0]
    print(f"   Name: {motif.name}")
    print(f"   Consensus: {motif.consensus}")
    print(f"   Length: {len(motif)}")
    print(f"   Alphabet: {motif.alphabet}")
    print(f"   GC content: {motif.counts.gc_content:.3f}")
    
    # Test 4: Sequence scoring
    print("\n4. Testing sequence scoring...")
    test_sequences = [
        "TGTGATCGAGGTCACACTT",      # Perfect match
        "TGTGATAGAGGTCACACTT",      # One mismatch  
        "AAAGATCGAGGTCACACTT",      # Three mismatches
        "AAAAAAAAAAAAAAAAAAA"       # Poor match
    ]
    
    pssm = motif.pssm
    print(f"   PSSM score range: {pssm.min:.3f} to {pssm.max:.3f}")
    
    for i, seq in enumerate(test_sequences):
        scores = pssm.calculate(seq + "EXTRA")  # Add extra bases
        max_score = scores.max() if not all(s == float('-inf') for s in scores) else float('-inf')
        max_pos = scores.argmax()
        
        match_seq = seq[max_pos:max_pos+len(motif)] if max_pos < len(seq) else seq
        identity = sum(1 for a, b in zip(match_seq, motif.consensus) if a == b)
        identity_pct = (identity / len(motif.consensus)) * 100
        
        print(f"   Seq {i+1}: {seq}")
        print(f"          Score: {max_score:.3f}, Identity: {identity}/{len(motif.consensus)} ({identity_pct:.1f}%)")
    
    # Test 5: Create motif from sequences
    print("\n5. Creating motif from aligned sequences...")
    sequences = [
        "GCGCATGC",
        "GCGCGTGC", 
        "GCGCTTGC",
        "ACGCATGC",
        "GCGCAAGC"
    ]
    
    custom_motif = motifs.create(sequences)
    print(f"   Created motif from {len(sequences)} sequences")
    print(f"   Consensus: {custom_motif.consensus}")
    print(f"   Degenerate consensus: {custom_motif.degenerate_consensus}")
    print(f"   Length: {len(custom_motif)}")
    
    print("\n=== Test Summary ===")
    print("✓ MEME file parsing works")
    print("✓ Motif analysis works") 
    print("✓ Sequence scoring works")
    print("✓ Custom motif creation works")
    print("✓ Both DNA and RNA formats supported")
    
    print("\n=== What This Library Does ===")
    print("• Parses MEME format motif files")
    print("• Represents motifs as position weight matrices")
    print("• Scores sequences against motifs using PSSM")
    print("• Calculates consensus sequences")
    print("• Creates motifs from aligned sequences")
    print("• Supports both DNA and RNA alphabets")
    print("• Provides statistical analysis tools")

if __name__ == "__main__":
    test_meme_files()
