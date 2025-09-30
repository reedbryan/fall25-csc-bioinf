#!/usr/bin/env python3
"""
Test the matrix functionality with our motif analysis code.
"""

import sys
import os
# Add the code directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code'))

import numpy as np
from matrix import FrequencyPositionMatrix, PositionWeightMatrix

def test_matrix_functionality():
    """Test basic matrix operations."""
    print("Testing matrix functionality...")
    
    # Create a simple frequency matrix for DNA motif
    alphabet = "ACGT"
    counts = {
        'A': [2, 1, 0, 1],
        'C': [1, 2, 3, 0], 
        'G': [0, 1, 1, 2],
        'T': [1, 0, 0, 1]
    }
    
    # Test FrequencyPositionMatrix
    print("Creating FrequencyPositionMatrix...")
    fpm = FrequencyPositionMatrix(alphabet, counts)
    print(f"Motif length: {fpm.length}")
    print(f"Alphabet: {fpm.alphabet}")
    print(f"Consensus: {fpm.consensus}")
    print(f"GC content: {fpm.gc_content:.2f}")
    
    # Test PositionWeightMatrix
    print("\nCreating PositionWeightMatrix...")
    pwm = fpm.normalize()
    print("PWM created successfully")
    
    # Test PSSM (this will use our Python fallback)
    print("\nCreating PSSM...")
    try:
        pssm = pwm.log_odds()
        print("PSSM created successfully")
        print(f"Max score: {pssm.max:.2f}")
        print(f"Min score: {pssm.min:.2f}")
        
        # Test sequence scoring
        test_seq = "ACGTACGT"
        print(f"\nTesting sequence scoring with: {test_seq}")
        scores = pssm.calculate(test_seq)
        print(f"Scores: {scores}")
        
    except Exception as e:
        print(f"PSSM test failed: {e}")
    
    print("\nMatrix functionality test completed!")

if __name__ == "__main__":
    test_matrix_functionality()