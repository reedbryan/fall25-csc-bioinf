#!/usr/bin/env python3
"""
Test the PSSM distribution functionality (which uses thresholds module)
"""

import sys
import os
# Add the code directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code'))

import numpy as np
from matrix import FrequencyPositionMatrix

def test_pssm_distribution():
    """Test PSSM distribution calculation."""
    print("Testing PSSM distribution functionality...")
    
    # Create a simple frequency matrix for DNA motif
    alphabet = "ACGT"
    counts = {
        'A': [5, 2, 1, 3],
        'C': [2, 5, 6, 1], 
        'G': [1, 2, 2, 5],
        'T': [2, 1, 1, 1]
    }
    
    # Create matrices
    fpm = FrequencyPositionMatrix(alphabet, counts)
    pwm = fpm.normalize()
    pssm = pwm.log_odds()
    
    print(f"Created PSSM with length: {pssm.length}")
    
    # Test distribution calculation
    try:
        distribution = pssm.distribution(precision=100)
        print("✓ PSSM distribution calculated successfully")
        print(f"Distribution type: {type(distribution)}")
        
        # Test threshold calculations
        fpr_threshold = distribution.threshold_fpr(0.01)
        print(f"✓ FPR threshold (1%): {fpr_threshold:.3f}")
        
        fnr_threshold = distribution.threshold_fnr(0.01)
        print(f"✓ FNR threshold (1%): {fnr_threshold:.3f}")
        
        balanced_threshold = distribution.threshold_balanced()
        print(f"✓ Balanced threshold: {balanced_threshold:.3f}")
        
    except Exception as e:
        print(f"✗ Distribution test failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("\nPSSM distribution test completed!")

if __name__ == "__main__":
    test_pssm_distribution()