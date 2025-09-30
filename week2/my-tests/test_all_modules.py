#!/usr/bin/env python3
"""
Test all modules working together
"""

import sys
import os
# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def test_all_modules():
    print("Testing all module imports...")
    
    # Test individual modules
    try:
        import matrix
        print("✓ matrix module imported successfully")
    except ImportError as e:
        print(f"✗ matrix import failed: {e}")
    
    try:
        import thresholds
        print("✓ thresholds module imported successfully")
    except ImportError as e:
        print(f"✗ thresholds import failed: {e}")
    
    try:
        import minimal
        print("✓ minimal module imported successfully")
    except ImportError as e:
        print(f"✗ minimal import failed: {e}")
    
    # Test the main module
    try:
        import __init__ as motifs
        print("✓ main motifs module imported successfully")
    except ImportError as e:
        print(f"✗ main motifs import failed: {e}")
    
    print("\nAll module tests completed!")

if __name__ == "__main__":
    test_all_modules()