#!/bin/bash

echo "=== Testing Global Alignment Algorithm ==="
echo ""

cd codon
echo "Running Codon implementation..."
codon run global_alignment.codon

echo ""
echo "Test complete!"