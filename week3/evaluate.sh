#!/bin/bash

echo "=== Week 3 Phylogenetic Algorithm Evaluation ==="
echo ""

# Create temporary files for Python and Codon versions
cp test.py test_python.py
cp test.py test_codon.py

# Set Python version flag
sed -i 's/__codon__ = False/__codon__ = False/' test_python.py

# Set Codon version flag  
sed -i 's/__codon__ = False/__codon__ = True/' test_codon.py

# Run Python tests and capture timing
echo "Running Python tests..."
python_output=$(python test_python.py 2>&1)
python_time=$(echo "$python_output" | grep "Total Python runtime:" | sed 's/.*: //' | sed 's/ms//')

# Run Codon tests and capture timing
echo "Running Codon tests..."
codon_output=$(~/.codon/bin/codon run test_codon.py 2>&1)
codon_time=$(echo "$codon_output" | grep "Total Codon runtime:" | sed 's/.*: //' | sed 's/ms//')

# Clean up temporary files
rm test_python.py test_codon.py

# Display results in requested format
echo ""
echo "Language    Runtime"
echo "-------------------"
echo "python      ${python_time}ms"
echo "codon       ${codon_time}ms"

echo ""
echo "Evaluation complete!"