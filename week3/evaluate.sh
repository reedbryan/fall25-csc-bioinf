#!/bin/bash

echo "=== Week 3 Phylogenetic Algorithm Evaluation ==="
echo ""

# Run Python tests and capture timing
echo "Running Python tests..."
python_output=$(python test_python_clean.py 2>&1)
echo "Python output: $python_output"

# Extract runtime - look for any line containing "runtime" and numbers
python_time=$(echo "$python_output" | grep -i "Total.*runtime:" | sed 's/.*runtime: //' | sed 's/ms.*//')
if [ -z "$python_time" ]; then
    # Check if it's a biotite availability error
    if echo "$python_output" | grep -q "biotite" && echo "$python_output" | grep -q "not"; then
        python_time="BIOTITE_MISSING"
    else
        python_time="ERROR"
    fi
fi

# Run Codon tests and capture timing  
echo "Running Codon tests..."
codon_output=$(codon run test_codon_clean.py 2>&1)
echo "Codon output: $codon_output"

# Extract runtime - look for any line containing "runtime" and numbers
codon_time=$(echo "$codon_output" | grep -i "Total.*runtime:" | sed 's/.*runtime: //' | sed 's/ms.*//')
if [ -z "$codon_time" ]; then
    codon_time="ERROR"
fi

# Display results in requested format
echo ""
echo "Language    Runtime"
echo "-------------------"
if [ "$python_time" = "ERROR" ]; then
    echo "python      FAILED"
elif [ "$python_time" = "BIOTITE_MISSING" ]; then
    echo "python      BIOTITE_NOT_INSTALLED"
else
    echo "python      ${python_time}ms"
fi

if [ "$codon_time" = "ERROR" ]; then
    echo "codon       FAILED"  
else
    echo "codon       ${codon_time}ms"
fi

echo ""
echo "Evaluation complete!"