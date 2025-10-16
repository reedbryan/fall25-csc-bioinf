#!/bin/bash

echo "=== Week 4 Global Alignment Performance Evaluation ==="
echo ""

# Test data files
TEST_FILES=(
    "data/q1.fa data/t1.fa"
    "data/MT-human.fa data/MT-orang.fa"
)

TEST_NAMES=(
    "q1_vs_t1"
    "human_vs_orang"
)

echo "Data         Language    Runtime"
echo "-----------------------------------"

# Run tests for each data pair
for i in "${!TEST_FILES[@]}"; do
    files=(${TEST_FILES[$i]})
    file1=${files[0]}
    file2=${files[1]}
    test_name=${TEST_NAMES[$i]}
    
    # Check if files exist
    if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
        echo "${test_name}    python      FILE_NOT_FOUND"
        echo "${test_name}    codon       FILE_NOT_FOUND"
        continue
    fi
    
    # Run Python version
    echo "Running Python test for ${test_name}..."
    python_output=$(cd python && python global_alignment.py "../${file1}" "../${file2}" 2>&1)
    python_time=$(echo "$python_output" | grep "Python runtime:" | sed 's/.*: //' | sed 's/ms.*//')
    
    if [ -z "$python_time" ]; then
        python_time="ERROR"
    fi
    
    # Run Codon version  
    echo "Running Codon test for ${test_name}..."
    codon_output=$(cd codon && codon run global_alignment.codon "../${file1}" "../${file2}" 2>&1)
    codon_time=$(echo "$codon_output" | grep "Codon runtime:" | sed 's/.*: //' | sed 's/ms.*//')
    
    if [ -z "$codon_time" ]; then
        echo "Codon error output: $codon_output"
        codon_time="ERROR"
    fi
    
    # Display results
    printf "%-12s python      %sms\n" "${test_name}" "${python_time}"
    printf "%-12s codon       %sms\n" "${test_name}" "${codon_time}"
done

echo ""
echo "Evaluation complete!"