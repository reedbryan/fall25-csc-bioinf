#!/bin/bash

echo "=== Week 4 Alignment Performance Evaluation ==="
echo ""

# Test data files
TEST_FILES=(
    "data/q1.fa data/t1.fa"
    "data/MT-human.fa data/MT-orang.fa"
)

TEST_NAMES=(
    "q1"
    "mt_human"
)

# Alignment methods
METHODS=("global" "local")

echo "Method            Language    Runtime"
echo "--------------------------------------"

# Run tests for each method and data pair
for method in "${METHODS[@]}"; do
    for i in "${!TEST_FILES[@]}"; do
        files=(${TEST_FILES[$i]})
        file1=${files[0]}
        file2=${files[1]}
        test_name=${TEST_NAMES[$i]}
        
        # Check if files exist
        if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
            printf "%-17s python      FILE_NOT_FOUND\n" "${method}-${test_name}"
            printf "%-17s codon       FILE_NOT_FOUND\n" "${method}-${test_name}"
            continue
        fi
        
        # Run Python version
        python_output=$(cd python && python ${method}_alignment.py "../${file1}" "../${file2}" 2>&1)
        python_time=$(echo "$python_output" | grep "Python runtime:" | sed 's/.*: //' | sed 's/ms.*//')
        
        if [ -z "$python_time" ]; then
            echo "Python error output: $python_output"
            python_time="ERROR"
        fi
        
        # Run Codon version  
        codon_output=$(cd codon && codon run ${method}_alignment.codon "../${file1}" "../${file2}" 2>&1)
        codon_time=$(echo "$codon_output" | grep "Codon runtime:" | sed 's/.*: //' | sed 's/ms.*//')
        
        if [ -z "$codon_time" ]; then
            echo "Codon error output: $codon_output"
            codon_time="ERROR"
        fi
        
        # Display results
        printf "%-17s python      %sms\n" "${method}-${test_name}" "${python_time}"
        printf "%-17s codon       %sms\n" "${method}-${test_name}" "${codon_time}"
    done
done

echo ""
echo "Evaluation complete!"