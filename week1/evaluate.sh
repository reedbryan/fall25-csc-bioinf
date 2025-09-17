#!/bin/bash

# Change to the code directory first
cd "$(dirname "$0")/code"

# Set dataset name
# Print header
echo -e "Dataset\tLanguage\tRuntime\t\tN50"
echo "-------------------------------------------------------------------------------------------------------"

# Function to format runtime as MM:SS:CC
format_runtime() {
    local total_seconds=$1
    local minutes=$((total_seconds / 60))
    local seconds=$((total_seconds % 60))
    local centiseconds=00
    if [[ $total_seconds == *.* ]]; then
        local decimal_part=$(echo $total_seconds | awk -F. '{print $2}')
        centiseconds=$(printf "%02d" $(echo "0.$decimal_part * 100" | bc | awk -F. '{print $1}'))
        total_seconds=$(echo $total_seconds | awk -F. '{print $1}')
        minutes=$((total_seconds / 60))
        seconds=$((total_seconds % 60))
    fi
    printf "%02d:%02d:%02d" $minutes $seconds $centiseconds
}

for DATASET in data1 data2 data3 data4; do

    # Run Python version and measure runtime
    echo -n -e "${DATASET}\tpython\t\t"
    python_start_time=$(date +%s)
    python_output=$(python3 main.py ${DATASET} 2>&1)
    python_exit_code=$?
    python_end_time=$(date +%s)
    python_runtime=$((python_end_time - python_start_time))
    python_runtime_formatted=$(format_runtime $python_runtime)

    # Extract N50 only if the command was successful
    if [ $python_exit_code -eq 0 ]; then
        python_n50=$(echo "${python_output}" | grep "N50:" | awk '{print $2}')
        echo -e "${python_runtime_formatted}\t${python_n50}"
    else
        echo -e "${python_runtime_formatted}\tERROR"
        echo "Python Error:"
        echo "${python_output}"
    fi

    echo

    # Run Codon version with release optimizations
    echo -n -e "${DATASET}\tcodon\t\t"
    codon_start_time=$(date +%s)
    codon_output=$(~/.codon/bin/codon run -release main_co.codon ${DATASET} 2>&1)
    codon_exit_code=$?
    codon_end_time=$(date +%s)
    codon_runtime=$((codon_end_time - codon_start_time))
    codon_runtime_formatted=$(format_runtime $codon_runtime)

    # Extract N50 only if the command was successful
    if [ $codon_exit_code -eq 0 ]; then
        codon_n50=$(echo "${codon_output}" | grep "N50:" | awk '{print $2}')
        echo -e "${codon_runtime_formatted}\t${codon_n50}"
    else
        echo -e "${codon_runtime_formatted}\tERROR"
        echo "Codon Error:"
        echo "${codon_output}"
    fi

    echo
done
codon_start_time=$(date +%s)
codon_output=$(~/.codon/bin/codon run -release main_co.codon ${DATASET} 2>&1)
codon_exit_code=$?
codon_end_time=$(date +%s)
codon_runtime=$((codon_end_time - codon_start_time))
codon_runtime_formatted=$(format_runtime $codon_runtime)

# Extract N50 only if the command was successful
if [ $codon_exit_code -eq 0 ]; then
    codon_n50=$(echo "${codon_output}" | grep "N50:" | awk '{print $2}')
    echo -e "${codon_runtime_formatted}\t${codon_n50}"
else
    echo -e "${codon_runtime_formatted}\tERROR"
    echo "Codon Error:"
    echo "${codon_output}"
fi