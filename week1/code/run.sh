#!/bin/bash

# Check if an argument was provided
if [ -z "$1" ]
then
    echo "Usage: ./run.sh <data_directory>"
    echo "Example: ./run.sh data1"
    exit 1
fi

echo "====================="
echo "Starting script execution"
echo "Working directory: $(pwd)"
echo "Data directory: $1"

echo "Files to be compiled:"
ls -la dbg_new.codon utils_co.codon n50_co.codon main_co.codon

echo "Checking if Codon is installed:"
which ~/.codon/bin/codon || echo "Codon not found"

echo "====================="
echo "Running Codon program..."
~/.codon/bin/codon run -plugin seq dbg_new.codon utils_co.codon n50_co.codon main_co.codon "$1"
RESULT=$?
echo "Command exit status: $RESULT"

if [ $RESULT -eq 0 ]; then
    echo "Program completed successfully"
    echo "Checking if output file was created:"
    ls -la "../data/$1/contig.fasta" || echo "Output file not found"
else
    echo "Program failed with exit code $RESULT"
fi
echo "====================="
