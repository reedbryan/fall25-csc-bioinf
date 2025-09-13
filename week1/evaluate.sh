#!/bin/bash

echo "Running main_co.codon standalone..."
# Change to the code directory first
cd "$(dirname "$0")/code"
# Run the codon file
~/.codon/bin/codon run main_co.codon data1

echo "Exit status: $?"