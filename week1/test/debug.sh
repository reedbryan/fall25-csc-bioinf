#!/bin/bash

echo "Checking data directory structure..."
ls -la ../data/

if [ -d "../data/data1" ]; then
  echo "data1 directory exists"
  ls -la ../data/data1
  
  if [ -f "../data/data1/short_1.fasta" ]; then
    echo "short_1.fasta exists"
    head -n 5 ../data/data1/short_1.fasta
  else
    echo "short_1.fasta not found!"
  fi
  
  if [ -f "../data/data1/short_2.fasta" ]; then
    echo "short_2.fasta exists"
    head -n 5 ../data/data1/short_2.fasta
  else
    echo "short_2.fasta not found!"
  fi
  
  if [ -f "../data/data1/long.fasta" ]; then
    echo "long.fasta exists"
    head -n 5 ../data/data1/long.fasta
  else
    echo "long.fasta not found!"
  fi
  
else
  echo "data1 directory not found!"
fi

echo "Running the Codon program with debug output..."
~/.codon/bin/codon run -plugin seq dbg_co.codon utils_co.codon n50_co.codon main_co.codon data1
