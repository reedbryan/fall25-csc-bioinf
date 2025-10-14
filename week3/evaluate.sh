#!/bin/bash

echo "=== Week 3 Phylogenetic Algorithm Evaluation ==="
echo ""

# Create temporary files for Python and Codon versions
cp test.py test_python.py
cp test.py test_codon.py

# Set Python version flag (leave as False)
# No change needed for test_python.py

# Remove python_imports.py to prevent Codon from trying to compile it
if [ -f python_imports.py ]; then
    mv python_imports.py python_imports.py.bak
fi

# Remove old codon_imports.py if it exists (we use codon_imports.codon now)
if [ -f codon_imports.py ]; then
    rm codon_imports.py
fi

# Set Codon version flag  
sed -i 's/__codon__ = False/__codon__ = True/' test_codon.py

# Run Python tests and capture timing
echo "Running Python tests..."
python_output=$(python test_python.py 2>&1)
echo "Python output: $python_output"

# Extract runtime - look for any line containing "runtime" and numbers
python_time=$(echo "$python_output" | grep -i "Total.*runtime:" | sed 's/.*runtime: //' | sed 's/ms.*//')
if [ -z "$python_time" ]; then
    # Check if it's a biotite availability error
    if echo "$python_output" | grep -q "biotite package not available"; then
        python_time="BIOTITE_MISSING"
    else
        python_time="ERROR"
    fi
fi

# Run Codon tests and capture timing  
echo "Running Codon tests..."
codon_output=$(codon run test_codon.py 2>&1)
echo "Codon output: $codon_output"

# Extract runtime - look for any line containing "runtime" and numbers
codon_time=$(echo "$codon_output" | grep -i "Total.*runtime:" | sed 's/.*runtime: //' | sed 's/ms.*//')
if [ -z "$codon_time" ]; then
    codon_time="ERROR"
fi

# Clean up temporary files
rm test_python.py test_codon.py

# Restore python_imports.py if it was backed up
if [ -f python_imports.py.bak ]; then
    mv python_imports.py.bak python_imports.py
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