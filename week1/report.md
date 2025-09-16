# 427 Deliverable 1
### Reed Bryan

## Step 1: Producing the configs
- Cloned the repo to my local machine from the provided link
- Insalled Codon using the provided commands (im on mac so no prob there)
- Ran the program:
    - chmod +x main.py && python main.py data1
    - it created a config.fasta file
    - tried this for all datasets minus data4 (all good)

## Repo setup
- Set up my own repo (for submission) and created the file structure to match the specifications
- Set up the yml file and tried to test it with a push (did not work)
    - after some troubleshooting I discovered that the yml file could not be indented with tabs, so I changed it to spaces
    - pushed again and got the green checkmark
- Moved all the nessesary files from the original repo to mine
- Ran the program again in my repo to ensure it still works (it does)

## Calculating N50
- Consulted Copilot on how to calculate and display the N50 value after getting the contigs
- Created a new n50.py file (called by main.py) to calculate and display the n50 value

## Converting to Codon
- Consulted Copilot for Codon-specific syntax and static typing requirements
- Created Codon versions of all core files
    - Asked Copilot to fill out these files with codon code that mimics the functionality of the python code
    - Addressed some Codon type system issues (explicit type annotations for dictionaries, avoiding `None` type for integers)
- Debug some stuff with Copilot to ensure all file I/O and string handling was compatible with Codon

## Automation
- Added runtime tracking to `main_co.codon` to measure performance
- Used the `-release` flag with Codon for optimized execution
- Added runtime trackingn to `main.py` to measure python performance
- Compared runtime and N50 results between Python and Codon implementations using an automated evaluation script
- Noticed a significant speedup with Codon while maintaining identical N50 output