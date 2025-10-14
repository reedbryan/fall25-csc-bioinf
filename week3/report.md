# 427 Deliverable 3
### Reed Bryan

## Step 1: Understanding the Assignment Requirements
- Reviewed the assignment specification requiring porting biotite's phylogenetic algorithms to Codon
- Identified the need for implementing UPGMA and neighbor-joining tree construction algorithms
- Set up test requirements:
  - `@test` decorators on all test functions
  - `if __codon__:` conditional imports for environment detection
  - Performance comparison between Python and Codon implementations

## Step 2: Analyzing the Biotite Phylo Module
- Studied biotite.sequence.phylo module to understand the phylogenetic algorithms
- Identified key components needed:
  - `upgma()` function for hierarchical clustering
  - `neighbor_joining()` function for tree construction
  - Tree data structures for representing phylogenetic relationships
- Created test cases using distance matrices to validate algorithm behavior

## Step 3: Porting to Codon
- Created biotite_codon package with Codon implementations:
  - `tree.codon` - TreeNode and Tree classes
  - `upgma.codon` - UPGMA clustering algorithm
  - `neighbor_joining.codon` - Neighbor-joining algorithm
  - `__init__.codon` - Package initialization
- Encountered challenges with Codon's type system and had to simplify some data structures

## Step 4: Resolving Compatibility Issues
- Fixed various issues with Codon compilation:
  - Simplified tree node relationships to avoid circular references
  - Removed complex type annotations that caused problems
  - Created separate test files for Python and Codon to avoid import conflicts
- Developed evaluation script to automate performance comparisons

## Step 5: Testing and Validation
- Created comprehensive test suite with three main tests:
  - Distance calculation validation
  - UPGMA algorithm verification
  - Neighbor-joining algorithm testing
- Implemented automated timing comparison showing runtime differences between implementations
- Successfully achieved working versions in both Python and Codon

## Key Moments
- Performance comparisons can be automated effectively with shell scripts
- Sometimes simplifying the project structure is necessary to achieve compatibility across different environments