NOTES for later:
- when importing the four og files, _pwm.c was sunk because we didnt have a way to compile the C file at runtime, solution: have AI rewrite the functionality in pure python.

- instead of trying to port the biopython test files, just download the data relevent to the 4 files we are porting and write my own test for them. this way we avoid port 5k lines and random unused dependencies