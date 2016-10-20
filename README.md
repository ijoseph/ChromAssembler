# ChromAssembler
## Install
Assembler module only depends on `Python 2.7`. 

Suggested packages, which are required for visualization module: [`graphviz`](https://pypi.python.org/pypi/graphviz), [`seaborn`](https://pypi.python.org/pypi/seaborn).

## Usage
`python Assemble.py [-h] --fragments FRAGMENTS [-o OUTPUT]`


## Overview
Assembles a chromosomal sequence from a `FASTA` text file contianing fragments using a de Bruijin Graph assembler. 

Input assumptions: 

1. Fragments overlap (share sequence) with at least one other fragment
2. Sharing region is ≥ ½ the length of each fragment
3. ∃ unique way to reconstruct entire chromsome from input sequences by aligning read 
4. Fragments of length ≤ 1000


