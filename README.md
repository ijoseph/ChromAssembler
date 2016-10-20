# ChromAssembler
## Install
Assembler module only depends on `Python 2.7`. 

Suggested packages, which are required for visualization module: [`graphviz`](https://pypi.python.org/pypi/graphviz), [`seaborn`](https://pypi.python.org/pypi/seaborn).

## Usage
`python Assemble.py [-h] --fragments FRAGMENTS [-o OUTPUT]`


## Overview
Assembles a chromosomal sequence from a `FASTA` text file contianing fragments using a de Bruijin Graph assembler adapted from [teaching materials of Dr. Ben Langmead](http://www.langmead-lab.org/teaching-materials/). 
### Algorithm

1. Make de Bruijn graph by chopping fragments into k-mers
  -k-mer length chosen by  ![choosing k](http://i.imgur.com/xgdp5kv.png)
2. Make Eulerian by:
    - Collapsing multi-edges into one by only considerng each k-mer once
    - Removing transitively-inferable edges
        - Remove edges skipping one or more nodes

3. Do Eulerian walk to get assembly prediction


### Input assumptions: 

1. Fragments overlap (share sequence) with at least one other fragment
2. Sharing region is ≥ ½ the length of each fragment
3. ∃ unique way to reconstruct entire chromsome from input sequences by aligning read 
4. Fragments of length ≤ 1000



## Results
See https://github.com/ijoseph/ChromAssembler/blob/master/Visualizations.ipynb for output and visualizations
