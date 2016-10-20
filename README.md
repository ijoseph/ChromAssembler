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
    - k-mer length chosen by <img src="http://i.imgur.com/xgdp5kv.png" width="200">
2. Make Eulerian by:
    - Collapsing multi-edges into one by only considerng each k-mer once


3. Find Eulerian path and do Eulerian walk to get assembly prediction

#### Possible Issues
May be impossible to make an Eulerian graph, at which point the assembler complains. This would generally be true if there are errors in reads, which I believe is prevented by the assumption at-least-half overalap (see below) and the setting of `k`. If I had more time, I would attempt a more formal proof of the above, and if I failed to do so, I would: 
- Removing transitively-inferable edges in the de Bruijn Graph (edges skipping one or more nodes) (if these were the cause of non-Eulerian nature)

        

### Input assumptions: 

1. Fragments overlap (share sequence) with at least one other fragment
2. Sharing region is ≥ ½ the length of each fragment
3. ∃ unique way to reconstruct entire chromsome from input sequences by aligning reads
4. Fragments of length ≤ 1000



## Results
See https://github.com/ijoseph/ChromAssembler/blob/master/Visualizations.ipynb for output and visualizations
