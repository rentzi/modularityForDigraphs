modularityForDigraphs: A package implementing a modularity measure for directed graphs (graphs) 
=================================
directedModularity.py finds communities or modules in directed networks (digraphs) and evaluates their modularity index. The algorithm which optimizes a modularity function and makes continuous bisections until no further improvement of the modularity function is possible, is from the following paper:
Leicht, Elizabeth A., and Mark EJ Newman. "Community structure in directed networks." Physical review letters 100.11 (2008): 118703.

I have also implemented the undirected version in ( https://github.com/rentzi/netRewireAnalyze ). The paper for this is the following:
Newman, M. E. (2006). Modularity and community structure in networks. Proceedings of the national academy of sciences, 103(23), 8577-8582.

Prerequisites
-------------

- Python 3+
- numpy

Installation
------------
the package ``modularityForDigraphs`` is the folder that will be created once you clone this repository. You can just copy the directedModularity.py to use the functions in it. The Jupyter notebook's purpose is to show examples of how the functions are used. Alternatively to run the functions of the package include the path with the package folder. For example

```
import sys
sys.path.append('the directory path')
import modularityForDigraphs
```


Getting Started
------------


Each of the functions is commented and relatively easy to understand. It is  recommendeded that you go through the jupyter notebook in the repository to understand how to get the modularity values and indices of the clusters. 

Author
------

Ilias Rentzeperis 
