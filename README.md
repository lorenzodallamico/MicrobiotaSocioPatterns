# Association Between Oral Microbiota and Close-Range Proximity in a Primary School

In this repository, we share the material for reproducing the results presented in *Dall'Amico et al, Association Between Oral Microbiota and Close-Range Proximity in a Primary School*. In this work, we investigated the association between close-range proximity of children in a primary school (with [SocioPatterns](http://www.sociopatterns.org/) proximity sensors) and the composition of their oral microbiota. The content is organized as follows:
* `Data/ContactNetwork.csv`. This is the aggregated contact network among the students. It features three columns: `pid`, `pid2`, and `weight`. The first two (`pid` and `pid2`) contain random identifiers of the children composed of $4$ letters. The third column contains the cumulative interaction time measured between `pid` and `pid2`. 
* `Data/microbiota.csv`: this file is indexed by the $4$-digit identifier of the children (that can be mapped to the proximity network) and its columns are the microbiota detected by the 16S analysis. Each entry is either $0$ or $1$ depending on whether or not that taxon was measured in any of the measurements.
* `Analysis.ipynb`: this file contains the Python codes to run our experiments.
* `src/utils.py`: this file contains some functions used in the `Analysis` notebook.



If you use these data please reference the following article: [https://www.biorxiv.org/content/10.1101/2024.12.27.628096v1](https://www.biorxiv.org/content/10.1101/2024.12.27.628096v1)

```@cite
@article {dallamico2024association,
 author = {Dall'Amico, Lorenzo and Bai, Xiangning and Weltzien, Sandra Marie and Rayner, Simon and Paolotti, Daniela and Budin Ljosne, Isabelle Sylvie and Matussek, Andreas and Furberg, Anne-Sofie and Cattuto, Ciro and Sivert Nielsen, Christopher},
 title = {Association Between Oral Microbiota and Close-Range Proximity in a Primary School},
 year = {2024},
 doi = {10.1101/2024.12.27.628096},
 journal = {bioRxiv}}
```
