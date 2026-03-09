# Microbiota and Social Contacts

This repository contains the code to reproduce the analyses presented in the paper
*Association Between Oral Microbiota and Close-Range Proximity in a Primary School*.

The study investigates whether face-to-face contact patterns — measured via
wearable proximity sensors in a Norwegian school cohort — are associated with
the composition of participants' gut microbiota.

---

## Repository structure

```

Analysis.ipynb          # Main analysis notebook
Data/
├── ContactNetwork.csv      # Weighted contact network (edge list)
├── microbiota_bool.csv     # Binary taxon presence/absence per individual
└── microbiota.csv          # Raw taxon abundance table (multiple time points)
src/
├── distances.py            # Microbiota distance/similarity functions
└── analysis_functions.py   # Higher-level analysis functions
```

## Data

| File | Description |
|------|-------------|
| `ContactNetwork.csv` | Edge list with columns `pid`, `pid2`, `weight` (cumulative contact duration in units of 10 s) |
| `microbiota_bool.csv` | Binary presence/absence table; rows = individuals, columns = semicolon-delimited ASV taxonomy strings |
| `microbiota.csv` | Full abundance table; first four columns are metadata (`Person ID`, `Sampling-date`, `Day code`, one more); remaining columns are ASV abundances |


## Citation

````
@article {dallamico2024association,
	author = {Dall'Amico, Lorenzo and Bai, Xiangning and Weltzien, Sandra Marie and Rayner, Simon and Paolotti, Daniela and Budin Ljosne, Isabelle Sylvie and Matussek, Andreas and Furberg, Anne-Sofie and Cattuto, Ciro and Sivert Nielsen, Christopher},
	title = {Association Between Oral Microbiota and Close-Range Proximity in a Primary School},
	year = {2024},
	doi = {10.1101/2024.12.27.628096},
	journal = {bioRxiv}
}
````