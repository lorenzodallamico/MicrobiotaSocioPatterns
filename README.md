# Social contacts and microbiome sharing (Molde, Norway)

Codes to reproduce the results of *Association Between Oral Microbiota and Close-Range Proximity in a Primary School*.

The study investigates whether face-to-face contact patterns — measured via wearable proximity sensors in a Norwegian school cohort — are associated with the composition of participants' oral microbiota.

> **NOTE**: these codes should be run in `Python 3.9`

---

## Repository structure

```
Analysis.ipynb          # Main notebook: loads the data, runs every analysis, and renders the figures inline
Data/
├── ASV_cleaned.csv          # Microbiota abundance table (multiple time points)
├── table__tax_by_ASV.xlsx   # ASV taxonomy
├── df_net.csv               # Weighted contact network (edge list)
└── molde__nwk_trees.tar-1/  # Rooted phylogenetic tree, used for UniFrac
src/
├── LoadData.py              # Load(): assembles the contact network, taxonomy, tree and microbiota table
└── analysis_functions.py    # Filtering, distance, statistical test and plotting-support functions
```

## Data

| File | Description |
|------|-------------|
| `ASV_cleaned.csv` | Abundance table; rows indexed by (pseudo-anonymized) `Person ID` and `Day code`, remaining columns are ASV abundances |
| `table__tax_by_ASV.xlsx` | ASV taxonomy table |
| `df_net.csv` | Edge list with columns `pid`, `pid2`, `weight` (cumulative contact duration), merged across the two deployments |
| `molde__nwk_trees.tar-1/` | Rooted phylogenetic tree (`6-rooted-tree.nwk`) |

Participant IDs (`Person ID` in `ASV_cleaned.csv`, `pid`/`pid2` in `df_net.csv`) are pseudo-anonymized hexadecimal strings (e.g. `0x168e`) and must match between the two files for the same person.

## Citation

```
@article {dallamico2024association,
	author = {Dall'Amico, Lorenzo and Bai, Xiangning and Weltzien, Sandra Marie and Rayner, Simon and Paolotti, Daniela and Budin Ljosne, Isabelle Sylvie and Matussek, Andreas and Furberg, Anne-Sofie and Cattuto, Ciro and Sivert Nielsen, Christopher},
	title = {Association Between Oral Microbiota and Close-Range Proximity in a Primary School},
	year = {2024},
	doi = {10.1101/2024.12.27.628096},
	journal = {bioRxiv}
}
```
