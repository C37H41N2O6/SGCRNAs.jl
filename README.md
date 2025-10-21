# SGCRNAs [![Build Status](https://github.com/C37H41N2O6/SGCRNAs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/C37H41N2O6/SGCRNAs.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage Status](https://coveralls.io/repos/github/C37H41N2O6/SGCRNAs.jl/badge.svg?branch=main)](https://coveralls.io/github/C37H41N2O6/SGCRNAs.jl?branch=main) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://C37H41N2O6.github.io/SGCRNAs.jl/)

## About SGCRNAs
[Excerpt from the paper's abstract]  
Weighted Gene Co-expression Network Analysis (WGCNA) is among the most widely employed methods in bioinformatics. WGCNA enables the identification of gene clusters (modules) exhibiting correlated expression patterns, the association of these modules with traits, and the exploration of candidate biomarker genes by focusing on hub genes within the modules. WGCNA has been successfully applied in diverse biological contexts. However, conventional algorithms manifest three principal limitations: the assumption of scale-free topology, the requirement for parameter tuning, and the neglect of regression line slopes. These limitations are addressed by SGCRNA (Spectral-clustering Generalised Correlation Regression Network Analysis).

## Installation
Install the package using Julia's package manager:
```julia
using Pkg; Pkg.add("SGCRNAs")
```

## Usage
Here's a basic example:
```julia
using CSV, DataFrames
using SGCRNAs

# Load Data
Data = CSV.read("test_data/normalizedCounts_coding.tsv", header=1, comment="#", delim='\t', DataFrame);

# Pre-proccessing
CorData, GradData = CGM(Data.Symbol, Matrix(Data[:,5:end]), mode=:FTEST);

# Clustering
clust, pos, edge_data = SpectralClustering(CorData, GradData);

# Draw Network
d = 1; k = maximum(clust[d]);
nw, new_pos, cnctdf, new_clust, score = SetNetwork(edge_data, clust[d], pos, il=collect(1:k));
DrawNetwork("Result/Fig/AllNetWork-0.5.png", nw, new_pos, cnctdf, new_clust, k, node_scores=score, edge_mode=:ALL, edge_threshold=0.5);

# module-phenomenon correlation
Phen = CSV.read("SraRunTable.csv", header=1, comment="#", delim=',', DataFrame);
sort!(Phen);
Data = innerjoin(Data, DataFrame(Symbol=names(edge_data)), on=:Symbol, order=:right);
CorPhenMod(Data[:,5:end], Phen[:,[2,3,5]], new_clust, "Result/Fig/CorPhenMod.png");
```

## Documentation
Full documentation is available at: https://C37H41N2O6.github.io/SGCRNAs.jl/

## License
This project is licensed under the MIT License. See LICENSE for details.

## Related Publication
If you use this package in your research, please cite the following paper:  
T. Osone, et al., SGCRNA: Spectral Clustering-Guided Co-Expression Network Analysis Without Scale-Free Constraints for Multi-Omic Data, bioRxiv, 2025.
https://www.biorxiv.org/cgi/content/short/2025.04.27.650628v1

