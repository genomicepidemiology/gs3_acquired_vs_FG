# Geographics and bacterial networks shape the global urban sewage resistome
This contains the code for analyzing the antimicrobial resistance genes and bacteria presented in the paper.

## Overview
### The [ARGprofiler pipeline](https://github.com/genomicepidemiology/ARGprofiler)
The [ARGprofiler](https://github.com/genomicepidemiology/ARGprofiler) pipeline was used to take the raw sequencing reads for each metagenomic sample, trim and quality check them with [fastp](https://github.com/OpenGene/fastp) 0.23.2, and map them against the [PanRes](https://zenodo.org/records/8055116) database and the [mOTUs](https://zenodo.org/records/4923187) database with [KMA](https://bitbucket.org/genomicepidemiology/kma/src) 1.4.12a.

More details on the ARGprofiler pipeline, please refer to the paper:
> Martiny, H. M., Pyrounakis, N., Petersen, T. N., Lukjančenko, O., Aarestrup, F. M., Clausen, P. T., & Munk, P. (2024). ARGprofiler—a pipeline for large-scale analysis of antimicrobial resistance genes and their flanking regions in metagenomic datasets. Bioinformatics, btae086. https://doi.org/10.1093/bioinformatics/btae086

### Software requirements
This repository contains a mix of R and Python code. 

We used USEARCH v.11.0.66 was used to cluster PanRes reference sequences to 98% nucleotide identity:
```{bash}
usearch -cluster_fast panres_genes.fa -id 0.98 -query_cov 0.98 -centroids panres.nr98.fa -uc panres.uc98.tsv
```

Python 3.12 was used to analyze the compositional data in the Jupyter notebook [Diversities_Abundances.ipynb](notebooks/Diversities_Abundances). The following Python packages were used:
```
pycodamath
matplotlib
seaborn
geopandas
matplotlib_venn
scipy
```

R 4.3.2 was used to analyze distance-decays in the notebook [Distance_Decays.ipynb](notebooks/Distance_Decays.ipynb). The following R packages were used:
```
tidyverse
ggpubr
treemapify
vegan
umap
ggrepel
rnaturalearth
rnaturalearthdata
sf
ggpmisc
mgcv
```

In the network analysis notebooks ([1_network_construction.Rmd](notebooks/1_network_construction.Rmd), [2_network_plot.Rmd](notebooks/2_network_plot.Rmd), [3_communities_plot.Rmd](notebooks/3_communities_plot.Rmd), [4_plot_community_human_composition.Rmd](notebooks/4_plot_community_human_composition.Rmd), and (5_make_summary_tables.Rmd)[notebooks/5_make_summary_tables.Rmd]), the following R packages was used:
```
tidyverse
mgnet
ggraph
tidygraph
```
