# pcms
Phylogenetic Covariance Matrix Sparisifcation

## Datasets

### Greengenes reference database and phylogeny

We make use of the [Greengenes database](https://ftp.microbio.me/greengenes_release/) [1-2], the most up-to-date version of which can be found under `/greengenes_release/current` at the link provided.
As an example, the `gg_13_8` dataset can be downloaded as follows:

```bash
wget https://ftp.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
tar -xzf gg_13_8_otus.tar.gz
```

### Guerrero Negro microbial mat samples

Some notebooks use the [Guerrero Negro microbial mat datasat](https://www.ncbi.nlm.nih.gov/nuccore/?term=JN427016%3AJN539989%5BAccession%5D) collected and analyzed by Harris, Caporaso, et al. [3].
In particular, the sequences available on GenBank are used to pick OTUs clustered against the Greengenes database as described in [3].
These data can be downloaded with Entrez Direct (see [install instructions](https://www.ncbi.nlm.nih.gov/books/NBK179288/)) by the following command:

The script `gg_13_8_97ref_97clust_JN427016_JN539989.sh` is provided to help with obtaining OTU counts against the 97% Greengenes tree. 
The script has the following dependencies:

1. The environment variable $DATA (see below). The script will output into `$DATA/guerrero_negro` by default.

1. QIIME2 Amplicon [4] (see [install instructions](https://library.qiime2.org/quickstart/amplicon)).

1. A reference set for chimera checking, which may be obtained with the following commands:

    ```bash
    wget https://drive5.com/uchime/gold.fa
    vsearch --derep_fulllength gold.fa \
            --output gold_nodup.fa \
            --sizeout \
            --relabel RDP_
    ```

1. RAxML-NG ([install instructions](https://anaconda.org/bioconda/raxml-ng) and [GitHub](https://github.com/amkozlov/raxml-ng)) for computing model statistics and EPA-NG ([install instructioins](https://anaconda.org/bioconda/epa-ng) and [GitHub](https://github.com/pierrebarbera/epa-ng)) for frament insertion.

### Data directory structure

The notebooks in this respository are written assuming the following directory structure:

```
$DATA/
├── greengenes/
│   ├── gg_13_5_otus/
│   │   └── <contents of gg_13_5.tar.gz>
│   ├── gg_13_8_otus/
│   │   └── <contents of gg_13_8.tar.gz>
│   └── ...
└── greengenes2/
    ├── gg_22_10_otus/
    |   ├── 2022.10.phylogeny.id.nwk
    │   └── ...
    └── ...
```

## References

1. McDonald, D. et al. Greengenes2 unifies microbial data in a single reference tree. Nat Biotechnol 42, 715–718 (2024). [https://doi.org/10.1038/s41587-023-01845-1](https://doi.org/10.1038/s41587-023-01845-1)
2. McDonald, D. et al. An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea. ISME J 6, 610–618 (2012). [https://doi.org/10.1038/ismej.2011.139](https://doi.org/10.1038/ismej.2011.139)
3. Kirk Harris, J. et al. Phylogenetic stratigraphy in the Guerrero Negro hypersaline microbial mat. The ISME Journal 7, 50–60 (2013). [https://doi.org/10.1038/ismej.2012.79](https://doi.org/10.1038/ismej.2012.79).
4. Bolyen E. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9.
5. Svihla, S. and Lladser, M. E. Sparsification of Phylogenetic Covariance Matrices
of Critical Beta-splitting Random Trees. In preparation.
6. Kozlov, A. M., Darriba, D., Flouri, T., Morel, B. & Stamatakis, A. RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics 35, 4453–4455 (2019). [https://doi.org/10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)
7. Barbera, P. et al. EPA-ng: Massively Parallel Evolutionary Placement of Genetic Sequences. Systematic Biology 68, 365–369 (2019). [https://doi.org/10.1093/sysbio/syy054](https://doi.org/10.1093/sysbio/syy054)
