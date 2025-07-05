# pcms
Phylogenetic Covariance Matrix Sparisifcation

## Datasets

### Greengenes reference database and phylogeny

We make use of the [Greengenes database](https://ftp.microbio.me/greengenes_release/) [1-2], the most up-to-date version of which can be found under `/greengenes_release/current` at the link provided.
As an example, the `gg_13_8` dataset can be downloaded as follows:

```
$ wget https://ftp.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
$ tar -xzf gg_13_8_otus.tar.gz
```

### Guerrero Negro microbial mat samples

Other notebooks use the [Guerrero Negro microbial mat datasat](https://www.ncbi.nlm.nih.gov/nuccore/?term=JN427016%3AJN539989%5BAccession%5D) collected and analyzed by Harris, Caporaso, et al. [3].
In particular, the sequences available on GenBank are used to pick OTUs clustered against the Greengenes database as described in [3].
These data can be downloaded at the link provided by selecting `Send to > Gene features` above the search results and creating a file with the "FASTA Nucleotide" format.

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
4. Svihla, S. and Lladser, M. E. Sparsification of Phylogenetic Covariance Matrices
of Critical Beta-splitting Random Trees. In preparation.
5. Gorman, E. & Lladser, M. E. Sparsification of large ultrametric matrices: insights into the microbial Tree of Life. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 479, 20220847 (2023). [https://doi.org/10.1098/rspa.2022.0847](https://doi.org/10.1098/rspa.2022.0847)
