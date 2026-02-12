# pcms
Phylogenetic Covariance Matrix Sparisifcation

## Installation

The package and its dependencies can then be installed as follows:

1. The sparsification algorithm requires the Intel Math Kernel Library (MKL). See installation instructions [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux&linux-install=apt).

1. The permutation test requires the `mvhg` [package](https://github.com/spsvihla/mvhg), which should be installed into `pcms/lib/`. First clone the `pcms` repository

    ```bash
    git clone https://github.com/spsvihla/pcms.git
    ```

    To avoid issues, it is recommended to install this package in a new virtual environment inside `pcms/`.

    ```bash
    cd pcms/
    python -m venv .venv
    ```

    Now install the `mvhg` package in `pcms/lib/`.

    ```bash
    mkdir -p lib/ && cd lib/
    git clone git@github.com:spsvihla/mvhg.git && cd mvhg/
    ./build.sh && ./build.sh -- clean
    ```

1. Finally, to install the `pcms` package itself, run

    ```bash
    ./build.sh && ./build.sh --clean
    ```

    The remaining Python dependencies are included in `pyptoject.toml`. 

1. The Python dependencies installed by default are only those required to run the package itself. If you wish to run notebooks, you must install aditional Python dependencies with the command 

    ```bash
    pip install ".[notebooks]"
    ```

## Datasets

### Greengenes reference database and phylogeny

We make use of the [Greengenes database](https://ftp.microbio.me/greengenes_release/) [1-2], the most up-to-date version of which can be found under `/greengenes_release/current` at the link provided.
As an example, the `gg_13_8` dataset can be downloaded as follows:

```
wget https://ftp.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
tar -xzf gg_13_8_otus.tar.gz
```

### Guerrero Negro microbial mat samples

Other notebooks use the [Guerrero Negro microbial mat datasat](https://www.ncbi.nlm.nih.gov/nuccore/?term=JN427016%3AJN539989%5BAccession%5D) collected and analyzed by Harris, Caporaso, et al. [3].
In particular, the sequences available on GenBank are used to pick OTUs clustered against the Greengenes database as described in [3].
These data can be downloaded at the link provided by selecting `Send to > Gene features` above the search results and creating a file with the "FASTA Nucleotide" format.

**Note:** The above analysis is for the full-length Sanger sequences.
The corresponding table for the 454 partial-length sequences described in the paper may be found on [Qiita](https://qiita.ucsd.edu/study/description/1200#).

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

## Other Information

We have provided a [Mathematica notebook](https://drive.google.com/file/d/1DCYi6A4dRPyst1LvJSTARhLdZLP85fAd/view?usp=sharing) used in the original publication of this work.

## References

1. McDonald, D. et al. Greengenes2 unifies microbial data in a single reference tree. Nat Biotechnol 42, 715–718 (2024). [https://doi.org/10.1038/s41587-023-01845-1](https://doi.org/10.1038/s41587-023-01845-1)
2. McDonald, D. et al. An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea. ISME J 6, 610–618 (2012). [https://doi.org/10.1038/ismej.2011.139](https://doi.org/10.1038/ismej.2011.139)
3. Kirk Harris, J. et al. Phylogenetic stratigraphy in the Guerrero Negro hypersaline microbial mat. The ISME Journal 7, 50–60 (2013). [https://doi.org/10.1038/ismej.2012.79](https://doi.org/10.1038/ismej.2012.79).
4. Svihla, S. and Lladser, M. E. Sparsification of Phylogenetic Covariance Matrices
of Critical Beta-splitting Random Trees. In preparation.
5. Gorman, E. & Lladser, M. E. Sparsification of large ultrametric matrices: insights into the microbial Tree of Life. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 479, 20220847 (2023). [https://doi.org/10.1098/rspa.2022.0847](https://doi.org/10.1098/rspa.2022.0847)
