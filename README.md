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

#### Downloading the dataset

Other notebooks use the [Guerrero Negro microbial mat datasat](https://www.ncbi.nlm.nih.gov/nuccore/?term=JN427016%3AJN539989%5BAccession%5D) collected and analyzed by Harris, Caporaso, et al. [3].
In particular, the sequences available on GenBank are used to pick OTUs clustered against the Greengenes database as described in [3].
These data can be downloaded with Entrez Direct (see [install instructions](https://www.ncbi.nlm.nih.gov/books/NBK179288/)) by the following command:

```bash
export DATASET_DIR=$DATA/guerrero_negro \
&& mkdir -p "$DATASET_DIR" \
&& ./build_accession_list.sh \
    -a "JN427016:JN539989" \
    -o "$DATASET_DIR" \
&& ./fetch_query_seqs.sh \
    -i "$DATASET_DIR/accessions_JN427016_JN539989.txt" \
    -o "$DATASET_DIR/query_seqs_JN427016_JN539989.fasta"
```

#### Installing QIIME2

We matched the sequences to the 97% Greengenes reference phylogeny using QIIME2 [4]. See the [quickstart guide](https://library.qiime2.org/quickstart) for QIIME2 Amplicon for installation instructions.

#### Obtaining OTU counts

To cluster your query sequences against a reference database and generate OTU counts, run the following commands:

```bash
export REF_SEQS=$DATA/greengenes/gg_13_8_otus/rep_set/97_otus.fasta \
&& ./import_query_seqs.sh \
    -f "$DATASET_DIR/query_seqs_JN427016_JN539989.fasta" \
    -o "$DATASET_DIR/query_seqs_JN427016_JN539989.qza" \
&& ./cluster_query_seqs.sh \
    -r "$REF_SEQS" \
    -q "$DATASET_DIR/query_seqs_JN427016_JN539989.qza" \
    -o "$DATASET_DIR"
```

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