<br>

<p align="center">
<strong>HiCDOC has moved to: <a href="https://github.com/mzytnicki/HiCDOC">mzytnicki/HiCDOC</a></strong>
</p>

<br><hr><br>

HiCDOC normalizes intrachromosomal Hi-C matrices and detects A/B compartments
with multiple replicates using unsupervised learning. This repository is a proof of concept.

<br>

## Prerequisites

The following dependencies are required:

- R >= 3.5
- Python >= 3.6
- MultiHiCcompare<sup>[[publication][multihiccompare-publication]][[installation][multihiccompare-installation]]</sup>

    ```bash
    R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
    R -e 'BiocManager::install("multiHiCcompare")'
    ```

- gcMapExplorer<sup>[[publication][gcmapexplorer-publication]][[installation][gcmapexplorer-installation]]</sup>

    ```bash
    pip3 install Cython gcMapExplorer
    ```

    On macOS, you might need this preliminary (with homebrew):
    ```bash
    brew install cmake gcc
    export CC=/usr/local/bin/gcc-9
    ```

- argparse, numpy, sklearn, statsmodels, matplotlib, plotly

    ```bash
    R -e 'install.packages("argparse")'
    pip3 install argparse numpy sklearn statsmodels matplotlib plotly
    ```
<br>

## Usage

HiCDOC is segmented in a collection of scripts, which can be used to construct a
custom pipeline.

The default pipeline can be run with `hicdoc.sh`. It consists of 5 steps:
1. Normalize technical biases with `normalize_cyclic_loess.r`
2. Normalize biological biases with `normalize_knight_ruiz.py`
3. Normalize distance effect with `normalize_distance_rnr_combined.py`
4. Detect compartments and compute measures with `detect_constrained_k_means.py`
5. Plot compartment changes and measures with `plot_compartment_changes.py`

<br>

### Input format

Each script accepts a tab-separated multi-replicate sparse matrix with a header
line and optional comment lines.

    # Optional comment lines
    # ...
    chromosome    position 1    position 2    replicate 1.1    replicate 1.2    replicate 2.1    ...
    3             1500000       7500000       145              184              72               ...
    ...

The interaction proportions between `position 1` and `position 2` are reported
in each replicate column, named `replicate <condition.replicate>`. There is no
limit to the number of replicates and conditions.

<br>

### Scripts and arguments

###### `hicdoc.sh`

    ./hicdoc.sh
      -i <file>                                  Input matrix file
      -d <directory>                             Output directory

Run the default pipeline on the input matrix. The matrix interactions will be
normalized by `normalize_cyclic_loess.r`, `normalize_knight_ruiz.py` and
`normalize_distance_rnr_combined.py`, then A/B compartments will be detected
with `detect_constrained_k_means.py` and plotted alongside various measures with
`plot_compartment_changes.py`.

<br>

###### `join_replicates.py`

    ./join_replicates.py
      -i <file> ...                              Input matrix files
      -o <file>                                  Output matrix file
      --replicates <condition.replicate> ...     Replicate column names
                                                 in the same order as the input files
      [--inputs-have-headers]                    Add if the input matrices have a header line
      [--comments "<comment line>" ...]          Comment lines to add to the top of the output file

Join single-replicate matrices (chromosome, position 1, position 2, interaction)
into one multi-replicate sparse matrix.

<br>

###### `normalize_cyclic_loess.r`

    ./normalize_cyclic_loess.r
      -i <file>                                  Input matrix file
      -o <file>                                  Output matrix file

Normalize technical biases (sequencing depth, restriction enzyme, etc.) with a
cyclic loess<sup>[[publication][multihiccompare-publication]][[implementation][cyclic-loess-implementation]]</sup>.

The cyclic loess method constructs a MD plot (difference ~ genomic distance) for
each replicate pair, then recursively corrects their interaction proportions,
such that the mean difference between each replicate pair reaches zero at each
genomic distance.

Low-proportions interaction vectors are NOT filtered before normalization.

<br>

###### `normalize_knight_ruiz.py`

    ./normalize_knight_ruiz.py
      -i <file>                                  Input matrix file
      -o <file>                                  Output matrix file

Normalize biological biases (GC content, repeated sequences, etc.) with the
Knight-Ruiz algorithm<sup>[[publication][knight-ruiz-publication]][[implementation][knight-ruiz-implementation]]</sup>.

The Knight-Ruiz matrix balancing algorithm normalizes coverage by transforming
each replicate matrix into a doubly stochastic matrix.

A filter is applied before normalization, removing low-proportions interaction
vectors whose number of zeros exceeds the 99th percentile of the distribution of
zeros per interaction vector.

<br>

###### `normalize_distance_rnr_combined.py`

    ./normalize_distance_rnr_combined.py
      -i <file>                                  Input matrix file
      -o <file>                                  Output matrix file
      [--expected <file>]                        Output "expected" interaction proportions
                                                 at each genomic distance

Normalize distance effect (linear proximity affecting interaction proportions)
with a combined radius-neighbors regression<sup>[[implementation][rnr-implementation]]</sup>.

First, an interactions ~ distance plot is constructed from all the combined
replicates. Then a radius-neighbors regression is applied to estimate the
"expected" interaction proportion at each genomic distance. Finally, each
"observed" interaction proportion is divided by its "expected" value for its
genomic distance.

A filter is applied to ignore empty interaction vectors before estimating
"expected" interactions.

<br>

###### `normalize_distance_rnr_individual.py`

    ./normalize_distance_rnr_individual.py
      -i <file>                                  Input matrix file
      -o <file>                                  Output matrix file
      [--expected <file>]                        Output "expected" interaction proportions
                                                 at each genomic distance

Normalize distance effect with an individual radius-neighbors
regression<sup>[[implementation][rnr-implementation]]</sup> for each replicate.

First, an interactions ~ distance plot is constructed for each individual
replicate. Then a radius-neighbors regression is applied to estimate the
"expected" interaction proportion at each genomic distance for each replicate.
Finally, each "observed" interaction proportion is divided by its "expected"
value for its genomic distance.

A filter is applied to ignore empty interaction vectors before estimating
"expected" interactions.

<br>

###### `normalize_distance_mean_individual.py`

    ./normalize_distance_mean_individual.py
      -i <file>                                  Input matrix file
      -o <file>                                  Output matrix file
      [--expected <file>]                        Output "expected" interaction proportions
                                                 at each genomic distance

Normalize distance effect with an individual interaction mean estimation<sup>[[implementation][interaction-mean-implementation]]</sup>
for each replicate.

First, an interactions ~ distance plot is constructed for each individual
replicate. Then the mean of interactions is computed to estimate the "expected"
interaction proportion at each genomic distance for each replicate. Finally,
each "observed" interaction proportion is divided by its "expected" value for
its genomic distance.

A filter is applied to ignore empty interaction vectors before estimating
"expected" interactions.

<br>

###### `normalize_vectors_min_max.py`

    ./normalize_vectors_min_max.py
      -i <file>                                  Input matrix file
      -o <file>                                  Output matrix file

Min-max scale each interaction vector to [0, 1].

<br>

###### `detect_constrained_k_means.py`

    ./detect_constrained_k_means.py
      -i <file>                                  Input matrix file
      -o <file>                                  Output compartments file
      [-k <n>]                                   Number of compartments to detect
                                                 Default: 2
      [--distances <file> ...]                   Output distances to centroids
                                                 One file per compartment
      [--concordance <file>]                     Output concordance file
      [--silhouette <file>]                      Output Silhouette coefficient file

Detect compartments using constrained k-means<sup>[[publication][constrained-k-means-publication]][[implementation][constrained-k-means-implementation]]</sup>.
The algorithm applies a compartment label to each genomic position based on
interaction vectors differences, with the additional constraint that, for each
given condition, different replicates of the same genomic position must belong to
the same compartment.

The original algorithm operates somewhat naively. At each iteration, it selects
a compartment for the first replicate it encounters, then forces that
compartment onto subsequent replicates. For HiCDOC, the implementation has been
modified. At each iteration, the algorithm selects a compartment that fits best
for the majority of replicates.

To determine the strength of membership of each genomic position replicate to its
compartment, a Silhouette
coefficient<sup>[[implementation][silhouette-implementation]]</sup> is computed,
as well as a concordance value. The concordance is a confidence measure designed
specifically for HiCDOC. Each genomic position replicate is given a concordance
between -1 and 1. The value is positive or negative depending on which
compartment fits best, and its distance from 0 indicates the strength of
membership. Concordance is the ratio of distance to each centroid, normalized by
the distance between the two centroids, scaled to the [-1, 1] interval:

<p align="center"><img src="https://user-images.githubusercontent.com/7478535/59969298-65d3ab00-954a-11e9-8b0f-30a0139ab08a.png"/></p>

<br>

###### `plot_matrix.py`

    ./plot_matrix.py
      -i <file>                                  Input matrix file
      -p <prefix>                                Output figure prefix
      [--measure <file>]                         Input measure file
      [--name <name>]                            Name of the measure to write on the figure

Plot a matrix, with an optional measure (concordance, silhouette or distance).
One figure will be created per chromosome and replicate. Each figure is saved to
`prefix_resolution_chromosome_replicate.png`.

<br>

###### `plot_ma.py`

    ./plot_ma.py
      -i <raw file> <normalized file>            Input raw and normalized matrix file
      -p <prefix>                                Output figure prefix

Create MA plots (difference ~ average) for each pair of replicates. One figure
will be created per chromosome. Each figure is saved to
`prefix_resolution_chromosome.png`.

<br>

###### `plot_expected.py`

    ./plot_expected.py
      -i <file>                                  Input "expected" interaction proportions file
      -p <prefix>                                Output figure prefix

Create an interactions ~ distance plot with the "expected" interaction
proportions. One figure will be created per chromosome. Each figure is saved to
`prefix_resolution_chromosome.png`.

<br>

###### `plot_compartment_changes.py`

    ./plot_compartment_changes.py
      -i <file>                                  Input compartments file
      -p <prefix>                                Output figure prefix
      [--distances <file> ...]                   Input distances to centroids
                                                 One file per compartment
      [--concordance <file>]                     Input concordance file
      [--silhouette <file>]                      Input Silhouette coefficient file

Plot compartment changes. One figure will be created per chromosome and
condition pair. Each figure is saved to
`prefix_resolution_chromosome_conditions.png`.

<br>

###### `plot_concordance_changes.py`

    ./plot_concordance_changes.py
      -i <compartments file> <concordance file>  Input compartments and concordance files
      -p <prefix>                                Output figure prefix

Plot concordance changes. One figure will be created per chromosome and
condition pair. Each figure is saved to
`prefix_resolution_chromosome_conditions.png`.

<br>

## References

Philip A. Knight, Daniel Ruiz, A fast algorithm for matrix balancing, _IMA
Journal of Numerical Analysis_, Volume 33, Issue 3, July 2013, Pages 1029–1047,
https://doi.org/10.1093/imanum/drs019

Rajendra Kumar, Haitham Sobhy, Per Stenberg, Ludvig Lizana, Genome contact map
explorer: a platform for the comparison, interactive visualization and analysis
of genome contact maps, _Nucleic Acids Research_, Volume 45, Issue 17, 29
September 2017, Page e152, https://doi.org/10.1093/nar/gkx644

John C Stansfield, Kellen G Cresswell, Mikhail G Dozmorov, multiHiCcompare:
joint normalization and comparative analysis of complex Hi-C experiments,
_Bioinformatics_, 2019, https://doi.org/10.1093/bioinformatics/btz048

Kiri Wagstaff, Claire Cardie, Seth Rogers, Stefan Schrödl, Constrained K-means
Clustering with Background Knowledge, _Proceedings of 18th International
Conference on Machine Learning_, 2001, Pages 577-584,
https://pdfs.semanticscholar.org/0bac/ca0993a3f51649a6bb8dbb093fc8d8481ad4.pdf

[multihiccompare-publication]: https://doi.org/10.1093/bioinformatics/btz048
[multihiccompare-installation]: https://bioconductor.org/packages/release/bioc/html/multiHiCcompare.html
[gcmapexplorer-publication]: https://doi.org/10.1093/nar/gkx644
[gcmapexplorer-installation]: https://gcmapexplorer.readthedocs.io/en/latest/install.html
[orca-installation]: https://github.com/plotly/orca#installation
[cyclic-loess-implementation]: https://bioconductor.org/packages/release/bioc/vignettes/multiHiCcompare/inst/doc/multiHiCcompare.html#cyclic-loess-normalization
[knight-ruiz-publication]: https://doi.org/10.1093/imanum/drs019
[knight-ruiz-implementation]: https://gcmapexplorer.readthedocs.io/en/latest/commands/normKR.html
[rnr-implementation]: https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.RadiusNeighborsRegressor.html
[interaction-mean-implementation]: https://gcmapexplorer.readthedocs.io/en/latest/commands/normMCFS.html
[constrained-k-means-publication]: https://pdfs.semanticscholar.org/0bac/ca0993a3f51649a6bb8dbb093fc8d8481ad4.pdf
[constrained-k-means-implementation]: https://github.com/Behrouz-Babaki/COP-Kmeans
[silhouette-implementation]: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_samples.html
