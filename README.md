# BarPhy

**Bar**code **Phy**logenetics

<div align="center">
<p align="center">
<img src="https://github.com/PeterMulhair/BarcodePlacer/blob/master/data/barcode_logo.png" width="500" height="250">
</p>
</div>

This script is required when an expected species ID does not match the result species/genus ID from a sequence similarity search against the [BOLD](https://www.boldsystems.org/index.php) database. It requires a merged barcode sequence as input, along with an expected species ID (**EXID**) and a result species/genus ID (**REXID**). Using the BOLD API it pulls down barcode data for EXID and REXID, creates a multiple sequence alignment, and constructed a phylogenetic tree. The output is a pdf image of the tree, where the identification of the barcode query can be confidently assigned in a phylogenetic context. 

## Requirements

* python3
* Mafft - [install](https://mafft.cbrc.jp/alignment/software/source.html)
* IQTree - [install](http://www.iqtree.org/doc/Quickstart)
* R v3 - [install](https://cran.r-project.org/doc/manuals/r-release/R-admin.html)
* Python modules; `pandas, os, glob, argparse, subprocess, ete3, Bio, joblib`. If you get an error for any of these, install using `pip` before running script.
* R libraries; `ggtree, getopt`

Clone repository locally using `git clone https://github.com/PeterMulhair/BarPhy.git`

## Usage

BarPhy can be run in two ways:

1. With an excel sheet or csv file of barcode results as input (see the excel sheet in [data/](https://github.com/PeterMulhair/BarcodePlacer/tree/main/data))

  - `python barcode_queries.py --barcode Barcoding_results.xlsx`

2. With a fasta file containing a barcode sequence named as EXID and a query species/genus to search against

  - `python barcode_queries.py --query EXID.fasta REXID`

The script also requires certain directories. To run the `--query` version, place your fasta files in a directory called `queries/`

**Output**

The output folder consists of a number of files including raw fasta, MSA, and tree output files.
The tree image file, ending in .pdf, is what you want to check to see where your barcode query fits in the tree. 

## Examples

Using the fasta file from `queries/` you can ID the barcode sequences using the `--query` version (the `--barcode` version can be run using the excel sheet in `data/`)


Then run the script:

```

$ python barcode_queries.py --query queries/Melangyna_labiatarum.fasta Melangyna_compositarum

OR to search against the REXID genus rather than species:

$ python barcode_queries.py --query queries/Melangyna_labiatarum.fasta Melangyna_

```

**NOTE** for `--query` version the input is required as the species bionomial name, or genus, separated by an underscore eg. Drosophila_melanogaster or Drosophila_

**Output**

For `--query` runs, the query species will be labelled with `_query` in the output tree image (labelled with `_DToL` in `--barcode` version). 

<div align="center">
<p align="center">
<img src="https://github.com/PeterMulhair/BarcodePlacer/blob/master/data/example_tree.png" width="700" height="550">
</p>
</div>


**To do**

- Add R error output line
- Label query tip in tree output image
