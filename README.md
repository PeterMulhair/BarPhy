# BarcodePlacer

<img src="https://github.com/PeterMulhair/BarcodePlacer/blob/master/data/barcode_logo.png" width="500" height="250">


This script is required when an expected species ID does not match the result species/genus ID from a sequence similarity search agains the [BOLD](https://www.boldsystems.org/index.php) database. It requires a merged barcode sequence as input, along with an expected species ID (EXID) and a result species/genus ID (REXID). Using the BOLD API it pulls down barcode data for EXID and REXID, creates a multiple sequence alignment, and constructed a phylogenetic tree. The output is a pdf image of the tree, where the identification of the barcode query can be confidently assigned in a phylogenetic context. 

## Requirements

* python3
* Mafft - [install](https://mafft.cbrc.jp/alignment/software/source.html)
* IQTree - [install](http://www.iqtree.org/doc/Quickstart)

## Usage

BarcodePlacer can be run in two ways:

1. With an excel sheet of barcode results as input

  - `python barcode_queries.py --barcode Barcoding_results.xlsx`

2. With a fasta file with a barcode sequence named as EXID and a query species/genus to search against

  - `python barcode_queries.py --query EXID.fasta REXID`

The script also requires certain directories. Both require a dir called `output/`. To run the `--query` version, place your fasta files in a directory called `queries/`

## Examples

Using the fasta files from `queries/` you can ID the barcode sequences using the `--query` version (the `--barcode` version can be run using the excel sheet in `data/`)

First make output dir with `mkdir output`

Then run script
eg. `python barcode_queries.py --query queries/Melangyna_labiatarum.fasta Melangyna_compositarum`
OR
eg. `python barcode_queries.py --query queries/Melangyna_labiatarum.fasta Melangyna_` to search against the REXID genus rather than species.


**NOTE** for `--query` version the input is required as the species bionomial name or genus separated by an underscore eg. Drosophila_melanogaster or Drosophila_
