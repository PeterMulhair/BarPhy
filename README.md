# BarcodePlacer
This script is required when an expected species ID does not match the result species/genus ID from a sequence similarity search agains the [BOLD](https://www.boldsystems.org/index.php) database. It requires a merged barcode sequence as input, along with an expected species ID (EXID) and a result species/genus ID (REXID). Using the BOLD API it pulls down barcode data for EXID and REXID, creates a multiple sequence alignment, and constructed a phylogenetic tree. The output is a pdf image of the tree, where the identification of the barcode query can be confidently assigned in a phylogenetic context. 

## Requirements

* python3
* BLAST+ suite - [install for command line](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
* MMseqs - [install](https://github.com/soedinglab/MMseqs2#installation)

## Usage

BarcodePlacer can be run in two ways:

    1. With an excel sheet of barcode results as input

    `python barcode_queries.py --barcode Barcoding_results.xlsx`

    2. With a fasta file with a barcode sequence named as EXID and a query species/genus to search against

    `python barcode_queries.py --query EXID.fasta REXID`