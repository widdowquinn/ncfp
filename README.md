# README.md - `ncbi_cds_from_protein`

This repository contains code for a script that identifies and writes the corresponding nucleotide sequences for each protein in an input multiple sequence file. It uses the NCBI accession or UniProt gene name (as appropriate) to identify source nucleotide sequences in the NCBI databases, and writes them to a file.

## Usage

Example command line:

```bash
ncbi_cds_from_protein tests/test_input/sequences/input_uniprot.fasta tests/test_output/script dev@null --uniprot -c test -v
```

