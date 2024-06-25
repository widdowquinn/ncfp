# CHANGES.md - `ncfp`

## v0.2.1a1

- add ruff.toml configuration
- move from CircleCI to GitHub Actions for CI
- if all other attempts to match the linked CDS fail, and there is only one CDS in the retrieved record, use that CDS as the candidate sequence (#46)
- fix badges (Zenodo/PyPI) in `README.md`
- make Stockholm domains work with NCBI sequence inputs (#27)
- update requirements to avoid deprecated `requests-cache` (#26)
- add `--use-protein-ids` option to use the same protein sequence seqID for the protein and the downloaded nucleotide sequence (#21, Emma Hobbs' first contribution)
- updated code to use `pathlib` internally instead of `os.path` operations
- fixed loglevel bug with logging to file (#25)
- update progress bar labels with step numbering (#24)
- update UniProt queries to use the new API (#35)
- add `--allow_alternative_start_codon` option (#34)
- warn, rather than halt, when a suspected IPG input sequence is found (#34)
- add extra information when any other than one coding sequence is found for an input protein (#34)
- add last ditch attempt to identify proteins in an NCBI nucleotide file using the gene name
- use the UniProt accession number if the GN field is not specific enough
- when UniProt records point to NCBI gene records, obtain the nucleotide accession from the gene record, and protein ID from the RefSeq cross-reference
- use the UniProt ORF gene name field for identifying CDS in larger records, to avoid issues with ambiguous GN fields
- adapt to bioservices 1.12 change that modifies UniProt return string format

## v0.2.0

- bump version and add tag for release

## v0.2.0-a1

- Update copyright notices
- Convert parsers to use `pathlib`
- Add requirements file for development tools
- Convert tests to use `pytest`, not `nose`
- Revise logging usage
- Change continuous integration from TravisCI to CircleCI
- Guess sequence origin rather than asking user to provide at CLI (allows mixed origin files)
- Update CLI parser and docs to reflect new sequence origin guessing

## v0.1.1

- Add installation instructions and other improvements to documentation
- Tidied codebase in some places (removing `print` statements, unused functions, etc.).
- Add CLI tests.
- Correct `ncfp` program name in help/usage text.

## v0.1.0

*- First release of `ncfp`
