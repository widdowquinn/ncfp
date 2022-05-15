# CHANGES.md - `ncfp`

## v0.2.1a1

- fix badges (Zenodo/PyPI) in `README.md`
- make Stockholm domains work with NCBI sequence inputs (#27)
- update requirements to avoid deprecated `requests-cache` (#26)

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