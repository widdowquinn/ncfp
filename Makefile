# Makefile
#
# This file is part of the ncfp package distribution
# (https://github.com/widdowquinn/ncfp)

# Set up all development dependencies in the current conda environment
setup_env:
	@conda install --file requirements-dev.txt --yes
	@conda install --file requirements.txt --yes
	@pip install -r requirements-pip.txt
	@pre-commit install
	@pip install -U -e .

# Run all tests and display coverage report in a browser
test:
	@pytest --cov-report=html --cov=ncbi_cds_from_protein -v tests/ && open htmlcov/index.html

# Build and display documentation
docs: clean_docs uml
	@cd docs && make html && open _build/html/index.html

uml:
	pyreverse -o pdf -p ncbi_cds_from_protein ncbi_cds_from_protein

# Clean up outputs
clean: clean_docs clean_tests clean_examples

clean_docs:
	@rm -rf docs/_build/html && \
	rm -rf classes_ncbi_cds_From_protein.pdf && \
	rm -rf packages_ncbi_cds_From_protein.pdf

clean_examples:
	@rm -rf tests/test_output/*

clean_tests:
	@rm -rf tests/test_output/*

# Run examples from documentation
examples:
	# NCBI no introns
	@ncfp tests/test_input/sequences/input_ncbi.fasta \
        tests/examples/ncbi dev@null.com -v
	# UniProt no introns
	@ncfp tests/test_input/sequences/input_uniprot.fasta \
        tests/examples/uniprot dev@null.com -v
	# UniProt/Stockholm no introns
	@ncfp -s tests/test_input/sequences/input_uniprot_stockholm.fasta \
        tests/examples/uniprot_stockholm dev@null.com -v
	# Human isoforms/intron-exon
	@ncfp tests/test_input/sequences/human.fasta \
        tests/examples/human dev@null.com -v
	# Logging
	@ncfp tests/test_input/sequences/human.fasta \
        tests/examples/logging dev@null.com \
        -l tests/examples/logging/human.log
	# Cache location
	@ncfp tests/test_input/sequences/human.fasta \
        tests/examples/caches dev@null.com \
        -d tests/examples/caches \
        -c ncfp_cache
	# Cache reuse
	@ncfp tests/test_input/sequences/human.fasta \
        tests/examples/caches1 dev@null.com \
        -d tests/examples/caches \
        -c ncfp_cache
	@ncfp tests/test_input/sequences/human.fasta \
        tests/examples/caches2 dev@null.com \
        -d tests/examples/caches \
        -c ncfp_cache \
        --filestem cached \
        --keepcache