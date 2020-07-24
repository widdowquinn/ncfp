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
clean_docs:
	@rm -rf docs/_build/html