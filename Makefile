PYTHON ?= .venv/bin/python
PERMUTATIONS ?= 9999
SEED ?= 70899
EXAMPLE_SEED ?=
export OPENBLAS_NUM_THREADS = 4
export OMP_NUM_THREADS = 4

.PHONY: fetch prepare limma permutation count-test change-test shuffled-example plots test all notebook execute-notebook
fetch:
	$(PYTHON) scripts/fetch_data.py
prepare:
	$(PYTHON) scripts/prepare_data.py
limma:
	bash scripts/run_r.sh scripts/figure1b_limma.R diet_code
	bash scripts/run_r.sh scripts/figure1b_limma.R diet_named
permutation:
	$(PYTHON) scripts/permutation.py --permutations $(PERMUTATIONS) --seed $(SEED)
count-test:
	$(PYTHON) scripts/count_test.py
change-test:
	$(PYTHON) scripts/change_test.py
shuffled-example:
	$(PYTHON) scripts/shuffled_example.py $(if $(EXAMPLE_SEED),--seed $(EXAMPLE_SEED))
plots: count-test change-test shuffled-example
	$(PYTHON) scripts/plot_results.py
test:
	$(PYTHON) -m unittest discover -s tests -v
notebook:
	$(PYTHON) -m jupyter lab notebooks/figure1b_results.ipynb
execute-notebook:
	$(PYTHON) scripts/execute_notebook.py
all:
	$(MAKE) fetch
	$(MAKE) prepare
	$(MAKE) test
	$(MAKE) limma
	$(MAKE) permutation
	$(MAKE) plots
