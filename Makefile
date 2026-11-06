PYTHON ?= .venv/bin/python
PERMUTATIONS ?= 9999
SEED ?= 70899
EXAMPLE_SEED ?=
SCENARIO ?= blocked_randomization
export OPENBLAS_NUM_THREADS = 4
export OMP_NUM_THREADS = 4

.PHONY: fetch prepare corrections randomization limma permutation count-test change-test shuffled-example plots test audit audit-deposited all notebook execute-notebook
fetch:
	$(PYTHON) scripts/fetch_data.py
prepare:
	$(PYTHON) scripts/prepare_data.py
	$(PYTHON) scripts/prepare_corrections.py
corrections:
	$(PYTHON) scripts/prepare_corrections.py
randomization:
	$(PYTHON) scripts/prepare_randomization.py
limma:
	bash scripts/run_r.sh scripts/figure1b_limma.R diet_code '' $(SCENARIO)
ifneq ($(SCENARIO),blocked_randomization)
	bash scripts/run_r.sh scripts/figure1b_limma.R diet_named '' $(SCENARIO)
endif
permutation:
	$(PYTHON) scripts/permutation.py --scenario $(SCENARIO) --permutations $(PERMUTATIONS) --seed $(SEED)
count-test:
	$(PYTHON) scripts/count_test.py --scenario $(SCENARIO)
change-test:
	$(PYTHON) scripts/change_test.py --scenario $(SCENARIO)
shuffled-example:
	$(PYTHON) scripts/shuffled_example.py --scenario $(SCENARIO) $(if $(EXAMPLE_SEED),--seed $(EXAMPLE_SEED))
plots: count-test change-test shuffled-example
	$(PYTHON) scripts/plot_results.py --scenario $(SCENARIO)
test:
	$(PYTHON) -m unittest discover -s tests -v
audit:
	$(PYTHON) scripts/audit_blocks.py
	$(PYTHON) scripts/audit_invariants.py --scenario blocked_randomization
	$(PYTHON) scripts/audit_results.py --scenario blocked_randomization --replay-only
audit-deposited:
	$(PYTHON) scripts/audit_results.py --replay
	$(PYTHON) scripts/audit_invariants.py --scenario deposited
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
