SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c

.PHONY: check figures paper

check:
	python3 -u code/exact_cones.py 1 2 3 4 --verify | tee paper/generated/exact_check.txt

figures:
	python3 -u code/generate_plot_data.py

paper: figures
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex
	cd paper && bibtex main
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex
