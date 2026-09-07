.PHONY: all tests figures paper docs guides clean

PYTHON := .venv/bin/python
PYTEST := .venv/bin/pytest

all: tests figures paper docs guides

tests:
	$(PYTEST) tests

figures:
	$(PYTHON) scripts/reproduce_all.py
	$(PYTHON) scripts/analyze_convergence.py

paper: docs/scientific_paper/paper_reimplementation.pdf
	cp docs/scientific_paper/paper_reimplementation.pdf paper_reimplementation.pdf

docs: docs/software_documentation/software_documentation.pdf
	cp docs/software_documentation/software_documentation.pdf software_documentation.pdf

guides: docs/guia_conceptual/guia_conceptual.pdf docs/solucion_numerica/solucion_numerica.pdf docs/arquitectura/arquitectura.pdf
	cp docs/guia_conceptual/guia_conceptual.pdf guia_conceptual.pdf
	cp docs/solucion_numerica/solucion_numerica.pdf solucion_numerica.pdf
	cp docs/arquitectura/arquitectura.pdf arquitectura.pdf

docs/scientific_paper/paper_reimplementation.pdf: docs/scientific_paper/paper_reimplementation.tex
	cd docs/scientific_paper && pdflatex -interaction=nonstopmode paper_reimplementation.tex
	cd docs/scientific_paper && pdflatex -interaction=nonstopmode paper_reimplementation.tex

docs/software_documentation/software_documentation.pdf: docs/software_documentation/software_documentation.tex
	cd docs/software_documentation && pdflatex -interaction=nonstopmode software_documentation.tex
	cd docs/software_documentation && pdflatex -interaction=nonstopmode software_documentation.tex

docs/guia_conceptual/guia_conceptual.pdf: docs/guia_conceptual/guia_conceptual.tex
	cd docs/guia_conceptual && pdflatex -interaction=nonstopmode guia_conceptual.tex
	cd docs/guia_conceptual && pdflatex -interaction=nonstopmode guia_conceptual.tex

docs/solucion_numerica/solucion_numerica.pdf: docs/solucion_numerica/solucion_numerica.tex
	cd docs/solucion_numerica && pdflatex -interaction=nonstopmode solucion_numerica.tex
	cd docs/solucion_numerica && pdflatex -interaction=nonstopmode solucion_numerica.tex

docs/arquitectura/arquitectura.pdf: docs/arquitectura/arquitectura.tex
	cd docs/arquitectura && pdflatex -interaction=nonstopmode arquitectura.tex
	cd docs/arquitectura && pdflatex -interaction=nonstopmode arquitectura.tex

clean:
	rm -f docs/scientific_paper/*.aux docs/scientific_paper/*.log docs/scientific_paper/*.out docs/scientific_paper/*.toc
	rm -f docs/software_documentation/*.aux docs/software_documentation/*.log docs/software_documentation/*.out docs/software_documentation/*.toc
	rm -f docs/guia_conceptual/*.aux docs/guia_conceptual/*.log docs/guia_conceptual/*.out docs/guia_conceptual/*.toc
	rm -f docs/solucion_numerica/*.aux docs/solucion_numerica/*.log docs/solucion_numerica/*.out docs/solucion_numerica/*.toc
	rm -f docs/arquitectura/*.aux docs/arquitectura/*.log docs/arquitectura/*.out docs/arquitectura/*.toc
