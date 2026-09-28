.PHONY: all reproduce tests lint results tmm pwe comparison figures documents paper docs guides pwe-docs \
        verify presentation presentation-pwe clean

PYTHON := .venv/bin/python
export MPLCONFIGDIR ?= /tmp/mpl

# Documentos LaTeX (docs/<carpeta>/<nombre>.tex).
PWE_DOCS   := docs/metodo_ondas_planas/metodo_ondas_planas.pdf docs/comparacion_tmm_pwe/comparacion_tmm_pwe.pdf
GUIDES     := docs/guia_conceptual/guia_conceptual.pdf docs/solucion_numerica/solucion_numerica.pdf \
              docs/arquitectura/arquitectura.pdf
PAPER      := docs/scientific_paper/paper_reimplementation.pdf
SOFTWARE   := docs/software_documentation/software_documentation.pdf
ALL_DOCS   := $(PWE_DOCS) $(GUIDES) $(PAPER) $(SOFTWARE)

# Entradas de cada documento: si cambian las figuras o tablas, el PDF se recompila.
TMM_FIGURES := $(wildcard results/tmm/figure_0*/output/*.pdf)
PWE_INPUTS  := $(wildcard results/pwe/*/output/*.pdf results/comparison/*/output/*.pdf \
                          results/comparison/tables/*.tex)

all: tests results documents verify

# Todo desde cero: resultados, documentos y verificación.
reproduce: results documents verify

tests:
	$(PYTHON) -m pytest tests

lint:
	$(PYTHON) -m ruff check src tests scripts presentacion/scripts
	$(PYTHON) -m mypy --cache-dir /tmp/mypy_cache

# Resultados por método: results/tmm, results/pwe y results/comparison.
results: tmm pwe comparison

tmm:
	$(PYTHON) scripts/tmm/run_all.py

pwe:
	$(PYTHON) scripts/pwe/run_all.py

comparison:
	$(PYTHON) scripts/comparison/run_all.py

# Compatibilidad: 'make figures' reproduce las figuras del paper con la TMM.
figures: tmm

documents: pwe-docs guides paper docs

pwe-docs: $(PWE_DOCS)
guides: $(GUIDES)
paper: $(PAPER)
docs: $(SOFTWARE)

$(PWE_DOCS): $(PWE_INPUTS) $(TMM_FIGURES)
$(GUIDES) $(PAPER) $(SOFTWARE): $(TMM_FIGURES)

docs/%.pdf: docs/%.tex
	cd $(dir $<) && for pass in 1 2; do \
	  pdflatex -interaction=nonstopmode -halt-on-error $(notdir $<) > /dev/null \
	  || { grep -a -A4 '^!' $(notdir $(basename $<)).log; exit 1; }; done

# Comprueba salidas, criterios TMM vs PWE y que los PDF estén al día y compilen limpios.
verify:
	$(PYTHON) scripts/verify_results.py

presentation:
	$(PYTHON) presentacion/scripts/generar_imagenes.py
	cd presentacion && pdflatex -interaction=nonstopmode presentacion.tex
	cd presentacion && pdflatex -interaction=nonstopmode presentacion.tex

# Presentación del método de ondas planas: usa results/pwe y results/comparison ya generados.
presentation-pwe:
	$(PYTHON) presentacion/scripts/generar_imagenes_pwe.py
	cd presentacion && pdflatex -interaction=nonstopmode presentacion_pwe.tex
	cd presentacion && pdflatex -interaction=nonstopmode presentacion_pwe.tex

clean:
	find docs presentacion -type f \( -name '*.aux' -o -name '*.log' -o -name '*.out' -o -name '*.toc' \
	  -o -name '*.nav' -o -name '*.snm' -o -name '*.vrb' -o -name '*.lof' -o -name '*.lot' \) -delete
