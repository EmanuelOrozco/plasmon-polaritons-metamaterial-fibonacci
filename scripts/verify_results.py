#!/usr/bin/env python3
"""Verificación automática de resultados y documentos.

Comprueba, sin recalcular nada:

1. que existan todas las salidas de ``results/{tmm,pwe,comparison}``;
2. los criterios de aceptación de la comparación TMM vs PWE (conteos F_{m-2},
   error de semitraza, errores de borde y de ancho, acuerdo banda/gap);
3. que cada tabla (``\\tabla{...}``) y figura (``\\includegraphics``) citada por
   los documentos LaTeX exista, que cada PDF sea más reciente que sus entradas y
   que su log no tenga errores, referencias indefinidas ni desbordes;
4. con ``--reference DIR``, que los ``summary.json`` y ``.npz`` coincidan con los
   de una corrida anterior guardada en DIR (se ignoran los tiempos).

Sale con código 1 si alguna comprobación falla.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"
FIGURES = [f"figure_0{i}" for i in range(1, 7)]

MAX_DELTA_R_IN_BANDS = {"figure_01": 5e-2, "figure_02": 5e-3, "figure_03": 2e-3, "figure_04": 2e-4}
"""Tolerancia de max|R_PWE − R_TMM| en bandas; la Fig. 1 llega a contraste 400:1 a 0.15 GHz."""
MIN_CLASSIFICATION_AGREEMENT = 0.999
MAX_EDGE_ERROR_GHZ = 1e-6
MAX_RELATIVE_WIDTH_ERROR = 1e-5
MAX_BOOK_RESIDUAL_INVERSE = 1e-3
TIMING_KEY = re.compile(r"seconds|timing|efficiency|ratio|outputs|^tables$|evaluations|^time")
DOCUMENTS = [
    "comparacion_tmm_pwe/comparacion_tmm_pwe",
    "metodo_ondas_planas/metodo_ondas_planas",
    "scientific_paper/paper_reimplementation",
    "software_documentation/software_documentation",
    "guia_conceptual/guia_conceptual",
    "solucion_numerica/solucion_numerica",
    "arquitectura/arquitectura",
]


class Report:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.passed = 0

    def check(self, ok: bool, label: str, detail: str = "") -> bool:
        if ok:
            self.passed += 1
            print(f"  ok     {label}")
        else:
            self.failures.append(f"{label}: {detail}")
            print(f"  FALLA  {label}  {detail}")
        return ok


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


# ------------------------------------------------------------------ 1. salidas
def expected_outputs() -> list[Path]:
    paths = []
    for fig in FIGURES:
        paths += [RESULTS / "tmm" / fig / "output" / f"{fig}.{ext}" for ext in ("pdf", "png", "svg")]
        paths += [RESULTS / "tmm" / fig / "data" / "summary.json"]
        paths += [RESULTS / "pwe" / fig / "output" / f"{fig}.{ext}" for ext in ("pdf", "png")]
        paths += [RESULTS / "pwe" / fig / "data" / "summary.json", RESULTS / "comparison" / fig / "data" / "summary.json"]
    paths.append(RESULTS / "tmm" / "convergence" / "convergence_figure04.json")
    for study in ("book_benchmark", "cell_fourier", "spectral_pollution", "convergence"):
        paths.append(RESULTS / "pwe" / study / "data" / "summary.json")
    paths += [RESULTS / "comparison" / "summary.json", RESULTS / "comparison" / "efficiency" / "data" / "summary.json",
              RESULTS / "comparison" / "efficiency" / "output" / "efficiency.pdf",
              RESULTS / "comparison" / "timing" / "output" / "timing.pdf"]
    return paths


def check_outputs(report: Report) -> None:
    print("\n1. Salidas")
    missing = [p for p in expected_outputs() if not p.exists()]
    report.check(not missing, f"{len(expected_outputs())} archivos esperados en results/",
                 ", ".join(str(p.relative_to(ROOT)) for p in missing[:5]))
    for fig in ("figure_01", "figure_02", "figure_03", "figure_04"):
        n = len(list((RESULTS / "pwe" / fig / "data").glob("*.npz")))
        report.check(n > 0, f"barridos PWE de {fig} ({n} archivos .npz)", "faltan datos para la comparación")


# ------------------------------------------------------------- 2. aceptación
def check_acceptance(report: Report) -> None:
    print("\n2. Criterios de aceptación TMM vs PWE")
    summary = load(RESULTS / "comparison" / "summary.json")
    for fig, data in summary["dispersion"].items():
        delta = data["worst_max_abs_delta_r_in_bands"]
        agree = data["worst_classification_agreement"]
        report.check(delta < MAX_DELTA_R_IN_BANDS[fig], f"{fig}: max|ΔR| en bandas = {delta:.1e}",
                     f"supera {MAX_DELTA_R_IN_BANDS[fig]:.0e}")
        report.check(agree >= MIN_CLASSIFICATION_AGREEMENT, f"{fig}: acuerdo banda/gap = {agree:.5f}",
                     f"menor que {MIN_CLASSIFICATION_AGREEMENT}")
    report.check(summary["mode_counts"]["all_equal"], "Figs. 2 y 4: conteo TMM = PWE = F_{m-2}")
    plasmon = summary["plasmon"]
    report.check(plasmon["counts_equal"], f"Figs. 5–6: conteo TMM = PWE en los {plasmon['n_cases']} casos")
    for fig in ("figure_05", "figure_06"):
        cases = load(RESULTS / "comparison" / fig / "data" / "summary.json")["cases"]
        flat = [c for group in cases.values() for c in group]
        bad = [(c["m"], round(c["theta"], 4)) for c in flat if not c["n_tmm"] == c["n_pwe"] == c["expected"]]
        report.check(not bad, f"{fig}: {len(flat)} casos siguen la regla F_{{m-2}}", f"fallan {bad[:5]}")
    report.check(plasmon["worst_edge_error_ghz"] < MAX_EDGE_ERROR_GHZ,
                 f"peor error de borde = {plasmon['worst_edge_error_ghz']:.1e} GHz", f"supera {MAX_EDGE_ERROR_GHZ:.0e}")
    report.check(plasmon["worst_relative_width_error"] < MAX_RELATIVE_WIDTH_ERROR,
                 f"peor error relativo de ancho = {plasmon['worst_relative_width_error']:.1e}",
                 f"supera {MAX_RELATIVE_WIDTH_ERROR:.0e}")
    numerics = summary["pwe_numerics"]
    report.check(numerics["book_residual_inverse"] < MAX_BOOK_RESIDUAL_INVERSE < numerics["book_residual_laurent"],
                 f"libro: residuo regla inversa {numerics['book_residual_inverse']:.1e} "
                 f"< Laurent {numerics['book_residual_laurent']:.2f}")
    report.check(all(int(v) > 0 for v in numerics["pollution"].values()),
                 "forma ω(k) con Drude: aparecen autovalores espurios en todos los N")
    tmm04 = load(RESULTS / "tmm" / "convergence" / "convergence_figure04.json")
    report.check(bool(tmm04), "convergencia TMM de la Fig. 4 disponible")
    eff = summary["efficiency"]
    low, high = eff["figure_ratio_range"]
    report.check(eff["recommended_method"] == "TMM" and low > 1,
                 f"tiempo: PWE/TMM por figura entre {low:.1e} y {high:.1e}; método recomendado {eff['recommended_method']}")


# -------------------------------------------------------------- 3. documentos
def graphics_paths(tex: str, doc_dir: Path) -> list[Path]:
    match = re.search(r"\\graphicspath\{((?:\{[^}]*\})+)\}", tex)
    dirs = re.findall(r"\{([^}]*)\}", match.group(1)) if match else []
    return [(doc_dir / d).resolve() for d in dirs] + [doc_dir]


def resolve_graphic(name: str, bases: list[Path]) -> Path | None:
    for base in bases:
        for suffix in ("", ".pdf", ".png", ".jpg"):
            candidate = base / f"{name}{suffix}"
            if candidate.is_file():
                return candidate
    return None


def check_documents(report: Report) -> None:
    print("\n3. Documentos")
    for doc in DOCUMENTS:
        tex_path = DOCS / f"{doc}.tex"
        pdf_path, log_path = tex_path.with_suffix(".pdf"), tex_path.with_suffix(".log")
        tex = re.sub(r"(?<!\\)%.*", "", tex_path.read_text(encoding="utf-8"))
        name = Path(doc).name
        inputs = [tex_path]
        tables = re.findall(r"\\tabla\{([^}]+)\}", tex)
        missing_tables = [t for t in tables if not (RESULTS / "comparison" / "tables" / f"{t}.tex").is_file()]
        inputs += [RESULTS / "comparison" / "tables" / f"{t}.tex" for t in tables if t not in missing_tables]
        bases = graphics_paths(tex, tex_path.parent)
        graphics = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex)
        resolved = {g: resolve_graphic(g, bases) for g in graphics}
        missing_graphics = [g for g, p in resolved.items() if p is None]
        inputs += [p for p in resolved.values() if p is not None]
        report.check(not missing_tables and not missing_graphics,
                     f"{name}: {len(tables)} tablas y {len(graphics)} figuras citadas existen",
                     f"faltan {missing_tables + missing_graphics}")
        if not pdf_path.exists():
            report.check(False, f"{name}: PDF compilado", "no existe; ejecute make")
            continue
        stale = [p for p in inputs if p.stat().st_mtime > pdf_path.stat().st_mtime + 1]
        report.check(not stale, f"{name}: PDF al día con sus entradas",
                     f"{len(stale)} entradas más nuevas, p. ej. {stale[0].relative_to(ROOT)}" if stale else "")
        if log_path.exists():
            log = log_path.read_text(encoding="latin-1")
            problems = {"errores": len(re.findall(r"^! ", log, re.M)),
                        "referencias indefinidas": len(re.findall(r"undefined", log)),
                        "desbordes > 5 pt": len([v for v in re.findall(r"Overfull \\hbox \(([\d.]+)pt", log)
                                                  if float(v) > 5])}
            report.check(not any(problems.values()), f"{name}: log de LaTeX limpio",
                         ", ".join(f"{v} {k}" for k, v in problems.items() if v))


# --------------------------------------------------------- 4. reproducibilidad
def compare_values(new, old, path: str, tol: float, diffs: list[str]) -> None:
    if isinstance(new, dict) and isinstance(old, dict):
        for key in new.keys() & old.keys():
            if not TIMING_KEY.search(str(key)):
                compare_values(new[key], old[key], f"{path}.{key}", tol, diffs)
    elif isinstance(new, list) and isinstance(old, list):
        if len(new) != len(old):
            diffs.append(f"{path}: longitud {len(new)} ≠ {len(old)}")
        for i, (a, b) in enumerate(zip(new, old)):
            compare_values(a, b, f"{path}[{i}]", tol, diffs)
    elif isinstance(new, (int, float)) and isinstance(old, (int, float)) and not isinstance(new, bool):
        if not np.isclose(new, old, rtol=tol, atol=1e-12, equal_nan=True):
            diffs.append(f"{path}: {new!r} ≠ {old!r}")
    elif new != old and not isinstance(new, str):
        diffs.append(f"{path}: {new!r} ≠ {old!r}")


def check_reference(report: Report, reference: Path, tol: float) -> None:
    print(f"\n4. Reproducibilidad frente a {reference}")
    old_files = sorted(p for p in reference.rglob("*") if p.suffix in (".json", ".npz"))
    for old_path in old_files:
        rel = old_path.relative_to(reference)
        new_path = ROOT / rel
        if "efficiency" in rel.parts or "timing" in rel.parts:
            continue
        if not new_path.exists():
            report.check(False, str(rel), "falta en la corrida nueva")
            continue
        diffs: list[str] = []
        if old_path.suffix == ".json":
            compare_values(load(new_path), load(old_path), rel.name, tol, diffs)
        else:
            with np.load(new_path) as new, np.load(old_path) as old:
                for key in set(new.files) & set(old.files) - {"seconds"}:
                    if new[key].shape != old[key].shape or not np.allclose(new[key], old[key], rtol=tol,
                                                                           atol=1e-12, equal_nan=True):
                        diffs.append(key)
        report.check(not diffs, f"{rel} coincide", "; ".join(diffs[:3]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference", type=Path, help="carpeta con una copia previa de results/")
    parser.add_argument("--rtol", type=float, default=1e-6, help="tolerancia relativa para --reference")
    parser.add_argument("--skip-docs", action="store_true", help="no revisar los documentos LaTeX")
    args = parser.parse_args()
    report = Report()
    check_outputs(report)
    check_acceptance(report)
    if not args.skip_docs:
        check_documents(report)
    if args.reference:
        check_reference(report, args.reference, args.rtol)
    print(f"\n{report.passed} comprobaciones correctas, {len(report.failures)} fallas")
    if report.failures:
        for failure in report.failures:
            print("  -", failure)
        sys.exit(1)


if __name__ == "__main__":
    main()
