#!/usr/bin/env python3
"""Resumen global TMM vs PWE: JSON, tablas LaTeX y figura de tiempos.

Las tablas en ``results/comparison/tables/*.tex`` se incluyen tal cual en el
anexo de comparación y en el paper del PWE, así las cifras de los documentos
salen siempre de los datos.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from _shared import PWE_RESULTS, ROOT
from _common import dump_json, load_json, relative, result_dirs
from fibonacci_photonics.plotting import apply_prb_style, save_figure

COMPARISON = ROOT / "results" / "comparison"
TABLES = COMPARISON / "tables"
THETA_NAMES = {0.0: "0", round(np.pi / 12, 6): r"\pi/12", round(np.pi / 6, 6): r"\pi/6", round(np.pi / 3, 6): r"\pi/3"}


def theta_name(theta: float) -> str:
    return THETA_NAMES.get(round(theta, 6), f"{np.degrees(theta):.0f}^\\circ")


def sci(value: float, digits: int = 1) -> str:
    if value == 0 or not np.isfinite(value):
        return "0" if value == 0 else "--"
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / 10**exponent
    return rf"${mantissa:.{digits}f}\times10^{{{exponent}}}$"


def write_table(name: str, header: str, rows: list[str], columns: str) -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    body = "\n".join(rows)
    (TABLES / f"{name}.tex").write_text(
        f"\\begin{{tabular}}{{{columns}}}\n\\toprule\n{header} \\\\\n\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}\n",
        encoding="utf-8",
    )


def dispersion_section() -> dict:
    rows, data = [], {}
    for figure_id in ("figure_01", "figure_02", "figure_03", "figure_04"):
        summary = load_json(COMPARISON / figure_id / "data" / "summary.json")
        data[figure_id] = summary
        for scan in sorted(summary["scans"], key=lambda s: (s["theta"], s["m"])):
            rows.append(
                f"{figure_id[-1]} & {scan['m']} & ${theta_name(scan['theta'])}$ & {scan['n_points']} & "
                f"{sci(scan['max_abs_delta_r_in_bands'])} & {scan['classification_agreement'] * 100:.2f}\\% & "
                f"{scan['tmm_seconds'] * 1e3:.1f} & {scan['pwe_seconds']:.1f} \\\\"
            )
    write_table("dispersion", r"Fig. & $m$ & $\theta$ & puntos & $\max|\Delta R|$ en bandas & acuerdo banda/gap"
                r" & TMM (ms) & PWE, pared (s)", rows, "rrrrrrrr")
    return data


def mode_count_section(dispersion: dict) -> dict:
    rows, out = [], []
    for figure_id in ("figure_02", "figure_04"):
        counts = dispersion[figure_id]["mode_counts"]
        for t, p in zip(counts["tmm"], counts["pwe"], strict=True):
            label = f"$\\theta={theta_name(np.pi / 3)}$" if figure_id == "figure_02" else f"$\\theta={t['theta'].replace('pi', chr(92) + 'pi')}$"
            expected = t.get("expected", {3: 1, 4: 2}.get(t["m"]))
            rows.append(f"{figure_id[-1]} & {t['m']} & {label} & {expected} & {t['detected']} & {p['detected']} \\\\")
            out.append({"figure": figure_id, "m": t["m"], "expected": expected, "tmm": t["detected"], "pwe": p["detected"]})
    write_table("mode_counts", r"Fig. & $m$ & ángulo & $F_{m-2}$ & TMM & PWE", rows, "rrrrrr")
    return {"rows": out, "all_equal": all(r["tmm"] == r["pwe"] == r["expected"] for r in out)}


def plasmon_section() -> dict:
    fig5 = load_json(COMPARISON / "figure_05" / "data" / "summary.json")
    fig6 = load_json(COMPARISON / "figure_06" / "data" / "summary.json")
    rows = []
    for label, cases in fig5["cases"].items():
        for c in cases:
            rows.append(
                f"${theta_name(c['theta'])}$ & {c['m']} & {c['expected']} & {c['n_tmm']} & {c['n_pwe']} & "
                f"{sci(c['max_abs_edge_error_ghz'])} & {sci(c['max_relative_width_error'])} & "
                f"{c['pwe_evaluations']} & {c['pwe_seconds']:.0f} \\\\"
            )
    write_table("figure_05", r"$\theta$ & $m$ & $F_{m-2}$ & TMM & PWE & $\max|\Delta\nu_{\rm borde}|$ (GHz)"
                r" & $\max|\Delta w|/w$ & evaluaciones & PWE (s)", rows, "rrrrrrrrr")
    rows = []
    for m, cases in sorted(fig6["cases"].items(), key=lambda kv: int(kv[0])):
        n_ok = sum(c["n_tmm"] == c["n_pwe"] == c["expected"] for c in cases)
        worst_edge = max(c["max_abs_edge_error_ghz"] for c in cases)
        worst_width = max(c["max_relative_width_error"] for c in cases)
        seconds = sum(c["pwe_seconds"] for c in cases)
        rows.append(f"{m} & {len(cases)} & {n_ok} & {sci(worst_edge)} & {sci(worst_width)} & {seconds:.0f} \\\\")
    write_table("figure_06", r"$m$ & ángulos & conteo $=F_{m-2}$ & $\max|\Delta\nu_{\rm borde}|$ (GHz)"
                r" & $\max|\Delta w|/w$ & PWE total (s)", rows, "rrrrrr")
    all_cases = [c for cases in fig5["cases"].values() for c in cases] + [c for cases in fig6["cases"].values() for c in cases]
    return {
        "worst_edge_error_ghz": max(c["max_abs_edge_error_ghz"] for c in all_cases),
        "worst_relative_width_error": max(c["max_relative_width_error"] for c in all_cases),
        "counts_equal": all(c["n_tmm"] == c["n_pwe"] for c in all_cases),
        "n_cases": len(all_cases),
    }


def pwe_numerics_section() -> dict:
    conv = load_json(PWE_RESULTS / "convergence" / "data" / "summary.json")
    book = load_json(PWE_RESULTS / "book_benchmark" / "data" / "summary.json")
    pollution = load_json(PWE_RESULTS / "spectral_pollution" / "data" / "summary.json")
    rows = [f"{r['m']} & {r['layers']} & {r['harmonics_per_layer']} & {r['n_plane_waves']} & "
            f"{r['tmm_seconds_per_frequency'] * 1e6:.1f} & {r['pwe_seconds_per_frequency']:.3f} & "
            f"{sci(r['pwe_seconds_per_frequency'] / r['tmm_seconds_per_frequency'], 0)} & "
            f"{sci(r['max_abs_delta_r'])} & {sci(r['max_abs_r_tmm'])} \\\\" for r in conv["cost"]]
    write_table("cost", r"$m$ & capas & armónicos/capa & $N$ & TMM ($\mu$s) & PWE (s) & PWE/TMM"
                r" & $\max|\Delta R|$ & $\max|R|$", rows, "rrrrrrrrr")
    rows = []
    for n in sorted({r["n_plane_waves"] for r in conv["book"]}):
        lau = next(r for r in conv["book"] if r["rule"] == "laurent" and r["n_plane_waves"] == n)
        inv = next(r for r in conv["book"] if r["rule"] == "inverse" and r["n_plane_waves"] == n)
        rows.append(f"{n} & {sci(lau['max_rel_error'])} & {sci(inv['max_rel_error'])} \\\\")
    write_table("book_convergence", r"$N$ & Laurent (libro) & regla inversa", rows, "rrr")
    rows = [f"{n} & {len(v)} & {pollution['spurious_counts'][n]} \\\\" for n, v in pollution["spectra_ghz"].items()]
    write_table("pollution", r"$N$ & autovalores en 0.5--3 GHz & espurios", rows, "rrr")
    return {"book_residual_laurent": book["max_abs_tmm_residual_laurent"],
            "book_residual_inverse": book["max_abs_tmm_residual_inverse"],
            "cost": conv["cost"], "pollution": pollution["spurious_counts"]}


def timing_figure() -> dict:
    """Tiempo de cada figura en un núcleo (de ``efficiency.py``), con la razón PWE/TMM."""
    rows = load_json(COMPARISON / "efficiency" / "data" / "summary.json")["figures"]["rows"]
    labels = [f"Fig. {r['figure'][-1]}" for r in rows]
    tmm = [r["tmm_seconds"] for r in rows]
    pwe = [r["pwe_single_core_seconds"] for r in rows]
    apply_prb_style()
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    x = np.arange(len(labels))
    ax.bar(x - 0.2, tmm, width=0.4, color="#2c4d8c", label="TMM")
    ax.bar(x + 0.2, pwe, width=0.4, color="#c23b22", label="PWE")
    for xi, t, p in zip(x, tmm, pwe, strict=True):
        exponent = int(np.floor(np.log10(p / t)))
        ax.annotate(rf"$\times{p / t / 10**exponent:.1f}{{\cdot}}10^{{{exponent}}}$", (xi + 0.2, p), ha="center",
                    va="bottom", fontsize=5.5, xytext=(0, 2), textcoords="offset points")
    ax.set_yscale("log")
    ax.set_ylim(top=max(pwe) * 20)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("tiempo de cómputo en un núcleo (s)")
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    dirs = result_dirs("comparison", "timing")
    paths = save_figure(fig, dirs["output"], "timing", formats=("pdf", "png"))
    return {"labels": labels, "tmm_seconds": tmm, "pwe_single_core_seconds": pwe,
            "outputs": {k: relative(v) for k, v in paths.items()}}


def efficiency_section() -> dict:
    eff = load_json(COMPARISON / "efficiency" / "data" / "summary.json")
    ratios = [r["ratio"] for r in eff["figures"]["rows"]]
    wp = eff["work_precision"]["rows"]
    tmm_best = min((r for r in wp if r["method"].startswith("TMM")), key=lambda r: r["seconds_per_frequency"])
    pwe_best = min((r for r in wp if r["method"] == "PWE"), key=lambda r: r["max_abs_error"])
    return {
        "pwe_time_exponent_in_n": eff["scaling"]["fit_pwe"]["alpha"],
        "tmm_product_time_exponent_in_layers": eff["scaling"]["fit_tmm_product"]["alpha"],
        "figure_ratio_range": [min(ratios), max(ratios)],
        "tmm_fastest": tmm_best,
        "pwe_most_accurate": pwe_best,
        "book": eff["book"],
        "recommended_method": "TMM",
    }


def main() -> None:
    dispersion = dispersion_section()
    summary = {
        "dispersion": {k: {"worst_max_abs_delta_r_in_bands": v["worst_max_abs_delta_r_in_bands"],
                           "worst_classification_agreement": v["worst_classification_agreement"]}
                       for k, v in dispersion.items()},
        "mode_counts": mode_count_section(dispersion),
        "plasmon": plasmon_section(),
        "pwe_numerics": pwe_numerics_section(),
        "timing": timing_figure(),
        "efficiency": efficiency_section(),
        "tables": sorted(relative(p) for p in TABLES.glob("*.tex")),
    }
    dump_json(COMPARISON / "summary.json", summary)
    print("Resumen escrito en", relative(COMPARISON / "summary.json"))
    print("  peor max|ΔR| (Figs. 1–4):", max(v["worst_max_abs_delta_r_in_bands"] for v in summary["dispersion"].values()))
    print("  peor error de borde (Figs. 5–6):", summary["plasmon"]["worst_edge_error_ghz"], "GHz")
    print("  conteos iguales:", summary["mode_counts"]["all_equal"], summary["plasmon"]["counts_equal"])


if __name__ == "__main__":
    main()
