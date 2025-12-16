#!/usr/bin/env python3
# analysis/unfiltered_mining/render_chosen_unicorns.py

from __future__ import annotations

import csv
import sys
from pathlib import Path

from indigo import Indigo
from indigo.renderer import IndigoRenderer

csv.field_size_limit(sys.maxsize)


# ----------------------------
# Paths
# ----------------------------

BASE_DIR = Path(__file__).resolve().parent
INPUT_CSV = BASE_DIR / "chosen_unicorns.csv"
OUT_BASE = BASE_DIR / "output"
OUT_DIR = OUT_BASE / "chosen_unicorn_svgs"
FAILURES_CSV = OUT_BASE / "chosen_unicorn_svgs_failures.csv"


# ----------------------------
# Indigo setup + rendering
# ----------------------------

def _make_indigo() -> tuple[Indigo, IndigoRenderer]:
    indigo = Indigo()
    renderer = IndigoRenderer(indigo)

    indigo.setOption("render-output-format", "svg")
    indigo.setOption("render-implicit-hydrogens-visible", False)
    indigo.setOption("embedding-uniqueness", "none")
    indigo.setOption("render-bond-length", 40)
    indigo.setOption("render-atom-ids-visible", False)
    indigo.setOption("render-aam-color", "0, 0, 0")

    return indigo, renderer


def render_indigo_reaction_svg(indigo: Indigo, renderer: IndigoRenderer, smirks: str, out_svg: Path) -> None:
    rxn = indigo.loadReactionSmarts(smirks)
    renderer.renderToFile(rxn, str(out_svg))


# ----------------------------
# IO
# ----------------------------

def _read_rows(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if "rxn_idx" not in (reader.fieldnames or []) or "smirks" not in (reader.fieldnames or []):
            raise RuntimeError(f"Expected columns rxn_idx, smirks in {path}")

        for row in reader:
            rxn_idx = (row.get("rxn_idx") or "").strip()
            smirks = (row.get("smirks") or "").strip()
            if not rxn_idx or not smirks:
                continue
            rows.append((rxn_idx, smirks))
    return rows


def main() -> None:
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Missing input: {INPUT_CSV}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_BASE.mkdir(parents=True, exist_ok=True)

    indigo, renderer = _make_indigo()
    rows = _read_rows(INPUT_CSV)

    n_ok = 0
    with FAILURES_CSV.open("w", newline="") as ff:
        writer = csv.DictWriter(ff, fieldnames=["rxn_idx", "error"])
        writer.writeheader()

        for rxn_idx, smirks in rows:
            out_svg = OUT_DIR / f"{rxn_idx}.svg"
            try:
                render_indigo_reaction_svg(indigo, renderer, smirks, out_svg)
                n_ok += 1
            except Exception as e:
                writer.writerow({"rxn_idx": rxn_idx, "error": str(e)})

    print(f"Read:  {INPUT_CSV}")
    print(f"Wrote: {OUT_DIR}  ({n_ok}/{len(rows)} SVGs)")
    print(f"Wrote: {FAILURES_CSV}")


if __name__ == "__main__":
    main()
