#!/usr/bin/env python3
"""
Figure 5 SVG grid generator

Generates:
  - structure_svgs/fig5_grid.svg
  - structure_svgs/fig5_<LABEL>.svg for each row
  - structure_svgs/fig5_failures.csv if any entries fail

Inputs:
  - analysis/evodex_predicted_products.csv (same directory as this script)

Run:
  python analysis/fig5_svg_generator.py
"""

from __future__ import annotations

import csv
import math
import os
import re
import sys
from dataclasses import dataclass
from typing import List, Optional, Tuple

# ----------------------------
# Tweakable layout parameters
# ----------------------------
N_COLS = 5

CELL_W = 220
CELL_H = 190

PAD = 12
LABEL_H = 26

LABEL_FONT_SIZE = 16
LABEL_FONT_FAMILY = "Arial, Helvetica, sans-serif"

BOX_STROKE_W = 2
BOX_DASH = "8,6"

# Derived
MOL_W = CELL_W - 2 * PAD
MOL_H = CELL_H - (2 * PAD + LABEL_H)

# ----------------------------
# Tweakable depiction params
# ----------------------------
REMOVE_EXPLICIT_H_FOR_DRAWING = True

BOND_LINE_W = 2.2
BASE_FONT_SIZE = 0.6  # RDKit uses a relative scale; this is a reasonable knob
PADDING_FRAC = 0.08   # RDKit drawer padding fraction in the mol canvas


@dataclass
class Entry:
    evodex_id: str
    predicted: str
    curation: bool
    label: str
    safe_label: str


def _script_dir() -> str:
    return os.path.dirname(os.path.abspath(__file__))


def _repo_root(script_dir: str) -> str:
    # script is analysis/fig5_svg_generator.py
    return os.path.dirname(script_dir)


def _out_dir(repo_root: str) -> str:
    return os.path.join(repo_root, "structure_svgs")


def _read_csv(csv_path: str) -> List[Tuple[str, str, str]]:
    rows: List[Tuple[str, str, str]] = []
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"EVODEX ID", "Predicted Products", "curation"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {sorted(missing)}")

        for r in reader:
            evodex_id = (r.get("EVODEX ID") or "").strip()
            predicted = (r.get("Predicted Products") or "").strip()
            curation_raw = (r.get("curation") or "").strip()
            rows.append((evodex_id, predicted, curation_raw))
    return rows


def _parse_bool(s: str) -> bool:
    # Deterministic, permissive parsing
    v = (s or "").strip().lower()
    if v in {"true", "t", "1", "yes", "y"}:
        return True
    if v in {"false", "f", "0", "no", "n", ""}:
        return False
    # Unrecognized values treated as False (dashed box), but caller can log
    return False


def _extract_label(evodex_id: str) -> str:
    # Expected: EVODEX.2-D369 -> D369
    if not evodex_id:
        raise ValueError("Empty EVODEX ID")
    parts = evodex_id.split("-")
    if len(parts) < 2 or not parts[-1]:
        raise ValueError(f"Unexpected EVODEX ID format: {evodex_id}")
    label = parts[-1]
    # We expect something like D### but do not enforce beyond non-empty
    return label


def _sanitize_label(label: str) -> str:
    # Allow only [A-Za-z0-9_], replace others with _
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", label)
    safe = safe.strip("_")
    return safe or "ID"


def _svg_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;")
            .replace("'", "&apos;")
    )


def _strip_svg_wrapper(svg: str) -> str:
    """
    Given an SVG string, return the inner contents (everything inside the outer <svg> ... </svg>).
    Also strips optional XML header.
    """
    s = svg.strip()
    s = re.sub(r"^\s*<\?xml[^>]*\?>\s*", "", s)
    m = re.search(r"<svg[^>]*>", s, flags=re.IGNORECASE)
    if not m:
        return s
    start = m.end()
    end = re.search(r"</svg>\s*$", s, flags=re.IGNORECASE)
    if not end:
        return s[start:]
    return s[start:end.start()]


def _remove_common_background_rect(inner: str) -> str:
    """
    RDKit sometimes inserts a background rect. We keep our own full-canvas white background,
    so remove obvious full-rect backgrounds from the mol fragment.
    """
    s = inner
    # Remove rects that look like "background" (very common RDKit pattern)
    s = re.sub(
        r'<rect[^>]*\bwidth="100%"\b[^>]*\bheight="100%"\b[^>]*/?>',
        "",
        s,
        flags=re.IGNORECASE,
    )
    return s


def _draw_mol_svg_fragment(predicted: str) -> str:
    """
    Returns an SVG fragment (no outer <svg>) for a molecule drawn into MOL_W x MOL_H.
    Raises an exception on failure.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D
    except Exception as e:
        raise RuntimeError(f"RDKit import failed: {e}")

    mol = Chem.MolFromSmiles(predicted)
    if mol is None:
        raise ValueError("RDKit could not parse molecule string as SMILES")

    if REMOVE_EXPLICIT_H_FOR_DRAWING:
        try:
            mol = Chem.RemoveHs(mol)
        except Exception:
            # Keep original if removal fails; do not hard-fail
            pass

    # Strip atom-map numbers and isotopes for clean publication depiction
    for a in mol.GetAtoms():
        try:
            a.SetAtomMapNum(0)
        except Exception:
            pass
        try:
            a.SetIsotope(0)
        except Exception:
            pass

    mol = Draw.PrepareMolForDrawing(mol, addChiralHs=False)

    drawer = rdMolDraw2D.MolDraw2DSVG(int(MOL_W), int(MOL_H))
    opts = drawer.drawOptions()

    # Monochrome, JACS-like linework
    try:
        opts.useBWAtomPalette()
    except Exception:
        # fallback: do nothing, most RDKit builds still allow monochrome via palette
        pass

    opts.bondLineWidth = float(BOND_LINE_W)
    opts.padding = float(PADDING_FRAC)
    opts.baseFontSize = float(BASE_FONT_SIZE)

    # Avoid RDKit painting a background; we provide backgrounds ourselves
    try:
        opts.clearBackground = False
    except Exception:
        pass

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    inner = _strip_svg_wrapper(svg)
    inner = _remove_common_background_rect(inner)
    return inner.strip()


def _cell_svg(
    entry: Entry,
    mol_inner_svg: str,
    draw_box: bool,
) -> str:
    """
    Build a full per-entry SVG with CELL_W x CELL_H, including white background,
    optional dashed border, label, and mol fragment placed in mol region.
    """
    label = _svg_escape(entry.label)

    box_svg = ""
    if draw_box:
        inset = BOX_STROKE_W / 2.0
        box_svg = (
            f'<rect x="{inset:.3f}" y="{inset:.3f}" '
            f'width="{CELL_W - BOX_STROKE_W:.3f}" height="{CELL_H - BOX_STROKE_W:.3f}" '
            f'fill="none" stroke="black" stroke-width="{BOX_STROKE_W}" '
            f'stroke-dasharray="{BOX_DASH}" />'
        )

    # Label positioning: keep it simple and robust
    label_x = PAD
    label_y = PAD + LABEL_FONT_SIZE  # baseline

    mol_group = ""
    if mol_inner_svg:
        mol_group = (
            f'<g transform="translate({PAD},{PAD + LABEL_H})">\n'
            f'{mol_inner_svg}\n'
            f"</g>"
        )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{CELL_W}" height="{CELL_H}" '
        f'viewBox="0 0 {CELL_W} {CELL_H}">\n'
        f'<rect x="0" y="0" width="{CELL_W}" height="{CELL_H}" fill="white" />\n'
        f"{box_svg}\n"
        f'<text x="{label_x}" y="{label_y}" font-family="{LABEL_FONT_FAMILY}" '
        f'font-size="{LABEL_FONT_SIZE}" fill="black">{label}</text>\n'
        f"{mol_group}\n"
        f"</svg>\n"
    )
    return svg


def _grid_svg(entries: List[Entry], mol_fragments: List[str]) -> str:
    n = len(entries)
    n_rows = int(math.ceil(n / float(N_COLS))) if n else 1

    width = N_COLS * CELL_W
    height = n_rows * CELL_H

    parts: List[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
    )
    parts.append(f'<rect x="0" y="0" width="{width}" height="{height}" fill="white" />\n')

    for i, (entry, mol_inner) in enumerate(zip(entries, mol_fragments)):
        row = i // N_COLS
        col = i % N_COLS
        x0 = col * CELL_W
        y0 = row * CELL_H

        draw_box = (entry.curation is False)

        box_svg = ""
        if draw_box:
            inset = BOX_STROKE_W / 2.0
            box_svg = (
                f'<rect x="{x0 + inset:.3f}" y="{y0 + inset:.3f}" '
                f'width="{CELL_W - BOX_STROKE_W:.3f}" height="{CELL_H - BOX_STROKE_W:.3f}" '
                f'fill="none" stroke="black" stroke-width="{BOX_STROKE_W}" '
                f'stroke-dasharray="{BOX_DASH}" />\n'
            )

        label = _svg_escape(entry.label)
        label_x = x0 + PAD
        label_y = y0 + PAD + LABEL_FONT_SIZE

        parts.append(box_svg)
        parts.append(
            f'<text x="{label_x}" y="{label_y}" font-family="{LABEL_FONT_FAMILY}" '
            f'font-size="{LABEL_FONT_SIZE}" fill="black">{label}</text>\n'
        )

        if mol_inner:
            parts.append(
                f'<g transform="translate({x0 + PAD},{y0 + PAD + LABEL_H})">\n'
                f"{mol_inner}\n"
                f"</g>\n"
            )

    parts.append("</svg>\n")
    return "".join(parts)


def main() -> int:
    script_dir = _script_dir()
    repo_root = _repo_root(script_dir)

    csv_path = os.path.join(script_dir, "evodex_predicted_products.csv")
    out_dir = _out_dir(repo_root)
    os.makedirs(out_dir, exist_ok=True)

    rows = _read_csv(csv_path)

    entries: List[Entry] = []
    failures: List[Tuple[str, str, str]] = []
    mol_fragments: List[str] = []

    # Build entries first (label parsing failures are logged but still included)
    for evodex_id, predicted, curation_raw in rows:
        curation = _parse_bool(curation_raw)

        label = ""
        safe_label = ""
        label_error: Optional[str] = None
        try:
            label = _extract_label(evodex_id)
            safe_label = _sanitize_label(label)
        except Exception as e:
            label_error = str(e)
            label = evodex_id if evodex_id else "ID"
            safe_label = _sanitize_label(label)

        entry = Entry(
            evodex_id=evodex_id,
            predicted=predicted,
            curation=curation,
            label=label,
            safe_label=safe_label,
        )
        entries.append(entry)

        if label_error:
            failures.append((evodex_id, predicted, f"label_parse_error: {label_error}"))

    # Draw and write per-entry svgs
    for entry in entries:
        mol_inner = ""
        depiction_error: Optional[str] = None
        if entry.predicted:
            try:
                mol_inner = _draw_mol_svg_fragment(entry.predicted)
            except Exception as e:
                depiction_error = str(e)
                mol_inner = ""
        else:
            depiction_error = "Empty Predicted Products field"

        if depiction_error:
            failures.append((entry.evodex_id, entry.predicted, depiction_error))

        mol_fragments.append(mol_inner)

        per_svg = _cell_svg(
            entry=entry,
            mol_inner_svg=mol_inner,
            draw_box=(entry.curation is False),
        )
        per_path = os.path.join(out_dir, f"fig5_{entry.safe_label}.svg")
        with open(per_path, "w", encoding="utf-8", newline="") as f:
            f.write(per_svg)

    # Write grid svg
    grid = _grid_svg(entries, mol_fragments)
    grid_path = os.path.join(out_dir, "fig5_grid.svg")
    with open(grid_path, "w", encoding="utf-8", newline="") as f:
        f.write(grid)

    # Write failures report if needed
    if failures:
        fail_path = os.path.join(out_dir, "fig5_failures.csv")
        with open(fail_path, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["EVODEX ID", "Predicted Products", "error"])
            for evodex_id, predicted, err in failures:
                w.writerow([evodex_id, predicted, err])

    print(f"Wrote: {grid_path}")
    print(f"Wrote per-entry SVGs: {out_dir}/fig5_<LABEL>.svg")
    if failures:
        print(f"Wrote failures: {os.path.join(out_dir, 'fig5_failures.csv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())