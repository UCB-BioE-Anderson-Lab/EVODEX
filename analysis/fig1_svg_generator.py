#!/usr/bin/env python3
"""
Compound SVG generator (Fig 5 style, constant bond length + constant font size)

Generates:
  - structure_svgs/compound_<NAME>.svg for each compound
  - structure_svgs/compound_failures.csv if any entries fail

Run:
  python analysis/compound_svg_generator.py
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from typing import List, Optional, Tuple

# ----------------------------
# Depiction params (style knobs)
# ----------------------------
REMOVE_EXPLICIT_H_FOR_DRAWING = True

BOND_LINE_W = 2.2
PADDING_FRAC = 0.08  # still applied, even with flexicanvas

# The two knobs that actually enforce constant visual scale across molecules:
FIXED_BOND_LENGTH_PX = 28  # constant bond length in pixels
FIXED_FONT_SIZE_PX = 16    # constant atom-label font size in pixels

# ----------------------------
# Compounds (names must be unique)
# ----------------------------
COMPOUNDS: List[Tuple[str, str]] = [
    ("propanol", "CCC[OH]"),
    ("propanal", "CCC=O"),
    ("methylether", "CCCOC"),
    ("propane", "CCC"),
    ("alpha_keto", "O=C(Cc1ccc(OC)cc1)C(=O)O"),
    ("alpha_hydroxy", "OC(Cc1ccc(OC)cc1)C(=O)O"),
    ("methylalphahydroxy", "COC(Cc1ccc(OC)cc1)C(=O)O"),
]


@dataclass
class Entry:
    name: str
    smiles: str
    safe_name: str


def _script_dir() -> str:
    return os.path.dirname(os.path.abspath(__file__))


def _repo_root(script_dir: str) -> str:
    # script is analysis/compound_svg_generator.py
    return os.path.dirname(script_dir)


def _out_dir(repo_root: str) -> str:
    return os.path.join(repo_root, "structure_svgs")


def _sanitize_label(label: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", label)
    safe = safe.strip("_")
    return safe or "ID"


def _strip_xml_header(svg: str) -> str:
    return re.sub(r"^\s*<\?xml[^>]*\?>\s*", "", svg.strip())


def _remove_common_background_rect(svg: str) -> str:
    # RDKit sometimes inserts a full-canvas rect background. Remove it.
    return re.sub(
        r'<rect[^>]*\bwidth="100%"\b[^>]*\bheight="100%"\b[^>]*/?>',
        "",
        svg,
        flags=re.IGNORECASE,
    )


def _prep_mol(smiles: str):
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        from rdkit.Chem import AllChem
    except Exception as e:
        raise RuntimeError(f"RDKit import failed: {e}")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse molecule string as SMILES")

    if REMOVE_EXPLICIT_H_FOR_DRAWING:
        try:
            mol = Chem.RemoveHs(mol)
        except Exception:
            pass

    # Clean depiction: no isotopes / atom maps
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

    # Ensure 2D coords exist
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass

    return mol


def _draw_mol_svg(smiles: str) -> str:
    """
    Draw molecule as an SVG with:
      - constant bond length (pixels)
      - constant font size (pixels)
      - monochrome palette
      - no forced bounding box sizing (flexicanvas)
    """
    try:
        from rdkit.Chem.Draw import rdMolDraw2D
    except Exception as e:
        raise RuntimeError(f"RDKit import failed: {e}")

    mol = _prep_mol(smiles)

    # Flexicanvas: let RDKit choose the canvas size needed for the fixed bond length.
    drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    opts = drawer.drawOptions()

    try:
        opts.useBWAtomPalette()
    except Exception:
        pass

    opts.bondLineWidth = float(BOND_LINE_W)
    opts.padding = float(PADDING_FRAC)

    # These are the key settings for constant appearance across different molecules:
    try:
        opts.fixedBondLength = float(FIXED_BOND_LENGTH_PX)
    except Exception:
        pass

    try:
        opts.fixedFontSize = int(FIXED_FONT_SIZE_PX)
    except Exception:
        pass

    # Do NOT set opts.baseFontSize when using fixedFontSize.

    try:
        opts.clearBackground = False
    except Exception:
        pass

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    svg = _strip_xml_header(svg)
    svg = _remove_common_background_rect(svg)
    return svg.strip() + "\n"


def main() -> int:
    script_dir = _script_dir()
    repo_root = _repo_root(script_dir)
    out_dir = _out_dir(repo_root)
    os.makedirs(out_dir, exist_ok=True)

    entries: List[Entry] = [
        Entry(name=name, smiles=smiles, safe_name=_sanitize_label(name))
        for name, smiles in COMPOUNDS
    ]

    failures: List[Tuple[str, str, str]] = []

    for entry in entries:
        out_path = os.path.join(out_dir, f"compound_{entry.safe_name}.svg")
        try:
            svg = _draw_mol_svg(entry.smiles)
            with open(out_path, "w", encoding="utf-8", newline="") as f:
                f.write(svg)
        except Exception as e:
            failures.append((entry.name, entry.smiles, str(e)))

    if failures:
        fail_path = os.path.join(out_dir, "compound_failures.csv")
        with open(fail_path, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["name", "smiles", "error"])
            for name, smiles, err in failures:
                w.writerow([name, smiles, err])

    print(f"Wrote per-entry SVGs: {out_dir}/compound_<NAME>.svg")
    if failures:
        print(f"Wrote failures: {os.path.join(out_dir, 'compound_failures.csv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())