#!/usr/bin/env python3
"""
Project the aryl_bromination EVODEX-E and RdChiral operators onto the
Hammett substrate set plus pyridine and indole, and draw a table that
shows all products.

Columns:
    Substituent | Substrate | EVODEX-E products | RdChiral products
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from rdchiral.main import rdchiralRun, rdchiralReactants, rdchiralReaction
from evodex.projection import project_operator


# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent

PROJECT_PATH = BASE_DIR / "out" / "2_ops" / "aryl_bromination.project.jsonl"

OUTPUT_DIR = BASE_DIR / "out" / "rdchiral_projection"
OUTPUT_SVG = OUTPUT_DIR / "aryl_bromination_dual_projection.svg"

# Hard coded EVODEX E operator (most active, updated mapping)
EVODEX_E_SMIRKS = (
    "[#1][c:1]1[c:2]([#1])[c:3]([#1])[c:4]([#1])[c:5]([#1])[c:6]1[N:7]"
    ">>"
    "[#1][c:1]1[c:2]([#1])[c:3](Br)[c:4]([#1])[c:5]([#1])[c:6]1[N:7]"
)

# RdChiral operator will be read from JSONL
OP_FIELDS = {
    "RdChiral": "rdchiral_smirks",
}

# Drawing geometry
ROW_HEIGHT = 130
ROW_GAP = 10
MARGIN_X = 20
MARGIN_TOP = 40
LABEL_COL_W = 80
SUBSTRATE_COL_W = 150
ARROW_COL_W = 40
HEADER_HEIGHT = 80
BACKGROUND_COLOR = "white"


# ---------------------------------------------------------------------------
# Helpers for SVG and RDKit drawings
# ---------------------------------------------------------------------------

def _strip_xml_header(svg_text: str) -> str:
    return re.sub(r"^\s*<\?xml[^>]*\?>\s*", "", svg_text)


def _svg_inner(svg_text: str) -> str:
    s = _strip_xml_header(svg_text).strip()
    m_open = re.search(r"<svg\b[^>]*>", s, flags=re.I)
    m_close = re.search(r"</svg\s*>", s, flags=re.I)

    if not (m_open and m_close and m_close.start() > m_open.end()):
        return s
    return s[m_open.end(): m_close.start()]


def mol_to_group(smiles: str, x: int, y: int, w: int, h: int) -> str:
    if not smiles:
        return error_panel(x, y, w, h, "No SMILES")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return error_panel(x, y, w, h, "Parse error")

    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    inner = _svg_inner(drawer.GetDrawingText())
    return f'<g transform="translate({x},{y})">{inner}</g>'


def error_panel(x: int, y: int, w: int, h: int, msg: str = "Error") -> str:
    return (
        f'<g transform="translate({x},{y})">'
        f'<rect x="0" y="0" width="{w}" height="{h}" '
        f'rx="8" ry="8" fill="#ffe6e6" stroke="#cc0000" stroke-width="1"/>'
        f'<text x="{w/2}" y="{h/2}" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#cc0000">'
        f"{msg}</text>"
        f"</g>"
    )


def text_label(x: int, y: int, text: str, anchor: str = "middle", size: int = 12) -> str:
    return (
        f'<text x="{x}" y="{y}" text-anchor="{anchor}" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="{size}">{text}</text>'
    )


# ---------------------------------------------------------------------------
# Data loading and operators
# ---------------------------------------------------------------------------

def load_aryl_bromination_records(path: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError as e:
                raise RuntimeError(
                    f"Malformed JSON on line {lineno} of {path}: {e}"
                ) from e
    if not rows:
        raise RuntimeError(f"No records loaded from {path}")
    return rows


def extract_operator(rows: List[Dict[str, Any]], field: str) -> str:
    values = {row[field] for row in rows if row.get(field)}
    if not values:
        raise RuntimeError(f"No {field} field found in project file")
    if len(values) > 1:
        raise RuntimeError(f"Expected a single operator in {field}, found {len(values)}")
    return next(iter(values))


def preprocess_for_display(smiles: str) -> str:
    # Replace At with H so the first row shows benzene
    return smiles.replace("[At]", "[H]")


def preprocess_for_reaction(smiles: str) -> str:
    return preprocess_for_display(smiles)


def build_substrate_panel(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    panel: List[Dict[str, Any]] = []

    for row in rows:
        label = row.get("label", f"row_{row.get('row', len(panel))}")
        smiles_raw = row.get("smiles_input") or row.get("smiles")
        if not smiles_raw:
            continue
        panel.append(
            {
                "label": label,
                "smiles_display": preprocess_for_display(smiles_raw),
                "smiles_react": preprocess_for_reaction(smiles_raw),
            }
        )

    panel.append(
        {
            "label": "pyridine",
            "smiles_display": "c1ccccn1",
            "smiles_react": "c1ccccn1",
        }
    )
    panel.append(
        {
            "label": "indole",
            "smiles_display": "c1ccc2[nH]ccc2c1",
            "smiles_react": "c1ccc2[nH]ccc2c1",
        }
    )

    return panel


# ---------------------------------------------------------------------------
# Projection
# ---------------------------------------------------------------------------

def run_rdchiral_projection(
    smirks: str,
    substrates: List[Dict[str, Any]],
    field_out: str,
) -> None:
    """Apply a RdChiral operator to each substrate."""
    rxn = rdchiralReaction(smirks)

    for rec in substrates:
        label = rec["label"]
        smiles_for_rxn = rec["smiles_react"]

        try:
            rct = rdchiralReactants(smiles_for_rxn)
            raw_products = rdchiralRun(rxn, rct)
        except Exception as exc:
            print(f"[WARN] RdChiral error for {field_out} on {label} ({smiles_for_rxn}): {exc}")
            products: List[str] = []
        else:
            flat: List[str] = []
            for p in raw_products or []:
                if isinstance(p, (list, tuple)):
                    if len(p) == 1:
                        flat.append(p[0])
                    else:
                        flat.append(".".join(p))
                else:
                    flat.append(p)
            seen = set()
            products = []
            for smi in flat:
                if smi not in seen:
                    seen.add(smi)
                    products.append(smi)

        rec[field_out] = products


def run_evodex_projection(
    smirks: str,
    substrates: List[Dict[str, Any]],
    field_out: str,
) -> None:
    """
    Apply EVODEX-E operator using evodex.projection.project_operator.

    project_operator already adds explicit Hs and handles permutations.
    """
    for rec in substrates:
        label = rec["label"]
        smiles_for_rxn = rec["smiles_react"]

        try:
            raw_sets = project_operator(smirks, smiles_for_rxn)
        except Exception as exc:
            print(f"[WARN] EVODEX-E error on {label} ({smiles_for_rxn}): {exc}")
            rec[field_out] = []
            continue

        # raw_sets is a list of strings, each a '.'-joined product set
        flat: set[str] = set()
        for s in raw_sets or []:
            for smi in s.split("."):
                if smi:
                    flat.add(smi)

        rec[field_out] = sorted(flat)


# ---------------------------------------------------------------------------
# SVG table construction
# ---------------------------------------------------------------------------

def build_svg_table(
    substrates: List[Dict[str, Any]],
    op_smirks: Dict[str, str],
) -> str:
    op_labels = ["EVODEX-E", "RdChiral"]

    # For sizing: find max number of products per operator across rows
    max_prods: Dict[str, int] = {lab: 0 for lab in op_labels}
    for rec in substrates:
        for lab in op_labels:
            n = len(rec.get(f"{lab}_products", []))
            if n > max_prods[lab]:
                max_prods[lab] = n
    for lab in op_labels:
        if max_prods[lab] < 1:
            max_prods[lab] = 1

    n_rows = len(substrates)

    # Operator column widths: number of slots times substrate width
    op_widths: Dict[str, int] = {
        lab: max_prods[lab] * SUBSTRATE_COL_W for lab in op_labels
    }

    total_w = (
        MARGIN_X * 2
        + LABEL_COL_W
        + SUBSTRATE_COL_W
        + ARROW_COL_W
        + sum(op_widths[lab] for lab in op_labels)
    )
    total_h = (
        MARGIN_TOP
        + HEADER_HEIGHT
        + n_rows * ROW_HEIGHT
        + (n_rows - 1) * ROW_GAP
        + MARGIN_TOP
    )

    parts: List[str] = []

    parts.append(
        f'<rect x="0" y="0" width="{int(total_w)}" height="{int(total_h)}" '
        f'fill="{BACKGROUND_COLOR}"/>'
    )

    title = "EVODEX-E vs RdChiral projection of aryl_bromination operators"
    title_x = total_w / 2
    title_y = MARGIN_TOP

    # Main title
    parts.append(
        text_label(int(title_x), int(title_y), title, anchor="middle", size=16)
    )

    # Put the two operator strings on separate lines, slightly smaller font
    evodex_line = f"EVODEX-E: {op_smirks['EVODEX-E']}"
    rdchiral_line = f"RdChiral: {op_smirks['RdChiral']}"

    parts.append(
        text_label(
            int(title_x),
            int(title_y + 22),
            evodex_line,
            anchor="middle",
            size=9,
        )
    )
    parts.append(
        text_label(
            int(title_x),
            int(title_y + 38),
            rdchiral_line,
            anchor="middle",
            size=9,
        )
    )

    # Now the column headings sit lower, with room above for the SMIRKS
    header_y = MARGIN_TOP + HEADER_HEIGHT - 10

    header_y = MARGIN_TOP + HEADER_HEIGHT - 10
    x_label = MARGIN_X + LABEL_COL_W / 2
    x_sub = MARGIN_X + LABEL_COL_W + SUBSTRATE_COL_W / 2
    # center of the small arrow column (no label)
    x_arrow = MARGIN_X + LABEL_COL_W + SUBSTRATE_COL_W + ARROW_COL_W / 2

    parts.append(text_label(int(x_label), header_y, "Substituent", anchor="middle", size=12))
    parts.append(text_label(int(x_sub), header_y, "Substrate", anchor="middle", size=12))
    # no label over the arrow; leave x_arrow unused here

    # Operator column headings
    op_x0 = MARGIN_X + LABEL_COL_W + SUBSTRATE_COL_W + ARROW_COL_W
    cur_x = op_x0
    for lab in op_labels:
        cx = cur_x + op_widths[lab] / 2
        parts.append(text_label(int(cx), header_y, lab, anchor="middle", size=12))
        cur_x += op_widths[lab]

    start_y = MARGIN_TOP + HEADER_HEIGHT

    for row_idx, rec in enumerate(substrates):
        y_top = start_y + row_idx * (ROW_HEIGHT + ROW_GAP)
        y_center = y_top + ROW_HEIGHT / 2

        label = rec["label"]
        smiles_disp = rec["smiles_display"]

        parts.append(
            text_label(
                int(MARGIN_X + LABEL_COL_W / 2),
                int(y_center),
                label,
                anchor="middle",
                size=12,
            )
        )

        sub_x = MARGIN_X + LABEL_COL_W
        parts.append(
            mol_to_group(
                smiles_disp,
                int(sub_x),
                int(y_top),
                int(SUBSTRATE_COL_W),
                int(ROW_HEIGHT),
            )
        )

        arrow_x = MARGIN_X + LABEL_COL_W + SUBSTRATE_COL_W
        arrow_y = y_center
        parts.append(
            f'<line x1="{arrow_x + 10}" y1="{arrow_y}" '
            f'x2="{arrow_x + ARROW_COL_W - 10}" y2="{arrow_y}" '
            f'stroke="black" stroke-width="2"/>'
        )
        parts.append(
            f'<polygon points="'
            f'{arrow_x + ARROW_COL_W - 10},{arrow_y} '
            f'{arrow_x + ARROW_COL_W - 18},{arrow_y - 4} '
            f'{arrow_x + ARROW_COL_W - 18},{arrow_y + 4}" '
            f'fill="black"/>'
        )

        # Operator panels
        cur_x = op_x0
        for lab in op_labels:
            prod_list = rec.get(f"{lab}_products", [])
            width = op_widths[lab]

            if prod_list:
                slot_w = SUBSTRATE_COL_W
                for j, smi in enumerate(prod_list):
                    x = cur_x + j * slot_w
                    parts.append(
                        mol_to_group(
                            smi,
                            int(x),
                            int(y_top),
                            int(slot_w),
                            int(ROW_HEIGHT),
                        )
                    )
            else:
                parts.append(
                    text_label(
                        int(cur_x + width / 2),
                        int(y_center),
                        "no match",
                        anchor="middle",
                        size=10,
                    )
                )

            cur_x += width

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{int(total_w)}" height="{int(total_h)}" '
        f'viewBox="0 0 {int(total_w)} {int(total_h)}">\n'
        + "\n".join(parts)
        + "\n</svg>\n"
    )
    return svg


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    if not PROJECT_PATH.exists():
        raise SystemExit(f"Missing project file: {PROJECT_PATH}")

    rows = load_aryl_bromination_records(PROJECT_PATH)

    rdchiral_smirks = extract_operator(rows, OP_FIELDS["RdChiral"])
    op_smirks = {"EVODEX-E": EVODEX_E_SMIRKS, "RdChiral": rdchiral_smirks}

    substrates = build_substrate_panel(rows)

    # EVODEX-E via evodex.projection (adds Hs internally)
    run_evodex_projection(EVODEX_E_SMIRKS, substrates, "EVODEX-E_products")
    # RdChiral via rdchiral
    run_rdchiral_projection(rdchiral_smirks, substrates, "RdChiral_products")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    svg = build_svg_table(substrates, op_smirks)

    with OUTPUT_SVG.open("w", encoding="utf-8") as f:
        f.write(svg)

    print(f"Wrote {OUTPUT_SVG}")


if __name__ == "__main__":
    main()