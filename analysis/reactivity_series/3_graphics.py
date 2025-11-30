#!/usr/bin/env python3
"""
Generate SVG summary tables for reactivity series data.

Inputs
------
- JSONL files in analysis/reactivity_series/out/2_ops
  Each line is a JSON record with fields including:
    - row: zero-based row index
    - data: dict with at least a "substrate" key and possibly "value"
    - evodex_*_smirks, rdchiral_smirks: SMIRKS / SMARTS strings

- CSV files in analysis/reactivity_series/data/tables
  Column C (zero-based index 2) is used as the "value" column when present.

Outputs
-------
- One SVG per input JSONL file in analysis/reactivity_series/out/3_graphics
- A JSONL error log in out/_logs/3_graphics.errors.jsonl

All layout constants and drawing behavior are fixed to preserve the original figures.
"""

import html
import json
import os
import re
import traceback
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem.Draw import rdMolDraw2D

# ----------------------------------------------------------------------
# Output and logging configuration
# ----------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
LOG_DIR = BASE_DIR / "out/_logs"
TABLES_DIR = BASE_DIR / "data/tables"


def log_error(
    fh,
    *,
    table_id,
    stage,
    error,
    row_idx=None,
    label=None,
    extra=None,
) -> None:
    """
    Write a single JSON error record to the given file handle.

    Parameters
    ----------
    fh
        Open text file handle for JSONL output.
    table_id
        Table identifier, usually derived from the JSONL filename.
    stage
        Name of the processing stage, for example "graphics_substrate".
    error
        Exception instance or message string.
    row_idx
        Optional row index associated with the error.
    label
        Optional human readable label for the row.
    extra
        Optional dict with additional context (for example a SMIRKS string).
    """
    rec: dict[str, object] = {
        "table_id": table_id,
        "stage": stage,
    }

    if isinstance(error, Exception):
        rec["error_type"] = type(error).__name__
        rec["error_message"] = str(error)
    else:
        rec["error"] = str(error)

    if row_idx is not None:
        rec["row"] = int(row_idx)
    if label is not None:
        rec["label"] = str(label)
    if extra:
        rec["extra"] = extra

    if isinstance(error, Exception):
        rec["traceback"] = traceback.format_exc()

    fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
    fh.flush()


# ----------------------------------------------------------------------
# Low level drawing helpers
# ----------------------------------------------------------------------

def _draw_mol_svg(pattern: str, w: int, h: int) -> str:
    """
    Render a substrate side SMARTS or SMILES pattern to an SVG string.

    The pattern is parsed as SMARTS first, then as SMILES if SMARTS parsing fails.
    An empty or whitespace only pattern raises ValueError.
    """
    if not pattern or pattern.strip() == "":
        raise ValueError("SMARTS/SMILES render error: empty pattern")

    mol = Chem.MolFromSmarts(pattern)
    if mol is None:
        mol = Chem.MolFromSmiles(pattern)
    if mol is None:
        raise ValueError(
            f"SMARTS/SMILES render error: could not parse pattern: {pattern}"
        )

    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _draw_rdchiral_mol_svg(pattern: str, w: int, h: int) -> str:
    """
    Render an RdChiral substrate side SMARTS pattern with explicit labels.

    Atom labels preserve:
    - aromaticity (for example "c" for aromatic carbon)
    - explicit hydrogen decorators (for example "cH", "nH", "CH")
    - atom map numbers (for example "c:5", "nH:4").
    """
    if not pattern or pattern.strip() == "":
        raise ValueError("SMARTS/SMILES render error: empty pattern (RdChiral)")

    mol = Chem.MolFromSmarts(pattern)
    if mol is None:
        mol = Chem.MolFromSmiles(pattern)
    if mol is None:
        raise ValueError(
            f"SMARTS/SMILES render error: could not parse pattern: {pattern}"
        )

    # Map number -> base symbol with explicit H decorator, for example "cH" or "nH".
    H_decorated: dict[int, str] = {}
    try:
        for match in re.finditer(r"\[([A-Za-z][a-z]?)H[^:]*:(\d+)\]", pattern):
            base_sym = match.group(1)
            mapnum = int(match.group(2))
            H_decorated[mapnum] = base_sym + "H"
    except Exception:
        H_decorated = {}

    # Work on a copy so we do not mutate any cached templates.
    mol = Chem.Mol(mol)

    for atom in mol.GetAtoms():
        atom_num = atom.GetAtomicNum()
        amap = atom.GetAtomMapNum()
        is_arom = atom.GetIsAromatic()

        if amap and amap in H_decorated:
            base = H_decorated[amap]
        elif atom_num == 6 and is_arom:
            base = "c"
        elif atom_num == 1:
            base = "H"
        else:
            base = atom.GetSymbol()

        label = base
        if amap:
            label = f"{base}:{amap}"

        atom.SetProp("atomLabel", label)

    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _draw_rxn_svg(smirks: str, w: int, h: int) -> str:
    """
    Render a full reaction SMIRKS or SMARTS string to an SVG string.
    """
    if not smirks or smirks.strip() == "":
        raise ValueError("Reaction render error: empty SMIRKS")

    rxn = Reactions.ReactionFromSmarts(smirks)
    if rxn is None:
        raise ValueError(
            f"Reaction render error: could not parse SMIRKS: {smirks}"
        )

    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.DrawReaction(rxn)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _strip_xml_header(svg_text: str) -> str:
    """
    Remove any XML declaration at the start of an SVG string.
    """
    return re.sub(r"^\s*<\?xml[^>]*\?>\s*", "", svg_text)


def _svg_as_group(svg_text: str, x: int, y: int) -> str:
    """
    Wrap an RDKit generated <svg> block in a translated <g> element.

    The inner content between <svg> and </svg> is extracted to avoid nested
    <svg> elements.
    """
    s = _strip_xml_header(svg_text).strip()
    m_open = re.search(r"<svg\b[^>]*>", s, flags=re.I)
    m_close = re.search(r"</svg\s*>", s, flags=re.I)

    if not (m_open and m_close and m_close.start() > m_open.end()):
        inner = s
    else:
        inner = s[m_open.end() : m_close.start()]

    return f'<g transform="translate({x},{y})">{inner}</g>'


def _error_panel(x: int, y: int, w: int, h: int, msg: str = "Error") -> str:
    """
    Draw a light red rounded rectangle with a centered error label.
    """
    return (
        f'<g transform="translate({x},{y})">'
        f'<rect x="0" y="0" width="{w}" height="{h}" '
        f'rx="8" ry="8" fill="#ffe6e6" stroke="#cc0000" stroke-width="1"/>'
        f'<text x="{w/2}" y="{h/2}" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#cc0000">'
        f"{html.escape(msg)}</text>"
        f"</g>"
    )


def _safe_mol_group(
    pattern: str | None,
    x: int,
    y: int,
    w: int,
    h: int,
    *,
    err_f=None,
    table_id: str | None = None,
    stage: str = "graphics_mol",
    row_idx: int | None = None,
    label: str | None = None,
    field: str | None = None,
) -> str:
    """
    Draw a molecule pattern in a translated group.

    On success, returns a <g> with drawn content.
    On failure or missing pattern, returns an error panel and logs the error
    if an error file handle is provided.
    """
    if not pattern:
        return _error_panel(x, y, w, h, "No pattern")

    try:
        if field == "rdchiral_smirks":
            svg = _draw_rdchiral_mol_svg(pattern, w, h)
        else:
            svg = _draw_mol_svg(pattern, w, h)
        return _svg_as_group(svg, x, y)
    except Exception as exc:
        if err_f and table_id:
            extra = {"pattern": pattern}
            if field:
                extra["field"] = field
            log_error(
                err_f,
                table_id=table_id,
                stage=stage,
                error=exc,
                row_idx=row_idx,
                label=label,
                extra=extra,
            )
        else:
            print(
                f"[3_graphics] _safe_mol_group error for pattern={pattern!r}: {exc}"
            )
        return _error_panel(x, y, w, h, "Parse error")


def _safe_rxn_group(
    smirks: str | None,
    x: int,
    y: int,
    w: int,
    h: int,
    *,
    err_f=None,
    table_id: str | None = None,
    stage: str = "graphics_rxn",
) -> str:
    """
    Draw a reaction SMIRKS string in a translated group.

    On success, returns a <g> with drawn content.
    On failure or missing SMIRKS, returns an error panel and logs the error
    if an error file handle is provided.
    """
    if not smirks:
        return _error_panel(x, y, w, h, "No SMIRKS")

    try:
        svg = _draw_rxn_svg(smirks, w, h)
        return _svg_as_group(svg, x, y)
    except Exception as exc:
        if err_f and table_id:
            log_error(
                err_f,
                table_id=table_id,
                stage=stage,
                error=exc,
                extra={"smirks": smirks},
            )
        else:
            print(
                f"[3_graphics] _safe_rxn_group error for smirks={smirks!r}: {exc}"
            )
        return _error_panel(x, y, w, h, "Parse error")


# ----------------------------------------------------------------------
# Table SVG builder
# ----------------------------------------------------------------------

print(
    "[3_graphics] OUTPUT_DIR:",
    os.path.abspath("analysis/reactivity_series/out/3_graphics"),
)


def _truncate(text: str | None, max_chars: int = 60) -> str:
    """
    Truncate a string to at most max_chars characters with an ellipsis.
    """
    if text is None:
        return ""
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 1] + "…"


def _build_table_svg(table_id: str, rows: list[dict], err_f=None) -> str:
    """
    Build a single SVG string for one table.

    Columns
    -------
    Substrate name, Substrate structure, Value, A, B, C, D, E, RdChiral.

    For each abstraction column (A through E, RdChiral) only the substrate
    side of the SMIRKS string is rendered.
    """
    # Layout constants
    margin_x = 20
    margin_y = 20
    name_col_w = 180
    struct_col_w = 150
    value_col_w = 100
    op_col_w = 150
    op_cols = ["A", "B", "C", "D", "E", "RdChiral"]
    n_ops = len(op_cols)

    cell_h = 120
    row_gap = 30

    header_reaction_h = 160
    header_gap = 80

    # Horizontal positions
    x_name = margin_x
    x_struct = x_name + name_col_w
    x_value = x_struct + struct_col_w
    x_ops = [x_value + value_col_w + i * op_col_w for i in range(n_ops)]

    # Total size
    n_rows = len(rows)
    table_height = n_rows * (cell_h + row_gap)
    total_w = x_ops[-1] + op_col_w + margin_x
    total_h = margin_y + header_reaction_h + header_gap + table_height + margin_y

    # Choose an exemplar reaction for the header (full SMIRKS)
    exemplar = None
    exemplar_keys = [
        "evodex_p_smirks",
        "rdchiral_smirks",
        "evodex_e_smirks",
        "evodex_d_smirks",
        "evodex_c_smirks",
        "evodex_b_smirks",
        "evodex_a_smirks",
    ]
    for key in exemplar_keys:
        for rec in rows:
            if rec.get(key):
                exemplar = rec[key]
                break
        if exemplar:
            break

    content_parts: list[str] = []

    # Title
    title = f"Table: {table_id}"
    content_parts.append(
        f'<text x="{total_w/2}" y="{margin_y}" text-anchor="middle" '
        f'font-size="20" font-family="Helvetica,Arial,sans-serif">'
        f"{html.escape(title)}</text>"
    )

    # Header reaction
    header_y = margin_y + 10
    if exemplar:
        content_parts.append(
            _safe_rxn_group(
                exemplar,
                int(margin_x),
                int(header_y + 10),
                int(total_w - 2 * margin_x),
                int(header_reaction_h - 20),
                err_f=err_f,
                table_id=table_id,
                stage="graphics_header",
            )
        )

    # Caption for the exemplar SMIRKS
    if exemplar:
        caption = _truncate(exemplar)
        caption_x = total_w / 2
        caption_y = margin_y + header_reaction_h + 15
        content_parts.append(
            f'<text x="{caption_x}" y="{caption_y}" text-anchor="middle" '
            f'font-size="9" font-family="ui-monospace,Menlo,Consolas,monospace">'
            f"{html.escape(caption)}</text>"
        )

    # Column headers
    col_header_y = margin_y + header_reaction_h + 40
    content_parts.append(
        f'<text x="{x_name}" y="{col_header_y}" font-size="12" '
        f'font-weight="bold" font-family="Helvetica,Arial,sans-serif">Substrate</text>'
    )
    content_parts.append(
        f'<text x="{x_struct}" y="{col_header_y}" font-size="12" '
        f'font-weight="bold" font-family="Helvetica,Arial,sans-serif">Structure</text>'
    )
    content_parts.append(
        f'<text x="{x_value}" y="{col_header_y}" font-size="12" '
        f'font-weight="bold" font-family="Helvetica,Arial,sans-serif">Value</text>'
    )
    for op_label, x in zip(op_cols, x_ops):
        content_parts.append(
            f'<text x="{x}" y="{col_header_y}" font-size="12" '
            f'font-weight="bold" font-family="Helvetica,Arial,sans-serif">{op_label}</text>'
        )

    # Rows
    y0 = col_header_y + 20
    draw_h = cell_h - 26

    for idx, rec in enumerate(rows):
        row_y = y0 + idx * (cell_h + row_gap)
        label_y = row_y + cell_h / 2 + 4

        # Substrate label
        substrate_name = (
            rec.get("data", {}).get("substrate")
            or rec.get("substrate")
            or rec.get("label")
            or f"row {rec.get('row', idx)}"
        )
        content_parts.append(
            f'<text x="{x_name}" y="{label_y}" font-size="11" '
            f'font-family="Helvetica,Arial,sans-serif">'
            f"{html.escape(substrate_name)}</text>"
        )

        # Substrate structure pattern:
        # prefer raw input SMILES to avoid showing map numbers or explicit H.
        sub_pattern = (
            rec.get("smiles_input")
            or rec.get("smiles")
            or rec.get("mapped_substrate_smiles")
        )
        content_parts.append(
            _safe_mol_group(
                sub_pattern,
                int(x_struct),
                int(row_y),
                int(struct_col_w - 10),
                int(draw_h),
                err_f=err_f,
                table_id=table_id,
                stage="graphics_substrate",
                row_idx=rec.get("row", idx),
                label=substrate_name,
                field="substrate",
            )
        )

        # Reactivity value: use data["value"] when present, otherwise
        # first non "substrate" entry in the data dict.
        value = None
        data_dict = rec.get("data", {})
        if data_dict:
            if "value" in data_dict:
                value = data_dict["value"]
            else:
                for k, v in data_dict.items():
                    if k != "substrate":
                        value = v
                        break
        value_str = str(value) if value is not None else ""
        content_parts.append(
            f'<text x="{x_value}" y="{label_y}" font-size="12" '
            f'font-family="Helvetica,Arial,sans-serif">'
            f"{html.escape(value_str)}</text>"
        )

        # Abstractions: A..E and RdChiral, substrate side only.
        def left_side(field: str) -> str | None:
            smirks = rec.get(field)
            if not smirks:
                return None
            return smirks.split(">>", 1)[0]

        field_map = {
            "A": "evodex_a_smirks",
            "B": "evodex_b_smirks",
            "C": "evodex_c_smirks",
            "D": "evodex_d_smirks",
            "E": "evodex_e_smirks",
            "RdChiral": "rdchiral_smirks",
        }

        for op_label, x in zip(op_cols, x_ops):
            pattern = left_side(field_map[op_label])
            content_parts.append(
                _safe_mol_group(
                    pattern,
                    int(x),
                    int(row_y),
                    int(op_col_w - 10),
                    int(draw_h),
                    err_f=err_f,
                    table_id=table_id,
                    stage=f"graphics_{op_label}",
                    row_idx=rec.get("row", idx),
                    label=substrate_name,
                    field=field_map[op_label],
                )
            )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{int(total_w)}" height="{int(total_h)}" '
        f'viewBox="0 0 {int(total_w)} {int(total_h)}">\n'
        f'<rect x="0" y="0" width="{int(total_w)}" height="{int(total_h)}" fill="white"/>\n'
        + "\n".join(content_parts)
        + "\n</svg>\n"
    )
    return svg


# ----------------------------------------------------------------------
# Main: read JSONL and write SVGs
# ----------------------------------------------------------------------

def main() -> None:
    from glob import glob

    output_dir = os.path.abspath("analysis/reactivity_series/out/3_graphics")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    errors_path = LOG_DIR / "3_graphics.errors.jsonl"
    err_f = open(errors_path, "w", encoding="utf-8")

    input_dir = "analysis/reactivity_series/out/2_ops"
    paths = sorted(glob(os.path.join(input_dir, "*.jsonl")))
    if not paths:
        print(f"[3_graphics] No input files found in {input_dir}")
        return

    for path in paths:
        table_id = os.path.splitext(os.path.basename(path))[0]
        print(f"[3_graphics] Processing {table_id} from {path}")
        rows: list[dict] = []

        # Load JSONL rows from 2_ops
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rows.append(json.loads(line))

        # Attach values from the original CSV (column C) if available
        base_id = table_id.replace(".project", "")
        csv_path = TABLES_DIR / f"{base_id}.csv"
        value_by_row: dict[int, object] = {}

        if csv_path.exists():
            try:
                df = pd.read_csv(csv_path)
                if df.shape[1] >= 3:
                    for i in range(len(df)):
                        value_by_row[i] = df.iloc[i, 2]
            except Exception as exc:
                log_error(
                    err_f,
                    table_id=base_id,
                    stage="graphics_value_load",
                    error=exc,
                    extra={"csv_path": str(csv_path)},
                )

        # Push those values into each record data dict under key "value"
        if value_by_row:
            for rec in rows:
                idx = rec.get("row")
                if isinstance(idx, int) and idx in value_by_row:
                    v = value_by_row[idx]
                    data_dict = rec.get("data") or {}
                    if "data" not in rec:
                        rec["data"] = data_dict
                    if "value" not in data_dict:
                        data_dict["value"] = v

        print(f"[3_graphics]  - {len(rows)} rows")
        try:
            svg = _build_table_svg(table_id, rows, err_f=err_f)
        except Exception as exc:
            print(f"[3_graphics] ERROR building SVG for {table_id}")
            traceback.print_exc()
            log_error(
                err_f,
                table_id=table_id,
                stage="graphics_table",
                error=exc,
            )
            continue

        out_path = os.path.join(output_dir, f"{table_id}.svg")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(svg)
        print(f"[3_graphics]  - wrote {out_path}")


if __name__ == "__main__":
    main()