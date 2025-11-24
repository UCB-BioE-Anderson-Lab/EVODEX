import os
import html
import traceback

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdChemReactions as Reactions


# ------------------------------------------------------------
# Low-level drawing utilities (RDKit-only)
# ------------------------------------------------------------

def _draw_mol_svg(pattern: str, w: int, h: int) -> str:
    """
    Render a substrate-side pattern (SMARTS / mapped SMIRKS or plain SMILES)
    to an SVG string using RDKit only.

    This function first tries to parse the pattern as SMARTS (including
    rdchiral-style decorators) via RDKit's MolFromSmarts. If that fails,
    it falls back to parsing as plain SMILES. An empty or whitespace-only
    pattern raises a ValueError.
    """
    if not pattern or pattern.strip() == "":
        raise ValueError("SMARTS/SMILES render error: empty pattern")

    mol = Chem.MolFromSmarts(pattern)
    if mol is None:
        mol = Chem.MolFromSmiles(pattern)

    if mol is None:
        raise ValueError(f"SMARTS/SMILES render error: could not parse pattern: {pattern}")

    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _draw_rxn_svg(smirks: str, w: int, h: int) -> str:
    """
    Render a full reaction SMIRKS/SMARTS to SVG via RDKit.
    """
    if not smirks or smirks.strip() == "":
        raise ValueError("Reaction render error: empty SMIRKS")

    rxn = Reactions.ReactionFromSmarts(smirks)
    if rxn is None:
        raise ValueError(f"Reaction render error: could not parse SMIRKS: {smirks}")

    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.DrawReaction(rxn)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _strip_xml_header(svg_text: str) -> str:
    """
    Remove any XML declaration from the start of the SVG text.
    """
    import re
    return re.sub(r'^\s*<\?xml[^>]*\?>\s*', '', svg_text)


def _svg_as_group(svg_text: str, x: int, y: int) -> str:
    """
    Wrap an RDKit <svg>…</svg> block as a <g transform="translate(x,y)">…</g>,
    extracting the inner contents to avoid nested <svg> issues.
    """
    import re

    s = _strip_xml_header(svg_text).strip()
    m_open = re.search(r'<svg\b[^>]*>', s, flags=re.I)
    m_close = re.search(r'</svg\s*>', s, flags=re.I)

    if not (m_open and m_close and m_close.start() > m_open.end()):
        inner = s  # fallback
    else:
        inner = s[m_open.end():m_close.start()]

    return f'<g transform="translate({x},{y})">{inner}</g>'


def _error_panel(x: int, y: int, w: int, h: int, msg: str = "Error") -> str:
    """
    Simple light-red rounded panel with a centered error label.
    Only used when drawing fails.
    """
    return (
        f'<g transform="translate({x},{y})">'
        f'<rect x="0" y="0" width="{w}" height="{h}" '
        f'rx="8" ry="8" fill="#ffe6e6" stroke="#cc0000" stroke-width="1"/>'
        f'<text x="{w/2}" y="{h/2}" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#cc0000">'
        f'{html.escape(msg)}</text>'
        f'</g>'
    )


def _safe_mol_group(pattern: str | None, x: int, y: int, w: int, h: int) -> str:
    """
    Wrap a molecule drawing (pattern may be SMARTS/SMILES) in a <g> block.
    On failure, emits an error panel and prints a message to stdout.
    """
    if not pattern:
        return _error_panel(x, y, w, h, "No pattern")

    try:
        svg = _draw_mol_svg(pattern, w, h)
        return _svg_as_group(svg, x, y)
    except Exception as e:
        print(f"[5_graphics] _safe_mol_group error for pattern={pattern!r}: {e}")
        return _error_panel(x, y, w, h, "Parse error")


def _safe_rxn_group(smirks: str | None, x: int, y: int, w: int, h: int) -> str:
    """
    Wrap a reaction drawing in a <g> block.
    On failure, emits an error panel and prints a message to stdout.
    """
    if not smirks:
        return _error_panel(x, y, w, h, "No SMIRKS")

    try:
        svg = _draw_rxn_svg(smirks, w, h)
        return _svg_as_group(svg, x, y)
    except Exception as e:
        print(f"[5_graphics] _safe_rxn_group error for smirks={smirks!r}: {e}")
        return _error_panel(x, y, w, h, "Parse error")


# ------------------------------------------------------------
# Table SVG builder
# ------------------------------------------------------------

print(
    "[5_graphics] OUTPUT_DIR:",
    os.path.abspath("analysis/reactivity_series/out/5_graphics"),
)


def _build_table_svg(table_id: str, rows: list[dict]) -> str:
    print(f"[5_graphics] _build_table_svg called for {table_id} ({len(rows)} rows)")
    """
    Build a single SVG string for one table.

    Columns: Substrate name, Substrate structure, Value, A, B, C, D, E, RdChiral.
    For each abstraction column we render only the left hand side of the
    stored SMIRKS (substrate side).
    """

    def _truncate(text: str, max_chars: int = 60) -> str:
        if text is None:
            return ""
        return text if len(text) <= max_chars else text[: max_chars - 1] + "…"

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

    # Compute total size
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
            )
        )
    else:
        print(f"[5_graphics] No exemplar reaction found for {table_id}")

    # Caption of full exemplar string below header reaction
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
    for label, x in zip(op_cols, x_ops):
        content_parts.append(
            f'<text x="{x}" y="{col_header_y}" font-size="12" '
            f'font-weight="bold" font-family="Helvetica,Arial,sans-serif">{label}</text>'
        )

    # Rows
    y0 = col_header_y + 20
    draw_h = cell_h - 26

    for idx, rec in enumerate(rows):
        print(f"[5_graphics]   row {idx} (table {table_id})")
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

        # Substrate structure (SMILES/SMARTS)
        sub_pattern = (
            rec.get("mapped_substrate_smiles")
            or rec.get("smiles")
            or rec.get("smiles_input")
        )
        content_parts.append(
            _safe_mol_group(sub_pattern, int(x_struct), int(row_y), int(struct_col_w - 10), int(draw_h))
        )
        if sub_pattern:
            struct_cx = x_struct + (struct_col_w - 10) / 2
            caption = _truncate(sub_pattern)
            content_parts.append(
                f'<text x="{struct_cx}" y="{row_y + draw_h + 14}" font-size="8" '
                f'text-anchor="middle" font-family="ui-monospace,Menlo,Consolas,monospace">'
                f"{html.escape(caption)}</text>"
            )

        # Reactivity value: first non-substrate key in data
        value = None
        data_dict = rec.get("data", {})
        if data_dict:
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

        # Abstractions: A..E, RdChiral (substrate side only)
        def left_side(field: str) -> str | None:
            smirks = rec.get(field)
            if not smirks:
                if field == "rdchiral_smirks" and rec.get("rdchiral_error"):
                    print(
                        f"[5_graphics] RdChiral error for {table_id} row "
                        f"{rec.get('row')}: {rec.get('rdchiral_error')}"
                    )
                else:
                    print(
                        f"[5_graphics] Missing {field} for {table_id} row "
                        f"{rec.get('row')}"
                    )
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

        for label, x in zip(op_cols, x_ops):
            pattern = left_side(field_map[label])
            content_parts.append(
                _safe_mol_group(pattern, int(x), int(row_y), int(op_col_w - 10), int(draw_h))
            )
            if pattern:
                op_cx = x + (op_col_w - 10) / 2
                caption = _truncate(pattern)
                content_parts.append(
                    f'<text x="{op_cx}" y="{row_y + draw_h + 14}" font-size="8" '
                    f'text-anchor="middle" font-family="ui-monospace,Menlo,Consolas,monospace">'
                    f"{html.escape(caption)}</text>"
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


# ------------------------------------------------------------
# Main: read JSONL and write SVGs
# ------------------------------------------------------------

def main():
    import json
    from glob import glob

    OUTPUT_DIR = os.path.abspath("analysis/reactivity_series/out/5_graphics")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    input_dir = "analysis/reactivity_series/out/4_ops"
    paths = sorted(glob(os.path.join(input_dir, "*.jsonl")))
    if not paths:
        print(f"[5_graphics] No input files found in {input_dir}")
        return

    for path in paths:
        table_id = os.path.splitext(os.path.basename(path))[0]
        print(f"[5_graphics] Processing {table_id} from {path}")
        rows: list[dict] = []

        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rows.append(json.loads(line))

        print(f"[5_graphics]  - {len(rows)} rows")
        try:
            svg = _build_table_svg(table_id, rows)
        except Exception:
            print(f"[5_graphics] ERROR building SVG for {table_id}")
            traceback.print_exc()
            continue

        out_path = os.path.join(OUTPUT_DIR, f"{table_id}.svg")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(svg)
        print(f"[5_graphics]  - wrote {out_path}")


if __name__ == "__main__":
    main()