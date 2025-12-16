#!/usr/bin/env python3
# analysis/fig3_svg_generator.py

from pathlib import Path

from indigo import Indigo
from indigo.renderer import IndigoRenderer


# ----------------------------
# Inputs
# ----------------------------

FULL_SMIRKS_WITHH = "[C:1]([H:5])([H:6])([H:7])[C:2]([H:8])([H:9])[C:3]([H:10])([H:11])[O:4][H]>>[C:1]([H:5])([H:6])([H:7])[C:2]([H:8])([H:9])[C:3]([H:10])([H:11])[O:4]C(=O)C([H])([H])C([H])([H])C([H])([H])C([H])([H])([H])"

FULL_FOR_OPERATOR = FULL_SMIRKS_WITHH


# ----------------------------
# Paths
# ----------------------------

BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = (BASE_DIR / "../structure_svgs").resolve()
OUT_DIR.mkdir(parents=True, exist_ok=True)


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
    print(f"Rendered: {out_svg}")


# ----------------------------
# Operator generation: E vs Em
# ----------------------------

def make_E_and_Em(full_smirks: str) -> tuple[str, str]:
    from evodex.operators import extract_operator_by_abstraction

    E = extract_operator_by_abstraction(full_smirks, "E", matched=False)
    Em = extract_operator_by_abstraction(full_smirks, "E", matched=True)
    return E, Em


def main() -> None:
    indigo, renderer = _make_indigo()

    # Full reaction (same renderer as operators)
    render_indigo_reaction_svg(indigo, renderer, FULL_SMIRKS_WITHH, OUT_DIR / "fig3_full_indigo.svg")

    # Operators
    E, Em = make_E_and_Em(FULL_FOR_OPERATOR)

    (OUT_DIR / "fig3_E.smirks.txt").write_text(E + "\n", encoding="utf-8")
    (OUT_DIR / "fig3_Em.smirks.txt").write_text(Em + "\n", encoding="utf-8")

    render_indigo_reaction_svg(indigo, renderer, E,  OUT_DIR / "fig3_E_indigo.svg")
    render_indigo_reaction_svg(indigo, renderer, Em, OUT_DIR / "fig3_Em_indigo.svg")


if __name__ == "__main__":
    main()