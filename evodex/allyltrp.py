"""
AllylTrp – RDKit → SVG renderer for σ/π visualization.

- Big atom discs (π-capable = orange, others = light blue)
- σ bonds = opaque blue strokes with round caps and a darker outline (underlay)
- π system = translucent yellow strokes (aromatic rings as rounded polylines; isolated C=C as single strokes)
- Correct bond endpoints: each stroke shortened so the *outer edge* of the round cap passes through the atom center.

Dependencies:
  pip install rdkit-pypi svgwrite  # (optional) cairosvg for standalone PNG export

Notebook usage: set SMILES below and run the cell.
"""
from __future__ import annotations

import math
from typing import List, Tuple

import svgwrite
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import GetSymmSSSR
from rdkit.Chem.Draw import rdMolDraw2D


# =================== Styling & Geometry ===================
COLORS = {
    "atom_sp3":   "#A8D8E8",   # discs for non-π-capable centers
    "atom_sp2":   "#FFD27A",   # discs for π-capable centers
    # σ colors derived by scaling atom_sp3 toward black (single-axis RGB):
    "sigma_fill": "#688690",   # ≈ atom_sp3 * 0.62
    "sigma_edge": "#4C6168",   # ≈ atom_sp3 * 0.45 (used only if outline enabled)
    "pi_fill":    "#D99900",   # π strokes (aromatic + isolated C=C)
}

ALPHA_PI = 0.65  # Only π uses transparency

GEOM = {
    "margin": 30,                # px margin around drawing
    "atom_radius": 18,           # big atom discs (visual)
    "sigma_scale": 0.55,         # σ stroke thickness = 2 * atom_radius * sigma_scale (slimmer)
    "sigma_edge_px": 6,          # extra pixels added to σ width for dark outline underlay
    "sigma_gap_px": 5,           # visible gap between atom disc edge and σ outer cap edge
    "sigma_tip_inset_px": 8,           # extra pull-back for σ endpoints to avoid overlap
    "sigma_tip_len_px": 14,          # length of tapered σ tips
    "sigma_tip_taper_frac": 0.35,      # taper sharpness fraction (0–1; smaller = pointier)
    "sigma_outline": False,     # draw darker σ underlay (False removes outlines)
    "pi_scale": 1.30,            # π stroke thickness = 2 * atom_radius * pi_scale
    "rdkit_panel": True,         # show raw RDKit depiction panel on the right
    "rdkit_panel_frac": 0.38,    # fraction of total width reserved for RDKit panel
    "panel_gutter_px": 16,       # space between left render and RDKit panel
}

# =================== Chemistry helpers ===================

def is_pi_capable(atom: Chem.Atom) -> bool:
    """Atom considered π-capable if any incident bond is not single."""
    for b in atom.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE:
            return True
    return False


def aromatic_rings(mol: Chem.Mol) -> List[List[int]]:
    rings: List[List[int]] = []
    for r in GetSymmSSSR(mol):
        idxs = list(r)
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in idxs):
            rings.append(idxs)
    return rings

# =================== Geometry helpers ===================

def compute_canvas_map(coords: List[Tuple[float, float]], size: Tuple[int, int]=(950, 520), margin: int=30):
    W, H = size
    xs = [x for x, _ in coords]
    ys = [y for _, y in coords]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    sx = (W - 2 * margin) / max(maxx - minx, 1e-6)
    sy = (H - 2 * margin) / max(maxy - miny, 1e-6)
    scale = min(sx, sy)

    def to_px(p: Tuple[float, float]) -> Tuple[float, float]:
        x = margin + (p[0] - minx) * scale
        y = margin + (p[1] - miny) * scale
        return (x, H - y)  # SVG y-down

    return to_px


def ring_points_px(mol: Chem.Mol, ring_idxs: List[int], px_coords: List[Tuple[float, float]]):
    pts = [px_coords[i] for i in ring_idxs]
    pts.append(pts[0])  # close path
    return pts


def path_from_points(points: List[Tuple[float, float]]) -> str:
    return " ".join(("M" if k == 0 else "L") + f"{x:.2f},{y:.2f}" for k, (x, y) in enumerate(points))


def shorten_segment(p1: Tuple[float, float], p2: Tuple[float, float], cap_r: float):
    """Move each endpoint toward the other by cap_r so the outer cap edge passes through atom center."""
    x1, y1 = p1
    x2, y2 = p2
    dx, dy = (x2 - x1), (y2 - y1)
    L = math.hypot(dx, dy)
    if L < 1e-9 or 2 * cap_r >= L:
        # degenerate: collapse to midpoint
        mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
        return (mx, my), (mx, my)
    ux, uy = dx / L, dy / L
    return (x1 + ux * cap_r, y1 + uy * cap_r), (x2 - ux * cap_r, y2 - uy * cap_r)

def tangent_endpoints(pA, pB, cut_r, cap_r):
    """Return endpoints S,E along AB so a round cap of radius cap_r is tangent to a circle of radius `cut_r` (target outer radius) at each atom center.
    If the bond is too short, collapse to the midpoint.
    """
    x1, y1 = pA; x2, y2 = pB
    dx, dy = (x2 - x1), (y2 - y1)
    L = math.hypot(dx, dy)
    if L < 1e-9:
        return pA, pA
    ux, uy = dx / L, dy / L
    base = max(cut_r - cap_r, 0.0)
    if 2 * base >= L:  # degenerate
        mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
        return (mx, my), (mx, my)
    S = (x1 + ux * base, y1 + uy * base)
    E = (x2 - ux * base, y2 - uy * base)
    return S, E

def tapered_segment_points(s: Tuple[float, float], e: Tuple[float, float], width: float, tip_len: float, taper_frac: float=0.5, tip_hw_abs: float=None):
    """Return polygon points for a tapered bar from s to e.
    The central bar has width `width` and flat sides; each end is an isosceles triangle of length `tip_len`.
    Points are ordered to form a single closed polygon suitable for dwg.polygon().
    """
    x1, y1 = s; x2, y2 = e
    dx, dy = (x2 - x1), (y2 - y1)
    L = math.hypot(dx, dy)
    if L < 1e-9:
        return [s]
    ux, uy = dx / L, dy / L
    # perpendicular unit
    nx, ny = -uy, ux
    hw = width / 2.0
    if tip_hw_abs is not None:
        tip_hw = min(tip_hw_abs, hw)
    else:
        tip_hw = hw * taper_frac
    s1 = (x1 + ux * tip_len, y1 + uy * tip_len)
    e1 = (x2 - ux * tip_len, y2 - uy * tip_len)
    pts = [
        (s1[0] + nx * hw, s1[1] + ny * hw),
        (e1[0] + nx * hw, e1[1] + ny * hw),
        (x2 + nx * tip_hw, y2 + ny * tip_hw),
        (x2 - nx * tip_hw, y2 - ny * tip_hw),
        (e1[0] - nx * hw, e1[1] - ny * hw),
        (s1[0] - nx * hw, s1[1] - ny * hw),
        (x1 - nx * tip_hw, y1 - ny * tip_hw),
        (x1 + nx * tip_hw, y1 + ny * tip_hw),
    ]
    return [(round(px, 1), round(py, 1)) for (px, py) in pts]

# =================== Core renderer ===================

def draw_circles_and_bonds_svg(
    smiles: str,
    size: Tuple[int, int]=(1000, 600),
    out_path: str="molecule.svg",
    export_png: bool=False,
    png_path: str="molecule.png",
    png_dpi: int=384,
):
    """Render a vector illustration of σ/π systems.

    - Big atom discs on bottom (colored by π-capability).
    - Opaque σ bonds as single round-capped strokes with a darker outline underlay.
    - Translucent π system (aromatic rings as single rounded polylines; isolated C=C as strokes).
    - Writes an SVG; optionally exports a high-DPI PNG (requires cairosvg).
    """
    rdkit_group_markup = None  # will hold inline SVG for the right-hand RDKit panel
    m = Chem.MolFromSmiles(smiles)
    if not m:
        raise ValueError("Invalid SMILES")
    m = Chem.AddHs(m)
    AllChem.Compute2DCoords(m)

    conf = m.GetConformer()
    coords = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(m.GetNumAtoms())]

    # Split the canvas: draw our custom rendering on the left; RDKit panel on the right
    W, H = size
    if GEOM.get("rdkit_panel", True):
        right_frac = GEOM["rdkit_panel_frac"]
        gutter = GEOM["panel_gutter_px"]
        left_W = int(W * (1.0 - right_frac) - gutter)
        left_W = max(left_W, 200)
    else:
        left_W = W

    to_px = compute_canvas_map(coords, size=(left_W, H), margin=GEOM["margin"])
    px_coords = [to_px(p) for p in coords]

    # Sizes
    R_atom = GEOM["atom_radius"]
    sigma_w = 2.0 * R_atom * GEOM["sigma_scale"]
    pi_w    = 2.0 * R_atom * GEOM["pi_scale"]
    sigma_cap_r = sigma_w / 2.0
    pi_cap_r    = pi_w / 2.0

    # --- SVG doc ---
    dwg = svgwrite.Drawing(out_path, size=(W, H), profile="full")
    dwg.add(dwg.rect(insert=(0, 0), size=(W, H), fill="white"))

    # Layer 1: atom discs
    g_atoms = dwg.g(id="atoms")
    for i, atom in enumerate(m.GetAtoms()):
        cx, cy = px_coords[i]
        fill = COLORS["atom_sp2"] if is_pi_capable(atom) else COLORS["atom_sp3"]
        g_atoms.add(dwg.circle(center=(cx, cy), r=R_atom, fill=fill))
    dwg.add(g_atoms)

    # Layer 2: σ bonds (edge underlay then fill overlay), endpoints tangent to an inner atom circle
    g_sigma = dwg.g(id="sigma")
    edge_w = GEOM["sigma_edge_px"]
    target_outer_r = max(R_atom - GEOM["sigma_gap_px"] - GEOM["sigma_tip_inset_px"], 0)  # pull σ endpoints further back to avoid overlap
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        s, e = tangent_endpoints(px_coords[i], px_coords[j], target_outer_r, sigma_cap_r)
        tip_len = GEOM["sigma_tip_len_px"]
        tip_hw_consistent = (sigma_w * 0.5) * GEOM["sigma_tip_taper_frac"]
        if GEOM.get("sigma_outline", True):
            # outline underlay (slightly thicker) as a tapered polygon
            poly_edge_pts = tapered_segment_points(s, e, sigma_w + 2 * edge_w, tip_len, GEOM["sigma_tip_taper_frac"], tip_hw_abs=tip_hw_consistent)
            g_sigma.add(
                dwg.polygon(points=poly_edge_pts, fill=COLORS["sigma_edge"], shape_rendering="geometricPrecision")
            )
        # fill overlay as a tapered polygon
        poly_fill_pts = tapered_segment_points(s, e, sigma_w, tip_len, GEOM["sigma_tip_taper_frac"], tip_hw_abs=tip_hw_consistent)
        g_sigma.add(
            dwg.polygon(points=poly_fill_pts, fill=COLORS["sigma_fill"], shape_rendering="geometricPrecision")
        )
    dwg.add(g_sigma)

    # Layer 3: π system (top) — stroke-only with alpha
    g_pi = dwg.g(id="pi", opacity=ALPHA_PI)

    # Aromatic rings as single rounded polylines
    for ring in aromatic_rings(m):
        pts = ring_points_px(m, ring, px_coords)
        d = path_from_points(pts)
        g_pi.add(
            dwg.path(
                d=d,
                fill="none",
                stroke=COLORS["pi_fill"], stroke_width=pi_w,
                stroke_linecap="round", stroke_linejoin="round",
                shape_rendering="geometricPrecision", stroke_miterlimit=1,
            )
        )

    # Non-aromatic C=C bonds as single round-capped strokes (endpoints shortened)
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE and not b.GetIsAromatic():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            s, e = shorten_segment(px_coords[i], px_coords[j], pi_cap_r * 1.02)
            g_pi.add(
                dwg.line(
                    start=s, end=e,
                    stroke=COLORS["pi_fill"], stroke_width=pi_w,
                    stroke_linecap="round", stroke_linejoin="round",
                    shape_rendering="geometricPrecision", stroke_miterlimit=1,
                )
            )

    dwg.add(g_pi)

    # ----- RDKit reference panel (right) -----
    if GEOM.get("rdkit_panel", True):
        right_frac = GEOM["rdkit_panel_frac"]
        gutter = GEOM["panel_gutter_px"]
        left_W = int(W * (1.0 - right_frac) - gutter)
        left_W = max(left_W, 200)
        right_W = W - left_W - gutter
        panel_x = left_W + gutter

        # Panel background
        dwg.add(dwg.rect(insert=(panel_x, 0), size=(right_W, H), fill="white"))

        # Render RDKit depiction to inline SVG and embed as raw XML
        draw_w = max(right_W, 200)
        draw_h = max(H, 200)
        d2d = rdMolDraw2D.MolDraw2DSVG(draw_w, draw_h)
        rdMolDraw2D.PrepareAndDrawMolecule(d2d, m)
        d2d.FinishDrawing()
        svg_panel = d2d.GetDrawingText()
        # Strip XML header/doctype and outer <svg> wrapper, then translate into the right panel
        lines = [ln for ln in svg_panel.splitlines() if not ln.startswith("<?xml") and "DOCTYPE" not in ln]
        svg_panel_clean = "\n".join(lines)
        start = svg_panel_clean.find("<svg")
        if start != -1:
            gt = svg_panel_clean.find(">", start)
            inner = svg_panel_clean[gt+1: svg_panel_clean.rfind("</svg>")]
        else:
            inner = svg_panel_clean
        rdkit_group_markup = f'<g transform="translate({panel_x},0)">{inner}</g>'

        # Optional divider line
        dwg.add(
            dwg.line(start=(panel_x - gutter * 0.5, 0), end=(panel_x - gutter * 0.5, H),
                     stroke="#cccccc", stroke_width=1)
        )

    dwg.save()

    # Post-process: inject RDKit inline SVG group at the end of root
    if GEOM.get("rdkit_panel", True) and rdkit_group_markup:
        try:
            with open(out_path, "r", encoding="utf-8") as f:
                svg_text = f.read()
            insert_at = svg_text.rfind("</svg>")
            if insert_at != -1:
                svg_text = svg_text[:insert_at] + rdkit_group_markup + svg_text[insert_at:]
                with open(out_path, "w", encoding="utf-8") as f:
                    f.write(svg_text)
        except Exception as e:
            print("Warning: failed to inline RDKit SVG panel:", e)

    # Optional PNG export
    if export_png:
        try:
            import cairosvg
            cairosvg.svg2png(url=out_path, write_to=png_path, dpi=png_dpi)
        except Exception as e:
            print("PNG export skipped (install cairosvg):", e)

    return out_path


# =================== Notebook invocation (KISS) ===================
# Set the SMILES you want to render and run this cell in Colab/Jupyter.
SMILES = "CC(=CCc1cccc2[nH]cc(C[C@H](N)C(=O)O)c12)C"  # <-- change as needed

svg_path = draw_circles_and_bonds_svg(
    SMILES,
    size=(1000, 600),
    out_path="molecule.svg",
)
print("Wrote:", svg_path)

# Optional: trigger a download in Colab (safe to leave on; ignored outside Colab)
try:
    from google.colab import files  # type: ignore
    files.download(svg_path)
except Exception:
    pass