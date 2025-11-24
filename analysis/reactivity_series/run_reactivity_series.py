from rdkit import Chem
from rdkit.Chem import Reactions
from typing import Optional, Tuple
from pathlib import Path

def apply_template_first(template: str, reactant_mol: Chem.Mol) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    """Apply the reaction and return the first valid product only."""
    try:
        rxn = Reactions.ReactionFromSmarts(template)
        if rxn is None:
            return None, "ReactionFromSmarts returned None"
        outcomes = rxn.RunReactants((reactant_mol,))
        if not outcomes:
            return None, "no_products"
        # Take the first outcome, first product molecule
        prod = Chem.Mol(outcomes[0][0])
        prod, err = sanitize_lenient(prod)
        if prod is None:
            repaired, rerr = fix_overvalent_nitrogens(Chem.Mol(outcomes[0][0]))
            if repaired is None:
                return None, f"product_{err}; repair={rerr}"
            prod = repaired
        # Strip isotopes for presentation
        prod = clear_isotopes(prod)
        return prod, None
    except Exception as e:
        return None, f"reaction_failed: {e}"

def process_table(table_id, table, out_dir: Path, verbose=False):
    error_rows = []
    results_rows = []
    tbl_out = {"errors": [], "substrates": []}
    for idx, row in table.iterrows():
        smiles_raw = row["substrate_smiles"]
        label = row.get("label", f"row_{idx}")
        template = row["template"]
        mol = Chem.MolFromSmiles(smiles_raw)
        if mol is None:
            error_rec = {"row": idx, "label": label, "substrate_smiles": smiles_raw, "stage": "input", "error": "invalid_smiles", "template": template}
            error_rows.append(error_rec)
            tbl_out["errors"].append(error_rec)
            if verbose:
                print(f"[ERROR] {table_id} row {idx} ({label}) invalid SMILES")
            continue

        # Apply reaction and take the first product only
        prod, perr = apply_template_first(template, mol)
        if prod is None:
            msg = perr or "unknown_reaction_error"
            error_rec = {"row": idx, "label": label, "substrate_smiles": smiles_raw, "stage": "reaction", "error": msg, "template": template}
            error_rows.append(error_rec)
            tbl_out["errors"].append(error_rec)
            if verbose and msg != "no_products":
                print(f"[WARN] {table_id} row {idx} ({label}) reaction: {msg}")
            continue

        # Render image and assemble JSON product
        substrate_canon = Chem.MolToSmiles(mol, kekuleSmiles=False)
        png_path = out_dir / table_id / f"{idx:03d}_{slug}.png"
        try:
            render_reaction_png(mol, prod, png_path)
        except Exception:
            try:
                render_reaction_png(mol, prod, png_path, sub_img_size=(200, 200))
            except Exception as e2:
                msg = f"{e2}"
                error_rec = {"row": idx, "label": label, "substrate_smiles": smiles_raw, "stage": "render", "error": msg, "template": template}
                error_rows.append(error_rec)
                tbl_out["errors"].append(error_rec)
                if verbose:
                    print(f"[ERROR] {table_id} row {idx} ({label}) render: {e2}")
                continue

        prod_smiles = Chem.MolToSmiles(prod, kekuleSmiles=False)
        prod_smiles_h = to_explicit_h_smiles(prod)

        sub_entry = {
            "row": idx,
            "label": label,
            "substrate_smiles_input": smiles_raw,
            "substrate_smiles": substrate_canon,
            "product": {
                "product_smiles": prod_smiles,
                "product_smiles_explicit_h": prod_smiles_h,
                "image": str(png_path.relative_to(out_dir)),
            },
        }

        # products.csv row per substrate
        results_rows.append(
            {
                "row": idx,
                "label": label,
                "substrate_smiles": substrate_canon,
                "product_smiles": prod_smiles,
                "product_smiles_explicit_h": prod_smiles_h,
                "image": str(png_path.relative_to(out_dir)),
            }
        )

        tbl_out["substrates"].append(sub_entry)

    return tbl_out, results_rows, error_rows
