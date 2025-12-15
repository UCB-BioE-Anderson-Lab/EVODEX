#!/usr/bin/env python3

import os
import pandas as pd
from indigo import Indigo
from indigo.renderer import IndigoRenderer

# Configure Indigo
indigo = Indigo()
renderer = IndigoRenderer(indigo)

indigo.setOption("render-output-format", "svg")
indigo.setOption("render-implicit-hydrogens-visible", False)
indigo.setOption("embedding-uniqueness", "none")
indigo.setOption("render-bond-length", 40)
indigo.setOption("render-atom-ids-visible", False)
indigo.setOption("render-aam-color", "0, 0, 0")
indigo.setOption("render-coloring", True)
indigo.setOption("render-stereo-style", "none")
indigo.setOption("render-label-mode", "hetero")


def generate_svg(smirks: str, filename: str, output_dir: str) -> None:
    try:
        rxn = indigo.loadReactionSmarts(smirks)
        rxn.aromatize()
        rxn.correctReactingCenters()
        renderer.renderToFile(rxn, os.path.join(output_dir, filename))
    except Exception as e:
        print(f"[!] Could not render {filename}: {e}")


def main() -> None:
    # Paths (updated to match the new top5 script output folder)
    csv_path = "analysis/top5_diverse_enzymes_operators/output/top5_enzyme_operators.csv"
    base_output_dir = "analysis/top5_diverse_enzymes_operators/output/operators_svg"
    os.makedirs(base_output_dir, exist_ok=True)

    # Load data
    df = pd.read_csv(csv_path)
    print("Columns:", df.columns.tolist())

    # Updated dedupe keys and operator tiers: Bm, Cm, Em
    df.drop_duplicates(
        subset=["organism", "protein_refs", "ec_num", "Bm_hash", "Cm_hash", "Em_hash"],
        inplace=True,
    )

    df["enzyme_id"] = df["organism"] + "_" + df["protein_refs"].astype(str) + "_" + df["ec_num"]

    # Melt operators
    melted_ops = pd.melt(
        df,
        id_vars=["enzyme_id"],
        value_vars=["Bm", "Cm", "Em"],
        var_name="operator_type",
        value_name="smirks",
    )

    # Melt hashes in the same row order, then attach
    melted_hashes = pd.melt(
        df,
        id_vars=["enzyme_id"],
        value_vars=["Bm_hash", "Cm_hash", "Em_hash"],
        var_name="hash_type",
        value_name="operator_hash",
    )

    melted_ops["operator_hash"] = melted_hashes["operator_hash"].values

    # Keep only rows with real operator + hash
    melted_ops = melted_ops.dropna(subset=["smirks", "operator_hash"])
    melted_ops = melted_ops[(melted_ops["smirks"].astype(str).str.len() > 0) & (melted_ops["operator_hash"].astype(str).str.len() > 0)]

    # De-duplicate per enzyme/operator/hash
    melted_ops.drop_duplicates(subset=["enzyme_id", "operator_type", "operator_hash"], inplace=True)

    # Deterministic ordering of enzymes for numbering
    top5_enzyme_ids = (
        df[["enzyme_id", "organism", "ec_num", "protein_refs"]]
        .drop_duplicates()
        .sort_values(["ec_num", "organism", "protein_refs", "enzyme_id"])
        ["enzyme_id"]
        .tolist()
    )[:5]

    enzyme_order = {eid: i + 1 for i, eid in enumerate(top5_enzyme_ids)}
    df_top5 = melted_ops[melted_ops["enzyme_id"].isin(top5_enzyme_ids)]

    # Render
    for (enzyme_id, operator_type), group in df_top5.groupby(["enzyme_id", "operator_type"]):
        enzyme_num = enzyme_order[enzyme_id]

        # Stable ordering within enzyme/operator_type: most frequent first, then hash
        counts = (
            group.groupby("operator_hash")
            .size()
            .reset_index(name="n")
            .sort_values(["n", "operator_hash"], ascending=[False, True])
        )
        ordered_hashes = counts["operator_hash"].tolist()

        # Build an ordered group view
        group_by_hash = {h: group[group["operator_hash"] == h].iloc[0] for h in ordered_hashes}

        output_dir = os.path.join(base_output_dir, operator_type)
        os.makedirs(output_dir, exist_ok=True)

        for idx, h in enumerate(ordered_hashes, start=1):
            row = group_by_hash[h]
            smirks = row["smirks"]

            # enzyme{N}_{TYPE}{k}.svg
            filename = f"enzyme{enzyme_num}_{operator_type}{idx}.svg"
            generate_svg(smirks, filename, output_dir)


if __name__ == "__main__":
    main()
