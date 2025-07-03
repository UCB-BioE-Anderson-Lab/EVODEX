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

def generate_svg(smirks, filename, output_dir):
    try:
        rxn = indigo.loadReactionSmarts(smirks)
        rxn.aromatize()
        rxn.correctReactingCenters()
        renderer.renderToFile(rxn, os.path.join(output_dir, filename))
    except Exception as e:
        print(f"[!] Could not render {filename}: {e}")

def main():
    # Paths
    csv_path = "data/processed/top5_enzyme_operators.csv"
    base_output_dir = "figures/operators_svg"
    os.makedirs(base_output_dir, exist_ok=True)

    # Load data
    df = pd.read_csv(csv_path)
    print("Columns:", df.columns.tolist())
    df.drop_duplicates(subset=["organism", "protein_refs", "ec_num", "Cm_hash", "Nm_hash", "Em_hash"], inplace=True)
    df["enzyme_id"] = df["organism"] + "_" + df["protein_refs"] + "_" + df["ec_num"]

    melted = pd.melt(df, 
                     id_vars=["enzyme_id"], 
                     value_vars=["Cm", "Nm", "Em"], 
                     var_name="operator_type", 
                     value_name="smirks")

    melted["operator_hash"] = pd.melt(df, 
                                      id_vars=["enzyme_id"], 
                                      value_vars=["Cm_hash", "Nm_hash", "Em_hash"]
                                     )["value"]
    
    melted.drop_duplicates(subset=["enzyme_id", "operator_type", "operator_hash"], inplace=True)

    top5_enzyme_ids = melted["enzyme_id"].unique()[:5]
    enzyme_order = {eid: i+1 for i, eid in enumerate(top5_enzyme_ids)}

    df_top5 = melted[melted["enzyme_id"].isin(top5_enzyme_ids)]

    for (enzyme_id, operator_type), group in df_top5.groupby(["enzyme_id", "operator_type"]):
        enzyme_num = enzyme_order[enzyme_id]
        for idx, row in enumerate(group.itertuples(index=False)):
            smirks = row.smirks
            output_dir = os.path.join(base_output_dir, operator_type)
            os.makedirs(output_dir, exist_ok=True)
            filename = f"enzyme{enzyme_num}_{operator_type}{idx+1}.svg"
            generate_svg(smirks, filename, output_dir)

if __name__ == "__main__":
    main()
