
#!/usr/bin/env python3

from pathlib import Path
import json
import glob

from .twinmap import twinmap as run_twinmap

BASE_DIR = Path(__file__).resolve().parent
IN_DIR = BASE_DIR / "out/1_project"
OUT_DIR = BASE_DIR / "out/2_map"
LOG_DIR = BASE_DIR / "out/_logs"

ONLY_TABLES = []  # e.g. ["trifluoroacetanilide_methanolysis"]


def stream_jsonl(path: Path):
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            yield json.loads(line)


def pick_substrate_smiles(row: dict):
    if "smiles" in row:
        return row["smiles"]
    if "substrate_smiles" in row:
        return row["substrate_smiles"]
    for k in row.keys():
        lk = str(k).lower()
        if "smiles" in lk and "product" not in lk:
            return row[k]
    return None


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    errors_path = LOG_DIR / "2_atommap.errors.jsonl"
    with errors_path.open("w", encoding="utf-8") as err_f:
        pattern = str(IN_DIR / "*.project.jsonl")
        input_files = sorted(glob.glob(pattern))

        for in_name in input_files:
            in_path = Path(in_name)
            table_id = in_path.name.replace(".project.jsonl", "")

            if ONLY_TABLES and table_id not in ONLY_TABLES:
                continue

            out_path = OUT_DIR / f"{table_id}.atommap.jsonl"
            with out_path.open("w", encoding="utf-8") as out_f:
                for row in stream_jsonl(in_path):
                    sub = pick_substrate_smiles(row)
                    prod = row.get("product_smiles")
                    template = row.get("reaction_template")

                    if not sub or not isinstance(sub, str):
                        err_rec = {
                            "table_id": table_id,
                            "row": row.get("row"),
                            "label": row.get("label"),
                            "stage": "input",
                            "error": "missing_or_invalid_substrate_smiles",
                        }
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        continue

                    if not prod or not isinstance(prod, str):
                        err_rec = {
                            "table_id": table_id,
                            "row": row.get("row"),
                            "label": row.get("label"),
                            "stage": "input",
                            "error": "missing_or_invalid_product_smiles",
                        }
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        continue

                    if not template or not isinstance(template, str):
                        err_rec = {
                            "table_id": table_id,
                            "row": row.get("row"),
                            "label": row.get("label"),
                            "stage": "input",
                            "error": "missing_or_invalid_reaction_template",
                        }
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        continue

                    try:
                        res = run_twinmap(sub, prod, template)
                        mapped_sub = res.get("mapped_substrate_smiles", sub)
                        mapped_prod = res.get("mapped_product_smiles", prod)
                    except Exception as e:
                        err_rec = {
                            "table_id": table_id,
                            "row": row.get("row"),
                            "label": row.get("label"),
                            "stage": "mapping",
                            "error": f"{e}",
                        }
                        err_f.write(json.dumps(err_rec, ensure_ascii=False) + "\n")
                        continue

                    rec = dict(row)
                    rec["mapped_substrate_smiles"] = mapped_sub
                    rec["mapped_product_smiles"] = mapped_prod
                    rec["mapping"] = res

                    out_f.write(json.dumps(rec, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
