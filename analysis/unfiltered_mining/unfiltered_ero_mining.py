import csv
import os
import sys
import time
import statistics
from collections import defaultdict

from evodex.operators import extract_operator_by_abstraction
from evodex.utils import reaction_hash

csv.field_size_limit(sys.maxsize)


def _ensure_dir_for_file(path: str) -> None:
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


def main():
    start_time = time.time()
    print("Unfiltered EVODEX-E mining started...")

    input_csv = "evodex/data/EVODEX-R.csv"
    out_csv = "data/processed/unfiltered_eros.csv"
    errors_csv = "data/processed/unfiltered_ero_errors.csv"
    report_txt = "data/processed/unfiltered_ero_report.txt"

    _ensure_dir_for_file(out_csv)
    _ensure_dir_for_file(errors_csv)
    _ensure_dir_for_file(report_txt)

    operator_map = defaultdict(lambda: {"smirks": None, "sources": set()})

    with open(errors_csv, "w", newline="") as errorfile:
        err_writer = csv.DictWriter(errorfile, fieldnames=["id", "smirks", "error_message"])
        err_writer.writeheader()

        with open(input_csv, "r", newline="") as infile:
            reader = csv.DictReader(infile)
            for i, row in enumerate(reader, start=1):
                if i % 100 == 0:
                    print(f"Processing EVODEX-R row {i} for operator extraction...")

                try:
                    op_smirks = extract_operator_by_abstraction(
                        row["smirks"],
                        abstraction="E",
                        matched=False,
                    )

                    if not op_smirks or op_smirks.startswith(">>") or op_smirks.endswith(">>"):
                        continue

                    op_hash = reaction_hash(op_smirks)
                    operator_map[op_hash]["smirks"] = op_smirks
                    operator_map[op_hash]["sources"].add(row["id"])

                except Exception as e:
                    err_writer.writerow(
                        {
                            "id": row.get("id", ""),
                            "smirks": row.get("smirks", ""),
                            "error_message": str(e),
                        }
                    )

    with open(out_csv, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=["id", "smirks", "sources"])
        writer.writeheader()
        for op_hash, data in operator_map.items():
            sources = sorted(data["sources"]) if data["sources"] else []
            writer.writerow(
                {
                    "id": op_hash,
                    "smirks": data["smirks"],
                    "sources": ",".join(sources),
                }
            )

    num_sources_list = [len(v["sources"]) for v in operator_map.values()]
    min_sources = min(num_sources_list) if num_sources_list else 0
    max_sources = max(num_sources_list) if num_sources_list else 0
    mean_sources = statistics.mean(num_sources_list) if num_sources_list else 0
    median_sources = statistics.median(num_sources_list) if num_sources_list else 0

    try:
        with open(errors_csv, "r", newline="") as errfile:
            error_count = sum(1 for _ in errfile) - 1
    except Exception:
        error_count = 0

    print("Statistics for extracted EVODEX-E operators:")
    print(f"  Number of operators: {len(operator_map)}")
    print(
        f"  Sources per operator: min={min_sources}, max={max_sources}, "
        f"mean={mean_sources:.2f}, median={median_sources}"
    )
    print(f"Number of extraction errors recorded: {error_count}")

    with open(report_txt, "w") as report_file:
        report_file.write("Unfiltered EVODEX-E Operator Extraction Report\n")
        report_file.write(f"Input CSV: {input_csv}\n")
        report_file.write(f"Output CSV: {out_csv}\n")
        report_file.write(f"Errors CSV: {errors_csv}\n")
        report_file.write(f"Number of operators: {len(operator_map)}\n")
        report_file.write("Sources per operator:\n")
        report_file.write(f"  min: {min_sources}\n")
        report_file.write(f"  max: {max_sources}\n")
        report_file.write(f"  mean: {mean_sources:.2f}\n")
        report_file.write(f"  median: {median_sources}\n")
        report_file.write(f"Extraction errors: {error_count}\n")

    elapsed_time = time.time() - start_time
    print(f"Unfiltered EVODEX-E mining completed in {elapsed_time:.2f} seconds.")


if __name__ == "__main__":
    main()
