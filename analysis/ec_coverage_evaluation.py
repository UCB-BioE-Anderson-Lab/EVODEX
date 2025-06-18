import pandas as pd
import random
from evodex.evaluation import match_operators
from evodex.evaluation import assign_evodex_F
import time
import socket
import platform
import datetime
import traceback

start_time = time.time()
start_timestamp = datetime.datetime.now().isoformat()

# Load prior results
input_path = "data/processed/ec_sampled_reactions.csv"
df = pd.read_csv(input_path)

from collections import defaultdict

total = len(df)
successful = 0
non_empty = 0
multiple = 0
match_counts = defaultdict(int)
matches_column = []
evaluation_times = []
error_messages = defaultdict(int)

f_assignments = []
f_successful = 0
f_empty = 0
f_counts = defaultdict(int)

for smiles in df['smiles']:
    start_eval = time.time()
    try:
        matches = match_operators(smiles)
        matches_column.append(matches)
        successful += 1
        if matches:
            non_empty += 1
            if len(matches) > 1:
                multiple += 1
            for m in matches:
                match_counts[m] += 1
    except Exception as e:
        matches_column.append([])
        error_msg = str(e).strip().split('\n')[0][:100]
        error_messages[error_msg] += 1
    evaluation_times.append(time.time() - start_eval)
    try:
        f_ids = assign_evodex_F(smiles)
        f_assignments.append(f_ids)
        if f_ids:
            f_successful += 1
            for f in f_ids:
                f_counts[f] += 1
        else:
            f_empty += 1
    except Exception as e:
        f_assignments.append([])

end_time = time.time()
end_timestamp = datetime.datetime.now().isoformat()
duration = end_time - start_time
avg_time = sum(evaluation_times) / len(evaluation_times) if evaluation_times else 0
max_time = max(evaluation_times) if evaluation_times else 0
min_time = min(evaluation_times) if evaluation_times else 0

df['evodex_matches'] = matches_column
df['evodex_F_ids'] = f_assignments

# Save updated CSV
output_path = "data/processed/ec_sampled_reactions_with_evodex.csv"
df.to_csv(output_path, index=False)

# Report and save stats
report_lines = [
    f"Evaluation Summary Report",
    f"===========================",
    f"Start time: {start_timestamp}",
    f"End time: {end_timestamp}",
    f"Total runtime: {duration:.2f} seconds",
    f"Host machine: {socket.gethostname()}",
    f"Platform: {platform.platform()}",
    "",
    f"Reaction Dataset:",
    f"- Total reactions in input: {total}",
    f"- Successful evaluations (no crash): {successful}",
    f"- Reactions with ≥1 operator match: {non_empty}",
    f"- Reactions with >1 operator match: {multiple}",
    f"- Reactions with 0 operator matches: {total - non_empty}",
    "",
    f"Timing Information (per reaction):",
    f"- Average evaluation time: {avg_time:.4f} s",
    f"- Fastest evaluation time: {min_time:.4f} s",
    f"- Slowest evaluation time: {max_time:.4f} s",
    "",
    f"Operator Match Frequency:",
]

for k, v in sorted(match_counts.items(), key=lambda x: -x[1]):
    report_lines.append(f"{k}: {v}")

report_lines.append("")
report_lines.append("Top Error Messages (if any):")
for msg, count in sorted(error_messages.items(), key=lambda x: -x[1])[:10]:
    report_lines.append(f"{count}x - {msg}")

report_lines.extend([
    "",
    "EVODEX-F Assignments:",
    f"- Successful F-ID assignments: {f_successful}",
    f"- Reactions with no F-ID match: {f_empty}",
    "",
    "EVODEX-F Frequency:",
])
for k, v in sorted(f_counts.items(), key=lambda x: -x[1]):
    report_lines.append(f"{k}: {v}")

report_path = "data/processed/evaluation_stats.txt"
with open(report_path, "w") as f:
    f.write("\n".join(report_lines))

print("\n".join(report_lines))
print(f"\nSaved evaluated reactions to: {output_path}")
print(f"Saved stats report to: {report_path}")
