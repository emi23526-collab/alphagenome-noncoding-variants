"""Count variants and variant x condition scoring tasks in every CSV of data/processed. No API calls."""
import glob

import pandas as pd

print("File | rows | unique keys | status counts")
total_ok_keys = set()
for f in sorted(glob.glob("data/processed/*.csv")):
    try:
        d = pd.read_csv(f)
    except Exception as e:
        print(f"{f}: could not read ({e})")
        continue
    line = f"{f} | rows={len(d)}"
    if "key" in d.columns:
        line += f" | unique keys={d['key'].nunique()}"
    if "status" in d.columns:
        line += " | " + ", ".join(f"{k}={v}" for k, v in d["status"].value_counts().items())
        if "key" in d.columns:
            total_ok_keys |= set(d.loc[d["status"] == "ok", "key"])
    print(line)
    print("    columns:", ", ".join(map(str, d.columns)))

print("\nUnique variants scored ok in all score files:", len(total_ok_keys))
print("That number x 2 conditions:", len(total_ok_keys) * 2)
print("\nCompare these numbers with 5,902 in your protocol.")