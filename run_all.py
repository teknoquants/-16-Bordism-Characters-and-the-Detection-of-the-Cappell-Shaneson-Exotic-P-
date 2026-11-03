#!/usr/bin/env python3
"""Run every validation script for the Problem 5 white paper and print a summary table.

Usage:  python run_all.py            (from this folder; needs Python 3.9+ and NumPy)
Exit status is 0 only if every script passes.
"""
import pathlib, re, subprocess, sys

HERE = pathlib.Path(__file__).resolve().parent
ORDER = ["1_1", "2_1", "2_2", "2_3", "2_4", "2_5", "2_6", "2_7",
         "3_1", "3_2", "3_3", "3_4", "3_5", "theorem_3_6",
         "4_1", "4_2", "4_3", "4_4", "4_5",
         "5_1", "5_2", "5_3", "5_4", "5_5", "5_6", "5_7",
         "A_1", "A_2", "A_3", "A_4", "A_5", "B_1"]

def path(key):
    name = key if key.startswith("theorem") else f"lemma_{key}"
    return HERE / "scripts" / f"{name}.py"

def main():
    verbose = "-v" in sys.argv
    rows, failures = [], 0
    for key in ORDER:
        p = path(key)
        r = subprocess.run([sys.executable, str(p)], capture_output=True, text=True)
        out = r.stdout + r.stderr
        if verbose:
            print(out)
        m = re.search(r": (\d+)/(\d+) checks passed", out)
        passed, total = (m.group(1), m.group(2)) if m else ("?", "?")
        notes = out.count("[NOTE]")
        status = "PASS" if r.returncode == 0 else "FAIL"
        failures += r.returncode != 0
        rows.append((p.stem, status, f"{passed}/{total}", notes))
    print(f"{'script':<16}{'status':<8}{'checks':<10}{'external/structural notes'}")
    print("-" * 60)
    for name, status, checks, notes in rows:
        print(f"{name:<16}{status:<8}{checks:<10}{notes}")
    print("-" * 60)
    print(f"{len(rows) - failures}/{len(rows)} scripts passed")
    sys.exit(1 if failures else 0)

if __name__ == "__main__":
    main()
