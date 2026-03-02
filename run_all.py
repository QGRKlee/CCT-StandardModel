#!/usr/bin/env python3
"""Master verification runner for all five theorems.

Usage:
  python run_all.py                  # Run all theorems
  python run_all.py --skip-partition # Skip Theorem B (slow, ~10-15 min)
  python run_all.py --theorem A      # Run only Theorem A
"""

import subprocess
import sys
import time
import re
import argparse

SCRIPTS = {
    "A": ("verify_hierarchy.py",    "Theorem A (Hierarchy)"),
    "B": ("verify_partition.py",    "Theorem B (Partition)"),
    "C": ("verify_gut_breaking.py", "Theorem C (GUT Breaking)"),
    "D": ("verify_triality.py",     "Theorem D (Triality)"),
    "E": ("verify_angles.py",       "Theorem E (Angles)"),
}

def run_script(script_path, label):
    """Run a verification script and parse results."""
    start = time.time()
    try:
        result = subprocess.run(
            [sys.executable, script_path],
            capture_output=True, text=True, timeout=1200
        )
        elapsed = time.time() - start
        stdout = result.stdout
        passes = len(re.findall(r'\[PASS\]', stdout))
        fails = len(re.findall(r'\[FAIL\]', stdout))
        return {
            "label": label,
            "passes": passes,
            "fails": fails,
            "elapsed": elapsed,
            "returncode": result.returncode,
            "stdout": stdout,
            "stderr": result.stderr,
        }
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        return {
            "label": label,
            "passes": 0,
            "fails": 1,
            "elapsed": elapsed,
            "returncode": -1,
            "stdout": "",
            "stderr": "TIMEOUT (>20 min)",
        }

def main():
    parser = argparse.ArgumentParser(description="Run machine verification suite")
    parser.add_argument("--skip-partition", action="store_true",
                        help="Skip Theorem B (slow, ~10-15 min)")
    parser.add_argument("--theorem", type=str, choices=list(SCRIPTS.keys()),
                        help="Run only a specific theorem")
    args = parser.parse_args()

    print("=" * 66)
    print("  MACHINE VERIFICATION SUITE")
    print("  Kinematic Spinors on Division-Algebra Root Systems")
    print("=" * 66)
    print()

    # Determine which scripts to run
    if args.theorem:
        to_run = {args.theorem: SCRIPTS[args.theorem]}
    elif args.skip_partition:
        to_run = {k: v for k, v in SCRIPTS.items() if k != "B"}
    else:
        to_run = SCRIPTS

    results = []
    total_start = time.time()

    for key in ["A", "B", "C", "D", "E"]:
        if key not in to_run:
            continue
        script, label = to_run[key]
        print(f"Running {script} ...")
        r = run_script(script, label)
        results.append(r)

        # Print output
        if r["stdout"]:
            for line in r["stdout"].strip().split("\n"):
                print(f"  {line}")
        if r["stderr"] and r["returncode"] != 0:
            print(f"  STDERR: {r['stderr'][:200]}")
        print()

    total_elapsed = time.time() - total_start

    # Summary table
    print("=" * 66)
    print("  VERIFICATION SUMMARY")
    print("=" * 66)
    total_pass = 0
    total_fail = 0
    all_ok = True
    for r in results:
        total_pass += r["passes"]
        total_fail += r["fails"]
        status = "PASS" if r["fails"] == 0 and r["returncode"] == 0 else "FAIL"
        if status == "FAIL":
            all_ok = False
        print(f"  {r['label']:30s} : {r['passes']:2d}/{r['passes']+r['fails']:2d} {status:4s}"
              f"   ({r['elapsed']:.1f}s)")
    print("-" * 66)
    print(f"  {'TOTAL':30s} : {total_pass:2d}/{total_pass+total_fail:2d}")
    print(f"  {'Time':30s} : {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")
    print()
    if all_ok:
        print("  *** ALL PASS ***")
    else:
        print("  *** FAILURES DETECTED ***")
    print("=" * 66)

    sys.exit(0 if all_ok else 1)

if __name__ == "__main__":
    main()
