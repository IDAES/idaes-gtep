#!/usr/bin/env python3
"""Generate a compact text report for a GTEP PH run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    config_path = Path(args.config)

    try:
        import yaml
    except ImportError as err:
        raise RuntimeError("PyYAML is required for this report script.") from err

    with config_path.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    summary = _read_json(run_dir / "final" / "ph_summary.json")
    final_nonants = _read_json(run_dir / "final" / "nonanticipative_values.json")
    scenario_solutions = _read_json(run_dir / "final" / "scenario_solutions.json")

    history_path = run_dir / "convergence" / "history.csv"
    history = _read_csv(history_path) if history_path.exists() else []

    weights = config["data"]["representative_weights"]
    len_reps = config["data"]["len_reps"]
    num_reps = config["data"]["num_reps"]

    represented_hours = sum(float(w) for w in weights) * float(len_reps)
    unweighted_rep_hours = float(num_reps) * float(len_reps)

    print("=" * 80)
    print("GTEP Progressive Hedging Report")
    print("=" * 80)
    print(f"Run directory: {run_dir}")
    print(f"Config:        {config_path}")
    print()

    print("PH status")
    print("-" * 80)
    print(f"Converged:                    {summary.get('converged')}")
    print(f"Final iteration:              {summary.get('final_iteration')}")
    print(f"State iteration:              {summary.get('state_iteration')}")
    print(f"Number of scenarios:          {summary.get('num_scenarios')}")
    print(f"Nonanticipative variables:    {summary.get('num_nonanticipative_variables')}")
    print(f"Expected objective:           {summary.get('expected_objective')}")
    print()

    residuals = summary.get("residuals", {})
    print("Final residuals")
    print("-" * 80)
    for key, value in residuals.items():
        print(f"{key:32s} {value}")
    print()

    print("Operational time represented")
    print("-" * 80)
    print(f"Representative periods:       {num_reps}")
    print(f"Representative length:        {len_reps} hours")
    print(f"Unweighted modeled hours:     {unweighted_rep_hours} hours")
    print(f"Representative weights:       {weights}")
    print(f"Weighted represented hours:   {represented_hours} hours")
    print(f"Weighted represented days:    {represented_hours / 24.0} days")
    print()

    print("Scenario objectives")
    print("-" * 80)
    for item in scenario_solutions.get("scenario_solutions", []):
        solver = item.get("solver", {})
        print(
            f"Scenario {item.get('scenario_id')}: "
            f"date={item.get('representative_date')}, "
            f"acceptable={item.get('acceptable')}, "
            f"objective={solver.get('objective')}, "
            f"gap={solver.get('gap')}, "
            f"termination={solver.get('termination_condition')}"
        )
    print()

    if history:
        print("Convergence history")
        print("-" * 80)
        for row in history:
            print(
                f"iter={row.get('iteration')}, "
                f"converged={row.get('converged')}, "
                f"expected_obj={row.get('expected_objective')}, "
                f"max_abs_resid={row.get('max_abs_residual')}, "
                f"weighted_l2={row.get('weighted_l2_residual')}"
            )
        print()

    nonant_records = final_nonants.get("nonanticipative_values", [])
    nonzero = [
        rec for rec in nonant_records
        if abs(float(rec.get("final_value", rec.get("xbar", 0.0)))) > 1.0e-8
    ]

    print("Final nonanticipative solution")
    print("-" * 80)
    print(f"Total records:                {len(nonant_records)}")
    print(f"Nonzero records:              {len(nonzero)}")
    print("First 25 nonzero records:")
    for rec in nonzero[:25]:
        print(
            f"  {rec.get('variable_id')} = "
            f"{rec.get('final_value', rec.get('xbar'))}"
        )


def _read_json(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _read_csv(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


if __name__ == "__main__":
    main()