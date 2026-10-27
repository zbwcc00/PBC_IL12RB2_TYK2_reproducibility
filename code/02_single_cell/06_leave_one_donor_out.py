from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--work-root",
        type=Path,
        default=REPO / "work" / "single_cell" / "hra008003",
    )
    args = parser.parse_args()
    work_root = args.work_root if args.work_root.is_absolute() else REPO / args.work_root
    statistics_root = work_root / "pseudobulk"
    input_path = statistics_root / "04_donor_target_logCPM.tsv"
    output_path = statistics_root / "07_target_leave_one_donor_out.tsv"

    with input_path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if (
            row["scheme"] == "primary_conservative"
            and row["cell_filter"] == "all_final_singlets"
            and row["minimum_cells"] == "20"
            and row["logCPM"] not in {"", "NA"}
        ):
            groups[(row["cell_type"], row["gene_symbol"])].append(row)

    output_rows: list[dict[str, object]] = []
    for (cell_type, gene_symbol), group_rows in sorted(groups.items()):
        conditions = {row["condition"] for row in group_rows}
        control_rows = [row for row in group_rows if row["condition"] == "control"]
        pbc_rows = [row for row in group_rows if row["condition"] == "PBC"]
        if conditions != {"control", "PBC"} or len(control_rows) < 3 or len(pbc_rows) < 3:
            continue

        control_values = [float(row["logCPM"]) for row in control_rows]
        pbc_values = [float(row["logCPM"]) for row in pbc_rows]
        full_effect = mean(pbc_values) - mean(control_values)
        leave_one_effects: list[tuple[str, str, float]] = []
        for omitted in group_rows:
            retained_control = [
                float(row["logCPM"])
                for row in group_rows
                if row["condition"] == "control" and row["donor_id"] != omitted["donor_id"]
            ]
            retained_pbc = [
                float(row["logCPM"])
                for row in group_rows
                if row["condition"] == "PBC" and row["donor_id"] != omitted["donor_id"]
            ]
            leave_one_effects.append(
                (omitted["donor_id"], omitted["condition"], mean(retained_pbc) - mean(retained_control))
            )

        same_direction = [
            effect == 0 if full_effect == 0 else math.copysign(1, effect) == math.copysign(1, full_effect)
            for _, _, effect in leave_one_effects
        ]
        for (omitted_donor, omitted_condition, effect), direction_matches in zip(
            leave_one_effects, same_direction
        ):
            output_rows.append(
                {
                    "cell_type": cell_type,
                    "gene_symbol": gene_symbol,
                    "n_control_donors": len(control_rows),
                    "n_PBC_donors": len(pbc_rows),
                    "full_effect_PBC_minus_control_logCPM": full_effect,
                    "omitted_donor": omitted_donor,
                    "omitted_condition": omitted_condition,
                    "leave_one_out_effect_logCPM": effect,
                    "same_direction_as_full": direction_matches,
                    "minimum_leave_one_out_effect": min(item[2] for item in leave_one_effects),
                    "maximum_leave_one_out_effect": max(item[2] for item in leave_one_effects),
                    "direction_preserved_fraction": sum(same_direction) / len(same_direction),
                }
            )

    fieldnames = list(output_rows[0])
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {len(output_rows)} leave-one-donor-out rows to {output_path}")


if __name__ == "__main__":
    main()
