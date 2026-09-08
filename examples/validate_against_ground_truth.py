from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd

from wolf_elispot import ELISpotAnalyzer


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets"
OUT = ROOT / "validation_output"
OUT.mkdir(exist_ok=True)

truth = pd.read_csv(ASSETS / "ELISpot_synthetic_ground_truth.csv")
analyzer = ELISpotAnalyzer()

mapping = {
    "Plate 1 IFN-gamma": ASSETS / "Plate_1_IFN-gamma_ELISpot.png",
    "Plate 2 IL-2": ASSETS / "Plate_2_IL-2_ELISpot.png",
}

frames = []
for truth_name, image_path in mapping.items():
    plate_id = truth_name.replace(" ", "_")
    results = analyzer.analyze_image(image_path, OUT / plate_id, plate_id=plate_id)
    observed = pd.DataFrame([r.__dict__ for r in results])
    observed["plate"] = truth_name
    frames.append(observed)

observed = pd.concat(frames, ignore_index=True)
merged = truth.merge(observed, on=["plate", "well"], suffixes=("_truth", "_observed"))

merged["absolute_error"] = (merged["spot_count"] - merged["synthetic_spot_count"]).abs()
merged["relative_error"] = merged["absolute_error"] / merged["synthetic_spot_count"].clip(lower=1)

summary = pd.DataFrame({
    "metric": [
        "MAE",
        "Median absolute error",
        "Mean relative error",
        "Pearson correlation"
    ],
    "value": [
        merged["absolute_error"].mean(),
        merged["absolute_error"].median(),
        merged["relative_error"].mean(),
        merged[["synthetic_spot_count", "spot_count"]].corr().iloc[0,1],
    ],
})

merged.to_csv(OUT / "well_level_validation.csv", index=False)
summary.to_csv(OUT / "validation_summary.csv", index=False)

print(summary.to_string(index=False))
