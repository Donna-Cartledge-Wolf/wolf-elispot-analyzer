from __future__ import annotations

import argparse
from pathlib import Path
from .analyzer import ELISpotAnalyzer, AnalysisConfig


def build_parser():
    p = argparse.ArgumentParser(
        prog="wolf-elispot",
        description="Automatic 96-well ELISpot image analysis."
    )
    p.add_argument("input", help="Plate image or directory of plate images.")
    p.add_argument("-o", "--output", default="elispot_output", help="Output directory.")
    p.add_argument("--plate-id", default=None, help="Optional identifier for a single plate.")
    p.add_argument("--threshold-sigma", type=float, default=0.20,
                   help="Spot sensitivity. Lower = more permissive.")
    p.add_argument("--min-area", type=int, default=1,
                   help="Minimum thresholded component area in pixels.")
    p.add_argument("--max-area", type=int, default=180,
                   help="Maximum thresholded component area in pixels.")
    return p


def main():
    args = build_parser().parse_args()

    config = AnalysisConfig(
        threshold_sigma=args.threshold_sigma,
        min_spot_area_px=args.min_area,
        max_spot_area_px=args.max_area,
    )
    analyzer = ELISpotAnalyzer(config)

    input_path = Path(args.input)
    output_path = Path(args.output)

    if input_path.is_dir():
        batches = analyzer.analyze_folder(input_path, output_path)
        print(f"Analyzed {len(batches)} plate images.")
    else:
        results = analyzer.analyze_image(
            input_path, output_path, plate_id=args.plate_id
        )
        print(f"Analyzed {len(results)} wells.")

    print(f"Results written to: {output_path.resolve()}")


if __name__ == "__main__":
    main()
