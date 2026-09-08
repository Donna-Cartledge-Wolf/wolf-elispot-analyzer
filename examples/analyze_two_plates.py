from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from wolf_elispot import ELISpotAnalyzer

analyzer = ELISpotAnalyzer()

analyzer.analyze_image(
    image_path="assets/Plate_1_IFN-gamma_ELISpot.png",
    output_dir="results/plate_1",
    plate_id="Plate_1_IFN_gamma",
)

analyzer.analyze_image(
    image_path="assets/Plate_2_IL-2_ELISpot.png",
    output_dir="results/plate_2",
    plate_id="Plate_2_IL_2",
)
