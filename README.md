<p align="center">
  <img src="Wolf_Analytics_Logo.png" alt="Wolf Scientific Data Analysis & Automation" width="350">
</p>

# Wolf ELISpot Analyzer

**Automated image analysis for 96-well ELISpot assay plates**

Wolf ELISpot Analyzer is a Python module for loading an ELISpot plate image, locating the plate, mapping the standard 8 × 12 well layout, identifying discrete spot-forming units within each well, applying basic image-quality checks, and exporting quantitative per-well results.

The project was developed as part of **Wolf Analytics**, a scientific data-analysis and automation portfolio focused on practical assay workflows, reproducible analysis, and decision-support tools.

> **Project status:** portfolio / research prototype. The included synthetic plates have known spot counts and are intended for algorithm development and validation. Results from real assay images should be independently validated before use in regulated, diagnostic, or release-testing workflows.

## Key Results

- Automated mapping and analysis of all 96 wells
- Batch processing of plate images
- Per-well spot counts, signal metrics, QC flags, and heatmaps
- Synthetic ground-truth benchmark across 192 wells
- Pearson correlation between generated and detected spot counts: **r = 0.922**

---

## Why this project matters

Traditional ELISpot analysis can require specialized reader software or manual review. This project explores a transparent, reproducible Python workflow that converts a plate image into structured well-level data.

The analyzer is designed to answer a simple question:

**Can a standard plate image be transformed automatically into auditable, well-by-well spot measurements using open scientific Python tools?**

The initial module performs:

- automatic image loading;
- detection of the plate boundary;
- perspective correction;
- mapping of all 96 wells;
- local background estimation;
- dark-spot segmentation;
- local-intensity-peak spot counting for improved handling of crowded wells;
- connected-component area measurements;
- spot-area and integrated-darkness measurements;
- basic well-level QC flags;
- annotated-image export;
- spot-count heatmap export;
- single-image or batch-directory processing;
- CSV and JSON output; and
- validation against synthetic ground-truth data.

---

## Scientific context

ELISpot detects secreted analytes at the single-cell level. After analyte capture and development, individual secreting cells produce localized membrane spots. A useful image-analysis workflow therefore needs to distinguish true localized signal from membrane background, diffuse staining, artifacts, plate edges, and neighboring wells.

This module treats every well as an independent image-analysis region. For each well it:

1. estimates local background;
2. calculates a darkness signal relative to that background;
3. applies a robust threshold;
4. identifies local intensity peaks as candidate ELISpot events;
5. independently segments connected signal regions for area measurements; and
6. returns spot count plus supporting intensity, area, and QC metrics.

The workflow is deliberately transparent so that the scientific assumptions can be inspected and modified.

---

## Plate geometry

The demonstration plates use a standard 96-well assay-plate footprint and an 8 × 12 well layout. The synthetic images are saved at approximately **127.8 mm × 85.5 mm** with a **9.0 mm well-to-well pitch**, matching the intended ANSI/SLAS-style geometry closely enough for controlled algorithm testing.

The software does not assume that every uploaded photograph has identical pixel dimensions. It first attempts to identify and rectify the outer plate boundary, then maps the standardized well geometry onto the rectified image.

---

## Repository structure

```text
wolf_elispot_project/
├── wolf_elispot/
│   ├── __init__.py
│   ├── analyzer.py
│   └── cli.py
├── assets/
│   ├── Plate_1_IFN-gamma_ELISpot.png
│   ├── Plate_2_IL-2_ELISpot.png
│   └── ELISpot_synthetic_ground_truth.csv
├── examples/
│   ├── analyze_two_plates.py
│   └── validate_against_ground_truth.py
├── tests/
│   └── test_geometry.py
├── main.py
├── pyproject.toml
├── requirements.txt
└── README.md
```

---

## Installation

Clone the repository and install the required packages:

```bash
pip install -r requirements.txt
```

Or install the project as a local package:

```bash
pip install -e .
```

---

## Basic use

### From Python

```python
from wolf_elispot import ELISpotAnalyzer

analyzer = ELISpotAnalyzer()

results = analyzer.analyze_image(
    image_path="my_elispot_plate.png",
    output_dir="results/my_plate",
    plate_id="experiment_001",
)
```

### From the command line

```bash
wolf-elispot my_elispot_plate.png -o results/my_plate
```

A directory of plate images can also be processed as a batch:

```bash
wolf-elispot path/to/plate_images -o results/batch_run
```

You can tune spot-detection sensitivity:

```bash
wolf-elispot my_elispot_plate.png \
    -o results/my_plate \
    --threshold-sigma 2.0 \
    --min-area 2 \
    --max-area 120
```

---

## Outputs

For every analyzed plate the module generates:

### 1. Well-level CSV

Example fields:

| Field | Meaning |
|---|---|
| `well` | Well identifier such as A1 or H12 |
| `spot_count` | Number of accepted connected spot objects |
| `spot_area_px` | Total segmented spot area |
| `mean_spot_area_px` | Mean connected-component area |
| `median_spot_area_px` | Median connected-component area |
| `integrated_darkness` | Cumulative local darkness signal |
| `mean_background_intensity` | Median/mean well background estimate |
| `well_qc` | PASS or a warning flag |

### 2. Rectified plate image

The perspective-corrected plate used for analysis.

### 3. Annotated plate image

An overlay showing the 96 mapped wells and detected spot signal.

### 4. Spot-count heatmap

A plate-format visualization of detected events per well.

### 5. Analysis metadata JSON

Stores the detected plate corners and the exact analysis parameters used, supporting reproducibility and auditability.

---

## Synthetic validation set

Two synthetic 96-well ELISpot plates are included:

- **Plate 1 – IFN-γ ELISpot**
- **Plate 2 – IL-2 ELISpot**

A companion CSV contains the known number of synthetic spots generated in every well.

This makes it possible to test the analyzer objectively rather than judging performance only by eye.

Run:

```bash
python examples/validate_against_ground_truth.py
```

The validation script compares detected counts with the known synthetic counts and reports:

- mean absolute error;
- median absolute error;
- mean relative error; and
- Pearson correlation.

This is an important design feature of the project: the algorithm can be evaluated quantitatively and improved iteratively.

---

## Baseline validation results

The current v0.1 algorithm was tested against the two included synthetic plates (**192 wells total**).

- **Pearson correlation between generated and detected spot counts:** r = 0.922
- **Mean absolute error:** `24.39` spots/well
- **Median absolute error:** `14.0` spots/well

![Synthetic validation scatter](assets/validation_scatter.png)

Across the synthetic benchmark, detected counts tracked generated spot counts strongly overall (Pearson r = 0.922). Absolute enumeration became progressively less accurate in crowded wells, where overlapping signals produced fewer individually resolvable intensity peaks. The current prototype therefore performs better for relative signal discrimination than for absolute spot counting at high densities. Further development will focus on separating touching/overlapping spots and improving high-density quantitation. See [`VALIDATION.md`](VALIDATION.md) for additional interpretation and the planned improvement path.

---

## Current analysis approach

### Plate detection

The program searches for a large approximately rectangular contour with an aspect ratio near the standard assay-plate footprint. If a reliable quadrilateral is found, the image is perspective-corrected.

A full-image fallback is used when the plate already fills the image.

### Well mapping

After rectification, well centers are generated from:

- first-well position;
- 9-mm pitch;
- eight rows; and
- twelve columns.

This avoids requiring a separate circle-detection result for every well.

### Background correction

A Gaussian-smoothed image is used as a local background estimate. Dark signal is calculated relative to the smooth background rather than by applying a single global grayscale threshold.

### Spot thresholding

A robust threshold is calculated from the median and median absolute deviation of the local darkness signal within each well.

### Spot identification

Local intensity maxima in the background-corrected darkness surface are used as candidate spot centers. This is intentionally separate from connected-component area measurement because neighboring ELISpot signals can touch or merge in dense wells. Connected components remain useful for total segmented area and artifact/QC measurements.

---

## Quality-control philosophy

The module intentionally keeps QC separate from the final biological interpretation.

Current image-level/well-level flags include:

- low-contrast wells; and
- potentially saturated or excessively dark wells.

Future versions can add:

- replicate precision;
- positive-control acceptance;
- negative-control limits;
- maximum allowable background;
- TNTC / too-numerous-to-count classification;
- plate-edge artifacts;
- merged-spot detection; and
- experiment-specific acceptance criteria.

---

## Important limitations

This is an **image-analysis research prototype**, not a validated commercial ELISpot reader.

Real ELISpot images vary substantially with:

- scanner or camera optics;
- illumination;
- membrane color;
- spot color;
- focus;
- spot morphology;
- plate manufacturer;
- image compression;
- diffuse background;
- merged spots; and
- assay development chemistry.

For real experimental use, the algorithm should be validated against representative plate images and, ideally, against manually reviewed or instrument-generated reference counts.

The current algorithm may under-count touching or highly confluent spots and may over-count textured background when thresholds are too permissive.

---

## Planned improvements

Potential next steps include:

- automatic tuning from negative-control wells;
- watershed separation of touching spots;
- morphology-based rejection of irregular artifacts;
- color-space analysis for chromogenic substrates;
- replicate grouping from a sample map;
- positive/negative control identification;
- SFU normalization per input cell number;
- plate heatmaps;
- concentration- or stimulation-response visualization;
- interactive review of flagged wells;
- batch processing of multiple plate images;
- HTML/PDF analysis reports; and
- machine-learning-assisted spot classification.

---

## Portfolio significance

This project demonstrates an end-to-end scientific automation workflow combining:

**ELISpot · image analysis · assay QC · computer vision · Python · OpenCV · data automation · validation · quantitative error analysis · scientific decision support**

The emphasis is not only on producing a spot count, but on building a workflow that is reproducible, testable, inspectable, and suitable for iterative scientific improvement.

---

## Author / project

**Wolf Analytics**  
Scientific Data Analysis & Automation  
March 2026 – Present

Synthetic data are used for portfolio demonstration and algorithm testing. No client, patient, or proprietary experimental data are included.
