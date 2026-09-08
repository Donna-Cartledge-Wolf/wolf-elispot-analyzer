from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Sequence, Tuple, Dict, List
import csv
import json

import cv2
import numpy as np


ROWS = "ABCDEFGH"
COLS = list(range(1, 13))


@dataclass
class WellResult:
    plate_id: str
    well: str
    row: str
    column: int
    spot_count: int
    segmented_spot_area_px: int
    mean_component_area_px: float
    median_component_area_px: float
    integrated_darkness: float
    median_background_intensity: float
    well_qc: str


@dataclass
class AnalysisConfig:
    # Standardized working canvas.
    warp_width: int = 1278
    warp_height: int = 855

    # Fallback ANSI/SLAS-style geometry at ~10 px/mm.
    first_well_x: float = 144.0
    first_well_y: float = 112.0
    pitch_x: float = 90.0
    pitch_y: float = 90.0
    well_radius: float = 29.0

    # Automatic well-lattice detection.
    use_hough_well_detection: bool = True
    hough_dp: float = 1.2
    hough_min_dist: float = 55.0
    hough_param1: float = 100.0
    hough_param2: float = 28.0
    hough_min_radius: int = 23
    hough_max_radius: int = 40
    lattice_cluster_tolerance_px: float = 24.0

    # Spot detection.
    background_blur_sigma: float = 5.0
    threshold_sigma: float = 0.20
    peak_kernel_size: int = 3
    minimum_darkness: float = 1.0
    min_spot_area_px: int = 1
    max_spot_area_px: int = 180
    border_margin_px: int = 2

    # QC heuristics.
    min_valid_well_std: float = 1.5
    max_fraction_dark: float = 0.45


class ELISpotAnalyzer:
    """
    Automated 96-well ELISpot image analyzer.

    Main workflow
    -------------
    1. Load an image.
    2. Detect the assay plate boundary and correct perspective.
    3. Detect the 8 x 12 well lattice automatically when possible.
    4. Estimate local background within each well.
    5. Detect local dark-intensity peaks as candidate ELISpot events.
    6. Measure thresholded spot area and integrated darkness.
    7. Apply simple image-QC flags.
    8. Export per-well CSV, metadata JSON, rectified image,
       annotated image, and a plate heatmap.

    The included synthetic images have known ground-truth spot counts.
    Real ELISpot imaging systems vary, so the algorithm must be validated
    on representative images before scientific, regulated, diagnostic,
    or release-testing use.
    """

    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()

    @staticmethod
    def load_image(path: str | Path) -> np.ndarray:
        image = cv2.imread(str(path), cv2.IMREAD_COLOR)
        if image is None:
            raise FileNotFoundError(f"Could not read image: {path}")
        return image

    @staticmethod
    def _order_points(pts: np.ndarray) -> np.ndarray:
        pts = np.asarray(pts, dtype=np.float32)
        s = pts.sum(axis=1)
        d = np.diff(pts, axis=1).ravel()
        ordered = np.zeros((4, 2), dtype=np.float32)
        ordered[0] = pts[np.argmin(s)]
        ordered[2] = pts[np.argmax(s)]
        ordered[1] = pts[np.argmin(d)]
        ordered[3] = pts[np.argmax(d)]
        return ordered

    def detect_plate_quad(self, image: np.ndarray) -> np.ndarray:
        """
        Detect the outer plate rectangle.

        Large quadrilateral contours with an aspect ratio close to the
        standard microplate footprint are preferred. When detection is
        uncertain, the full image is used.
        """
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        blur = cv2.GaussianBlur(gray, (5, 5), 0)
        edges = cv2.Canny(blur, 40, 120)
        edges = cv2.dilate(edges, np.ones((5, 5), np.uint8), iterations=1)

        contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        h, w = gray.shape
        image_area = h * w
        candidates = []

        for contour in contours:
            area = cv2.contourArea(contour)
            if area < 0.25 * image_area:
                continue

            peri = cv2.arcLength(contour, True)
            approx = cv2.approxPolyDP(contour, 0.02 * peri, True)
            if len(approx) != 4:
                continue

            rect = cv2.minAreaRect(contour)
            rw, rh = rect[1]
            if rw <= 0 or rh <= 0:
                continue

            ratio = max(rw, rh) / min(rw, rh)
            ratio_score = abs(ratio - 1.495)
            area_score = -area
            candidates.append((ratio_score, area_score, approx.reshape(4, 2)))

        if candidates:
            candidates.sort(key=lambda x: (x[0], x[1]))
            best = candidates[0]
            # Avoid aggressive cropping when the detected rectangle is already
            # essentially the whole image.
            quad = self._order_points(best[2])
            quad_area = cv2.contourArea(quad.astype(np.float32))
            if quad_area / image_area > 0.80:
                return np.array(
                    [[0, 0], [w - 1, 0], [w - 1, h - 1], [0, h - 1]],
                    dtype=np.float32,
                )
            return quad

        return np.array(
            [[0, 0], [w - 1, 0], [w - 1, h - 1], [0, h - 1]],
            dtype=np.float32,
        )

    def rectify_plate(self, image: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        quad = self.detect_plate_quad(image)
        dst = np.array(
            [
                [0, 0],
                [self.config.warp_width - 1, 0],
                [self.config.warp_width - 1, self.config.warp_height - 1],
                [0, self.config.warp_height - 1],
            ],
            dtype=np.float32,
        )
        matrix = cv2.getPerspectiveTransform(quad.astype(np.float32), dst)
        warped = cv2.warpPerspective(
            image, matrix, (self.config.warp_width, self.config.warp_height)
        )
        return warped, quad

    @staticmethod
    def _cluster_1d(values: np.ndarray, tolerance: float) -> List[Tuple[float, int]]:
        if len(values) == 0:
            return []
        groups: List[List[float]] = []
        for value in sorted(float(v) for v in values):
            if not groups or value - float(np.mean(groups[-1])) > tolerance:
                groups.append([value])
            else:
                groups[-1].append(value)
        return [(float(np.mean(g)), len(g)) for g in groups]

    def detect_well_lattice(
        self, plate_image: np.ndarray
    ) -> Tuple[Dict[str, Tuple[float, float]], Dict]:
        """
        Detect circular wells and infer the 12 x 8 lattice.

        The Hough detections are clustered separately along x and y.
        True column positions tend to have many circle detections across
        rows, while true row positions have many detections across columns.
        """
        fallback = self.well_centers_fallback()
        diagnostics = {
            "method": "fallback_standard_geometry",
            "detected_circle_count": 0,
            "detected_x_positions": [],
            "detected_y_positions": [],
        }

        if not self.config.use_hough_well_detection:
            return fallback, diagnostics

        gray = cv2.cvtColor(plate_image, cv2.COLOR_BGR2GRAY)
        gray = cv2.GaussianBlur(gray, (5, 5), 1.2)

        circles = cv2.HoughCircles(
            gray,
            cv2.HOUGH_GRADIENT,
            dp=self.config.hough_dp,
            minDist=self.config.hough_min_dist,
            param1=self.config.hough_param1,
            param2=self.config.hough_param2,
            minRadius=self.config.hough_min_radius,
            maxRadius=self.config.hough_max_radius,
        )

        if circles is None:
            return fallback, diagnostics

        circles = np.asarray(circles[0], dtype=float)
        diagnostics["detected_circle_count"] = int(len(circles))

        x_clusters = self._cluster_1d(
            circles[:, 0], self.config.lattice_cluster_tolerance_px
        )
        y_clusters = self._cluster_1d(
            circles[:, 1], self.config.lattice_cluster_tolerance_px
        )

        # A true column should usually be detected in several of 8 rows;
        # a true row should usually be detected in several of 12 columns.
        x_candidates = [x for x, n in x_clusters if n >= 4]
        y_candidates = [y for y, n in y_clusters if n >= 6]

        # If extra candidates remain, favor those with stronger support.
        if len(x_candidates) != 12:
            supported = sorted(x_clusters, key=lambda t: (-t[1], t[0]))[:12]
            x_candidates = sorted(x for x, _ in supported)
        else:
            x_candidates = sorted(x_candidates)

        if len(y_candidates) != 8:
            supported = sorted(y_clusters, key=lambda t: (-t[1], t[0]))[:8]
            y_candidates = sorted(y for y, _ in supported)
        else:
            y_candidates = sorted(y_candidates)

        if len(x_candidates) != 12 or len(y_candidates) != 8:
            return fallback, diagnostics

        # Reject implausible lattices.
        x_pitch = np.diff(x_candidates)
        y_pitch = np.diff(y_candidates)
        if (
            np.min(x_pitch) < 50
            or np.max(x_pitch) > 130
            or np.min(y_pitch) < 50
            or np.max(y_pitch) > 130
        ):
            return fallback, diagnostics

        centers: Dict[str, Tuple[float, float]] = {}
        for r, row in enumerate(ROWS):
            for c, col in enumerate(COLS):
                centers[f"{row}{col}"] = (
                    float(x_candidates[c]),
                    float(y_candidates[r]),
                )

        diagnostics.update(
            {
                "method": "automatic_hough_lattice",
                "detected_x_positions": [round(x, 3) for x in x_candidates],
                "detected_y_positions": [round(y, 3) for y in y_candidates],
                "median_pitch_x_px": float(np.median(x_pitch)),
                "median_pitch_y_px": float(np.median(y_pitch)),
            }
        )
        return centers, diagnostics

    def well_centers_fallback(self) -> Dict[str, Tuple[float, float]]:
        centers = {}
        for r, row in enumerate(ROWS):
            for c, col in enumerate(COLS):
                x = self.config.first_well_x + c * self.config.pitch_x
                y = self.config.first_well_y + r * self.config.pitch_y
                centers[f"{row}{col}"] = (x, y)
        return centers

    def well_centers(self) -> Dict[str, Tuple[float, float]]:
        """Public accessor for the fallback standard geometry."""
        return self.well_centers_fallback()

    def _well_mask(
        self, shape: Tuple[int, int], center: Tuple[float, float]
    ) -> np.ndarray:
        mask = np.zeros(shape, dtype=np.uint8)
        cx, cy = int(round(center[0])), int(round(center[1]))
        radius = int(round(self.config.well_radius - self.config.border_margin_px))
        cv2.circle(mask, (cx, cy), radius, 255, -1)
        return mask

    def _segment_spots(
        self, gray: np.ndarray, well_mask: np.ndarray
    ) -> Tuple[np.ndarray, List[int], np.ndarray, float, float, str]:
        """
        Detect spot candidates by local intensity maxima in a background-
        corrected darkness image.

        Connected components are retained as an independent area metric.
        Peak counting is used for spot_count because dense ELISpot wells
        can contain touching signals that merge into a single component.
        """
        pixels = gray[well_mask > 0]
        if pixels.size == 0:
            return (
                np.zeros_like(gray),
                [],
                np.empty((0, 2), dtype=int),
                0.0,
                0.0,
                "FAIL_EMPTY",
            )

        median_bg = float(np.median(pixels))
        well_std = float(np.std(pixels))

        smooth = cv2.GaussianBlur(
            gray.astype(np.float32),
            (0, 0),
            sigmaX=self.config.background_blur_sigma,
            sigmaY=self.config.background_blur_sigma,
        )
        darkness = smooth - gray.astype(np.float32)
        dark_pixels = darkness[well_mask > 0]

        med = float(np.median(dark_pixels))
        mad = float(np.median(np.abs(dark_pixels - med)))
        robust_sigma = max(0.5, 1.4826 * mad)
        threshold = max(
            self.config.minimum_darkness,
            med + self.config.threshold_sigma * robust_sigma,
        )

        binary = ((darkness > threshold) & (well_mask > 0)).astype(np.uint8) * 255

        n_labels, labels, stats, _ = cv2.connectedComponentsWithStats(binary, 8)
        keep = np.zeros_like(binary)
        component_areas: List[int] = []
        for label in range(1, n_labels):
            area = int(stats[label, cv2.CC_STAT_AREA])
            if self.config.min_spot_area_px <= area <= self.config.max_spot_area_px:
                keep[labels == label] = 255
                component_areas.append(area)

        # Local maxima in the darkness surface = candidate spot centers.
        k = max(3, int(self.config.peak_kernel_size))
        if k % 2 == 0:
            k += 1
        local_max = cv2.dilate(darkness, np.ones((k, k), np.float32))
        peak_mask = (
            (darkness >= local_max - 1e-6)
            & (darkness > threshold)
            & (well_mask > 0)
        ).astype(np.uint8)

        # Collapse flat maxima/plateaus to one event each.
        n_peaks, peak_labels, _, peak_centroids = cv2.connectedComponentsWithStats(
            peak_mask, 8
        )
        peak_xy = []
        for label in range(1, n_peaks):
            x, y = peak_centroids[label]
            peak_xy.append((int(round(x)), int(round(y))))
        peak_xy = np.asarray(peak_xy, dtype=int) if peak_xy else np.empty((0, 2), dtype=int)

        fraction_dark = float(np.count_nonzero(keep)) / max(
            1, np.count_nonzero(well_mask)
        )

        qc = "PASS"
        if well_std < self.config.min_valid_well_std:
            qc = "WARN_LOW_CONTRAST"
        if fraction_dark > self.config.max_fraction_dark:
            qc = "WARN_SATURATED"

        integrated_darkness = float(np.sum(darkness[keep > 0]))
        return (
            keep,
            component_areas,
            peak_xy,
            median_bg,
            integrated_darkness,
            qc,
        )

    def analyze_image(
        self,
        image_path: str | Path,
        output_dir: str | Path,
        plate_id: Optional[str] = None,
    ) -> List[WellResult]:
        image_path = Path(image_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        plate_id = plate_id or image_path.stem

        image = self.load_image(image_path)
        warped, detected_quad = self.rectify_plate(image)
        gray = cv2.cvtColor(warped, cv2.COLOR_BGR2GRAY)

        centers, lattice_diagnostics = self.detect_well_lattice(warped)
        overlay = warped.copy()
        results: List[WellResult] = []

        for r, row in enumerate(ROWS):
            for c, col in enumerate(COLS):
                well = f"{row}{col}"
                center = centers[well]
                mask = self._well_mask(gray.shape, center)

                (
                    spot_mask,
                    component_areas,
                    peak_xy,
                    median_bg,
                    integrated_darkness,
                    qc,
                ) = self._segment_spots(gray, mask)

                # Red point = detected spot center.
                for x, y in peak_xy:
                    cv2.circle(overlay, (int(x), int(y)), 1, (0, 0, 255), -1)

                cx, cy = int(round(center[0])), int(round(center[1]))
                cv2.circle(
                    overlay,
                    (cx, cy),
                    int(round(self.config.well_radius)),
                    (0, 180, 0) if qc == "PASS" else (0, 165, 255),
                    1,
                )

                results.append(
                    WellResult(
                        plate_id=plate_id,
                        well=well,
                        row=row,
                        column=col,
                        spot_count=int(len(peak_xy)),
                        segmented_spot_area_px=int(sum(component_areas)),
                        mean_component_area_px=float(np.mean(component_areas))
                        if component_areas
                        else 0.0,
                        median_component_area_px=float(np.median(component_areas))
                        if component_areas
                        else 0.0,
                        integrated_darkness=integrated_darkness,
                        median_background_intensity=median_bg,
                        well_qc=qc,
                    )
                )

        self._write_csv(results, output_dir / f"{plate_id}_well_results.csv")
        cv2.imwrite(str(output_dir / f"{plate_id}_rectified.png"), warped)
        cv2.imwrite(str(output_dir / f"{plate_id}_annotated.png"), overlay)
        self._write_heatmap(results, output_dir / f"{plate_id}_spot_count_heatmap.png")

        metadata = {
            "input_image": str(image_path),
            "plate_id": plate_id,
            "detected_plate_quad_xy": detected_quad.tolist(),
            "well_lattice": lattice_diagnostics,
            "config": asdict(self.config),
            "n_wells": 96,
        }
        with open(
            output_dir / f"{plate_id}_analysis_metadata.json",
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(metadata, f, indent=2)

        return results

    def analyze_folder(
        self,
        input_dir: str | Path,
        output_dir: str | Path,
        extensions: Sequence[str] = (".png", ".jpg", ".jpeg", ".tif", ".tiff"),
    ) -> Dict[str, List[WellResult]]:
        """Analyze every supported plate image in a directory."""
        input_dir = Path(input_dir)
        output_dir = Path(output_dir)
        all_results: Dict[str, List[WellResult]] = {}

        image_paths = sorted(
            p for p in input_dir.iterdir() if p.suffix.lower() in set(extensions)
        )
        for image_path in image_paths:
            plate_out = output_dir / image_path.stem
            all_results[image_path.stem] = self.analyze_image(
                image_path, plate_out, plate_id=image_path.stem
            )
        return all_results

    @staticmethod
    def _write_csv(results: Sequence[WellResult], path: Path) -> None:
        if not results:
            return
        rows = [asdict(r) for r in results]
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def _write_heatmap(results: Sequence[WellResult], path: Path) -> None:
        """
        Create a dependency-light plate heatmap with OpenCV.
        Values are scaled within each plate for visualization.
        """
        matrix = np.zeros((8, 12), dtype=float)
        for result in results:
            r = ROWS.index(result.row)
            c = result.column - 1
            matrix[r, c] = result.spot_count

        max_value = max(1.0, float(np.max(matrix)))
        cell_w, cell_h = 82, 72
        left, top = 65, 45
        canvas = np.full((top + 8 * cell_h + 55, left + 12 * cell_w + 25, 3), 255, dtype=np.uint8)

        normalized = np.clip(matrix / max_value * 255, 0, 255).astype(np.uint8)
        heat = cv2.applyColorMap(normalized, cv2.COLORMAP_VIRIDIS)

        for r in range(8):
            for c in range(12):
                x1 = left + c * cell_w
                y1 = top + r * cell_h
                x2 = x1 + cell_w - 3
                y2 = y1 + cell_h - 3
                color = tuple(int(v) for v in heat[r, c])
                cv2.rectangle(canvas, (x1, y1), (x2, y2), color, -1)
                value = str(int(matrix[r, c]))
                cv2.putText(
                    canvas, value, (x1 + 8, y1 + 43),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.55, (255, 255, 255), 1, cv2.LINE_AA
                )

        for c in range(12):
            cv2.putText(
                canvas, str(c + 1), (left + c * cell_w + 28, 29),
                cv2.FONT_HERSHEY_SIMPLEX, 0.55, (30, 30, 30), 1, cv2.LINE_AA
            )
        for r, row in enumerate(ROWS):
            cv2.putText(
                canvas, row, (24, top + r * cell_h + 43),
                cv2.FONT_HERSHEY_SIMPLEX, 0.65, (30, 30, 30), 1, cv2.LINE_AA
            )

        cv2.putText(
            canvas,
            "Detected ELISpot events per well",
            (left, canvas.shape[0] - 18),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.65,
            (30, 30, 30),
            1,
            cv2.LINE_AA,
        )
        cv2.imwrite(str(path), canvas)
