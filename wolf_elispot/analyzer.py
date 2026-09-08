from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Dict, List
import json
import cv2
import numpy as np
import pandas as pd

ROWS = "ABCDEFGH"
COLS = list(range(1, 13))

@dataclass
class WellResult:
    plate_id: str
    well: str
    row: str
    column: int
    spot_count: int
    resolved_spot_count: int
    signal_equivalent_count: float
    count_method: str
    segmented_spot_area_px: int
    mean_component_area_px: float
    median_component_area_px: float
    chromatic_signal: float
    integrated_darkness: float
    median_background_intensity: float
    well_qc: str

@dataclass
class AnalysisConfig:
    warp_width: int = 1278
    warp_height: int = 855

    first_well_x: float = 144.0
    first_well_y: float = 112.0
    pitch_x: float = 90.0
    pitch_y: float = 90.0
    well_radius: float = 29.0

    use_hough_well_detection: bool = True
    hough_dp: float = 1.2
    hough_min_dist: float = 55.0
    hough_param1: float = 100.0
    hough_param2: float = 28.0
    hough_min_radius: int = 23
    hough_max_radius: int = 40
    lattice_cluster_tolerance_px: float = 24.0

    # v0.2 sparse-well peak detection
    spot_peak_kernel_size: int = 3
    minimum_spot_saturation: float = 15.0
    saturation_above_background: float = 10.0
    chromatic_background_percentile: float = 5.0
    chromatic_background_offset: float = 3.0
    resolved_count_limit: int = 20

    # plate-adaptive crowded-well calibration
    calibration_min_resolved: int = 2
    calibration_max_resolved: int = 15
    fallback_signal_per_spot: float = 1180.0

    # overlap/crowding correction; developed on synthetic Plate 1
    crowding_capacity: float = 800.0

    min_component_area_px: int = 1
    max_component_area_px: int = 180
    border_margin_px: int = 3
    min_valid_well_std: float = 1.5
    max_fraction_signal: float = 0.55

class ELISpotAnalyzer:
    """Wolf ELISpot Analyzer v0.2 hybrid spot counter."""

    def __init__(self, config: Optional[AnalysisConfig] = None):
        self.config = config or AnalysisConfig()

    @staticmethod
    def load_image(path):
        image = cv2.imread(str(path), cv2.IMREAD_COLOR)
        if image is None:
            raise FileNotFoundError(f"Could not read image: {path}")
        return image

    @staticmethod
    def _order_points(pts):
        pts = np.asarray(pts, dtype=np.float32)
        s = pts.sum(axis=1)
        d = np.diff(pts, axis=1).ravel()
        ordered = np.zeros((4, 2), dtype=np.float32)
        ordered[0] = pts[np.argmin(s)]
        ordered[2] = pts[np.argmax(s)]
        ordered[1] = pts[np.argmin(d)]
        ordered[3] = pts[np.argmax(d)]
        return ordered

    def detect_plate_quad(self, image):
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
            candidates.append((abs(ratio - 1.495), -area, approx.reshape(4, 2)))

        if candidates:
            candidates.sort(key=lambda x: (x[0], x[1]))
            q = self._order_points(candidates[0][2])
            if cv2.contourArea(q.astype(np.float32)) / image_area > 0.80:
                return np.array([[0, 0], [w - 1, 0], [w - 1, h - 1], [0, h - 1]], dtype=np.float32)
            return q

        return np.array([[0, 0], [w - 1, 0], [w - 1, h - 1], [0, h - 1]], dtype=np.float32)

    def rectify_plate(self, image):
        quad = self.detect_plate_quad(image)
        dst = np.array([
            [0, 0],
            [self.config.warp_width - 1, 0],
            [self.config.warp_width - 1, self.config.warp_height - 1],
            [0, self.config.warp_height - 1],
        ], dtype=np.float32)
        matrix = cv2.getPerspectiveTransform(quad.astype(np.float32), dst)
        warped = cv2.warpPerspective(image, matrix, (self.config.warp_width, self.config.warp_height))
        return warped, quad

    @staticmethod
    def _cluster_1d(values, tolerance):
        if len(values) == 0:
            return []
        groups = []
        for value in sorted(float(v) for v in values):
            if not groups or value - float(np.mean(groups[-1])) > tolerance:
                groups.append([value])
            else:
                groups[-1].append(value)
        return [(float(np.mean(g)), len(g)) for g in groups]

    def well_centers_fallback(self):
        centers = {}
        for r, row in enumerate(ROWS):
            for c, col in enumerate(COLS):
                centers[f"{row}{col}"] = (
                    self.config.first_well_x + c * self.config.pitch_x,
                    self.config.first_well_y + r * self.config.pitch_y,
                )
        return centers

    def detect_well_lattice(self, plate_image):
        fallback = self.well_centers_fallback()
        diag = {"method": "fallback_standard_geometry", "detected_circle_count": 0}

        if not self.config.use_hough_well_detection:
            return fallback, diag

        gray = cv2.cvtColor(plate_image, cv2.COLOR_BGR2GRAY)
        gray = cv2.GaussianBlur(gray, (5, 5), 1.2)
        circles = cv2.HoughCircles(
            gray, cv2.HOUGH_GRADIENT,
            dp=self.config.hough_dp,
            minDist=self.config.hough_min_dist,
            param1=self.config.hough_param1,
            param2=self.config.hough_param2,
            minRadius=self.config.hough_min_radius,
            maxRadius=self.config.hough_max_radius,
        )
        if circles is None:
            return fallback, diag

        circles = np.asarray(circles[0], dtype=float)
        diag["detected_circle_count"] = int(len(circles))
        xc = self._cluster_1d(circles[:, 0], self.config.lattice_cluster_tolerance_px)
        yc = self._cluster_1d(circles[:, 1], self.config.lattice_cluster_tolerance_px)

        xs = [x for x, n in xc if n >= 4]
        ys = [y for y, n in yc if n >= 6]

        if len(xs) != 12:
            xs = sorted(x for x, _ in sorted(xc, key=lambda t: (-t[1], t[0]))[:12])
        else:
            xs = sorted(xs)
        if len(ys) != 8:
            ys = sorted(y for y, _ in sorted(yc, key=lambda t: (-t[1], t[0]))[:8])
        else:
            ys = sorted(ys)

        if len(xs) != 12 or len(ys) != 8:
            return fallback, diag

        if np.min(np.diff(xs)) < 50 or np.max(np.diff(xs)) > 130 or np.min(np.diff(ys)) < 50 or np.max(np.diff(ys)) > 130:
            return fallback, diag

        centers = {}
        for r, row in enumerate(ROWS):
            for c, col in enumerate(COLS):
                centers[f"{row}{col}"] = (float(xs[c]), float(ys[r]))

        diag.update({
            "method": "automatic_hough_lattice",
            "median_pitch_x_px": float(np.median(np.diff(xs))),
            "median_pitch_y_px": float(np.median(np.diff(ys))),
        })
        return centers, diag

    def _raw_well_features(self, image, center):
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY).astype(np.float32)
        hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
        sat = hsv[:, :, 1].astype(np.float32)

        H, W = gray.shape
        Y, X = np.ogrid[:H, :W]
        radius = self.config.well_radius - self.config.border_margin_px
        mask = ((X - center[0]) ** 2 + (Y - center[1]) ** 2 <= radius ** 2)

        s = sat[mask]
        g = gray[mask]
        sat_bg = float(np.percentile(s, self.config.chromatic_background_percentile))
        gray_bg = float(np.percentile(g, 80))
        well_std = float(np.std(g))

        chromatic_signal = float(
            np.clip(s - (sat_bg + self.config.chromatic_background_offset), 0, None).sum()
        )

        k = max(3, int(self.config.spot_peak_kernel_size))
        if k % 2 == 0:
            k += 1
        local_max = cv2.dilate(sat, np.ones((k, k), np.float32))
        peak_threshold = max(
            self.config.minimum_spot_saturation,
            sat_bg + self.config.saturation_above_background,
        )
        peak_mask = ((sat >= local_max - 1e-6) & (sat > peak_threshold) & mask).astype(np.uint8)

        n, labels, stats, centroids = cv2.connectedComponentsWithStats(peak_mask, 8)
        peaks = []
        for lab in range(1, n):
            x, y = centroids[lab]
            peaks.append((int(round(x)), int(round(y))))
        peaks = np.asarray(peaks, dtype=int) if peaks else np.empty((0, 2), dtype=int)

        binary = ((sat > peak_threshold) & mask).astype(np.uint8) * 255
        n2, lab2, stats2, _ = cv2.connectedComponentsWithStats(binary, 8)
        areas = []
        for lab in range(1, n2):
            area = int(stats2[lab, cv2.CC_STAT_AREA])
            if self.config.min_component_area_px <= area <= self.config.max_component_area_px:
                areas.append(area)

        darkness = np.clip(gray_bg - gray, 0, None)
        fraction_signal = float(np.count_nonzero(binary)) / max(1, np.count_nonzero(mask))

        qc = "PASS"
        if well_std < self.config.min_valid_well_std:
            qc = "WARN_LOW_CONTRAST"
        if fraction_signal > self.config.max_fraction_signal:
            qc = "WARN_CROWDED"

        return {
            "resolved": int(len(peaks)),
            "peaks": peaks,
            "chromatic_signal": chromatic_signal,
            "areas": areas,
            "integrated_darkness": float(darkness[mask].sum()),
            "median_background_intensity": float(np.median(g)),
            "qc": qc,
        }

    def _estimate_signal_per_spot(self, raw):
        ratios = []
        for d in raw.values():
            n = d["resolved"]
            if self.config.calibration_min_resolved <= n <= self.config.calibration_max_resolved and d["chromatic_signal"] > 0:
                ratios.append(d["chromatic_signal"] / n)

        if len(ratios) >= 5:
            return float(np.median(ratios)), len(ratios), "plate_adaptive"

        return float(self.config.fallback_signal_per_spot), len(ratios), "fallback_default"

    def _crowding_correct(self, equivalent):
        K = float(self.config.crowding_capacity)
        frac = min(max(equivalent / K, 0.0), 0.95)
        return float(-K * np.log(1.0 - frac))

    def analyze_image(self, image_path, output_dir, plate_id=None):
        image_path = Path(image_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        plate_id = plate_id or image_path.stem

        image = self.load_image(image_path)
        warped, quad = self.rectify_plate(image)
        centers, lattice = self.detect_well_lattice(warped)

        raw = {}
        for row in ROWS:
            for col in COLS:
                well = f"{row}{col}"
                raw[well] = self._raw_well_features(warped, centers[well])

        signal_per_spot, n_cal, cal_method = self._estimate_signal_per_spot(raw)

        overlay = warped.copy()
        results = []

        for row in ROWS:
            for col in COLS:
                well = f"{row}{col}"
                d = raw[well]
                resolved = d["resolved"]
                equivalent = d["chromatic_signal"] / max(signal_per_spot, 1e-9)

                if resolved <= self.config.resolved_count_limit:
                    final_count = resolved
                    method = "resolved_peaks"
                else:
                    final_count = int(round(self._crowding_correct(equivalent)))
                    method = "signal_equivalent_crowding_corrected"

                for x, y in d["peaks"]:
                    cv2.circle(overlay, (int(x), int(y)), 1, (0, 0, 255), -1)

                cx, cy = map(lambda z: int(round(z)), centers[well])
                cv2.circle(overlay, (cx, cy), int(round(self.config.well_radius)), (0, 180, 0), 1)

                areas = d["areas"]
                results.append(WellResult(
                    plate_id=plate_id,
                    well=well,
                    row=row,
                    column=col,
                    spot_count=final_count,
                    resolved_spot_count=resolved,
                    signal_equivalent_count=float(equivalent),
                    count_method=method,
                    segmented_spot_area_px=int(sum(areas)),
                    mean_component_area_px=float(np.mean(areas)) if areas else 0.0,
                    median_component_area_px=float(np.median(areas)) if areas else 0.0,
                    chromatic_signal=float(d["chromatic_signal"]),
                    integrated_darkness=float(d["integrated_darkness"]),
                    median_background_intensity=float(d["median_background_intensity"]),
                    well_qc=d["qc"],
                ))

        pd.DataFrame([asdict(r) for r in results]).to_csv(
            output_dir / f"{plate_id}_well_results.csv", index=False
        )
        cv2.imwrite(str(output_dir / f"{plate_id}_annotated.png"), overlay)
        cv2.imwrite(str(output_dir / f"{plate_id}_rectified.png"), warped)

        metadata = {
            "version": "0.2.0",
            "input_image": str(image_path),
            "plate_id": plate_id,
            "signal_per_spot": signal_per_spot,
            "calibration_wells_used": n_cal,
            "calibration_method": cal_method,
            "well_lattice": lattice,
            "config": asdict(self.config),
        }
        with open(output_dir / f"{plate_id}_analysis_metadata.json", "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)

        return results
