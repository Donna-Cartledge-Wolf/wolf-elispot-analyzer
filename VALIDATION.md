# Validation Notes

## Synthetic ground-truth benchmark

Version 0.1 was tested against two synthetic 96-well ELISpot plates containing **192 total wells** with known generated spot counts.

Baseline performance:

- **Pearson correlation:** 0.922
- **Mean absolute error (MAE):** 24.39 spots/well
- **Median absolute error:** 14.0 spots/well

![Synthetic validation scatter](assets/validation_scatter.png)

## Interpretation

The current local-intensity-peak algorithm captures the **relative biological signal pattern well**, as shown by the strong correlation between generated and detected spot counts.

Absolute spot counts are less accurate in densely populated wells. This is expected for the current prototype because overlapping or partially merged signals can produce fewer resolvable local intensity maxima than the number of originally generated spots.

This limitation is scientifically useful: it identifies the next algorithm-development target rather than masking the problem with hand-tuned counts.

## Planned improvement path

1. Add multiscale blob detection.
2. Add watershed/deconvolution for touching spots.
3. Use negative-control wells to estimate plate-specific detection thresholds.
4. Add TNTC / saturated-well classification.
5. Compare automatic results with manually reviewed reference counts.
6. Validate on real reader-generated ELISpot images from multiple imaging conditions.

The synthetic benchmark is intended for algorithm development only and does not constitute analytical-method validation.
