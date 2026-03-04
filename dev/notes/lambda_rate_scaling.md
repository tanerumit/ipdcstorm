## Lambda rate scaling calibration

- Scope: annual event-rate calibration only.
- No windfield physics, event severity generation, storm filters, gate definitions, or at-site TS/HUR classification logic changed.
- For each site and validation class (`TS34plus`, `HUR`), the model now derives a deterministic multiplicative scaler from the existing rate-check table.
- Target lambda is `lambda_ref * expected_ratio`, where `expected_ratio` is the existing gate-approximation adjustment already used by validation.
- Scalers are clamped to `[0.25, 4.0]` to avoid pathological rate inflation/deflation.
- If `lambda_ref` is missing, the scaler defaults to `1.0` with `scale_status = "no_ref"`.
- Runtime simulation converts the calibrated `TS34plus` and `HUR` targets back into the internal `TS` and `HUR` lambdas:
  - `lambda_total_adj = (lambda_TS + lambda_HUR) * scale_TS34plus`
  - `lambda_HUR_adj = lambda_HUR * scale_HUR`
  - `lambda_TS_adj = max(0, lambda_total_adj - lambda_HUR_adj)`
- If the adjusted hurricane lambda would exceed the adjusted total lambda, it is capped at the adjusted total so the implied TS lambda remains non-negative.
