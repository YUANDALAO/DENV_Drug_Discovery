# Reward Shaping Methodology for QSAR-Guided Molecular Generation

## Problem Statement
When using QSAR models (trained on biased datasets with mean pIC50 ≈ 7.35) as reward functions in reinforcement learning-based molecular generation (REINVENT4), the generated molecules converge prematurely to the training set mean due to:
1. **Flat reward landscape:** Raw pIC50 differences (7.5→8.0, Δ=0.05) are too small after geometric averaging
2. **Regression to mean:** Random Forest predictions for novel structures default to ~7.5
3. **Insufficient gradient:** Model sees no incentive to explore beyond "good enough" (pIC50≈7.5)

## Solution: Non-linear Reward Transformation

### Sigmoid Transform
$$
f(x) = \frac{1}{1 + e^{-k(x - m)}}
$$
where:
- $x$ = predicted pIC50
- $m$ = midpoint $= \frac{\text{low} + \text{high}}{2}$
- $k$ = steepness parameter
- $\text{low}, \text{high}$ = target pIC50 range

### Two-Stage Strategy
| Stage | low | high | k   | Purpose |
|-------|-----|------|-----|---------|
| 1     | 7.0 | 8.5  | 2.0 | Encourage exploration from medium→high activity |
| 2     | 7.5 | 9.0  | 2.5 | Focus on exceptional activity (pIC50>8.5) |

### Numerical Example
| pIC50 | Raw Score | Stage 1 Transform | Stage 2 Transform | Gradient Amplification |
|-------|-----------|-------------------|-------------------|------------------------|
| 7.0   | 0.700     | 0.119             | 0.018             | -                      |
| 7.5   | 0.750     | 0.500             | 0.215             | -                      |
| 8.0   | 0.800     | 0.731             | 0.500             | 4.6× (Stage 1)         |
| 8.5   | 0.850     | 0.881             | 0.785             | 5.7× (Stage 2)         |
| 9.0   | 0.900     | 0.952             | 0.920             | -                      |

**Impact:** The score difference for pIC50 7.5→8.0 increases from 0.05 to 0.23 (Stage 1), making optimization 4.6× more sensitive.

## Synergistic Hyperparameter Tuning

### Learning Rate Reduction
- **Before:** `rate = 0.0001` → Overshoots in steep gradient regions
- **After:** `rate = 0.00003` → Stable exploration with transformed rewards
- **Rationale:** Lower LR + higher gradient = controlled progress toward pIC50>8.5

### Sigma Reduction (DAP Algorithm)
- **Before:** `sigma = 120` → Reward ≈ linear function of score
- **After:** `sigma = 60` → Reward ≈ exponential, highly sensitive to high scores
- **Formula:** $\text{Reward} = e^{\text{score}/\sigma}$
- **Effect:** Amplifies preference for pIC50>8.0 molecules

### Synergy Equation
$$
\text{Effective Gradient} = \underbrace{\text{Sigmoid Stretch}}_{\times 4.6} \times \underbrace{\text{Low } \sigma \text{ Boost}}_{\times 2} \times \underbrace{\text{Stable LR}}_{\text{controlled}}
$$

## Implementation (REINVENT4 TOML)
```toml
[learning_strategy]
type = "dap"
sigma = 60              # Amplify high-score preference
rate = 0.00003          # Stable updates with steep gradients

[[stage.scoring.component]]
[stage.scoring.component.QSARScorer]
[[stage.scoring.component.QSARScorer.endpoint]]
name = "DENV_Activity"
weight = 0.80
params.model_path = "models/random_forest_tuned_model.joblib"

transform.type = "sigmoid"
transform.low = 7.0     # Lower bound (low activity)
transform.high = 8.5    # Upper bound (high activity)
transform.k = 2.0       # Steepness (higher = sharper transition)