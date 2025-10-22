# DoubleML Survival Analysis Integration

## Project Overview
Extending survival analysis and competing risk methods from `/original` into the DoubleML R package.

## Environment
- Working Directory: `/Users/xushenbo/Downloads/Github/DoubleML`
- Branch: `survival_outcome_integration`
- Target Package: `doubleml-for-r-main/`

## Key Constraints
- **ML Framework**: SuperLearner (treatment), survSuperLearner (survival/censoring) - mlr3 NOT ALLOWED
- **No placeholders**: Always use real ML, never dummy values like `matrix(0.5, n, m)`
- **Data Generation**: Use EXACT DGP from `@original/simulation/survival` - do NOT create new/custom DGP
- **Do NOT create code files unless explicitly requested**
- **Preserve original implementations**: Keep function signatures, logic flow, and mathematical formulas unchanged

## Development Priorities
- Test that extended code actually works - functionality over form
- Follow existing patterns exactly, don't be imaginative or creative
- Use real ML implementations in tests (no placeholders)
- Check `@original/simulation` for test parameters (sample size, admin.cens, taus, etc.)

## Implementation Status
**✅ Survival Analysis (Complete):**
- `DoubleMLSurvival` (base class)
- `DoubleMLSurvivalProb` (survival probability)
- `DoubleMLRMST` (restricted mean survival time)
- Inference infrastructure: `inference_point.R`, `inference_curve.R`
- ML nuisance functions: `survival_nuisance.R` (combined S and G estimation)

**❌ Competing Risks (Pending):**
- `DoubleMLCompeting` (base class - placeholder)
- `DoubleMLCIFJ` (cumulative incidence - placeholder)
- `DoubleMLRMTLJ` (restricted mean time lost - placeholder)

## Notes
- All tests run code under `/Users/xushenbo/Downloads/Github/DoubleML/doubleml-for-r-main`
- Do not modify existing files in `original/`
- Maintain compatibility with existing DoubleML API