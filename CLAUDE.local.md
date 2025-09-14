# DoubleML Survival Analysis Integration

## Project Overview
Integrating survival analysis and competing risk methods from `/original` into the DoubleML R package as a contributor.

## Development Environment
- Working Directory: `/Users/xushenbo/Downloads/Github/DoubleML`
- Branch: `survival_outcome_integration`
- Target Package: `doubleml-for-r-main/`

## Implementation Plan

### File Structure
```
doubleml-for-r-main/
├── R/
│   ├── double_ml_survival.R        # Base survival class
│   ├── double_ml_survival_prob.R   # Survival probability estimation
│   ├── double_ml_rmst.R           # Restricted mean survival time
│   ├── double_ml_competing.R      # Base competing risks class
│   ├── double_ml_cifj.R           # Cumulative incidence function
│   ├── double_ml_rmtlj.R          # Restricted mean time lost
│   └── helper_survival.R          # Survival-specific helpers
├── tests/testthat/
│   ├── test-double_ml_survival_prob.R
│   ├── test-double_ml_rmst.R
│   ├── test-double_ml_cifj.R
│   └── test-double_ml_rmtlj.R
└── inst/examples/
    ├── survival_mlr3_example.R
    └── survival_superlearner_example.R
```

### Key Design Decisions
1. **Dual Pipeline Support**: Users can choose between mlr3proba and survSuperLearner
2. **Preserve Original Methods**: Keep statistical methods from `original/estimators/` unchanged
3. **R6 Class Structure**: Adapt to DoubleML's class hierarchy while minimizing changes
4. **No Score Functions**: Maintain original influence function calculations without conversion

### ML Framework Options
- **mlr3 framework**: Uses mlr3 for treatment model, mlr3proba for survival/censoring models
- **SuperLearner framework**: Uses SuperLearner for treatment model, survSuperLearner for survival/censoring models

### Source Files to Migrate
- `original/estimators/double_ml_survival.R`: Contains `double_ml_survival_prob()`, `double_ml_rmst()`
- `original/estimators/double_ml_competing.R`: Contains `double_ml_cifj()`, `double_ml_rmtlj()`

## Development Commands

## Notes
- Do not modify existing files in doubleml-for-r-main/, consider adding
- Do not modify existing files in original/
- Maintain compatibility with existing DoubleML API
- Focus on class structure adaptation, not algorithm changes
- **All tests should run code under /Users/xushenbo/Downloads/Github/DoubleML/doubleml-for-r-main**

## Contributing Policy

From `doubleml-for-r-main/CONTRIBUTING.md`:

### Code Requirements
- **Code passes `R CMD check`** and all unit tests: `devtools::check()`
- **Unit tests mandatory** for enhancements/new features: `devtools::test()`

### New Model Class Contributions
- **Abstract base class `DoubleML`** implements core functionality via linear Neyman orthogonal scores
- **Only need to specify**:
  - All nuisance functions for estimation (regressions/classifications)
  - Score components for Neyman orthogonal score function
- **All other functionality inherited** automatically from base class

### Documentation
- **API docs**: Generated with **roxygen2** and **pkgdown**
- **Build docs**: `devtools::document()`
- **Build website**: `pkgdown::build_site()`
- **User guide/examples**: Separate repo [doubleml-docs](https://github.com/DoubleML/doubleml-docs)