# Missing-modality feasibility (PART 11)

## Historical mistake (do not repeat)

Zero-filling cell-wise proteins and calling that “missing modality” is invalid
for totalVI 1.3.3 because the API does not provide arbitrary cell-level protein
missingness masks with correct likelihood semantics.

## totalVI 1.3.3

Prior project audit (`scripts/09_validate_missing_modality_api.py` / PHASE9 logs):
**no honest cell-wise missing-protein likelihood mask** for arbitrary patterns.

Therefore:

- **Cell-wise missing protein evaluation: UNSUPPORTED** for totalVI 1.3.3.
- Do not fake it.

## Valid options

A. Panel-level protein missingness **only if** API genuinely masks likelihood — not confirmed as general cell-wise.
B. **Antibody-panel dropout** (remove antibodies / random 25–50% of panel) — feasible, interpretable, not cell-wise missingness.
C. A modern model with explicit missing-modality support — not yet selected/fairly configured.

## Antibody-panel dropout status

**Designed but not yet executed** in this revision round (requires GPU training
on RTX 4060 with frozen clean models). Script scaffold to be added under
`corruption_shift/` / `controls/` before any run.

## Cell-wise missing protein

**Documented as unsupported** under the current totalVI pin.
