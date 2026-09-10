# Reviewer contingency experiments (NOT currently required)

All items below are **optional post-review contingencies**. None are required for the frozen supervisor package.

## A. Baseline → LPS / CD3_CD28 condition shift (Lawlor)
- **Question:** Does reliability generalize from donor shift to stimulation-condition shift?
- **Compute:** Roughly similar to a reduced Lawlor grid (estimate: tens of GPU-hours depending on seed/fold subset)
- **Value:** High for breadth of “distribution shift”
- **Reviewer concern answered:** PBMC Baseline-only / unused stimulation conditions
- **Status:** HIGH-VALUE OPTIONAL — **NOT CURRENTLY REQUIRED**

## B. 39-ADT sensitivity (Lawlor full antibody panel)
- **Question:** Do conclusions change when totalVI sees the broader Lawlor protein panel?
- **Compute:** Moderate (retrain totalVI arm; scVI RNA arm largely unchanged)
- **Value:** Low–moderate for the matched-panel external claim; useful sensitivity note
- **Reviewer concern answered:** 12/14 protein overlap limitation
- **Status:** LOW-VALUE OPTIONAL — **NOT CURRENTLY REQUIRED**

## C. External source-composition perturbation
- **Question:** Does uncertainty predict fragility under source-subset perturbation on Lawlor?
- **Compute:** High (subset × seed grid on external data)
- **Value:** Medium; closes an explicitly “not testable externally” gap
- **Reviewer concern answered:** Source-composition fragility not externally replicated
- **Status:** OPTIONAL — **NOT CURRENTLY REQUIRED**
