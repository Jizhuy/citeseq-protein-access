# Lawlor donor macro-F1 primary definition

Chosen **before** computing totalVI−scVI donor deltas:

**Primary = `present_class_macro_f1` (PRESENT-CLASS)**

Reason: classes absent from a donor are not observed biological events for that
donor; assigning them F1=0 under a fixed ontology is a scoring convention that
can dilute or distort donor-level biology. Present-class macro-F1 averages only
over classes with >0 labeled support in that donor.

Sensitivity = `fixed_ontology_macro_f1` over the predeclared 7 eligible classes
with `zero_division=0` for both models.
