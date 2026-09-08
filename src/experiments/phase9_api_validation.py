"""PHASE 9 API validation: how scvi-tools 1.3.3 totalVI represents missing protein.

This is not the missing-modality experiment. It only documents what the installed
API actually supports so PHASE 9 does not encode missingness as ordinary zeros.
"""

from __future__ import annotations

from datetime import datetime, timezone
from traceback import format_exc
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scvi import REGISTRY_KEYS
from scvi.data.fields import ProteinObsmField
from scvi.model import TOTALVI

from src.utils.io import resolve_path, save_json

def _short_error(text: Any, limit: int = 240) -> str:
    if not text:
        return ""
    first = str(text).splitlines()[0]
    return first if len(first) <= limit else first[: limit - 3] + "..."


def _tiny_adata(n_obs: int = 40, n_genes: int = 20, n_proteins: int = 4, seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    counts = rng.poisson(5, size=(n_obs, n_genes)).astype(np.float32)
    protein = rng.poisson(8, size=(n_obs, n_proteins)).astype(np.float32)
    protein[:, 0] += 3
    adata = ad.AnnData(X=counts.copy())
    adata.layers["counts"] = counts.copy()
    adata.obs_names = [f"c{i}" for i in range(n_obs)]
    adata.var_names = [f"g{j}" for j in range(n_genes)]
    adata.obs["batch"] = np.array(["A"] * (n_obs // 2) + ["B"] * (n_obs - n_obs // 2))
    adata.obsm["protein_expression"] = pd.DataFrame(
        protein, index=adata.obs_names, columns=[f"p{k}" for k in range(n_proteins)]
    )
    return adata


def _protein_loss_split(model: TOTALVI, hide: np.ndarray) -> dict[str, Any]:
    """Measure protein reconstruction loss on hidden vs observed synthetic cells."""
    loader = model._make_data_loader(
        adata=model.adata,
        indices=np.arange(model.adata.n_obs),
        batch_size=model.adata.n_obs,
        shuffle=False,
    )
    tensors = next(iter(loader))
    inference_inputs = model.module._get_inference_input(tensors)
    inference_outputs = model.module.inference(**inference_inputs)
    generative_inputs = model.module._get_generative_input(tensors, inference_outputs)
    generative_outputs = model.module.generative(**generative_inputs)
    y = tensors[REGISTRY_KEYS.PROTEIN_EXP_KEY]
    panel_index = tensors[model.module.panel_key]
    if model.module.protein_batch_mask is not None:
        pro_batch_mask_minibatch = torch.zeros_like(y)
        for batch_id in torch.unique(panel_index):
            batch_rows = (panel_index == batch_id).reshape(-1)
            pro_batch_mask_minibatch[batch_rows] = torch.tensor(
                model.module.protein_batch_mask[str(int(batch_id.item()))].astype(np.float32),
                device=y.device,
            )
    else:
        pro_batch_mask_minibatch = None
    _, protein_loss = model.module.get_reconstruction_loss(
        tensors[REGISTRY_KEYS.X_KEY],
        y,
        generative_outputs["px_"],
        generative_outputs["py_"],
        pro_batch_mask_minibatch,
        generative_outputs["per_batch_efficiency"],
    )
    protein_loss = protein_loss.detach().cpu().numpy()
    y_np = y.detach().cpu().numpy()
    hide = np.asarray(hide, dtype=bool)
    return {
        "mean_protein_reconst_loss_hidden_cells": float(protein_loss[hide].mean()),
        "mean_protein_reconst_loss_observed_cells": float(protein_loss[~hide].mean()),
        "protein_loss_exactly_zero_for_all_hidden_cells": bool(np.allclose(protein_loss[hide], 0.0)),
        "encoder_protein_input_is_all_zero_for_hidden_cells": bool(np.allclose(y_np[hide], 0.0)),
        "n_hidden_cells": int(hide.sum()),
        "n_observed_cells": int((~hide).sum()),
    }


def _inspect_setup(
    adata: ad.AnnData,
    panel_key: str | None = None,
    hide: np.ndarray | None = None,
) -> dict[str, Any]:
    TOTALVI.setup_anndata(
        adata,
        protein_expression_obsm_key="protein_expression",
        batch_key="batch",
        layer="counts",
        panel_key=panel_key,
    )
    model = TOTALVI(adata, n_latent=4, n_hidden=16, n_layers_encoder=1, n_layers_decoder=1)
    registry = model.adata_manager.get_state_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY)
    mask = registry.get(ProteinObsmField.PROTEIN_BATCH_MASK)
    if mask is not None:
        mask = {k: np.asarray(v).astype(bool).tolist() for k, v in dict(mask).items()}
    y = np.asarray(adata.obsm["protein_expression"], dtype=float)
    out: dict[str, Any] = {
        "setup_ok": True,
        "n_nan_in_protein": int(np.isnan(y).sum()) if np.issubdtype(y.dtype, np.floating) else 0,
        "protein_batch_mask": mask,
        "use_adversarial_classifier": bool(model._use_adversarial_classifier),
        "encoder_concatenates_protein_counts": True,
        "note": (
            "TOTALVAE._regular_inference always concatenates RNA and protein "
            "into the encoder: encoder_input = cat(x_, y_). The batch/panel mask "
            "only zeros protein reconstruction loss, not encoder input."
        ),
    }
    if hide is not None:
        out["forward_pass"] = _protein_loss_split(model, hide)
    try:
        latent = np.asarray(model.get_latent_representation())
        out["latent_finite"] = bool(np.isfinite(latent).all())
        out["latent_shape"] = list(latent.shape)
    except Exception as exc:  # noqa: BLE001
        out["latent_finite"] = False
        out["latent_error"] = f"{type(exc).__name__}: {exc}"
    return out


def run_api_validation() -> dict[str, Any]:
    hide = np.zeros(40, dtype=bool)
    hide[:5] = True
    hide[20:25] = True

    tests: dict[str, Any] = {}

    # A. NaN cell-level missingness
    adata_nan = _tiny_adata()
    prot = adata_nan.obsm["protein_expression"].copy()
    prot.loc[hide, :] = np.nan
    adata_nan.obsm["protein_expression"] = prot
    try:
        tests["A_nan_cellwise"] = _inspect_setup(adata_nan, hide=hide)
        latent_ok = tests["A_nan_cellwise"].get("latent_finite")
        tests["A_nan_cellwise"]["interpreted_as"] = (
            "setup did not hard-fail, but NaN proteins are not a documented missingness encoding"
            if latent_ok
            else "setup did not hard-fail; latent extraction/training is not valid with NaN proteins"
        )
    except Exception as exc:  # noqa: BLE001
        tests["A_nan_cellwise"] = {
            "setup_ok": False,
            "error": f"{type(exc).__name__}: {str(exc).splitlines()[0][:300]}",
            "traceback": format_exc(limit=8),
            "interpreted_as": "not supported; NaNs are not a valid totalVI protein missingness encoding",
        }

    # B. Zero-fill cells inside existing batches (no panel)
    adata_zero = _tiny_adata()
    prot = adata_zero.obsm["protein_expression"].copy()
    prot.loc[hide, :] = 0
    adata_zero.obsm["protein_expression"] = prot
    try:
        tests["B_zero_fill_same_batch"] = _inspect_setup(adata_zero, hide=hide)
        mask_b = tests["B_zero_fill_same_batch"]["protein_batch_mask"]
        hidden_loss_nonzero = not tests["B_zero_fill_same_batch"]["forward_pass"][
            "protein_loss_exactly_zero_for_all_hidden_cells"
        ]
        tests["B_zero_fill_same_batch"]["interpreted_as"] = (
            "observed zero counts (protein reconstruction loss is applied to hidden cells)"
            if mask_b is None and hidden_loss_nonzero
            else f"batch mask={mask_b}; hidden_loss_zero={not hidden_loss_nonzero}"
        )
    except Exception as exc:  # noqa: BLE001
        tests["B_zero_fill_same_batch"] = {
            "setup_ok": False,
            "error": f"{type(exc).__name__}: {str(exc).splitlines()[0][:300]}",
            "traceback": format_exc(limit=8),
        }

    # C. Dedicated protein-missing panel (documented unmatched-panel mechanism)
    adata_panel = _tiny_adata()
    prot = adata_panel.obsm["protein_expression"].copy()
    prot.loc[hide, :] = 0
    adata_panel.obsm["protein_expression"] = prot
    adata_panel.obs["protein_panel"] = np.where(hide, "unmeasured", adata_panel.obs["batch"].astype(str))
    try:
        tests["C_zero_fill_with_panel_key"] = _inspect_setup(
            adata_panel, panel_key="protein_panel", hide=hide
        )
        tests["C_zero_fill_with_panel_key"]["interpreted_as"] = (
            "panel-level missing proteins if protein_batch_mask is False for the unmeasured panel; "
            "encoder still receives zeros"
        )
        tests["C_zero_fill_with_panel_key"]["panel_counts"] = (
            adata_panel.obs["protein_panel"].value_counts().to_dict()
        )
    except Exception as exc:  # noqa: BLE001
        tests["C_zero_fill_with_panel_key"] = {
            "setup_ok": False,
            "error": f"{type(exc).__name__}: {str(exc).splitlines()[0][:300]}",
            "traceback": format_exc(limit=8),
        }

    b_ok = tests.get("B_zero_fill_same_batch", {})
    summary = {
        "scvi_tools_version": __import__("scvi").__version__,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "native_missingness_grain": "batch or protein panel, not cell",
        "native_definition": (
            "A protein is missing in a batch/panel iff it is all zeros in every cell of that "
            "batch/panel (scvi.data.fields._protein.ProteinFieldMixin._get_batch_mask_protein_data)."
        ),
        "cell_level_nan_supported": bool(tests["A_nan_cellwise"].get("setup_ok")),
        "cell_level_zero_fill_without_panel_is_observed_zero": (
            b_ok.get("setup_ok") is True and b_ok.get("protein_batch_mask") is None
        ),
        "encoder_always_sees_protein_tensor": True,
        "tests": tests,
        "conclusion": (
            "scvi-tools 1.3.3 totalVI does not natively support cell-level missing protein "
            "profiles inside a CITE-seq batch. Zero-filling selected cells is interpreted as "
            "observed zeros (PHASE 8), not missing modality. NaNs are not a supported encoding. "
            "The only documented missingness mechanism is unmatched protein panels / batches."
        ),
        "safest_supported_alternative": (
            "Assign missing-protein cells a dedicated protein_panel category via TOTALVI.setup_anndata "
            "panel_key, and set those cells' protein matrix to all zeros so the panel is registered as "
            "unmeasured. This masks protein reconstruction loss for that panel. It is still not true "
            "cell-level missingness: the encoder concatenates those zeros, empirical priors for the "
            "unmeasured panel are random, and an adversarial classifier is enabled. Do not run this "
            "mapping without explicit approval."
        ),
        "phase9_experiment_started": False,
    }
    return summary


def main() -> None:
    result = run_api_validation()
    logs = resolve_path("results/logs")
    save_json(result, logs / "phase9_api_validation.json")
    md = logs / "phase9_api_validation.md"
    tests = result["tests"]
    md.write_text(
        f"""# PHASE 9 API validation (scvi-tools {result['scvi_tools_version']})

Date (UTC): {result['timestamp_utc']}

## Native missingness mechanism

{result['native_definition']}

Grain: **{result['native_missingness_grain']}**.

The protein encoder always concatenates RNA and protein (`encoder_input = cat(x_, y_)`).
The batch/panel mask only multiplies **protein reconstruction loss**. It does not drop protein
from the encoder.

## Synthetic tests (40 cells, 4 proteins, hide 10 cells)

### A. Cell-wise NaN
- setup_ok: {tests['A_nan_cellwise'].get('setup_ok')}
- latent_finite: {tests['A_nan_cellwise'].get('latent_finite')}
- forward_pass: {tests['A_nan_cellwise'].get('forward_pass')}
- result: {_short_error(tests['A_nan_cellwise'].get('error') or tests['A_nan_cellwise'].get('latent_error') or tests['A_nan_cellwise'].get('interpreted_as'))}

### B. Cell-wise zeros in the original batches (no panel_key)
- protein_batch_mask: {tests['B_zero_fill_same_batch'].get('protein_batch_mask')}
- forward_pass: {tests['B_zero_fill_same_batch'].get('forward_pass')}
- interpreted as: **{tests['B_zero_fill_same_batch'].get('interpreted_as')}**

### C. Cell-wise zeros plus a dedicated `panel_key` category
- protein_batch_mask: {tests['C_zero_fill_with_panel_key'].get('protein_batch_mask')}
- adversarial classifier: {tests['C_zero_fill_with_panel_key'].get('use_adversarial_classifier')}
- forward_pass: {tests['C_zero_fill_with_panel_key'].get('forward_pass')}
- interpreted as: {tests['C_zero_fill_with_panel_key'].get('interpreted_as')}

## Conclusion

{result['conclusion']}

PHASE 9 training was **not started**.

## Safest supported alternative (not executed)

{result['safest_supported_alternative']}
""",
        encoding="utf-8",
    )
    print(md.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
