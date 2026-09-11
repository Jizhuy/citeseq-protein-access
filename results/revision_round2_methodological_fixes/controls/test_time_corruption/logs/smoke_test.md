# Smoke test PASS

```json
{
  "direction_B": {
    "n_rows": 3,
    "meta": {
      "direction": "direction_B",
      "seed": 0,
      "source_batch": "PBMC5k",
      "target_batch": "PBMC10k",
      "n_eligible_classes": 21,
      "audit": {
        "source_n": 3994,
        "target_n": 6855,
        "n_genes": 15792,
        "n_proteins": 14
      },
      "models": {
        "totalvi": {
          "clean_macro_f1": 0.6157476799755611,
          "clean_vs_historical": "CONSISTENT",
          "param_hash": "314bb1fdf30eced58ab514dacd55bcb5f3fba522b716ed93640b3e2a11871dd9"
        },
        "scvi_matched": {
          "clean_macro_f1": 0.6402484093156036,
          "clean_vs_historical": "CONSISTENT",
          "param_hash": "63357b96f2873db63190b2bc9267c07951b3a95a89382c55cace89a1204d2e8d",
          "scvi_invariance_max_abs": 0.0
        }
      }
    },
    "integrity_ok": true,
    "scvi_invariance_max_abs": 0.0,
    "totalvi_f1_clean": 0.6157476799755611,
    "totalvi_f1_0.25": 0.6111421377055952
  },
  "direction_A": {
    "n_rows": 3,
    "meta": {
      "direction": "direction_A",
      "seed": 0,
      "source_batch": "PBMC10k",
      "target_batch": "PBMC5k",
      "n_eligible_classes": 22,
      "audit": {
        "source_n": 6855,
        "target_n": 3994,
        "n_genes": 15792,
        "n_proteins": 14
      },
      "models": {
        "totalvi": {
          "clean_macro_f1": 0.6357672182574826,
          "clean_vs_historical": "CONSISTENT",
          "param_hash": "4cd435ce6dd574ebd6a49f481433c3cf55931f344d303a447e1869f1f915505c"
        },
        "scvi_matched": {
          "clean_macro_f1": 0.6510217830627849,
          "clean_vs_historical": "CONSISTENT",
          "param_hash": "04174e0df6ab678447c38514c89a967e8ad3c89fa42937cc8a76f7ed80e98de1",
          "scvi_invariance_max_abs": 0.0
        }
      }
    },
    "integrity_ok": true,
    "scvi_invariance_max_abs": 0.0,
    "totalvi_f1_clean": 0.6357672182574826,
    "totalvi_f1_0.25": 0.6364947490663152
  }
}
```
