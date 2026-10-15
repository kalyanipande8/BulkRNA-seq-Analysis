"""Real differential expression: DESeq2-style median-of-ratios
normalization followed by a per-gene t-test with a combined p-value +
effect-size significance rule -- validated against known ground-truth DE
genes."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .simulate_rnaseq import TRUE_DE_GENES

GENE_COLS_PREFIX = "GENE_"


def _gene_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith(GENE_COLS_PREFIX)]


def median_of_ratios_normalize(df: pd.DataFrame) -> pd.DataFrame:
    genes = _gene_cols(df)
    counts = df[genes].replace(0, np.nan)
    log_geo_mean = np.log(counts).mean(axis=0)
    valid = log_geo_mean[np.isfinite(log_geo_mean)].index
    size_factors = np.exp((np.log(counts[valid]).sub(log_geo_mean[valid], axis=1)).median(axis=1))

    out = df.copy()
    out[genes] = df[genes].div(size_factors, axis=0)
    return out


def run_differential_expression(normalized: pd.DataFrame, alpha: float = 0.05, log2fc_threshold: float = 0.3) -> pd.DataFrame:
    genes = _gene_cols(normalized)
    diseased = normalized[normalized["condition"] == "disease"]
    control = normalized[normalized["condition"] == "control"]

    rows = []
    for gene in genes:
        log2fc = np.log2(diseased[gene].mean() + 1) - np.log2(control[gene].mean() + 1)
        t_stat, p_value = stats.ttest_ind(diseased[gene], control[gene])
        significant = bool(p_value < alpha and abs(log2fc) >= log2fc_threshold)
        rows.append({"gene": gene, "log2fc": round(float(log2fc), 4), "p_value": round(float(p_value), 6),
                      "significant": significant})

    return pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)


def validate_de_recovery(de_results: pd.DataFrame) -> dict:
    discovered = set(de_results.loc[de_results["significant"], "gene"])
    true_de = set(TRUE_DE_GENES)
    return {"n_discovered": len(discovered), "n_true_recovered": len(discovered & true_de),
            "n_true_total": len(true_de), "recovery_rate_pct": round(100 * len(discovered & true_de) / len(true_de), 2)}


if __name__ == "__main__":
    from .simulate_rnaseq import simulate

    df = simulate()
    normalized = median_of_ratios_normalize(df)
    de = run_differential_expression(normalized)
    print(de.head(10))
    print(validate_de_recovery(de))
