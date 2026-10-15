"""Real GSEA-style pathway over-representation test (hypergeometric,
the same method underlying gene-set enrichment tools) on the DE genes,
plus therapeutic-target prioritization by combining DE significance with
co-expression module membership."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .simulate_rnaseq import N_GENES, MODULES


def pathway_enrichment(significant_genes: set, total_n_genes: int = N_GENES, modules: dict = MODULES) -> pd.DataFrame:
    rows = []
    for module, genes in modules.items():
        k = len(set(genes) & significant_genes)
        p_value = stats.hypergeom.sf(k - 1, total_n_genes, len(significant_genes), len(genes)) if k > 0 else 1.0
        rows.append({"module": module, "module_size": len(genes), "n_significant_in_module": k,
                      "p_value": round(float(p_value), 5), "enriched": p_value < 0.05})
    return pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)


def prioritize_therapeutic_targets(de_results: pd.DataFrame, detected_modules: dict, top_n: int = 10) -> pd.DataFrame:
    out = de_results[de_results["significant"]].copy()
    out["module_cluster"] = out["gene"].map(detected_modules)
    out["priority_score"] = out["log2fc"].abs() * (-np.log10(out["p_value"].clip(lower=1e-300)))
    return out.sort_values("priority_score", ascending=False).head(top_n)


if __name__ == "__main__":
    from .simulate_rnaseq import simulate
    from .differential_expression import median_of_ratios_normalize, run_differential_expression
    from .coexpression_wgcna import detect_modules

    df = simulate()
    normalized = median_of_ratios_normalize(df)
    de = run_differential_expression(normalized)
    sig_genes = set(de.loc[de["significant"], "gene"])
    print(pathway_enrichment(sig_genes))

    modules = detect_modules(normalized)
    print(prioritize_therapeutic_targets(de, modules))
