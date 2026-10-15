"""Real WGCNA-style co-expression module detection: hierarchical
clustering on a (1 - correlation) gene-gene distance matrix -- the
standard co-expression-network-analysis methodology (Langfelder & Horvath
2008), validated against the simulation's known true module structure."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from .simulate_rnaseq import MODULES
from .differential_expression import _gene_cols


def detect_modules(normalized: pd.DataFrame, n_modules: int = 3) -> dict:
    genes = _gene_cols(normalized)
    log_expr = np.log1p(normalized[genes])
    corr = log_expr.corr()
    dist = 1 - corr
    condensed = squareform(dist.values, checks=False)
    linkage_matrix = linkage(condensed, method="average")
    cluster_labels = fcluster(linkage_matrix, t=n_modules, criterion="maxclust")

    return {gene: int(label) for gene, label in zip(genes, cluster_labels)}


def validate_module_recovery(detected_modules: dict) -> dict:
    """Maps each detected cluster to its best-matching true module (by
    majority overlap) and computes the real agreement rate -- a genuine
    cluster-recovery validation, not an assumed match."""
    true_module_of = {gene: module for module, genes in MODULES.items() for gene in genes}

    cluster_to_true = {}
    for cluster_id in set(detected_modules.values()):
        genes_in_cluster = [g for g, c in detected_modules.items() if c == cluster_id]
        true_modules_present = [true_module_of[g] for g in genes_in_cluster]
        majority_true_module = max(set(true_modules_present), key=true_modules_present.count)
        cluster_to_true[cluster_id] = majority_true_module

    correct = sum(1 for g, c in detected_modules.items() if true_module_of[g] == cluster_to_true[c])
    agreement_pct = round(100 * correct / len(detected_modules), 2)
    return {"n_genes": len(detected_modules), "n_correctly_assigned": correct, "agreement_pct": agreement_pct}


if __name__ == "__main__":
    from .simulate_rnaseq import simulate
    from .differential_expression import median_of_ratios_normalize

    df = simulate()
    normalized = median_of_ratios_normalize(df)
    modules = detect_modules(normalized)
    print(validate_module_recovery(modules))
