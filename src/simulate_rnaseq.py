"""Simulate bulk RNA-seq count data for 150 samples (disease vs. control)
across 80 genes, organized into 3 real co-expression modules (genes
within a module are correlated via a shared latent factor), with 15
genes in Module 1 genuinely differentially expressed in disease --
the ground truth for DE, WGCNA-style co-expression, and pathway
enrichment analysis to recover."""
from __future__ import annotations

import numpy as np
import pandas as pd

N_SAMPLES = 150
N_GENES = 80
SEED = 251

MODULES = {
    "Module_1": [f"GENE_{i:02d}" for i in range(0, 25)],
    "Module_2": [f"GENE_{i:02d}" for i in range(25, 50)],
    "Module_3": [f"GENE_{i:02d}" for i in range(50, 80)],
}
TRUE_DE_GENES = MODULES["Module_1"][:15]


def simulate(n: int = N_SAMPLES, seed: int = SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    condition = rng.choice(["disease", "control"], n, p=[0.5, 0.5])
    gene_names = [f"GENE_{i:02d}" for i in range(N_GENES)]

    # a per-sample latent factor per module drives genuine co-expression
    # (genes within a module move together across samples)
    module_factors = {m: rng.normal(0, 1, n) for m in MODULES}

    base_rate = rng.lognormal(3, 0.6, N_GENES)
    counts = np.zeros((n, N_GENES))
    for g_idx, gene in enumerate(gene_names):
        module = next(m for m, genes in MODULES.items() if gene in genes)
        factor = module_factors[module]
        de_effect = 1.4 if gene in TRUE_DE_GENES else 1.0
        de_effect_applied = np.where(condition == "disease", de_effect, 1.0)
        lam = base_rate[g_idx] * np.exp(0.3 * factor) * de_effect_applied
        counts[:, g_idx] = rng.poisson(lam.clip(1, None))

    df = pd.DataFrame(counts, columns=gene_names)
    df.insert(0, "sample_id", [f"S-{i:03d}" for i in range(n)])
    df.insert(1, "condition", condition)
    return df


if __name__ == "__main__":
    df = simulate()
    print(df.shape)
