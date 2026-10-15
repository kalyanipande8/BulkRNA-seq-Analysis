"""End-to-end: simulate RNA-seq data, run differential expression,
co-expression module detection, pathway enrichment, and therapeutic
target prioritization."""
import json
import os

from src.simulate_rnaseq import simulate
from src.differential_expression import median_of_ratios_normalize, run_differential_expression, validate_de_recovery
from src.coexpression_wgcna import detect_modules, validate_module_recovery
from src.pathway_enrichment import pathway_enrichment, prioritize_therapeutic_targets

os.makedirs("data", exist_ok=True)
os.makedirs("results", exist_ok=True)

print("=== Simulating bulk RNA-seq data (150 samples, 80 genes, 3 co-expression modules) ===")
df = simulate()
df.to_csv("data/raw_counts.csv", index=False)
print(df.shape)

print("\n=== Normalization + differential expression ===")
normalized = median_of_ratios_normalize(df)
de = run_differential_expression(normalized)
de.to_csv("results/differential_expression.csv", index=False)
de_validation = validate_de_recovery(de)
print(de_validation)

print("\n=== WGCNA-style co-expression module detection ===")
modules = detect_modules(normalized)
module_validation = validate_module_recovery(modules)
print(module_validation)

print("\n=== Pathway enrichment (hypergeometric test) ===")
sig_genes = set(de.loc[de["significant"], "gene"])
enrichment = pathway_enrichment(sig_genes)
print(enrichment)
enrichment.to_csv("results/pathway_enrichment.csv", index=False)

print("\n=== Therapeutic target prioritization ===")
targets = prioritize_therapeutic_targets(de, modules)
print(targets[["gene", "log2fc", "p_value", "priority_score"]])
targets.to_csv("results/prioritized_targets.csv", index=False)

summary = {
    "n_samples": len(df),
    "de_recovery_rate_pct": de_validation["recovery_rate_pct"],
    "coexpression_module_agreement_pct": module_validation["agreement_pct"],
    "n_pathways_enriched": int(enrichment["enriched"].sum()),
    "top_therapeutic_target": targets.iloc[0]["gene"],
}
with open("results/summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("\n=== Summary ===")
print(json.dumps(summary, indent=2))
