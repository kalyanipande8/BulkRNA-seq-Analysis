import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

from src.simulate_rnaseq import simulate
from src.differential_expression import median_of_ratios_normalize, run_differential_expression, validate_de_recovery
from src.coexpression_wgcna import detect_modules, validate_module_recovery
from src.pathway_enrichment import pathway_enrichment, prioritize_therapeutic_targets


@pytest.fixture(scope="module")
def normalized():
    df = simulate()
    return median_of_ratios_normalize(df)


@pytest.fixture(scope="module")
def de_results(normalized):
    return run_differential_expression(normalized)


def test_de_recovers_most_true_genes(de_results):
    result = validate_de_recovery(de_results)
    assert result["recovery_rate_pct"] >= 70.0


def test_de_does_not_flag_every_gene(de_results):
    assert de_results["significant"].sum() < len(de_results)


def test_coexpression_modules_recovered_with_high_agreement(normalized):
    modules = detect_modules(normalized)
    result = validate_module_recovery(modules)
    assert result["agreement_pct"] >= 90.0


def test_pathway_enrichment_flags_only_the_true_module(de_results):
    sig_genes = set(de_results.loc[de_results["significant"], "gene"])
    enrichment = pathway_enrichment(sig_genes)
    enriched_modules = set(enrichment.loc[enrichment["enriched"], "module"])
    assert enriched_modules == {"Module_1"}


def test_prioritized_targets_are_all_significant(normalized, de_results):
    modules = detect_modules(normalized)
    targets = prioritize_therapeutic_targets(de_results, modules)
    assert targets["significant"].all()
