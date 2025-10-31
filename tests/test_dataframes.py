from __future__ import annotations

from pandas.testing import assert_frame_equal

from pygorilla import (
    GeneHit,
    GORecord,
    GorillaResult,
    RunMetadata,
    go_terms_to_denormalized_dataframe,
    go_terms_to_normalized_dataframe,
)


def _sample_records() -> list[GORecord]:
    return [
        GORecord(
            go_id="GO:0000001",
            description="mitochondrion inheritance",
            p_value=1e-5,
            fdr_q_value=2e-4,
            enrichment=5.0,
            counts=(1000, 50, 25, 12),
            genes=[
                GeneHit(symbol="GENE1", description="first gene"),
                GeneHit(symbol="GENE2", description="second gene"),
            ],
        ),
        GORecord(
            go_id="GO:0000002",
            description="mitochondrial genome maintenance",
            p_value=None,
            fdr_q_value=None,
            enrichment=None,
            counts=(950, 40, 20, 0),
            genes=[],
        ),
    ]


def test_go_terms_to_normalized_dataframe() -> None:
    records = _sample_records()
    dataframe = go_terms_to_normalized_dataframe(records)

    assert dataframe.columns.tolist() == [
        "go_id",
        "description",
        "p_value",
        "fdr_q_value",
        "enrichment",
        "population_total",
        "go_term_total",
        "target_total",
        "overlap_total",
        "gene_count",
        "gene_symbols",
        "gene_descriptions",
    ]
    assert len(dataframe) == 2
    row = dataframe.iloc[0]
    assert row["go_id"] == "GO:0000001"
    assert row["gene_count"] == 2
    assert row["gene_symbols"] == ["GENE1", "GENE2"]
    assert row["gene_descriptions"] == ["first gene", "second gene"]
    empty_row = dataframe.iloc[1]
    assert empty_row["go_id"] == "GO:0000002"
    assert empty_row["gene_count"] == 0


def test_go_terms_to_denormalized_dataframe() -> None:
    records = _sample_records()
    dataframe = go_terms_to_denormalized_dataframe(records)

    assert dataframe.columns.tolist() == [
        "go_id",
        "description",
        "gene_symbol",
        "gene_description",
        "p_value",
        "fdr_q_value",
        "enrichment",
        "population_total",
        "go_term_total",
        "target_total",
        "overlap_total",
        "gene_index",
    ]
    assert len(dataframe) == 2
    first, second = dataframe.iloc[0], dataframe.iloc[1]
    assert first["gene_symbol"] == "GENE1"
    assert second["gene_symbol"] == "GENE2"
    assert first["go_id"] == second["go_id"] == "GO:0000001"
    assert first["gene_index"] == 0
    assert second["gene_index"] == 1


def test_gorilla_result_dataframe_helpers_match_free_functions() -> None:
    records = _sample_records()
    result = GorillaResult(
        run_id="run-123",
        result_url="http://example.com",
        go_terms=records,
        metadata=RunMetadata(),
        raw_html="<html></html>",
    )

    expected_normalized = go_terms_to_normalized_dataframe(records)
    expected_denormalized = go_terms_to_denormalized_dataframe(records)

    assert_frame_equal(result.to_go_terms_dataframe(), expected_normalized)
    assert_frame_equal(result.to_go_terms_dataframe(normalized=False), expected_denormalized)
