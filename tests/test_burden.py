#!/usr/bin/env python3
"""pytest suite for calculate_fp_burden.py v5.3.7"""

import sys
import os
import pandas as pd
import numpy as np
import pytest

# Add scripts/ to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from calculate_fp_burden import (
    load_happy, load_manifest, compute_burden, aggregate_summary,
    validate_file_exists, format_median_range, to_csv_table,
    FileAccessError, DataValidationError, ManifestError, AggregationError
)


class TestFileValidation:
    """Tests for input file validation."""

    def test_missing_file(self):
        with pytest.raises(FileAccessError):
            validate_file_exists("/nonexistent/file.csv", "test")

    def test_valid_file(self, tmp_path):
        f = tmp_path / "test.csv"
        f.write_text("a,b\n1,2\n")
        p = validate_file_exists(str(f), "test")
        assert p == f

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.csv"
        f.write_text("")
        with pytest.raises(FileAccessError):
            validate_file_exists(str(f), "test")


class TestHappyLoading:
    """Tests for hap.py summary file loading."""

    def test_load_valid_happy(self):
        df = load_happy("tests/fixtures/minimal_baseline.csv")
        assert len(df) == 2
        assert "METRIC.Recall" in df.columns
        assert df["METRIC.Recall"].iloc[0] == 0.985

    def test_load_missing_column(self, tmp_path):
        f = tmp_path / "bad.csv"
        f.write_text("Type,Subset\nSNV,AllAutosomes\n")
        with pytest.raises(DataValidationError):
            load_happy(str(f))

    def test_load_negative_fp(self, tmp_path):
        f = tmp_path / "bad_fp.csv"
        f.write_text(
            "Type,Subtype,Subset,Filter,Genotype,"
            "METRIC.Recall,METRIC.Precision,QUERY.FP,QUERY.TP,TRUTH.TOTAL,Subset.IS_CONF.Size\n"
            "SNV,*,AllAutosomes,ALL,*,0.985,0.995,-1,395000,400000,2500000000\n"
        )
        with pytest.raises(DataValidationError):
            load_happy(str(f))


class TestManifestLoading:
    """Tests for manifest file loading."""

    def test_load_valid_manifest(self, tmp_path):
        f = tmp_path / "manifest.csv"
        f.write_text(
            "sample,caller,condition,baseline_path,consensus_path\n"
            "HG001,bcftools,test,tests/fixtures/minimal_baseline.csv,tests/fixtures/minimal_consensus.csv\n"
        )
        df = load_manifest(str(f))
        assert len(df) == 1
        assert df["sample"].iloc[0] == "HG001"

    def test_load_missing_column(self, tmp_path):
        f = tmp_path / "bad_manifest.csv"
        f.write_text("sample,caller\nHG001,bcftools\n")
        with pytest.raises(ManifestError):
            load_manifest(str(f))


class TestComputeBurden:
    """Tests for core burden computation."""

    def test_recall_gain_pp(self):
        df_base = load_happy("tests/fixtures/minimal_baseline.csv")
        df_cons = load_happy("tests/fixtures/minimal_consensus.csv")
        result, _ = compute_burden(df_base, df_cons, "HG001", "bcftools", "test")

        # SNV recall: 0.988 - 0.985 = 0.003 -> 0.3 pp
        snv_recall = result[result["Type"] == "SNV"]["Recall_gain_pp"].iloc[0]
        assert abs(snv_recall - 0.3) < 1e-10

        # INDEL recall: 0.730 - 0.680 = 0.050 -> 5.0 pp
        indel_recall = result[result["Type"] == "INDEL"]["Recall_gain_pp"].iloc[0]
        assert abs(indel_recall - 5.0) < 1e-10

    def test_precision_gain_pp(self):
        df_base = load_happy("tests/fixtures/minimal_baseline.csv")
        df_cons = load_happy("tests/fixtures/minimal_consensus.csv")
        result, _ = compute_burden(df_base, df_cons, "HG001", "bcftools", "test")

        # SNV precision: 0.993 - 0.995 = -0.002 -> -0.2 pp
        snv_prec = result[result["Type"] == "SNV"]["Precision_gain_pp"].iloc[0]
        assert abs(snv_prec - (-0.2)) < 1e-10

    def test_additional_tp_fp(self):
        df_base = load_happy("tests/fixtures/minimal_baseline.csv")
        df_cons = load_happy("tests/fixtures/minimal_consensus.csv")
        result, _ = compute_burden(df_base, df_cons, "HG001", "bcftools", "test")

        snv = result[result["Type"] == "SNV"].iloc[0]
        assert snv["Additional_TP"] == 200   # 395200 - 395000
        assert snv["Additional_FP"] == 1500  # 3500 - 2000
        assert snv["Net_benefit"] == -1300   # 200 - 1500

    def test_pp_drift_check(self):
        df_base = load_happy("tests/fixtures/minimal_baseline.csv")
        df_cons = load_happy("tests/fixtures/minimal_consensus.csv")
        result, _ = compute_burden(df_base, df_cons, "HG001", "bcftools", "test")

        for _, row in result.iterrows():
            assert abs(row["Recall_gain_pp"] - row["Recall_gain_pp_check"]) < 1e-10
            assert abs(row["Precision_gain_pp"] - row["Precision_gain_pp_check"]) < 1e-10

    def test_tp_fp_ratio_edge_cases(self):
        """Test Laplace-adjusted ratio handling."""
        df_base = load_happy("tests/fixtures/minimal_baseline.csv")
        df_cons = load_happy("tests/fixtures/minimal_consensus.csv")
        result, _ = compute_burden(df_base, df_cons, "HG001", "bcftools", "test")

        # All rows should have defined adjusted ratios
        assert result["TP_FP_ratio_adj"].notna().all()
        # Net benefit should always be defined
        assert result["Net_benefit"].notna().all()

    def test_f1_computation(self):
        df_base = load_happy("tests/fixtures/minimal_baseline.csv")
        df_cons = load_happy("tests/fixtures/minimal_consensus.csv")
        result, _ = compute_burden(df_base, df_cons, "HG001", "bcftools", "test")

        # F1 should be between 0 and 1 (or NaN if undefined)
        assert (result["F1_base"] >= 0).all() or result["F1_base"].isna().any()
        assert (result["F1_cons"] >= 0).all() or result["F1_cons"].isna().any()


class TestAggregation:
    """Tests for cross-sample aggregation."""

    def test_aggregate_two_samples(self):
        df_base = load_happy("tests/fixtures/minimal_baseline.csv")
        df_cons = load_happy("tests/fixtures/minimal_consensus.csv")

        raw1, _ = compute_burden(df_base, df_cons, "HG001", "bcftools", "test")
        raw2, _ = compute_burden(df_base, df_cons, "HG002", "bcftools", "test")

        df_raw = pd.concat([raw1, raw2], ignore_index=True)
        summary = aggregate_summary(df_raw)

        # Should have 2 rows (SNV and INDEL)
        assert len(summary) == 2

        # Median should equal the single value (both samples identical)
        snv = summary[summary["Type"] == "SNV"].iloc[0]
        assert snv["Recall_gain_pp_median"] == 0.3
        assert snv["Recall_gain_pp_min"] == 0.3
        assert snv["Recall_gain_pp_max"] == 0.3

    def test_empty_group_error(self):
        df = pd.DataFrame({"A": [], "Metric": []})
        with pytest.raises(AggregationError):
            aggregate_summary(df, group_cols=["A"])


class TestFormatting:
    """Tests for output formatting."""

    def test_format_median_range(self):
        assert format_median_range(4.75, 2.27, 7.11) == "4.75 [2.27-7.11]"
        assert format_median_range(100, 90, 110) == "100 [90-110]"
        assert format_median_range(0.5, 0.3, 0.7) == "0.50 [0.30-0.70]"

    def test_format_single_sample(self):
        # When min == max == median, show just the value
        assert format_median_range(3.14, 3.14, 3.14) == "3.14"

    def test_csv_table_columns(self):
        summary = pd.DataFrame({
            "Caller": ["bcftools"],
            "Recall_gain_pp_formatted": ["0.30 [0.25-0.35]"],
            "Precision_gain_pp_formatted": ["-0.20 [-0.25--0.15]"],
            "Additional_TP_formatted": ["200 [150-250]"],
            "Additional_FP_formatted": ["1500 [1200-1800]"],
            "TP_FP_ratio_formatted": ["0.13 [0.10-0.17]"],
            "Net_benefit_formatted": ["-1300 [-1600--1000]"],
            "Additional_TP_per_Mbp_formatted": ["0.08 [0.06-0.10]"],
            "Additional_FP_per_Mbp_formatted": ["0.60 [0.48-0.72]"]
        })
        csv_df = to_csv_table(summary)

        assert "Recall_Gain_pp" in csv_df.columns
        assert "Precision_Gain_pp" in csv_df.columns
        assert "Additional_TP" in csv_df.columns
        assert "Net_Benefit" in csv_df.columns
        assert "TP_FP_Ratio" in csv_df.columns
        assert "Recall_gain_pp_formatted" not in csv_df.columns


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
