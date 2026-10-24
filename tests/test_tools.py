import pandas as pd
import pysam
import pytest

import dseq.tools  # noqa: F401 (registers the .dseq accessor)
from dseq.tools import _to_df


def test_to_df(df):
    assert df.shape == (56424, 6)
    assert df["reference position"].dtype == "Int64"
    assert df["read position"].dtype == "Int64"
    assert df["reference position"].between(74846395, 74846545).all()


@pytest.mark.parametrize(
    "method,shape",
    [
        ("get_coverage", (151, 3)),
        ("get_counts_per_mutation", (5,)),
        ("get_counts_per_nt_mutation", (212, 4)),
        ("get_counts_per_codon_mutation", (330, 4)),
        ("get_counts_per_aa_mutation", (330, 6)),
        ("get_frequencies_per_nt_mutation", (212, 6)),
    ],
)
def test_accessor_shapes(df, method, shape):
    assert getattr(df.dseq, method)().shape == shape


def test_indels(df):
    assert df.query("`CIGAR`.str.contains('I') or `CIGAR`.str.contains('D')").shape == (2699, 6)


def test_to_df_soft_clipped():
    """Read bases are taken from the full query sequence, so soft-clipped bases do not shift them."""
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr", "LN": 10}]})
    read = pysam.AlignedSegment(header)
    read.query_name = "r1"
    read.query_sequence = "ggACGT"
    read.reference_start = 0
    read.cigarstring = "2S4M"
    read.set_tag("MD", "4")
    df1 = _to_df(read).dropna(subset=["reference position"])
    assert "".join(df1["read"].str.upper()) == "ACGT"


def test_get_map_aggfunc():
    df0 = pd.DataFrame(
        {
            "reference position": [1, 1, 2],
            "reference": ["A", "A", "C"],
            "read": ["G", "G", "T"],
            "depth": [2, 5, 3],
        }
    )
    assert df0.dseq.get_submap(value="depth", aggfunc="sum").data.loc["G", "A"] == 7
    assert df0.dseq.get_submap(value="depth", aggfunc="max").data.loc["G", "A"] == 5
