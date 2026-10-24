from pathlib import Path

import pytest

ROOT_DIR_PATH = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def df():
    """Reads aligned to a target region of the example bam file (examples/inputs/test.bam)."""
    from dseq.tools import to_df

    return to_df(
        bam_path=str(ROOT_DIR_PATH / "examples" / "inputs" / "test.bam"),
        contig="14",
        start=74846395,
        end=74846545,
    )
