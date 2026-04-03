import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

@pytest.fixture
def rbs_checker():
    checker = InternalRBSChecker()
    checker.initiate()
    return checker

def test_no_internal_rbs(rbs_checker):
    """Test sequences that are safe and lack an internal RBS."""
    
    # Normal CDS with no strong Shine-Dalgarno sequence
    cds_safe = "ATGCTAGCTAGCTAGTAA"
    passed, motif = rbs_checker.run(cds_safe)
    assert passed is True
    assert motif is None

    # Sequence has a Shine-Dalgarno motif, but NO start codon nearby
    cds_no_start = "AGGAGGCCTTCCTTCCTT"
    passed, motif = rbs_checker.run(cds_no_start)
    assert passed is True
    assert motif is None

    # Sequence has a start codon, but NO Shine-Dalgarno motif upstream
    cds_no_sd = "CCTTCCTTCCTTATGCCC"
    passed, motif = rbs_checker.run(cds_no_sd)
    assert passed is True
    assert motif is None

def test_has_internal_rbs(rbs_checker):
    """Test sequences that mistakenly contain internal RBS components."""
    
    # Strong SD (AGGAGG) + 6bp spacer + Start (ATG)
    cds_strong_rbs = "CCTAGGAGGCCCCCCATGCCC"
    passed, motif = rbs_checker.run(cds_strong_rbs)
    assert passed is False
    assert motif == "AGGAGGCCCCCCATG"

    # Alt SD (GGAGG) + 12bp spacer + Alt Start (GTG)
    cds_alt_rbs = "CCGGAGGTTCCTTCCTTCCTGTGCC"
    passed, motif = rbs_checker.run(cds_alt_rbs)
    assert passed is False
    assert motif == "GGAGGTTCCTTCCTTCCTGTG"

    # Short SD (AGGA) + 4bp spacer + Alt Start (TTG)
    cds_short_rbs = "AAAGGACTACTTGC"
    passed, motif = rbs_checker.run(cds_short_rbs)
    assert passed is False
    assert motif == "AGGACTACTTG"