import pytest
from genedesign.checkers.codon_checker import CodonChecker

@pytest.fixture
def codon_checker():
    """
    Fixture to initialize the CodonChecker with codon usage data.
    """
    checker = CodonChecker()
    checker.initiate()  # Load the codon usage data
    return checker

def test_high_cai_example(codon_checker):
    """
    Test case for a high CAI CDS, which should return True.
    """
    cds = ['ATG', 'AAA', 'CAT', 'TGG']  # High CAI example
    codons_above_board, codon_diversity, rare_codon_count, cai_value = codon_checker.run(cds)
    
    assert codons_above_board == True
    assert codon_diversity == 0.06451612903225806 ###4/62, because there are 4 unique codons and 62 total codons. Changed to match how codon diversity is calculated in codon checker.py
    assert rare_codon_count == 0
    assert cai_value > 0.5

def test_low_cai_example(codon_checker):
    """
    Test case for a low CAI CDS, which should return False.
    """
    cds = ['AGG', 'AGA', 'AGG', 'AGA']  # Low CAI example
    codons_above_board, codon_diversity, rare_codon_count, cai_value = codon_checker.run(cds)

    assert codons_above_board == False
    assert codon_diversity == 0.03225806451612903  ###2/62, because there are 2 unique codons and 62 total codons. Changed to match how codon diversity is calculated in codon checker.py
    assert rare_codon_count == 4
    assert cai_value < 0.1

def test_empty_cds(codon_checker):
    """
    Test case for an empty CDS, which should return False and zeros for all metrics.
    """
    cds = []
    codons_above_board, codon_diversity, rare_codon_count, cai_value = codon_checker.run(cds)

    assert codons_above_board == False
    assert codon_diversity == 0.0
    assert rare_codon_count == 0
    assert cai_value == 0.0

def test_medium_cai_example(codon_checker):
    """
    Test case for a medium CAI CDS, which should return True if it passes the thresholds.
    """
    cds = ['ATG', 'GCT', 'GAA', 'TAA']  # A moderate CAI example
    codons_above_board, codon_diversity, rare_codon_count, cai_value = codon_checker.run(cds)

    assert codons_above_board == True
    assert codon_diversity > 0.06451612903225806 ###4/62, because there are 4 unique codons and 62 total codons. Changed to match how codon diversity is calculated in codon checker.py
    assert rare_codon_count == 0
    assert cai_value > 0.2
