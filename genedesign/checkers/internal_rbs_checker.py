import re

class InternalRBSChecker:
    """
    Description:
    Identifies potential internal ribosome binding sites within a sequence.
    Searches for a Shine-Dalgarno motif (AGGAGG, GGAGG, AGGA, GAGG) followed by a 
    spacer of 4-14 nucleotides, and terminating in a start codon (ATG, GTG, TTG).
    """

    def __init__(self):
        self.pattern = None

    def initiate(self) -> None:
        """
        Compiles the regex pattern for finding internal RBS sites.
        This allows the TranscriptDesigner to use .pattern.findall() for quick counting.
        """
        # Core Shine-Dalgarno motifs: AGGAGG, GGAGG, AGGA, GAGG
        # Spacer: 4 to 14 bases [ACGTacgt]{4,14}
        # Start Codons: ATG, GTG, TTG
        regex_string = r'((?:AGGAGG|GGAGG|AGGA|GAGG)[ACGTacgt]{4,14}(?:ATG|GTG|TTG))'
        self.pattern = re.compile(regex_string, re.IGNORECASE)

    def run(self, sequence: str) -> tuple[bool, str | None]:
        """
        Scans the sequence for internal ribosome binding sites.
        
        Returns:
            (True, None) if the sequence is safe (no internal RBS found).
            (False, problematic_sequence) if a potential internal RBS is detected.
        """
        match = self.pattern.search(sequence)
        
        if match:
            # Match found, return False and the exact sequence that triggered it
            return False, match.group(1).upper()
            
        return True, None

if __name__ == "__main__":
    checker = InternalRBSChecker()
    checker.initiate()
    
    # Example 1: Safe Sequence
    safe_seq = "ATGCCCCTTTTTAAA"
    print(f"Safe sequence check: {checker.run(safe_seq)}")
    
    # Example 2: Unintended Internal RBS (AGGAGG + 6bp spacer + ATG)
    dangerous_seq = "CCAGGAGGCCCCCCATGCCT"
    print(f"Dangerous sequence check: {checker.run(dangerous_seq)}")