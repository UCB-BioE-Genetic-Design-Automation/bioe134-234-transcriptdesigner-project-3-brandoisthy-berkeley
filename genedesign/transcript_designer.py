import csv
from collections import defaultdict
from itertools import product

from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.models.transcript import Transcript
from genedesign.rbs_chooser import RBSChooser


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a CDS using a 9 + 9 + 18 sliding window:
    - 9 bp locked preamble
    - 9 bp current decision zone
    - 18 bp look-ahead viability search
    """

    PREAMBLE_CODONS = 3 #9 NT window
    CURRENT_CODONS = 3 #9 NT window
    LOOKAHEAD_CODONS = 6 #18 NT window
    SIGMA70_WINDOW_CODONS = 12 #36 NT window
    PROMOTER_WINDOW_BP = 29
    PROMOTER_CONTEXT_CODONS = 9
    STOP_CODON = "TAA"

    def __init__(self):
        self.amino_acid_to_codons = {}
        self.codon_checker = CodonChecker()
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.promoter_checker = PromoterChecker()
        self.rbsChooser = None
        self._design_cache = {}
        self._lookahead_cache = {}
        self._promoter_check_cache = {}

    def initiate(self) -> None:
        """
        Initializes the codon table, checkers, and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.codon_checker.initiate()
        self.forbidden_checker.initiate()
        self.promoter_checker.initiate()

        self.amino_acid_to_codons = self._load_synonymous_codons()

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Designs a CDS for the peptide using the 9 + 9 + 18 local search.
        """
        if self.rbsChooser is None:
            raise RuntimeError("TranscriptDesigner must be initiated before run().")

        self._design_cache = {}
        self._lookahead_cache = {}
        self._promoter_check_cache = {}

        designed_codons = self._design_from(peptide, 0, tuple(), frozenset())
        if designed_codons is None:
            raise ValueError("Unable to design a CDS that satisfies the sliding-window constraints.")

        codons = designed_codons + [self.STOP_CODON]
        cds = ''.join(codons)
        selected_rbs = self._select_rbs(cds, ignores)
        return Transcript(selected_rbs, peptide, codons)

    def _select_rbs(self, cds: str, ignores: set):
        candidate_ignores = set(ignores)
        fallback_rbs = None

        while True:
            try:
                candidate_rbs = self.rbsChooser.run(cds, candidate_ignores)
            except Exception:
                if fallback_rbs is not None:
                    return fallback_rbs
                raise

            if fallback_rbs is None:
                fallback_rbs = candidate_rbs

            transcript_dna = candidate_rbs.utr.upper() + cds
            forbidden_safe, _ = self.forbidden_checker.run(transcript_dna)
            if self._is_promoter_safe(transcript_dna) and forbidden_safe:
                return candidate_rbs

            candidate_ignores.add(candidate_rbs)

    def _load_synonymous_codons(self) -> dict[str, list[str]]:
        amino_acid_to_codons = defaultdict(list)
        codon_usage_file = "genedesign/data/codon_usage.txt"

        with open(codon_usage_file, "r") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if len(row) < 3:
                    continue
                codon = row[0].strip().upper()
                amino_acid = row[1].strip()
                if amino_acid == "*":
                    continue
                usage = float(row[2].strip())
                amino_acid_to_codons[amino_acid].append((codon, usage))

        ranked_codons = {}
        for amino_acid, codons in amino_acid_to_codons.items():
            ranked_codons[amino_acid] = [
                codon for codon, _ in sorted(codons, key=lambda item: item[1], reverse=True)
            ]
        return ranked_codons

    def _design_from(
        self,
        peptide: str,
        position: int,
        locked_codons: tuple[str, ...],
        used_codons: frozenset[str],
    ) -> list[str] | None:
        if position >= len(peptide):
            return []

        key = (position, locked_codons[-self.PROMOTER_CONTEXT_CODONS :], used_codons)
        if key in self._design_cache:
            cached = self._design_cache[key]
            return list(cached) if cached is not None else None

        current_aas = peptide[position : position + self.CURRENT_CODONS]
        preamble = list(locked_codons[-self.PREAMBLE_CODONS :])
        promoter_context = tuple(locked_codons[-self.PROMOTER_CONTEXT_CODONS :])

        for current_codons in self._enumerate_chunk_codons(current_aas, used_codons):
            if not self._passes_current_suite(preamble, current_codons):
                continue

            lookahead_start = position + len(current_aas)
            lookahead_aas = peptide[lookahead_start : lookahead_start + self.LOOKAHEAD_CODONS]
            if not self._passes_lookahead_suite(promoter_context, current_codons, lookahead_aas):
                continue

            next_locked_codons = locked_codons + tuple(current_codons)
            next_used_codons = used_codons | frozenset(current_codons)
            suffix = self._design_from(peptide, lookahead_start, next_locked_codons, next_used_codons)
            if suffix is not None:
                solution = tuple(current_codons + suffix)
                self._design_cache[key] = solution
                return list(solution)

        self._design_cache[key] = None
        return None

    def _enumerate_chunk_codons(self, amino_acids: str, used_codons: frozenset[str]) -> list[list[str]]:
        if not amino_acids:
            return [[]]

        options = [self.amino_acid_to_codons[aa] for aa in amino_acids]
        ranked_chunks = []
        for candidate in product(*options):
            candidate_set = set(candidate)
            new_unique_codons = len(candidate_set - used_codons)
            intra_chunk_diversity = len(candidate_set)
            frequency_score = sum(self.codon_checker.codon_frequencies.get(codon, 0.0) for codon in candidate)
            ranked_chunks.append(((new_unique_codons, intra_chunk_diversity, frequency_score), list(candidate)))
        ranked_chunks.sort(key=lambda item: item[0], reverse=True)
        return [candidate for _, candidate in ranked_chunks]

    def _passes_current_suite(self, preamble: list[str], current_codons: list[str]) -> bool:
        combined_codons = preamble + current_codons
        combined_dna = ''.join(combined_codons)

        sequence_allowed, _ = self.forbidden_checker.run(combined_dna)
        if not sequence_allowed:
            return False

        return all(codon not in self.codon_checker.rare_codons for codon in current_codons)

    def _passes_lookahead_suite(
        self,
        promoter_context: tuple[str, ...],
        current_codons: list[str],
        lookahead_aas: str,
    ) -> bool:
        fixed_prefix = promoter_context + tuple(current_codons)

        if not lookahead_aas:
            return self._passes_terminal_stop_suite(fixed_prefix)

        cache_key = (''.join(fixed_prefix), lookahead_aas)
        if cache_key not in self._lookahead_cache:
            self._lookahead_cache[cache_key] = self._search_lookahead_path(fixed_prefix, lookahead_aas, tuple())
        return self._lookahead_cache[cache_key]

    def _search_lookahead_path(
        self,
        fixed_prefix: tuple[str, ...],
        lookahead_aas: str,
        candidate_codons: tuple[str, ...],
    ) -> bool:
        if len(candidate_codons) == len(lookahead_aas):
            return self._evaluate_full_lookahead(fixed_prefix, candidate_codons)

        amino_acid = lookahead_aas[len(candidate_codons)]
        for codon in self.amino_acid_to_codons[amino_acid]:
            if codon in self.codon_checker.rare_codons:
                continue

            next_candidate = candidate_codons + (codon,)
            if not self._passes_partial_window_checks(fixed_prefix, next_candidate):
                continue
            if self._search_lookahead_path(fixed_prefix, lookahead_aas, next_candidate):
                return True

        return False

    def _passes_partial_window_checks(
        self,
        fixed_prefix: tuple[str, ...],
        candidate_codons: tuple[str, ...],
    ) -> bool:
        combined_dna = ''.join(fixed_prefix + candidate_codons)

        sequence_allowed, _ = self.forbidden_checker.run(combined_dna)
        if not sequence_allowed:
            return False

        return True

    def _evaluate_full_lookahead(
        self,
        fixed_prefix: tuple[str, ...],
        lookahead_codons: tuple[str, ...],
    ) -> bool:
        combined_dna = ''.join(fixed_prefix + lookahead_codons)
        sequence_allowed, _ = self.forbidden_checker.run(combined_dna)
        if not sequence_allowed:
            return False

        return (
            self._passes_sigma70_suite(fixed_prefix, lookahead_codons)
            and
            all(codon not in self.codon_checker.rare_codons for codon in lookahead_codons)
            and self._passes_lookahead_codon_usage(list(lookahead_codons))
        )

    def _passes_lookahead_codon_usage(self, lookahead_codons: list[str]) -> bool:
        if not lookahead_codons:
            return True

        _, _, rare_codon_count, _ = self.codon_checker.run(lookahead_codons)
        return rare_codon_count == 0

    def _passes_terminal_stop_suite(self, codons: tuple[str, ...]) -> bool:
        terminal_dna = ''.join(codons) + self.STOP_CODON
        sequence_allowed, _ = self.forbidden_checker.run(terminal_dna)
        if not sequence_allowed:
            return False

        return self._passes_promoter_check(codons + (self.STOP_CODON,))

    def _passes_sigma70_suite(
        self,
        fixed_prefix: tuple[str, ...],
        lookahead_codons: tuple[str, ...],
    ) -> bool:
        sigma70_dna = ''.join(fixed_prefix + lookahead_codons)
        current_start = max(0, len(fixed_prefix) - self.CURRENT_CODONS) * 3
        current_end = current_start + (self.CURRENT_CODONS * 3)

        if len(sigma70_dna) < self.PROMOTER_WINDOW_BP:
            return True

        if len(sigma70_dna) <= self.SIGMA70_WINDOW_CODONS * 3:
            return self._is_promoter_safe(sigma70_dna)

        max_window_start = len(sigma70_dna) - (self.SIGMA70_WINDOW_CODONS * 3)
        for start in range(max_window_start + 1):
            end = start + (self.SIGMA70_WINDOW_CODONS * 3)
            overlaps_current = start < current_end and end > current_start
            if not overlaps_current:
                continue

            if not self._is_promoter_safe(sigma70_dna[start:end]):
                return False

        return True

    def _is_promoter_safe(self, dna: str) -> bool:
        if len(dna) < self.PROMOTER_WINDOW_BP:
            return True

        cached = self._promoter_check_cache.get(dna)
        if cached is not None:
            return cached

        promoter_safe, _ = self.promoter_checker.run(dna)
        self._promoter_check_cache[dna] = promoter_safe
        return promoter_safe

    def _passes_promoter_check(self, codons: tuple[str, ...]) -> bool:
        dna = ''.join(codons)
        return self._is_promoter_safe(dna)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"

    designer = TranscriptDesigner()
    designer.initiate()

    transcript = designer.run(peptide, set())
    print(transcript)