import itertools
import re
from collections import Counter
from genedesign.models.transcript import Transcript
from genedesign.rbs_chooser import RBSChooser
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.seq_utils.hairpin_counter import hairpin_counter
from genedesign.seq_utils.reverse_complement import reverse_complement
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker  # ADDED

class TranscriptDesigner:
    """
    Constructs a DNA CDS from a protein sequence using a sliding window algorithm.
    Features robust lookahead evaluation and targeted regional error correction 
    to resolve structural and functional constraints.
    """

    # Standard genetic code mapping
    AMINO_ACID_MAP = {
        'A': ['GCG', 'GCC', 'GCA', 'GCT'], 'C': ['TGC', 'TGT'], 'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATT', 'ATA'], 'K': ['AAA', 'AAG'],
        'L': ['CTG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA'], 'M': ['ATG'], 'N': ['AAC', 'AAT'],
        'P': ['CCG', 'CCT', 'CCC', 'CCA'], 'Q': ['CAG', 'CAA'], 
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'T': ['ACC', 'ACT', 'ACA', 'ACG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 
        'W': ['TGG'], 'Y': ['TAC', 'TAT'], 'STOP': ['TAA', 'TAG', 'TGA']
    }

    def __init__(self):
        # Tools
        self.rbs_selector = None
        self.banned_seq_scanner = None
        self.promoter_scanner = None
        self.usage_analyzer = None
        self.rbs_scanner = None  # ADDED
        
        # Data dictionaries
        self.synonymous_codons = {}
        self.cai_scores = {}
        self.primary_codon = {}
        self.restricted_motifs = [] 

    def initiate(self) -> None:
        """Initializes external checkers and precomputes codon statistics."""
        self.rbs_selector = RBSChooser()
        self.rbs_selector.initiate()

        self.banned_seq_scanner = ForbiddenSequenceChecker()
        self.banned_seq_scanner.initiate()

        self.promoter_scanner = PromoterChecker()
        self.promoter_scanner.initiate()

        self.usage_analyzer = CodonChecker()
        self.usage_analyzer.initiate()
        
        self.rbs_scanner = InternalRBSChecker()  # ADDED
        self.rbs_scanner.initiate()              # ADDED

        # Build comprehensive list of restricted sequences including reverse complements
        motif_set = set(self.banned_seq_scanner.forbidden)
        for motif in self.banned_seq_scanner.forbidden:
            motif_set.add(reverse_complement(motif))
        self.restricted_motifs = list(motif_set)

        # Precompute optimal CAI mappings
        for aa, codon_list in self.AMINO_ACID_MAP.items():
            max_usage = max([self.usage_analyzer.codon_frequencies.get(c, 0.001) for c in codon_list])
            approved_codons = []
            
            for c in codon_list:
                usage = self.usage_analyzer.codon_frequencies.get(c, 0.001)
                calculated_cai = usage / max_usage
                self.cai_scores[c] = calculated_cai
                
                if c not in self.usage_analyzer.rare_codons:
                    approved_codons.append((c, calculated_cai))
            
            approved_codons.sort(key=lambda item: item[1], reverse=True)
            self.synonymous_codons[aa] = [c[0] for c in approved_codons]
            
            if not self.synonymous_codons[aa]:
                self.synonymous_codons[aa].append(codon_list[0])
                
            self.primary_codon[aa] = self.synonymous_codons[aa][0]

    def run(self, protein_seq: str, rbs_exclusions: set) -> Transcript:
        """Main orchestrator for generating the optimized transcript."""
        # Ensure proper termination by treating the sequence as a list of tokens
        clean_protein = protein_seq.strip().replace('\n', '')
        if clean_protein.endswith('*'):
            clean_protein = clean_protein[:-1]
        elif clean_protein.endswith('STOP'):
            clean_protein = clean_protein[:-4]
            
        protein_tokens = list(clean_protein)
        protein_tokens.append('STOP')
        
        rbs_exclusions = rbs_exclusions or set()
        chosen_rbs = self.rbs_selector.run("ATG", rbs_exclusions) 
        full_5_prime_utr = chosen_rbs.utr.upper()

        constructed_cds = []
        global_codon_counts = Counter()

        cache_promoter = {}
        cache_hairpin = {}

        # Process tokens in chunks of 3
        for i in range(0, len(protein_tokens), 3):
            active_aa_chunk = protein_tokens[i : i + 3]
            lookahead_aa_chunk = protein_tokens[i + 3 : i + 20] 

            options_per_aa = [self.synonymous_codons.get(token, [self.primary_codon.get(token, "NNN")]) for token in active_aa_chunk]
            search_space = list(itertools.product(*options_per_aa))

            # Establish local context for boundaries
            cds_upstream_chunk = "".join(constructed_cds)[-50:]
            context_upstream = (full_5_prime_utr + cds_upstream_chunk)[-50:] 
            context_downstream = "".join([self.primary_codon.get(token, "NNN") for token in lookahead_aa_chunk])
            
            # Baseline metrics
            baseline_banned = sum(context_upstream.count(m) for m in self.restricted_motifs)

            if context_upstream not in cache_hairpin:
                cache_hairpin[context_upstream] = hairpin_counter(context_upstream, 3, 4, 9)[0]
            baseline_hp_up = cache_hairpin[context_upstream]

            if context_downstream not in cache_hairpin:
                cache_hairpin[context_downstream] = hairpin_counter(context_downstream, 3, 4, 9)[0]
            baseline_hp_down = cache_hairpin[context_downstream]
            
            # ADDED: Baseline metrics for Internal RBS
            baseline_rbs_up = len(self.rbs_scanner.pattern.findall(cds_upstream_chunk))
            baseline_rbs_down = len(self.rbs_scanner.pattern.findall(context_downstream))

            optimal_weight = -float('inf')
            optimal_selection = search_space[0]

            # Evaluate search space
            for candidate_codons in search_space:
                fitness = self._calculate_chunk_fitness(
                    candidate_codons, context_upstream, cds_upstream_chunk, context_downstream, 
                    baseline_banned, baseline_hp_up, baseline_hp_down, baseline_rbs_up, baseline_rbs_down,
                    cache_promoter, cache_hairpin, global_codon_counts
                )

                if fitness > optimal_weight:
                    optimal_weight = fitness
                    optimal_selection = candidate_codons

            constructed_cds.extend(optimal_selection)
            for c in optimal_selection: 
                global_codon_counts[c] += 1

        # Final targeted cleanup
        constructed_cds = self._iterative_error_correction(constructed_cds, protein_tokens, chosen_rbs)

        # Return the clean protein string without the artificial token at the end
        return Transcript(chosen_rbs, clean_protein, constructed_cds)

    def _calculate_chunk_fitness(
        self, candidate_codons, context_upstream, cds_upstream_chunk, context_downstream, 
        baseline_banned, baseline_hp_up, baseline_hp_down, baseline_rbs_up, baseline_rbs_down,
        cache_promoter, cache_hairpin, global_codon_counts
    ) -> float:
        """Evaluates a sliding window candidate against all biological constraints."""
        window_seq = "".join(candidate_codons)
        eval_region = context_upstream + window_seq + context_downstream
        eval_region_upper = eval_region.upper()
        fitness_score = 0
        
        # Penalty 1: Restricted Motifs
        current_banned = sum(eval_region_upper.count(m) for m in self.restricted_motifs)
        added_banned = current_banned - baseline_banned
        if added_banned > 0:
            fitness_score -= added_banned * 10000000

        # Penalty 2: Internal Promoters
        promoter_eval_string = context_upstream[-28:] + window_seq + context_downstream[:28]
        if promoter_eval_string not in cache_promoter:
            found_promoter = not self.promoter_scanner.run(promoter_eval_string)[0]
            cache_promoter[promoter_eval_string] = 1 if found_promoter else 0
            
        if cache_promoter[promoter_eval_string] > 0:
            fitness_score -= 500000

        # Penalty 3: Secondary Structure (Hairpins)
        if eval_region_upper not in cache_hairpin:
            cache_hairpin[eval_region_upper] = hairpin_counter(eval_region_upper, 3, 4, 9)[0]
        total_hairpins = cache_hairpin[eval_region_upper]
        
        added_hairpins = total_hairpins - baseline_hp_up - baseline_hp_down
        if added_hairpins > 0:
            fitness_score -= added_hairpins * 500000
            
        local_junction = (context_upstream + window_seq)[-50:]
        if local_junction not in cache_hairpin:
            cache_hairpin[local_junction] = hairpin_counter(local_junction, 3, 4, 9)[0]
        junction_hairpins = cache_hairpin[local_junction]
        
        if junction_hairpins > 1:
            fitness_score -= 2000000
        elif junction_hairpins == 1:
            fitness_score -= 20000

        # ADDED: Penalty 4: Internal RBS check (Scanning CDS only)
        test_cds_local = cds_upstream_chunk + window_seq + context_downstream
        test_rbs_count = len(self.rbs_scanner.pattern.findall(test_cds_local))
        added_rbs = test_rbs_count - baseline_rbs_up - baseline_rbs_down
        if added_rbs > 0:
            fitness_score -= added_rbs * 5000000

        # Reward: CAI & Diversity Mapping
        cumulative_cai = sum(self.cai_scores.get(c, 0.01) for c in candidate_codons)
        diversity_tax = sum((global_codon_counts[c] ** 3) for c in candidate_codons) 
        fitness_score += (cumulative_cai * 10) - (diversity_tax * 5)

        return fitness_score

    def _iterative_error_correction(self, current_cds, protein_tokens, rbs):
        """Scans the completed sequence and surgically patches specific violations."""
        max_iterations = 12 
        for _ in range(max_iterations):
            full_transcript = rbs.utr.upper() + "".join(current_cds)
            is_clean = True

            has_banned, banned_fragment = self.banned_seq_scanner.run(full_transcript)
            if has_banned and banned_fragment:
                current_cds = self._patch_local_region(banned_fragment, current_cds, protein_tokens, rbs, "forbidden")
                is_clean = False
                continue

            passed_promoter, promoter_fragment = self.promoter_scanner.run(full_transcript)
            if not passed_promoter and promoter_fragment:
                current_cds = self._patch_local_region(promoter_fragment, current_cds, protein_tokens, rbs, "promoter")
                is_clean = False
                continue

            passed_hp, hp_fragment = hairpin_checker(full_transcript)
            if not passed_hp and hp_fragment:
                current_cds = self._patch_local_region(hp_fragment, current_cds, protein_tokens, rbs, "hairpin")
                is_clean = False
                continue

            # ADDED: RBS Sweep strictly scanning the CDS
            cds_only = "".join(current_cds)
            passed_rbs, rbs_fragment = self.rbs_scanner.run(cds_only)
            if not passed_rbs and rbs_fragment:
                current_cds = self._patch_local_region(rbs_fragment, current_cds, protein_tokens, rbs, "rbs")
                is_clean = False
                continue

            if is_clean:
                break
                
        return current_cds

    def _patch_local_region(self, violation_data, current_cds, protein_tokens, rbs, category):
        """Identifies the coordinates of a violation and reroutes the sequence locally."""
        if category == "hairpin":
            match_data = re.search(r'([ACGTacgt]+)\(([ACGTacgt]+)\)([ACGTacgt]+)', violation_data)
            if match_data:
                pure_violation_str = (match_data.group(1) + match_data.group(2) + match_data.group(3)).upper()
            else:
                pure_violation_str = re.sub(r'[^ACGTacgt]', '', violation_data).upper()
        else:
            pure_violation_str = re.sub(r'[^ACGTacgt]', '', violation_data).upper()
            
        # ADDED: Safe coordinate branching for RBS
        if category == "rbs":
            search_seq = "".join(current_cds)
            anchor_idx = search_seq.find(pure_violation_str)
            if anchor_idx == -1: 
                anchor_idx = search_seq.find(reverse_complement(pure_violation_str))
            if anchor_idx == -1: 
                return current_cds 
            target_nucleotide = anchor_idx + 3 
            radius = 3 
        else:
            search_seq = rbs.utr.upper() + "".join(current_cds)
            anchor_idx = search_seq.find(pure_violation_str)
            if anchor_idx == -1: 
                anchor_idx = search_seq.find(reverse_complement(pure_violation_str))
            if anchor_idx == -1: 
                return current_cds 
            utr_offset = len(rbs.utr)
            
            # Configure search radius based on violation type
            if category == "promoter":
                target_nucleotide = anchor_idx - utr_offset + 22 
                radius = 2 
            elif category == "hairpin":
                target_nucleotide = anchor_idx - utr_offset + (len(pure_violation_str) // 2)
                radius = 5 
            else:
                target_nucleotide = anchor_idx - utr_offset + (len(pure_violation_str) // 2)
                radius = 2

        center_codon = max(0, target_nucleotide // 3)
        idx_start = max(0, center_codon - radius)
        idx_end = min(len(current_cds), center_codon + radius + 1)
        
        # target_aas handles the token list smoothly regardless of character length
        target_aas = protein_tokens[idx_start:idx_end]
        
        # Contextual Boundary Resolution
        compiled_cds_upstream = "".join(current_cds[:idx_start])
        cds_upstream_chunk = compiled_cds_upstream[-50:] # Added for RBS boundary checks
        physical_upstream = (rbs.utr.upper() + compiled_cds_upstream)[-50:]
        physical_downstream = "".join(current_cds[idx_end : idx_end + 16])
        
        # Attempt 1: Standard restricted search (High CAI bias)
        standard_options = [self.synonymous_codons.get(token, [self.primary_codon.get(token, "NNN")])[:2] for token in target_aas]
        standard_combos = list(itertools.product(*standard_options))
        top_fitness, best_solution = self._evaluate_repair_candidates(
            standard_combos, physical_upstream, cds_upstream_chunk, physical_downstream, pure_violation_str
        )
        
        # Attempt 2: Expanded Search Space
        if top_fitness < -1000000:
            expanded_options = []
            for i, token in enumerate(target_aas):
                # Provide max flexibility at the violation epicenter
                if len(target_aas)//2 - 2 <= i <= len(target_aas)//2 + 2:
                    expanded_options.append(self.synonymous_codons.get(token, [self.primary_codon.get(token, "NNN")])[:4])
                else:
                    expanded_options.append(self.synonymous_codons.get(token, [self.primary_codon.get(token, "NNN")])[:1])
            
            expanded_combos = list(itertools.product(*expanded_options))
            fallback_fitness, fallback_solution = self._evaluate_repair_candidates(
                expanded_combos, physical_upstream, cds_upstream_chunk, physical_downstream, pure_violation_str
            )
            
            if fallback_fitness > top_fitness:
                best_solution = fallback_solution
                
        current_cds[idx_start:idx_end] = list(best_solution)
        return current_cds

    def _evaluate_repair_candidates(self, permutations, physical_upstream, cds_upstream_chunk, physical_downstream, pure_violation_str):
        """Scores potential patch sequences for the targeted region."""
        highest_score = -float('inf')
        winning_sequence = permutations[0] if permutations else []
        
        for candidate in permutations:
            simulated_span = physical_upstream + "".join(candidate) + physical_downstream
            simulated_span_upper = simulated_span.upper()
            temp_score = 0
            
            if pure_violation_str in simulated_span_upper or reverse_complement(pure_violation_str) in simulated_span_upper:
                temp_score -= 10000000
            
            banned_hits = sum(simulated_span_upper.count(m) for m in self.restricted_motifs)
            temp_score -= banned_hits * 100000
            
            promoter_hits = 1 if not self.promoter_scanner.run(simulated_span)[0] else 0
            temp_score -= promoter_hits * 50000
            
            passed_hp, _ = hairpin_checker(simulated_span_upper)
            if not passed_hp:
                temp_score -= 2000000
            else:
                hp_count, _ = hairpin_counter(simulated_span, 3, 4, 9)
                temp_score -= hp_count * 20000
                
            # ADDED: RBS validation inside the patch sequence scoring
            test_cds_local = cds_upstream_chunk + "".join(candidate) + physical_downstream
            if not self.rbs_scanner.run(test_cds_local)[0]:
                temp_score -= 20000000
            
            cai_bonus = sum(self.cai_scores.get(c, 0.01) for c in candidate)
            temp_score += cai_bonus * 10
            
            if temp_score > highest_score:
                highest_score = temp_score
                winning_sequence = candidate
                
        return highest_score, winning_sequence
               
if __name__ == "__main__":
    test_peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    exclusions = set()
    final_transcript = designer.run(test_peptide, exclusions)
    
    print(final_transcript)