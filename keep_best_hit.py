import argparse
import os
import sys


class BlastHit:
    def __init__(self, my_hit):
        self.query_id = my_hit[0]
        self.subject_id = my_hit[1]
        self.percent_identity = float(my_hit[2])
        self.alignment_length = int(my_hit[3])
        self.mismatches = int(my_hit[4])
        self.gap_opens = int(my_hit[5])
        self.q_start = int(my_hit[6])
        self.q_end = int(my_hit[7])
        self.s_start = int(my_hit[8])
        self.s_end = int(my_hit[9])
        self.evalue = float(my_hit[10])
        self.bit_score = float(my_hit[11])
        self.q_len = int(my_hit[12])
        self.s_len = int(my_hit[13])
        self.q_coverage = round((self.q_end + 1 - self.q_start) / self.q_len * 100, 2)
        self.s_coverage = round((self.s_end + 1 - self.s_start) / self.s_len * 100, 2)

    def __str__(self):
        return f"{self.query_id}\t{self.subject_id}\t{self.percent_identity}\t{self.q_coverage}\t{self.s_coverage}\n"

    def get_min_coverage(self):
        return min(self.q_coverage, self.s_coverage)

def main():
    parser = argparse.ArgumentParser(
        description="Find best reciprocal matches between two sets of sequences."
    )
    parser.add_argument(
        "-i", "--input", type=str, help="Path to input tabulated blastp file (-outfmt '7 std qlen slen')"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Path to the output file for best reciprocal matches.",
    )
    args = parser.parse_args()

    input_path = args.input
    output_path = args.output

    if input_path is None or output_path is None:
        parser.print_help()
        sys.exit("\nMissing required arguments: input and/or output file paths.")
    elif not os.path.isfile(input_path):
        raise FileNotFoundError(f"Input file {input_path} does not exist.")

    best_hit_q: dict = {}
    best_hit_s: dict = {}
    best_score_q: dict = {}
    best_score_s: dict = {}
    PERC_IDENTITY_CUTOFF: float = 35.0
    PERC_COVERAGE_CUTOFF: float = 50.0

    with open(input_path, "r") as infile:
        for line in infile:
            my_line = line.strip().split("\t")
            if my_line[0].startswith("#"):
                continue

            hit = BlastHit(my_line)
            score = hit.get_min_coverage() + hit.percent_identity

            if score > best_score_q.get(hit.query_id, 0) and \
            score > best_score_s.get(hit.subject_id, 0) and score >= 110:
                best_score_q[hit.query_id] = score
                best_score_s[hit.subject_id] = score
                
                # Delete old best hits if new best hit is found:
                # Check first if either the query or subject already has a best hit
                # If so, remove the previous best hit from the opposite dictionary
                if hit.query_id in best_hit_q:
                    prev_best_hit = best_hit_q[hit.query_id]
                    if prev_best_hit.subject_id in best_hit_s:
                        del best_hit_s[prev_best_hit.subject_id]
                if hit.subject_id in best_hit_s:
                    prev_best_hit = best_hit_s[hit.subject_id]
                    if prev_best_hit.query_id in best_hit_q:
                        del best_hit_q[prev_best_hit.query_id]
                
                # Add new best hit
                best_hit_q[hit.query_id] = hit
                best_hit_s[hit.subject_id] = hit
    

            
               
    # Filtering hits based on best reciprocal matches stats (PERC_IDENTITY + PERC_COVERAGE)
    with open(output_path, "w") as outfile:
        for query_id in best_hit_q.keys():
            best_hit = best_hit_q[query_id]
            if best_hit.percent_identity >= PERC_IDENTITY_CUTOFF and \
            best_hit.get_min_coverage() >= PERC_COVERAGE_CUTOFF:
                outfile.write(str(best_hit))


if __name__ == "__main__":
    main()
