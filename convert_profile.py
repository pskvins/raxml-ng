#!/usr/bin/env python3
"""
Convert profile format from HHblits/HHsearch output to RAxML-NG profile format.

Input format (output_profile.tsv):
- Line: "Query profile of sequence {name}"
- Line: "     A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y"
- Following lines: space-separated probabilities (20 values per line)

Output format (RAxML-NG profile):
- Line: ">{sequence_name}"
- Following lines: "{consensus_aa}\t{prob_A}\t{prob_R}\t...\t{prob_V}" (tab-separated, 21 columns)
- Amino acid order: A R N D C Q E G H I L K M F P S T W Y V (libpll standard)

NOTE: RAxML-NG requires all sequences to have the same length (aligned).
If your input sequences have different lengths, you need to align them first.
"""

import sys
import argparse


INPUT_AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
OUTPUT_AA_ORDER = "ARNDCQEGHILKMFPSTWYV"

INPUT_TO_OUTPUT_MAP = {aa: OUTPUT_AA_ORDER.index(aa) for aa in INPUT_AA_ORDER}


def get_consensus(probs_dict):
    """Get the amino acid with the highest probability."""
    return max(probs_dict.keys(), key=lambda aa: probs_dict[aa])


def reorder_probs(input_probs):
    """
    Reorder probabilities from input order to output order.
    
    Args:
        input_probs: list of 20 probabilities in input order (ACDEFGHIKLMNPQRSTVWY)
    
    Returns:
        list of 20 probabilities in output order (ARNDCQEGHILKMFPSTWYV)
    """
    probs_dict = {aa: input_probs[i] for i, aa in enumerate(INPUT_AA_ORDER)}
    return [probs_dict[aa] for aa in OUTPUT_AA_ORDER]


def parse_input_file(filepath):
    """
    Parse the input profile file.
    
    Returns:
        list of (sequence_name, list of probability rows)
    """
    sequences = []
    current_name = None
    current_probs = []
    header_seen = False
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            
            if line.startswith("Query profile of sequence"):
                if current_name is not None:
                    sequences.append((current_name, current_probs))
                current_name = line.replace("Query profile of sequence", "").strip()
                current_probs = []
                header_seen = False
            elif line.strip().startswith("A") and "C" in line and "D" in line:
                header_seen = True
            elif header_seen and line.strip():
                parts = line.split()
                if len(parts) == 20:
                    try:
                        probs = [float(p) for p in parts]
                        current_probs.append(probs)
                    except ValueError:
                        pass
    
    if current_name is not None and current_probs:
        sequences.append((current_name, current_probs))
    
    return sequences


def check_alignment(sequences):
    """Check if all sequences have the same length."""
    if not sequences:
        return True, 0
    
    lengths = [(name, len(probs)) for name, probs in sequences]
    first_len = lengths[0][1]
    
    all_same = all(l == first_len for _, l in lengths)
    return all_same, lengths


def convert_profile(input_file, output_file, check_lengths=True):
    """
    Convert profile format from input to RAxML-NG format.
    """
    sequences = parse_input_file(input_file)
    
    if check_lengths:
        all_same, lengths = check_alignment(sequences)
        if not all_same:
            print("WARNING: Sequences have different lengths (not aligned)!")
            print("RAxML-NG requires all sequences to have the same length.")
            print("\nSequence lengths:")
            for name, length in lengths:
                print(f"  {name}: {length} positions")
            print("\nPlease align your sequences before using with RAxML-NG.")
            print("Use --no-check to skip this warning and convert anyway.\n")
    
    with open(output_file, 'w') as f:
        for seq_name, prob_rows in sequences:
            f.write(f">{seq_name}\n")
            
            for input_probs in prob_rows:
                output_probs = reorder_probs(input_probs)
                
                probs_dict = {aa: input_probs[i] for i, aa in enumerate(INPUT_AA_ORDER)}
                consensus = get_consensus(probs_dict)
                
                prob_str = "\t".join(f"{p:.6f}" for p in output_probs)
                f.write(f"{consensus}\t{prob_str}\n")
    
    print(f"Converted {len(sequences)} sequences from '{input_file}' to '{output_file}'")


def main():
    parser = argparse.ArgumentParser(
        description="Convert HHblits/HHsearch profile format to RAxML-NG profile format"
    )
    parser.add_argument("input", help="Input profile file (e.g., output_profile.tsv)")
    parser.add_argument("output", help="Output profile file for RAxML-NG")
    parser.add_argument("--no-check", action="store_true",
                        help="Skip alignment length check")
    parser.add_argument("--single", type=str, metavar="NAME",
                        help="Extract only a single sequence by name")
    
    args = parser.parse_args()
    
    if args.single:
        sequences = parse_input_file(args.input)
        matching = [(n, p) for n, p in sequences if n == args.single]
        if not matching:
            print(f"ERROR: Sequence '{args.single}' not found in input file.")
            print("Available sequences:")
            for name, _ in sequences:
                print(f"  {name}")
            sys.exit(1)
        
        with open(args.output, 'w') as f:
            for seq_name, prob_rows in matching:
                f.write(f">{seq_name}\n")
                for input_probs in prob_rows:
                    output_probs = reorder_probs(input_probs)
                    probs_dict = {aa: input_probs[i] for i, aa in enumerate(INPUT_AA_ORDER)}
                    consensus = get_consensus(probs_dict)
                    prob_str = "\t".join(f"{p:.6f}" for p in output_probs)
                    f.write(f"{consensus}\t{prob_str}\n")
        print(f"Extracted sequence '{args.single}' to '{args.output}'")
    else:
        convert_profile(args.input, args.output, check_lengths=not args.no_check)


if __name__ == "__main__":
    main()
