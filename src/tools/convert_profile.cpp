/*
 * Convert profile format from HHblits/HHsearch output to RAxML-NG profile format.
 *
 * Input format (output_profile.tsv):
 * - Line: "Query profile of sequence {name}"
 * - Line: "     A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y"
 * - Following lines: space-separated probabilities (20 values per line)
 *
 * Output format (RAxML-NG profile):
 * - Line: ">{sequence_name}"
 * - Following lines: "{consensus_aa}\t{prob_A}\t{prob_R}\t...\t{prob_V}" (tab-separated, 21 columns)
 * - Amino acid order: A R N D C Q E G H I L K M F P S T W Y V (libpll standard)
 *
 * NOTE: RAxML-NG requires all sequences to have the same length (aligned).
 * If your input sequences have different lengths, you need to align them first.
 *
 * Compile: g++ -O2 -o convert_profile convert_profile.cpp
 * Usage: ./convert_profile input.tsv output.txt [--no-check] [--single NAME]
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
#include <cstring>

const std::string INPUT_AA_ORDER = "ACDEFGHIKLMNPQRSTVWY";
const std::string OUTPUT_AA_ORDER = "ARNDCQEGHILKMFPSTWYV";
const size_t NUM_AA = 20;

struct Sequence {
    std::string name;
    std::vector<std::vector<double>> prob_rows;
};

std::map<char, size_t> build_aa_index_map(const std::string& order) {
    std::map<char, size_t> index_map;
    for (size_t i = 0; i < order.size(); ++i) {
        index_map[order[i]] = i;
    }
    return index_map;
}

std::vector<double> reorder_probs(const std::vector<double>& input_probs) {
    static const auto input_map = build_aa_index_map(INPUT_AA_ORDER);
    
    std::vector<double> output_probs(NUM_AA);
    for (size_t i = 0; i < NUM_AA; ++i) {
        char aa = OUTPUT_AA_ORDER[i];
        size_t input_idx = input_map.at(aa);
        output_probs[i] = input_probs[input_idx];
    }
    return output_probs;
}

char get_consensus(const std::vector<double>& input_probs) {
    size_t max_idx = 0;
    double max_prob = input_probs[0];
    for (size_t i = 1; i < input_probs.size(); ++i) {
        if (input_probs[i] > max_prob) {
            max_prob = input_probs[i];
            max_idx = i;
        }
    }
    return INPUT_AA_ORDER[max_idx];
}

std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

std::vector<std::string> split(const std::string& str) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<Sequence> parse_input_file(const std::string& filepath) {
    std::vector<Sequence> sequences;
    std::ifstream infile(filepath);
    
    if (!infile.is_open()) {
        std::cerr << "ERROR: Cannot open input file: " << filepath << std::endl;
        return sequences;
    }
    
    std::string line;
    Sequence current_seq;
    bool header_seen = false;
    bool has_current = false;
    
    const std::string prefix = "Query profile of sequence";
    
    while (std::getline(infile, line)) {
        if (line.find(prefix) == 0) {
            if (has_current && !current_seq.prob_rows.empty()) {
                sequences.push_back(std::move(current_seq));
            }
            current_seq = Sequence();
            current_seq.name = trim(line.substr(prefix.length()));
            header_seen = false;
            has_current = true;
        }
        else {
            std::string trimmed = trim(line);
            if (!trimmed.empty() && trimmed[0] == 'A' && 
                trimmed.find('C') != std::string::npos && 
                trimmed.find('D') != std::string::npos) {
                header_seen = true;
            }
            else if (header_seen && !trimmed.empty()) {
                std::vector<std::string> parts = split(trimmed);
                if (parts.size() == NUM_AA) {
                    std::vector<double> probs(NUM_AA);
                    bool valid = true;
                    for (size_t i = 0; i < NUM_AA; ++i) {
                        try {
                            probs[i] = std::stod(parts[i]);
                        } catch (...) {
                            valid = false;
                            break;
                        }
                    }
                    if (valid) {
                        current_seq.prob_rows.push_back(probs);
                    }
                }
            }
        }
    }
    
    if (has_current && !current_seq.prob_rows.empty()) {
        sequences.push_back(std::move(current_seq));
    }
    
    infile.close();
    return sequences;
}

bool check_alignment(const std::vector<Sequence>& sequences, bool print_warning) {
    if (sequences.empty()) return true;
    
    size_t first_len = sequences[0].prob_rows.size();
    bool all_same = true;
    
    for (const auto& seq : sequences) {
        if (seq.prob_rows.size() != first_len) {
            all_same = false;
            break;
        }
    }
    
    if (!all_same && print_warning) {
        std::cout << "WARNING: Sequences have different lengths (not aligned)!" << std::endl;
        std::cout << "RAxML-NG requires all sequences to have the same length." << std::endl;
        std::cout << "\nSequence lengths:" << std::endl;
        for (const auto& seq : sequences) {
            std::cout << "  " << seq.name << ": " << seq.prob_rows.size() << " positions" << std::endl;
        }
        std::cout << "\nPlease align your sequences before using with RAxML-NG." << std::endl;
        std::cout << "Use --no-check to skip this warning and convert anyway.\n" << std::endl;
    }
    
    return all_same;
}

void write_sequence(std::ofstream& outfile, const Sequence& seq) {
    outfile << ">" << seq.name << "\n";
    
    for (const auto& input_probs : seq.prob_rows) {
        std::vector<double> output_probs = reorder_probs(input_probs);
        char consensus = get_consensus(input_probs);
        
        outfile << consensus;
        outfile << std::fixed << std::setprecision(6);
        for (const auto& prob : output_probs) {
            outfile << "\t" << prob;
        }
        outfile << "\n";
    }
}

void convert_profile(const std::string& input_file, const std::string& output_file, bool check_lengths) {
    std::vector<Sequence> sequences = parse_input_file(input_file);
    
    if (sequences.empty()) {
        std::cerr << "ERROR: No sequences found in input file." << std::endl;
        return;
    }
    
    if (check_lengths) {
        check_alignment(sequences, true);
    }
    
    std::ofstream outfile(output_file);
    if (!outfile.is_open()) {
        std::cerr << "ERROR: Cannot open output file: " << output_file << std::endl;
        return;
    }
    
    for (const auto& seq : sequences) {
        write_sequence(outfile, seq);
    }
    
    outfile.close();
    std::cout << "Converted " << sequences.size() << " sequences from '" 
              << input_file << "' to '" << output_file << "'" << std::endl;
}

void extract_single(const std::string& input_file, const std::string& output_file, const std::string& seq_name) {
    std::vector<Sequence> sequences = parse_input_file(input_file);
    
    const Sequence* found = nullptr;
    for (const auto& seq : sequences) {
        if (seq.name == seq_name) {
            found = &seq;
            break;
        }
    }
    
    if (!found) {
        std::cerr << "ERROR: Sequence '" << seq_name << "' not found in input file." << std::endl;
        std::cerr << "Available sequences:" << std::endl;
        for (const auto& seq : sequences) {
            std::cerr << "  " << seq.name << std::endl;
        }
        return;
    }
    
    std::ofstream outfile(output_file);
    if (!outfile.is_open()) {
        std::cerr << "ERROR: Cannot open output file: " << output_file << std::endl;
        return;
    }
    
    write_sequence(outfile, *found);
    outfile.close();
    
    std::cout << "Extracted sequence '" << seq_name << "' to '" << output_file << "'" << std::endl;
}

void print_usage(const char* program) {
    std::cout << "Usage: " << program << " INPUT OUTPUT [OPTIONS]\n"
              << "\n"
              << "Convert HHblits/HHsearch profile format to RAxML-NG profile format.\n"
              << "\n"
              << "Arguments:\n"
              << "  INPUT          Input profile file (e.g., output_profile.tsv)\n"
              << "  OUTPUT         Output profile file for RAxML-NG\n"
              << "\n"
              << "Options:\n"
              << "  --no-check     Skip alignment length check\n"
              << "  --single NAME  Extract only a single sequence by name\n"
              << "  --help         Show this help message\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string input_file;
    std::string output_file;
    bool no_check = false;
    std::string single_name;
    
    int positional = 0;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        }
        else if (strcmp(argv[i], "--no-check") == 0) {
            no_check = true;
        }
        else if (strcmp(argv[i], "--single") == 0) {
            if (i + 1 < argc) {
                single_name = argv[++i];
            } else {
                std::cerr << "ERROR: --single requires a sequence name argument" << std::endl;
                return 1;
            }
        }
        else if (argv[i][0] != '-') {
            if (positional == 0) {
                input_file = argv[i];
                positional++;
            }
            else if (positional == 1) {
                output_file = argv[i];
                positional++;
            }
        }
        else {
            std::cerr << "ERROR: Unknown option: " << argv[i] << std::endl;
            return 1;
        }
    }
    
    if (input_file.empty() || output_file.empty()) {
        std::cerr << "ERROR: Input and output files are required." << std::endl;
        print_usage(argv[0]);
        return 1;
    }
    
    if (!single_name.empty()) {
        extract_single(input_file, output_file, single_name);
    }
    else {
        convert_profile(input_file, output_file, !no_check);
    }
    
    return 0;
}
