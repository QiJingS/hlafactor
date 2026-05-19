#include "header.h"

namespace {

struct AlignmentData {
    std::string locus;
    std::string reference_allele;
    std::string reference_aligned;
    std::unordered_map<std::string, std::string> allele_to_aligned;
    std::vector<int> positions;
    std::vector<bool> has_reference_position;
};

std::string trim(const std::string& s)
{
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");

    if (start == std::string::npos) {
        return "";
    }

    return s.substr(start, end - start + 1);
}

bool starts_with(const std::string& s, const std::string& prefix)
{
    return s.rfind(prefix, 0) == 0;
}

std::string strip_hla_prefix(const std::string& value)
{
    if (starts_with(value, "HLA_")) {
        return value.substr(4);
    }

    return value;
}

std::string normalize_locus_name(const std::string& locus)
{
    return strip_hla_prefix(locus);
}

std::string normalize_query_allele(const std::string& allele)
{
    std::string cleaned = strip_hla_prefix(allele);
    auto star_pos = cleaned.find('*');

    if (star_pos == std::string::npos) {
        return cleaned;
    }

    std::string locus = cleaned.substr(0, star_pos);
    std::string rest = cleaned.substr(star_pos + 1);
    std::stringstream ss(rest);
    std::string field1;
    std::string field2;

    std::getline(ss, field1, ':');
    std::getline(ss, field2, ':');

    if (field2.empty()) {
        return locus + "*" + field1;
    }

    return locus + "*" + field1 + ":" + field2;
}

std::string normalize_alignment_allele(const std::string& allele)
{
    return strip_hla_prefix(trim(allele));
}

std::string locus_from_alignment_path(const std::string& path)
{
    std::filesystem::path p(path);
    std::string filename = p.filename().string();
    auto pos = filename.find("_prot");

    if (pos == std::string::npos) {
        return filename;
    }

    return filename.substr(0, pos);
}

bool is_amino_acid_char(char c)
{
    return std::isalpha(static_cast<unsigned char>(c)) || c == 'X';
}

int parse_first_protein_position(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open alignment: " + filename);
    }

    std::string line;
    while (std::getline(in, line)) {
        std::string cleaned = trim(line);
        if (!starts_with(cleaned, "Prot")) {
            continue;
        }

        std::stringstream ss(cleaned);
        std::string token;
        while (ss >> token) {
            try {
                size_t pos = 0;
                int value = std::stoi(token, &pos);
                if (pos == token.size()) {
                    return value;
                }
            } catch (...) {
                continue;
            }
        }
    }

    throw std::runtime_error("Cannot find protein numbering header in: " + filename);
}

void build_reference_positions(AlignmentData& data, int first_position)
{
    int current_position = first_position;
    data.positions.clear();
    data.has_reference_position.clear();
    data.positions.reserve(data.reference_aligned.size());
    data.has_reference_position.reserve(data.reference_aligned.size());

    for (char c : data.reference_aligned) {
        if (is_amino_acid_char(c) || c == '*') {
            data.positions.push_back(current_position);
            data.has_reference_position.push_back(true);
            ++current_position;
            if (current_position == 0) {
                ++current_position;
            }
        } else {
            data.positions.push_back(0);
            data.has_reference_position.push_back(false);
        }
    }
}

std::string decode_alignment_char(char reference_char, char allele_char)
{
    if (allele_char == '-') {
        if (is_amino_acid_char(reference_char)) {
            return std::string(1, reference_char);
        }
        return "NA";
    }

    if (allele_char == '.') {
        return "DEL";
    }

    if (allele_char == '*') {
        return "NA";
    }

    if (allele_char == ' ' || allele_char == '\0') {
        return "NA";
    }

    if (is_amino_acid_char(allele_char)) {
        return std::string(1, allele_char);
    }

    return "NA";
}

std::string make_genotype_token(char reference_char, char a1)
{
    return decode_alignment_char(reference_char, a1) + "/.";
}

std::string make_genotype_token(char reference_char, char a1, char a2)
{
    return decode_alignment_char(reference_char, a1) + "/" +
           decode_alignment_char(reference_char, a2);
}

std::string resolve_allele_key(const AlignmentData& alignment,
                               const std::string& raw_allele)
{
    std::string normalized = normalize_query_allele(raw_allele);
    std::string best_match;
    size_t best_score = 0;

    for (const auto& entry : alignment.allele_to_aligned) {
        if (normalize_query_allele(entry.first) != normalized) {
            continue;
        }

        size_t score = 0;
        for (size_t i = 0; i < entry.second.size() &&
                           i < alignment.has_reference_position.size(); ++i) {
            if (!alignment.has_reference_position[i]) {
                continue;
            }

            char c = entry.second[i];
            if (c == '-' || c == '.' || c == 'X' || std::isalpha(static_cast<unsigned char>(c))) {
                ++score;
            }
        }

        if (best_match.empty() ||
            score > best_score ||
            (score == best_score && entry.first < best_match)) {
            best_match = entry.first;
            best_score = score;
        }
    }

    return best_match;
}

AlignmentData parse_imgt_alignment(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open alignment: " + filename);
    }

    AlignmentData data;
    data.locus = locus_from_alignment_path(filename);

    std::string line;
    bool block_reference_set = false;

    while (std::getline(in, line)) {
        if (trim(line).empty()) {
            block_reference_set = false;
            continue;
        }

        if (starts_with(trim(line), "#")) {
            continue;
        }

        if (line.find('*') == std::string::npos) {
            block_reference_set = false;
            continue;
        }

        std::stringstream ss(line);
        std::string allele;
        ss >> allele;
        allele = normalize_alignment_allele(allele);

        std::string segment_part;
        std::string aligned_segment;
        while (ss >> segment_part) {
            aligned_segment += segment_part;
        }

        if (aligned_segment.empty()) {
            continue;
        }

        data.allele_to_aligned[allele] += aligned_segment;

        if (!block_reference_set) {
            if (data.reference_allele.empty()) {
                data.reference_allele = allele;
            }

            if (allele == data.reference_allele) {
                data.reference_aligned += aligned_segment;
            }

            block_reference_set = true;
        }
    }

    if (data.reference_allele.empty() || data.reference_aligned.empty()) {
        throw std::runtime_error("Failed to parse reference sequence from: " + filename);
    }

    build_reference_positions(data, parse_first_protein_position(filename));
    return data;
}

std::vector<std::string> get_alignment_files()
{
    return {
        resource_path("IMGTHLA/alignments/A_prot.txt").string(),
        resource_path("IMGTHLA/alignments/B_prot.txt").string(),
        resource_path("IMGTHLA/alignments/C_prot.txt").string(),
        resource_path("IMGTHLA/alignments/DPA1_prot.txt").string(),
        resource_path("IMGTHLA/alignments/DPB1_prot.txt").string(),
        resource_path("IMGTHLA/alignments/DQA1_prot.txt").string(),
        resource_path("IMGTHLA/alignments/DQB1_prot.txt").string(),
        resource_path("IMGTHLA/alignments/DRA_prot.txt").string(),
        resource_path("IMGTHLA/alignments/DRB1_prot.txt").string(),
        resource_path("IMGTHLA/alignments/DRB3_prot.txt").string(),
        resource_path("IMGTHLA/alignments/DRB4_prot.txt").string(),
        resource_path("IMGTHLA/alignments/DRB5_prot.txt").string(),
        resource_path("IMGTHLA/alignments/E_prot.txt").string(),
        resource_path("IMGTHLA/alignments/F_prot.txt").string(),
        resource_path("IMGTHLA/alignments/G_prot.txt").string()
    };
}

} // namespace

struct aminoacid_genotype
{
    std::string locus;
    int position;
    std::string genotype;
};

struct ind_amiaci
{
    std::string sampleid;
    std::map<std::string, std::vector<aminoacid_genotype>> loci;

    std::vector<std::string> get_loci_names() const
    {
        std::vector<std::string> names;
        for (const auto& pair : loci) {
            names.push_back(pair.first);
        }
        return names;
    }
};

std::vector<ind_amiaci> parse_mapping_reference(std::vector<samplerecord>& samples)
{
    std::unordered_map<std::string, AlignmentData> alignments;
    for (const auto& file : get_alignment_files()) {
        AlignmentData parsed = parse_imgt_alignment(file);
        alignments.emplace(parsed.locus, std::move(parsed));
    }

    std::vector<ind_amiaci> result;
    result.reserve(samples.size());

    for (const auto& sample : samples) {
        ind_amiaci aa_ind;
        aa_ind.sampleid = sample.sampleid;

        for (const auto& locus_pair : sample.loci) {
            std::string locus = normalize_locus_name(locus_pair.first);
            auto alignment_it = alignments.find(locus);
            if (alignment_it == alignments.end()) {
                continue;
            }

            const AlignmentData& alignment = alignment_it->second;
            std::unordered_map<std::string, int> resolved_hardcalls;

            for (const auto& call : locus_pair.second) {
                if (call.alleles.empty() || call.dosages.empty()) {
                    continue;
                }

                double dosage = call.dosages[0];
                if (!std::isfinite(dosage) || dosage <= 0.0) {
                    continue;
                }

                int hardcall = static_cast<int>(std::lround(dosage));
                if (hardcall <= 0) {
                    continue;
                }

                std::string allele = resolve_allele_key(alignment, call.alleles[0]);
                if (allele.empty()) {
                    continue;
                }

                auto current = resolved_hardcalls.find(allele);
                if (current == resolved_hardcalls.end() || hardcall > current->second) {
                    resolved_hardcalls[allele] = hardcall;
                }
            }

            std::vector<std::string> allele_copies;
            for (const auto& [allele, hardcall] : resolved_hardcalls) {
                int copies = std::min(hardcall, 2);
                for (int i = 0; i < copies; ++i) {
                    allele_copies.push_back(allele);
                }
            }

            if (allele_copies.empty()) {
                continue;
            }

            if (allele_copies.size() > 2) {
                allele_copies.resize(2);
            }

            const std::string& seq1 = alignment.allele_to_aligned.at(allele_copies[0]);
            const std::string* seq2 = nullptr;
            size_t n = std::min({seq1.size(),
                                 alignment.reference_aligned.size(),
                                 alignment.positions.size(),
                                 alignment.has_reference_position.size()});

            if (allele_copies.size() >= 2) {
                seq2 = &alignment.allele_to_aligned.at(allele_copies[1]);
                n = std::min(n, seq2->size());
            }

            for (size_t i = 0; i < n; ++i) {
                if (!alignment.has_reference_position[i]) {
                    continue;
                }

                aminoacid_genotype aa;
                aa.locus = locus;
                aa.position = alignment.positions[i];
                aa.genotype = seq2
                    ? make_genotype_token(
                        alignment.reference_aligned[i],
                        seq1[i],
                        (*seq2)[i]
                    )
                    : make_genotype_token(
                        alignment.reference_aligned[i],
                        seq1[i]
                    );

                aa_ind.loci[locus].push_back(std::move(aa));
            }
        }

        result.push_back(std::move(aa_ind));
    }

    return result;
}

std::vector<ind_amiaci> map_imgt_amino_acids(int argc, char* argv[])
{
    auto samples = parse_file(argc, argv);
    return parse_mapping_reference(samples);
}

void print_amiaci_preview(const std::vector<ind_amiaci>& aa_samples)
{
    for (const auto& ind : aa_samples) {
        std::cout << "Sample: " << ind.sampleid << std::endl;

        for (const auto& locus_pair : ind.loci) {
            std::cout << "  Locus: " << locus_pair.first << std::endl;

            for (size_t i = 0;
                 i < std::min<size_t>(5, locus_pair.second.size());
                 ++i) {
                const auto& aa = locus_pair.second[i];
                std::cout << "    Pos " << aa.position
                          << " -> " << aa.genotype
                          << std::endl;
            }
        }
    }
}
