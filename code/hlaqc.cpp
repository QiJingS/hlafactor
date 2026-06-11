#include "header.h"

namespace fs = std::filesystem;
bool in_HLA_region(int pos, int hg) {
    int pos1 = pos + 1;
    if (hg == 38) {
        return (pos1 >= 25726063 && pos1 <= 33400644);
    } else if (hg == 19) {
        return (pos1 >= 25726291 && pos1 <= 33368421);
    }
    return false;
}

bool is_chr6(const std::string& chr) {
    return (chr == "6" || chr == "chr6");
}

bool id_starts_with_HLA_or_AA(bcf1_t* rec) {
    bcf_unpack(rec, BCF_UN_STR);
    if (!rec->d.id) return false;

    std::string id = rec->d.id;

    return (id.rfind("HLA", 0) == 0) ||
           (id.rfind("AA_B_82", 0) == 0) ||
           (id.rfind("AA_B_-21", 0) == 0) ||
           (id.rfind("AA_C_80", 0) == 0);
}
static std::string resolve_output_path(const std::string& output,
                                       const std::string& default_filename = "filtered.vcf.gz") {
    fs::path p(output);
    if (fs::exists(p) && fs::is_directory(p)) {
        return (p / default_filename).string();
    }
    if (!output.empty() && output.back() == '/') {
        return (p / default_filename).string();
    }
    return output;
}
static std::string choose_hts_mode(const std::string& output) {
    if (output.size() >= 4 && output.rfind(".bcf") == output.size() - 4) {
        return "wb";
    }
    if (output.size() >= 3 && output.rfind(".gz") == output.size() - 3) {
        return "wz";
    }
    return "w";
}

static bool calculate_af_from_gt(bcf_hdr_t* hdr, bcf1_t* rec, int n_samples, float& af) {
    if (n_samples <= 0) return false;

    int32_t* gt_arr = NULL;
    int ngt_arr = 0;
    int nret = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);

    if (nret <= 0) {
        free(gt_arr);
        return false;
    }

    int ploidy = nret / n_samples;
    if (ploidy <= 0) {
        free(gt_arr);
        return false;
    }

    int alt_count = 0;
    for (int i = 0; i < n_samples; ++i) {
        int base = i * ploidy;
        int alleles_to_read = std::min(ploidy, 2);

        for (int j = 0; j < alleles_to_read; ++j) {
            int32_t gt = gt_arr[base + j];
            if (gt == bcf_int32_vector_end) break;
            if (bcf_gt_is_missing(gt)) continue;

            int allele = bcf_gt_allele(gt);
            if (allele > 0) {
                ++alt_count;
            }
        }
    }

    af = static_cast<float>(alt_count) / static_cast<float>(n_samples * 2);
    free(gt_arr);
    return true;
}

static bool qc_is_missing_token(const std::string& token) {
    return token.empty() ||
           token == "."   ||
           token == "NA"  ||
           token == "na"  ||
           token == "Na"  ||
           token == "nA"  ||
           token == "NaN" ||
           token == "nan" ||
           token == "NAN";
}

static bool qc_parse_dosage(const std::string& token, double& value) {
    if (qc_is_missing_token(token)) return false;

    try {
        size_t pos = 0;
        value = std::stod(token, &pos);
        return pos == token.size() && std::isfinite(value);
    } catch (...) {
        return false;
    }
}

static bool qc_id_starts_with_HLA_or_AA(const std::string& id) {
    return (id.rfind("HLA", 0) == 0) ||
           (id.rfind("AA_B_82", 0) == 0) ||
           (id.rfind("AA_B_-21", 0) == 0) ||
           (id.rfind("AA_C_80", 0) == 0);
}

static bool passes_af_filter(double af, float af_min, float af_max,
                             bool use_af_min, bool use_af_max) {
    if (use_af_min && af < af_min) return false;
    if (use_af_max && af > af_max) return false;
    return true;
}

static bool prepare_text_output(const std::string& output, const std::string& default_filename,
                                std::ofstream& fout, std::string& outpath) {
    outpath = resolve_output_path(output, default_filename);
    fs::path outp(outpath);
    if (outp.has_parent_path()) {
        std::error_code ec;
        fs::create_directories(outp.parent_path(), ec);
        if (ec) {
            std::string error_msg = "!!ERROR: cannot create output directory: " + outp.parent_path().string() + ". " + ec.message();
            std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
            return false;
        }
    }

    fout.open(outpath);
    if (!fout.is_open()) {
        std::string error_msg = "!!ERROR: cannot open output: " + outpath;
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        return false;
    }
    return true;
}

static void print_filter_summary(int n_samples, int total_variants, int kept_variants,
                                 const std::string& outpath) {
    std::cout << "\nThanks for using the HLAfactor :)\n";
    std::cout << "*********************************\n";
    std::cout << "Samples (individuals): " << n_samples << "\n";
    std::cout << "Total variants read:  " << total_variants << "\n";
    std::cout << "Variants kept:        " << kept_variants << "\n";
    std::cout << "Variants removed:     " << (total_variants - kept_variants) << "\n";
    std::cout << "*********************************\n";
    std::cout << "\nFiltering complete, please find the output file at: " << outpath << "\n";
}

static void filter_txt_by_dosage_af(const std::string& input, const std::string& output,
                                    float af_min, float af_max,
                                    bool use_af_min, bool use_af_max) {
    std::ifstream fin(input);
    if (!fin.is_open()) {
        std::string error_msg = "!!ERROR: cannot open input file: " + input;
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        return;
    }

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        std::string error_msg = "!!ERROR: empty input file: " + input;
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        return;
    }

    std::istringstream header_stream(header_line);
    std::vector<std::string> header_tokens;
    std::string token;
    while (header_stream >> token) {
        header_tokens.push_back(token);
    }

    if (header_tokens.size() < 2) {
        std::string error_msg = "!!ERROR: malformed txt header: " + input;
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        return;
    }

    std::vector<std::vector<std::string>> rows;
    std::vector<double> dosage_sums(header_tokens.size(), 0.0);

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        std::istringstream line_stream(line);
        std::vector<std::string> row;
        while (line_stream >> token) {
            row.push_back(token);
        }
        if (row.empty()) continue;

        row.resize(header_tokens.size(), "NA");
        for (size_t i = 1; i < header_tokens.size(); ++i) {
            double dosage = 0.0;
            if (qc_parse_dosage(row[i], dosage)) {
                dosage_sums[i] += dosage;
            }
        }
        rows.push_back(std::move(row));
    }

    int n_samples = static_cast<int>(rows.size());
    double denom = static_cast<double>(n_samples * 2);
    std::vector<size_t> keep_indices;
    keep_indices.push_back(0);

    int total_variants = 0;
    for (size_t i = 1; i < header_tokens.size(); ++i) {
        if (!qc_id_starts_with_HLA_or_AA(header_tokens[i])) continue;
        ++total_variants;

        double af = denom > 0.0 ? dosage_sums[i] / denom : 0.0;
        if ((use_af_min || use_af_max) &&
            !passes_af_filter(af, af_min, af_max, use_af_min, use_af_max)) {
            continue;
        }
        keep_indices.push_back(i);
    }

    std::ofstream fout;
    std::string outpath;
    if (!prepare_text_output(output, "filtered.txt", fout, outpath)) return;

    for (size_t i = 0; i < keep_indices.size(); ++i) {
        if (i) fout << '\t';
        fout << header_tokens[keep_indices[i]];
    }
    fout << '\n';

    for (const auto& row : rows) {
        for (size_t i = 0; i < keep_indices.size(); ++i) {
            if (i) fout << '\t';
            fout << row[keep_indices[i]];
        }
        fout << '\n';
    }

    print_filter_summary(n_samples, total_variants, static_cast<int>(keep_indices.size()) - 1, outpath);
}

static void filter_dosage_by_dosage_af(int hg, const std::string& input, const std::string& output,
                                       float af_min, float af_max,
                                       bool use_af_min, bool use_af_max) {
    std::ifstream fin(input);
    if (!fin.is_open()) {
        std::string error_msg = "!!ERROR: cannot open input file: " + input;
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        return;
    }

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        std::string error_msg = "!!ERROR: empty input file: " + input;
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        return;
    }

    std::ofstream fout;
    std::string outpath;
    if (!prepare_text_output(output, "filtered.dosage", fout, outpath)) return;
    fout << header_line << '\n';

    std::istringstream header_stream(header_line);
    std::vector<std::string> header_tokens;
    std::string token;
    while (header_stream >> token) {
        header_tokens.push_back(token);
    }

    if (header_tokens.size() < 7) {
        std::string error_msg = "!!ERROR: malformed dosage header: " + input;
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        return;
    }

    int n_samples = static_cast<int>(header_tokens.size() - 6);
    double denom = static_cast<double>(n_samples * 2);
    int total_variants = 0;
    int kept_variants = 0;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        std::istringstream line_stream(line);
        std::vector<std::string> fields;
        while (line_stream >> token) {
            fields.push_back(token);
        }
        if (fields.size() < 6) continue;

        const std::string& chromosome = fields[0];
        const std::string& snpid = fields[1];
        int pos1 = 0;
        try {
            size_t pos = 0;
            pos1 = std::stoi(fields[3], &pos);
            if (pos != fields[3].size()) continue;
        } catch (...) {
            continue;
        }

        if (!(is_chr6(chromosome) && in_HLA_region(pos1 - 1, hg))) continue;
        if (!qc_id_starts_with_HLA_or_AA(snpid)) continue;

        ++total_variants;

        double dosage_sum = 0.0;
        for (size_t i = 6; i < fields.size(); ++i) {
            double dosage = 0.0;
            if (qc_parse_dosage(fields[i], dosage)) {
                dosage_sum += dosage;
            }
        }

        double af = denom > 0.0 ? dosage_sum / denom : 0.0;
        if ((use_af_min || use_af_max) &&
            !passes_af_filter(af, af_min, af_max, use_af_min, use_af_max)) {
            continue;
        }

        fout << line << '\n';
        ++kept_variants;
    }

    print_filter_summary(n_samples, total_variants, kept_variants, outpath);
}

void filter_HLA_region(int hg, const std::string& input, const std::string& output,
                       float af_min, float af_max,
                       bool use_af_min, bool use_af_max,
                       float r2_min, float r2_max,
                       bool use_r2_min, bool use_r2_max) {
    std::string outpath = resolve_output_path(output);
    std::string mode = choose_hts_mode(outpath);
    fs::path outp(outpath);
    if (outp.has_parent_path()) {
        std::error_code ec;
        fs::create_directories(outp.parent_path(), ec);
        if (ec) {
			std::string error_msg = "!!ERROR: cannot create output directory: " + outp.parent_path().string() + ". " + ec.message();
			std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
            return;
        }
    }
    htsFile* in = bcf_open(input.c_str(), "r");
    if (!in) {
		std::string error_msg = "!!ERROR: cannot open input file: " + input;
		std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        return;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(in);
    if (!hdr) {
		std::string error_msg = "!!ERROR: cannot read header from input file: " + input;
		std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        bcf_close(in);
        return;
    }
    int n_samples = bcf_hdr_nsamples(hdr);
    htsFile* out = bcf_open(outpath.c_str(), mode.c_str());
    if (!out) {
        std::string error_msg = "!!ERROR: cannot open output: " + outpath;
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        bcf_hdr_destroy(hdr);
        bcf_close(in);
        return;
    }
    if (bcf_hdr_write(out, hdr) != 0) {
        std::string error_msg = "!!ERROR: cannot write header\n";
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        bcf_hdr_destroy(hdr);
        bcf_close(in);
        bcf_close(out);
        return;
    }
    bcf1_t* rec = bcf_init();
    float* af_arr = NULL;
    int naf = 0;
    float* r2_arr = NULL;
    int nr2 = 0;

    int total_variants = 0;
    int kept_variants = 0;

    while (bcf_read(in, hdr, rec) == 0) {
        total_variants++;

        int rid = rec->rid;
        if (rid < 0) continue;

        const char* chrom = bcf_hdr_id2name(hdr, rid);
        int pos = rec->pos;

        if (!(chrom && is_chr6(chrom) && in_HLA_region(pos, hg))) continue;
        if (!id_starts_with_HLA_or_AA(rec)) continue;

        if (use_af_min || use_af_max) {
            float af = 0.0f;
            if (bcf_get_info_float(hdr, rec, "AF", &af_arr, &naf) > 0) {
                af = af_arr[0];
            } else if (!calculate_af_from_gt(hdr, rec, n_samples, af)) {
                continue;
            }
            if (use_af_min && af < af_min) continue;
            if (use_af_max && af > af_max) continue;
        }

        if (use_r2_min || use_r2_max) {
            if (bcf_get_info_float(hdr, rec, "R2", &r2_arr, &nr2) <= 0) continue;
            float r2 = r2_arr[0];
            if (use_r2_min && r2 < r2_min) continue;
            if (use_r2_max && r2 > r2_max) continue;
        }

        if (bcf_write(out, hdr, rec) != 0) {
            std::cerr << "Error writing record\n";
            break;
        }
        kept_variants++;
    }

    free(af_arr);
    free(r2_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(in);
    bcf_close(out);

    std::cout << "\nThanks for using the HLAfactor :)\n";
    std::cout << "*********************************\n";
    std::cout << "Samples (individuals): " << n_samples << "\n";
    std::cout << "Total variants read:  " << total_variants << "\n";
    std::cout << "Variants kept:        " << kept_variants << "\n";
    std::cout << "Variants removed:     " << (total_variants - kept_variants) << "\n";
    std::cout << "*********************************\n";
    std::cout << "\nFiltering complete, please find the output file at: " << outpath << "\n";
}

void function_hla_qc(int argc, char* argv[]) {
    std::string input;
    std::string output;  
    bool has_output = false;
    int hg = 0;

    float af_min = 0, af_max = 0;
    bool use_af_min = false, use_af_max = false;

    float r2_min = 0, r2_max = 0;
    bool use_r2_min = false, use_r2_max = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-o" && i + 1 < argc) {
            output = argv[++i];
            has_output = true;
        }
        else if (arg == "-hlarg38") {
            hg = 38;
        }
        else if (arg == "-hlarg37") {
            hg = 19;
        }
        else if (arg == "-i" && i + 1 < argc) {
            input = argv[++i];
        }
        else if (arg == "-af" && i + 1 < argc) {
            af_min = std::stof(argv[++i]);
            use_af_min = true;

            if (i + 1 < argc && argv[i + 1][0] != '-') {
                af_max = std::stof(argv[++i]);
                use_af_max = true;
            }
        }
        else if (arg == "-impr2" && i + 1 < argc) {
            r2_min = std::stof(argv[++i]);
            use_r2_min = true;

            if (i + 1 < argc && argv[i + 1][0] != '-') {
                r2_max = std::stof(argv[++i]);
                use_r2_max = true;
            }
        }
    }

    if (input.empty()) {
        std::string error_msg = "!!Error: no input file (-i)";
        std::cerr << error_msg << '\n';
        return;
    }

    if (!has_output || output.empty()) {
        output = "filtered.vcf";
    }

    if (hg == 0) {
        std::string error_msg = "!!Error: specify -hlarg37 or -hlarg38";
        std::cerr <<  error_msg << '\n';
        return;
    }

    fs::path input_path(input);
    std::string ext = input_path.extension().string();

    if (ext == ".txt") {
        if (use_r2_min || use_r2_max) {
            std::cerr << "!!ERROR: -impr2 requires imputation R2 information, but .txt input does not contain an R2 field.\n";
            return;
        }
        filter_txt_by_dosage_af(input, output, af_min, af_max, use_af_min, use_af_max);
    } else if (ext == ".dosage") {
        if (use_r2_min || use_r2_max) {
            std::cerr << "!!ERROR: -impr2 requires imputation R2 information, but .dosage input does not contain an R2 field.\n";
            return;
        }
        filter_dosage_by_dosage_af(hg, input, output, af_min, af_max, use_af_min, use_af_max);
    } else {
        filter_HLA_region(hg, input, output,
                          af_min, af_max, use_af_min, use_af_max,
                          r2_min, r2_max, use_r2_min, use_r2_max);
    }
}
// int main(int argc, char* argv[]) {
//     if (argc < 2) {
//         std::cout << "Usage:\n";
//         std::cout << "  ./hla_filter -hlarg38 -i input.vcf.gz -o out.vcf.gz\n";
//         std::cout << "Options:\n";
//         std::cout << "  -af x [y]      (AF >= x OR x <= AF <= y)\n";
//         std::cout << "  -impr2 x [y]   (R2 >= x OR x <= R2 <= y)\n";
//         return 0;
//     }
//
//     function_hla_qc(argc, argv);
//     return 0;
// }
