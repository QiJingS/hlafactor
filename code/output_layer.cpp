#include "header.h"
// #include "output_layer.cpp"

struct OutputSpec {
    std::filesystem::path dir;
    std::string base;
};

void log_info(const std::string& msg) {
    std::cout << "[INFO] " << msg << '\n';
}

void error_message(const std::string& msg) {
    std::cerr << num_decro_stars(msg) << '\n'
              << msg << '\n'
              << num_decro_stars(msg) << '\n';
}

OutputSpec get_output_spec(int argc, char* argv[], const std::string& default_base) {
    std::string raw;

    for (int i = 0; i < argc; ++i) {
        if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            raw = argv[i + 1];

            fs::path p(raw);
            if (p.extension() != ".txt") {
                error_message("Output file must have a .txt extension");
                exit(1);
            }
            break;
        }
    }

    if (raw.empty()) {
        std::filesystem::path dir = ".";
        std::filesystem::create_directories(dir);
        log_info("No output path provided, using default output directory: .");
        return {dir, default_base};
    }

    std::filesystem::path p(raw);
    const bool ends_with_sep = !raw.empty() && (raw.back() == '/' || raw.back() == '\\');

    if (ends_with_sep || (std::filesystem::exists(p) && std::filesystem::is_directory(p))) {
        std::filesystem::create_directories(p);
        log_info("Output directory detected: " + p.string());
        return {p, default_base};
    }

    std::filesystem::path dir = p.parent_path();
    if (dir.empty()) dir = ".";
    std::string base = p.stem().string();
    if (base.empty()) base = default_base;

    std::filesystem::create_directories(dir);

    log_info("Output file base name: " + base);
    log_info("Output directory: " + dir.string());
    return {dir, base};
}

std::filesystem::path make_outfile(const OutputSpec& out, const std::string& suffix) {
    auto filepath = out.dir / (out.base + suffix + ".txt");
    std::cout << "[INFO] Results will be saved to: " << filepath << '\n';
    return filepath;
}

void write_expression_outputs(
    const std::vector<alleles_individual_expression_A_C>& summary,
    const OutputSpec& out) {

    log_info("Writing expression outputs...");
    log_info("Input records: " + std::to_string(summary.size()));

    std::filesystem::path outA = make_outfile(out, "_A");
    std::filesystem::path outC = make_outfile(out, "_C");

    std::ofstream foutA(outA);
    std::ofstream foutC(outC);

    if (!foutA || !foutC) {
        throw std::runtime_error("Cannot open expression output files.");
    }

    size_t countA = 0, countC = 0;

    foutA << "SampleID\tlocus\talleles\tmiss_match\texpression_value\n";
    foutC << "SampleID\tlocus\talleles\tmiss_match\texpression_value\n";

    for (const auto& expr : summary) {
        if (expr.locus == "HLA_A") {
            foutA << expr.sampleid << '\t'
                  << expr.locus << '\t'
                  << expr.alleles << '\t'
                  << expr.missing_alleles << '\t'
                  << expr.expression_value << '\n';
            ++countA;
        } else if (expr.locus == "HLA_C") {
            foutC << expr.sampleid << '\t'
                  << expr.locus << '\t'
                  << expr.alleles << '\t'
                  << expr.missing_alleles << '\t'
                  << expr.expression_value << '\n';
            ++countC;
        }
    }

    std::cout << "[INFO] HLA_A rows written: " << countA << '\n';
    std::cout << "[INFO] HLA_C rows written: " << countC << '\n';
    log_info("Expression output done.");
}

void write_taps_outputs(const std::vector<tapasin_score_per_ind_saver>& taps_summary,
                        const OutputSpec& out) {

    log_info("Writing tapasin outputs...");
    log_info("Input records: " + std::to_string(taps_summary.size()));

    std::filesystem::path outA = make_outfile(out, "_A");
    std::filesystem::path outB = make_outfile(out, "_B");
    std::filesystem::path outC = make_outfile(out, "_C");
    std::filesystem::path outg = make_outfile(out, "_global");

    std::ofstream foutA(outA);
    std::ofstream foutB(outB);
    std::ofstream foutC(outC);
    std::ofstream foutg(outg);

    if (!foutA || !foutB || !foutC || !foutg) {
        throw std::runtime_error("Cannot open tapasin output files.");
    }

    foutA << "SampleID\tlocus\talleles\tmiss_match\ttapsA_value\n";
    foutB << "SampleID\tlocus\talleles\tmiss_match\ttapsB_value\n";
    foutC << "SampleID\tlocus\talleles\tmiss_match\ttapsC_value\n";
    foutg << "SampleID\ttapsA_value\ttapsB_value\ttapsC_value\ttaps_global_value\n";

    std::unordered_map<std::string, std::tuple<double, double, double>> sample_scores;
    std::vector<std::string> sample_order;

    size_t countA = 0, countB = 0, countC = 0;

    for (const auto& expr : taps_summary) {
        if (!sample_scores.count(expr.sample_id)) {
            sample_scores[expr.sample_id] = {
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()
            };
            sample_order.push_back(expr.sample_id);
        }

        auto& entry = sample_scores[expr.sample_id];

        if (expr.locus == "HLA_A") {
            foutA << expr.sample_id << '\t'
                  << expr.locus << '\t'
                  << expr.alleles << '\t'
                  << expr.missing_alleles << '\t'
                  << expr.scores << '\n';
            std::get<0>(entry) = expr.scores;
            ++countA;
        } else if (expr.locus == "HLA_B") {
            foutB << expr.sample_id << '\t'
                  << expr.locus << '\t'
                  << expr.alleles << '\t'
                  << expr.missing_alleles << '\t'
                  << expr.scores << '\n';
            std::get<1>(entry) = expr.scores;
            ++countB;
        } else if (expr.locus == "HLA_C") {
            foutC << expr.sample_id << '\t'
                  << expr.locus << '\t'
                  << expr.alleles << '\t'
                  << expr.missing_alleles << '\t'
                  << expr.scores << '\n';
            std::get<2>(entry) = expr.scores;
            ++countC;
        }
    }

    size_t countG = 0;
    for (const auto& sample_id : sample_order) {
        const auto& scores = sample_scores[sample_id];
        double a = std::get<0>(scores);
        double b = std::get<1>(scores);
        double c = std::get<2>(scores);

        double global_score = (!std::isnan(a) && !std::isnan(b) && !std::isnan(c))
                              ? (a + b + c)
                              : std::numeric_limits<double>::quiet_NaN();

        foutg << sample_id << '\t'
              << a << '\t'
              << b << '\t'
              << c << '\t'
              << global_score << '\n';
        ++countG;
    }

    std::cout << "[INFO] HLA_A rows written: " << countA << '\n';
    std::cout << "[INFO] HLA_B rows written: " << countB << '\n';
    std::cout << "[INFO] HLA_C rows written: " << countC << '\n';
    std::cout << "[INFO] Global rows written: " << countG << '\n';
    log_info("Tapasin output done.");
}

std::vector<std::string> get_sup_names(const std::string& folder_name) {
    const std::string path = resource_path(folder_name);
    std::ifstream fin(path);
    if (!fin.is_open()) {
        throw std::runtime_error("ERROR: could not open file: " + path);
    }

    std::string line;
    std::vector<std::string> supertype_names_saver;

    if (std::getline(fin, line)) {
        std::istringstream iss(line);
        std::string name;
        while (iss >> name) {
            supertype_names_saver.push_back(name);
        }
    }
    return supertype_names_saver;
}

void write_supertypes_app_results(const std::vector<ind_supertypes>& supertypes_res_const,
                                  const OutputSpec& out) {
    log_info("Writing supertypes output...");
    log_info("Input records: " + std::to_string(supertypes_res_const.size()));

    auto supertype_names_A = get_sup_names("superA.txt");
    auto supertype_names_B = get_sup_names("superB.txt");
    auto supertype_names_classII = get_sup_names("super_class_II.txt");

    std::vector<std::string> supertype_names_set;
    supertype_names_set.insert(supertype_names_set.end(), supertype_names_A.begin(), supertype_names_A.end());
    supertype_names_set.insert(supertype_names_set.end(), supertype_names_B.begin(), supertype_names_B.end());
    supertype_names_set.insert(supertype_names_set.end(), supertype_names_classII.begin(), supertype_names_classII.end());

    std::filesystem::path outg = make_outfile(out, "_supertypes");
    std::ofstream foutg(outg);

    if (!foutg) {
        throw std::runtime_error("Cannot open supertypes output file.");
    }

    std::string line_names = "SampleID\t";
    for (const auto& n : supertype_names_set) {
        line_names += n + '\t';
    }
    foutg << line_names << '\n';

    std::unordered_map<std::string, std::vector<double>> dosage_by_id;
    std::vector<std::string> sample_order;

    for (const auto& res : supertypes_res_const) {
        if (dosage_by_id.find(res.sampleid) == dosage_by_id.end()) {
            dosage_by_id[res.sampleid] = std::vector<double>(supertype_names_set.size(), 0.0);
            sample_order.push_back(res.sampleid);
        }

        auto& row = dosage_by_id[res.sampleid];
        for (size_t i = 0; i < res.supertypes.size(); ++i) {
            const auto& supertype_name = res.supertypes[i];
            auto it = std::find(supertype_names_set.begin(), supertype_names_set.end(), supertype_name);
            if (it != supertype_names_set.end() && i < res.dosages.size()) {
                size_t index = std::distance(supertype_names_set.begin(), it);
                row[index] += res.dosages[i];
            }
        }
    }

    for (const auto& sampleid : sample_order) {
        foutg << sampleid;
        const auto& row = dosage_by_id[sampleid];
        for (double v : row) {
            foutg << '\t' << v;
        }
        foutg << '\n';
    }

    std::cout << "[INFO] Supertypes written for samples: " << sample_order.size() << '\n';
    log_info("Supertypes output done.");
}

void write_hla_heterozygosity_results(const std::vector<ind_heter>& heter_res,
                                      const OutputSpec& out) {
    if (heter_res.empty()) {
        log_info("Heterozygosity input is empty, skip writing.");
        return;
    }

    log_info("Writing heterozygosity output...");
    log_info("Input records: " + std::to_string(heter_res.size()));

    std::filesystem::path outg = make_outfile(out, "_heterozygosity");
    std::ofstream foutg(outg);

    if (!foutg) {
        throw std::runtime_error("Cannot open heterozygosity output file.");
    }

    std::vector<std::string> locus_set;
    std::string sam1_id = heter_res[0].sampleid;

    for (const auto& res : heter_res) {
        if (res.sampleid != sam1_id) break;
        if (res.locus[0] == 'A') continue;
        locus_set.push_back(res.locus);
    }

    std::string lines_flow = "SampleID\t";
    for (const auto& locus : locus_set) {
        lines_flow += locus + '\t';
    }
    foutg << lines_flow << '\n';

    std::unordered_map<std::string, std::vector<std::string>> merged;
    std::vector<std::string> sample_order;
    std::unordered_map<std::string, size_t> locus_idx;

    for (size_t i = 0; i < locus_set.size(); ++i) {
        locus_idx[locus_set[i]] = i;
    }

    for (const auto& res : heter_res) {
        if (merged.find(res.sampleid) == merged.end()) {
            merged[res.sampleid] = std::vector<std::string>(locus_set.size(), "");
            sample_order.push_back(res.sampleid);
        }

        auto& row = merged[res.sampleid];
        auto it = locus_idx.find(res.locus);
        if (it != locus_idx.end()) {
            size_t pos = it->second;
            if (row[pos].empty()) {
                row[pos] = res.het_hom;
            }
        }
    }

    for (const auto& sampleid : sample_order) {
        foutg << sampleid;
        const auto& row = merged[sampleid];
        for (const auto& v : row) {
            foutg << '\t';
            if (!v.empty()) foutg << v;
        }
        foutg << '\n';
    }

    std::cout << "[INFO] Heterozygosity samples written: " << sample_order.size() << '\n';
    log_info("Heterozygosity output done.");
}

void write_HLA_variant_NK_results(const std::vector<sample_variants_NK>& NK_res, const OutputSpec& out) {
    log_info("Writing HLA variant NK output...");
    log_info("Input records: " + std::to_string(NK_res.size()));

    std::filesystem::path outg = make_outfile(out, "_HLA_variant_NK");
    std::ofstream foutg(outg);

    if (!foutg) {
        throw std::runtime_error("Cannot open HLA variant NK output file.");
    }

    foutg << "SampleID\tAA_B_-21\tAA_B_82\tAA_C_80\n";
    for (const auto& res : NK_res) {
        foutg << res.sample_id;
        for (size_t i = 0; i < res.genotypes.size(); ++i) {
            foutg << '\t' << res.genotypes[i];
        }
        foutg << '\n';
    }

    log_info("HLA variant NK output done.");
}

void write_hla_functional_heterozygosity_results(const std::vector<ind_heter>& heter_res,
                                                 const OutputSpec& out) {
    if (heter_res.empty()) {
        log_info("Functional heterozygosity input is empty, skip writing.");
        return;
    }

    log_info("Writing functional heterozygosity output...");
    log_info("Input records: " + std::to_string(heter_res.size()));

    std::filesystem::path outg = make_outfile(out, "_functional_zygosity");
    std::ofstream foutg(outg);
    if (!foutg) {
        throw std::runtime_error("Cannot open output file: " + outg.string());
    }

    std::unordered_map<std::string, std::unordered_map<std::string, ind_heter>> mp;

    for (const auto& res : heter_res) {
        if (res.locus == "HLA_A" || res.locus == "HLA_B" || res.locus == "HLA_C") {
            if (mp[res.sampleid].find(res.locus) == mp[res.sampleid].end()) {
                mp[res.sampleid][res.locus] = res;
            }
        }
    }

    foutg << "sampleid\tHLA_A_allele\tHLA_A_FH\tHLA_B_allele\tHLA_B_FH\tHLA_C_allele\tHLA_C_FH\n";

    size_t sample_count = 0;
    for (const auto& [sampleid, loci] : mp) {
        auto get_allele = [&](const std::string& loc) -> std::string {
            auto it = loci.find(loc);
            return (it != loci.end()) ? it->second.alleles : "NA";
        };

        auto get_fh = [&](const std::string& loc) -> std::string {
            auto it = loci.find(loc);
            if (it != loci.end()) return std::to_string(it->second.het_hom_fh);
            return "NA";
        };

        foutg << sampleid << '\t'
              << get_allele("HLA_A") << '\t' << get_fh("HLA_A") << '\t'
              << get_allele("HLA_B") << '\t' << get_fh("HLA_B") << '\t'
              << get_allele("HLA_C") << '\t' << get_fh("HLA_C") << '\n';
        ++sample_count;
    }

    std::cout << "[INFO] Functional heterozygosity samples written: " << sample_count << '\n';
    log_info("Functional heterozygosity output done.");
}