#include "header.h"
struct ind_alleles_protein_expression_A_C
{
    std::string locus;
    std::string alleles;
    double expression_value;
};

struct all_ind_alleles_protein_expression_A_C
{
    std::string sampleid;
    std::vector<ind_alleles_protein_expression_A_C> expressions;
};

struct alleles_individual_expression_A_C
{
    std::string sampleid;
    std::string locus;
    std::string alleles;
    std::string missing_alleles;
    double expression_value;
};

std::vector<ind_alleles_protein_expression_A_C> load_expression_A_C()
{
    auto path = resource_path("expression.txt");
    std::vector<ind_alleles_protein_expression_A_C> exprs;
    std::ifstream fin(path);

    if (!fin.is_open()) {
        throw std::runtime_error("!! ERROR! We need the expression data file at data_input/expression_A_C.txt");
    }

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string alleles;
        std::string expr_str;

        if (!(iss >> alleles >> expr_str)) {
            continue;  
        }
        if (alleles == "alleles" || alleles == "SNPID" || alleles == "locus") {
            continue;
        }
        if (alleles.find('*') == std::string::npos) {
            continue;
        }

        double expression_value;
        try {
            if (expr_str == "NA" || expr_str == "NaN" || expr_str == "nan") {
                continue;
            }
            expression_value = std::stod(expr_str);
        } catch (...) {
            continue;
        }

        if (!std::isfinite(expression_value)) {
            continue;
        }

        std::string locus = alleles.substr(0, alleles.find('*'));
        exprs.push_back({locus, alleles, expression_value});
    }

    return exprs;
}

std::vector<alleles_individual_expression_A_C> alleles_expression_A_C_summary(int argc, char* argv[])
{
    auto samples = parse_file(argc, argv);
    auto exp_ref = load_expression_A_C();

    std::vector<alleles_individual_expression_A_C> sam_expression_score_saver;
    std::unordered_map<std::string, double> expression_map;

    for (const auto& a : exp_ref) {
        expression_map[a.alleles] = a.expression_value;
    }

    std::unordered_set<std::string> expression_locus = {"HLA_A", "HLA_C"};

    for (const auto& s : samples) {
        for (const auto& [locus_name, allele_calls] : s.loci) {
            if (!expression_locus.count(locus_name)) {
                continue;
            }

            std::unordered_map<std::string, double> ind_temp_dic;

            for (const auto& call : allele_calls) {
                const size_t n = std::min(call.alleles.size(), call.dosages.size());

                for (size_t k = 0; k < n; ++k) {
                    const std::string& allele = call.alleles[k];
                    double dosage = call.dosages[k];

                    if (dosage == 0.0) continue;
                    if (!std::isfinite(dosage)) continue;
                    if (allele.find(':') != std::string::npos) continue;

                    ind_temp_dic[allele] += dosage;
                    ind_temp_dic[allele] += dosage;
                }
            }

            double ind_score = 0.0;
            int two_a_finder = 0;
            std::string final_alleles;
            std::string missing_alleles_exp;

            for (const auto& [key, value] : ind_temp_dic) {
                int v_int = static_cast<int>(std::lround(value));
                if (v_int > 0) {
                    final_alleles += key + " ";
                }

                if (expression_map.count(key)) {
                    two_a_finder += v_int;
                    ind_score += expression_map.at(key) * v_int;
                } else {
                    if (v_int == 1) {
                        missing_alleles_exp += key + " ";
                    }
                }
            }

            double final_score = (two_a_finder == 2 && !missing_alleles_exp.empty())
                ? ind_score
                : std::numeric_limits<double>::quiet_NaN();

            sam_expression_score_saver.push_back({
                s.sampleid,
                locus_name,
                final_alleles,
                missing_alleles_exp,
                final_score
            });
        }
    }

    return sam_expression_score_saver;
}