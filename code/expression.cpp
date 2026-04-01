#include "header.h"
// #include "input_files_parsers.cpp"
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

std::vector<ind_alleles_protein_expression_A_C> load_expression_A_C() {
    const std::string path = "../data_input/expression.txt";
    std::vector<ind_alleles_protein_expression_A_C> exprs;
    std::ifstream fin(path);
    if (!fin.is_open()) {
        throw std::runtime_error("!! ERROR! We need the expression data file at: " + path);
    }
    std::string line;
    while(std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string locus, alleles;
        double expression_value;
        iss >> alleles >> expression_value;
        locus = alleles.substr(0, alleles.find('*'));
        exprs.push_back({locus, alleles, expression_value});
    }
    return exprs;
}
// 
struct alleles_individual_expression_A_C
{
    std::string sampleid;
    std::string locus;
    std::string alleles;
    std::string missing_alleles;
    double expression_value;
};
std::vector<alleles_individual_expression_A_C> alleles_expression_A_C_summary(int argc, char* argv[]) {
    auto samples = parse_file(argc, argv); 
    auto exp_ref = load_expression_A_C();
    std::vector <alleles_individual_expression_A_C> sam_expression_score_saver;
    std::unordered_map<std::string, double> expression_map; 
    // save expression table to hashtable
    for(const auto &a : exp_ref){
        expression_map[a.alleles] = a.expression_value;
    }
    std::unordered_set<std::string> expression_locus = {"HLA_A","HLA_C"};
    for (const auto &s : samples) {
        for (const auto& [locus_name, allele_calls] : s.loci) {
            if (expression_locus.count(locus_name)) {
                std::unordered_map <std::string, double> ind_temp_dic;
                for (const auto& call : allele_calls) {
                    for (size_t k = 0; k < call.alleles.size(); ++k) {
                        if( call.dosages[k] != 0 & call.alleles[k].find(':') == std::string::npos){
                            ind_temp_dic[call.alleles[k]] = call.dosages[k];
                        }      
                    }
                }
                double ind_score = 0;
                int two_a_finder = 0;
                std::string final_alleles;
                std::string missing_alleles_exp;
                int v_int;
                for (const auto& [key, value] : ind_temp_dic) {
                    
                    v_int = std::round(value);
                    if(v_int > 0){
                        final_alleles  +=  key + " ";
                    }
                    if (expression_map.count(key)) {
                        two_a_finder += v_int;
                        if(v_int == 2){
                            final_alleles  +=  key + " ";
                        }
                        ind_score += expression_map.at(key) * v_int;
                    }else{
                        if(v_int > 0){
                            missing_alleles_exp += key + ' ';
                        }
                    }
                }    
                double final_score = (two_a_finder == 2) ? ind_score : std::numeric_limits<double>::quiet_NaN();
                // std::cout << final_alleles << '\t' << missing_alleles_exp << '\n';
                sam_expression_score_saver.push_back({s.sampleid, locus_name, final_alleles, missing_alleles_exp, final_score});
            }
        }
    }
    return sam_expression_score_saver;
}

