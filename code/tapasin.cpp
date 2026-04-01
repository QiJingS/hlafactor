#include "header.h"
// #include "input_files_parsers.cpp"
struct ind_alleles_tapasin_score
{
    std::string locus;
    std::string alleles;
    double tapasin_score;
};

std::vector<ind_alleles_tapasin_score> load_tapasin_A_B_C() {
    const std::string path = "../data_input/tapasin.txt";
    std::vector<ind_alleles_tapasin_score> taps;
    std::ifstream fin(path);
    if (!fin.is_open()) {
        throw std::runtime_error("!! ERROR! We need the tapasin data file at: " + path);
    }
    std::string line;
    while(std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string locus, alleles;
        double expression_value;
        iss >> alleles >> expression_value;
        locus = alleles.substr(0, alleles.find('*'));
        taps.push_back({locus, alleles, expression_value});
    }
    return taps;
};

struct tapasin_score_per_ind_saver{
    std::string sample_id;
    std::string locus;
    std::string alleles;
    std::string missing_alleles;
    double scores;
    
};

struct tapasin_ABC_G_final{
    std::string sample_id;
    std::string locus;
    std::string missing_alleles;
    double scores;
};

std::vector<tapasin_score_per_ind_saver> calculatoer_tapasin_A_B_C_globle(int argc, char* argv[]){
    auto orginal_c_table = load_tapasin_A_B_C();
    auto samples = parse_file(argc, argv); 
    std::vector<tapasin_score_per_ind_saver> tap_summary_saver;
    std::unordered_map<std::string, double> referece_tapasin_score;
    for(const auto &t:orginal_c_table){
        referece_tapasin_score[t.alleles] = t.tapasin_score;
    }
    std::unordered_set<std::string> tapasin_locus = {"HLA_A", "HLA_B", "HLA_C"};
    for (const auto &s : samples) {
        for (const auto& [locus_name, allele_calls] : s.loci) {
            if (tapasin_locus.count(locus_name)) {
                std::unordered_map <std::string, double> ind_temp_dic;
                for (const auto& call : allele_calls) {
                    for (size_t k = 0; k < call.alleles.size(); ++k) {
                        if( call.dosages[k] != 0 & call.alleles[k].find(':') != std::string::npos){
                            ind_temp_dic[call.alleles[k]] = call.dosages[k];
                        }      
                    }
                }
                double ind_score = 0;
                int two_a_finder = 0 ;
                std::string missing_alleles_tap;
                int v_int;
                std::string final_alleles;
                for (const auto& [key, value] : ind_temp_dic) {    
                    v_int = std::round(value);
                    if(v_int > 0){
                        final_alleles  +=  key + " ";
                    }
                    if (referece_tapasin_score.count(key)) {                            
                        ind_score += referece_tapasin_score.at(key) * v_int;
                        if(v_int == 2){
                            final_alleles += key + ' ';
                        }
                        two_a_finder += v_int;
                    }else{
                        if(v_int > 0){
                            missing_alleles_tap += key + ' ';
                        }
                    }
                }    
                double final_score = (two_a_finder == 2) ? std::log10(ind_score) : std::numeric_limits<double>::quiet_NaN();
                // std::cout << s.sampleid << '\t' << final_alleles << ' ' << ind_score << '\t' << two_a_finder <<'\t' << final_score << std::endl;
                tap_summary_saver.push_back({s.sampleid, locus_name, final_alleles, missing_alleles_tap, final_score});
            }
        }
    }
    return tap_summary_saver;
};

// int main(int argc, char* argv[]){
//     auto test = calculatoer_tapasin_A_B_C_globle(argc, argv);
//     for(const auto&t : test){
//         std::cout << t.sample_id << '\t' <<t.alleles << '\t' << t.locus << '\t' << t.missing_alleles << t.scores << std::endl;
//     }
//     return 0;
// }


