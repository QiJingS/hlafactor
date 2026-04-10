#include "header.h"
// #include "input_files_parsers.cpp"
struct ind_supertypes
{
    std::string sampleid;
    std::string locus;
    std::vector<std::string> supertypes;
    std::vector<double> dosages;
};
std::unordered_map<std::string, std::string> reference_supertypes(const std::string& folder_name) {
    auto path = resource_path(folder_name);
    // std::cout << "path = " << path << std::endl;
    std::ifstream fin(path);
    if (!fin.is_open()) {
        throw std::runtime_error("ERROR: could not open file at data_input/supertype.txt");
    }
    std::string line;
    std::vector<std::string> supertype_names_saver;
    std::unordered_map<std::string, std::string> supertypes_reference;
    if (std::getline(fin, line)) {
        std::istringstream iss(line);
        std::string name;
        while (iss >> name) {
            supertype_names_saver.push_back(name);
        }
    }
    int ind_temp;
    while (std::getline(fin, line)) {
        ind_temp = 0;
        std::istringstream iss(line);
        std::string name;
        while (iss >> name) { 
            supertypes_reference[name] = supertype_names_saver[ind_temp];
            ind_temp++;
        }
    }
    return supertypes_reference;
}

std::string getPrefix(const std::string& allele) {
    size_t pos = allele.find('*');
    return (pos != std::string::npos) ? allele.substr(0, pos) : allele;
};
std::vector <ind_supertypes> supertypes_mapping_data(int argc, char* argv[]){
    std::vector <ind_supertypes> outp;
    auto sup_A = reference_supertypes("superA.txt");
    auto sup_B = reference_supertypes("superB.txt");
    auto class_II = reference_supertypes("super_class_II.txt");

    auto combined = sup_A;

    combined.insert(sup_B.begin(), sup_B.end());
    combined.insert(class_II.begin(), class_II.end());
    
    std::unordered_set<std::string> HLA_locus;
    HLA_locus.insert({"HLA_A","HLA_B"});
    for(const auto&[key, value]: class_II){
        std::string prefix = getPrefix(key);
        if (HLA_locus.insert(prefix).second) {}
    }
    auto samples = parse_file(argc, argv); 
    
    int k = 0;

    for(const auto &s : samples){
        std::vector<std::string> temp_ind_user_saver;
        std::vector<std::string> temp_ind_super_saver;
        std::vector<double> temp_ind_dosages_saver;
        for (const auto& [locus_name, allele_calls] : s.loci) {
            if (HLA_locus.count(locus_name)) {
                for(const auto &a : allele_calls){   
                        if(a.dosages[k] != 0 && a.alleles[k].find(':') != std::string::npos){
                            temp_ind_user_saver.push_back(a.alleles[k]);
                            if(a.dosages[k] == 2){
                                temp_ind_dosages_saver.push_back(a.dosages[k]);
                                temp_ind_super_saver.push_back(combined[a.alleles[k]]);
                                continue;
                            }else{
                                auto it = combined.find(a.alleles[k]);
                                if (it == combined.end()) continue;
                                const std::string& super = it->second;
                                bool finder = true;
                                for(size_t t = 0; t < temp_ind_super_saver.size(); t++){                                                                    
                                    if(temp_ind_super_saver[t] == super){ 
                                        if(t < temp_ind_dosages_saver.size()){
                                            temp_ind_dosages_saver[t] += a.dosages[k];
                                        }
                                        finder = false;
                                    }
                                }
                                if(finder){
                                    temp_ind_dosages_saver.push_back(a.dosages[k]);
                                    temp_ind_super_saver.push_back(combined[a.alleles[k]]);
                                }
                            }    
                    }
                }
                outp.push_back({s.sampleid, locus_name, temp_ind_super_saver, temp_ind_dosages_saver});
                temp_ind_super_saver.clear();
                temp_ind_dosages_saver.clear();
            }
        }
    }
    return outp;
}

// int main( int argc, char* argv[]) {
//     auto test = supertypes_mapping_data(argc, argv);
//     for(auto &ttt :test){
//         std::cout << ttt.sampleid << '\t' << ttt.locus <<"\t" ;
//         for(size_t i = 0 ; i < ttt.supertypes.size(); i++){
//             std::cout << ttt.supertypes[i] << '\t' << ttt.dosages[i] <<"\t" ;
//         }
//         std::cout << '\n';
//     }
//     return 0;
// }