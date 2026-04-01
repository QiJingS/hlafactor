#include "header.h"
// #include "input_files_parsers.cpp"

struct PairKey {
    std::string a;
    std::string b;

    PairKey(const std::string& x, const std::string& y) {
        if (x < y) { a = x; b = y; }
        else       { a = y; b = x; }
    }

    bool operator==(const PairKey& other) const {
        return a == other.a && b == other.b;
    }
};
struct PairKeyHash {
    size_t operator()(const PairKey& k) const {
        size_t h1 = std::hash<std::string>{}(k.a);
        size_t h2 = std::hash<std::string>{}(k.b);
        return h1 ^ (h2 << 1);
    }
};
std::unordered_map<PairKey, double, PairKeyHash> read_reference_fh(const std::string& path) {
    std::unordered_map<PairKey, double, PairKeyHash> mapping;
    std::ifstream fin(path);
    if (!fin.is_open()) {
        throw std::runtime_error("ERROR: cannot open FH file: " + path);
    }
    std::string line;
    bool first_line = true;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string allele1, allele2;
        double fh;
        if (first_line) {
            first_line = false;
            if (!(iss >> allele1 >> allele2 >> fh)) {
                continue;
            }
        } else {
            if (!(iss >> allele1 >> allele2 >> fh)) {
                continue;
            }
        }
        mapping.emplace(PairKey(allele1, allele2), fh);
    }
    return mapping;
}


struct ind_heter
{
    std::string sampleid;
    std::string locus;
    std::string het_hom;
    std::string alleles;
    double het_hom_fh;
};

std::vector<ind_heter> map_het_samples(int arg, char* argv[], bool digits_choice) {
    std::vector<ind_heter> saver;
    auto fh_map_A = read_reference_fh("../data_input/heter_A.txt");
    auto fh_map_B = read_reference_fh("../data_input/heter_B.txt");
    auto fh_map_C = read_reference_fh("../data_input/heter_C.txt");

    auto samples = parse_file(arg, argv);

    for (const auto& s : samples) {
        for (const auto& [locus_name, allele_calls] : s.loci) {
            std::string het_hom_ind_temp = "NA";
            double fh_value = NAN;
            std::vector<std::string> temp_alleles;
            for (const auto& a : allele_calls) {
                for (size_t k = 0; k < a.alleles.size(); ++k) {
                    if (a.alleles[k].empty()) continue;
                    if (a.alleles[k][0] == 'A') continue;

                    bool match = digits_choice
                        ? (a.alleles[k].find(':') != std::string::npos)
                        : (a.alleles[k].find(':') == std::string::npos);

                    if (a.dosages[k] != 0 && match) {
                        int int_v = std::round(a.dosages[k]);
                        if (int_v == 2) {
                            het_hom_ind_temp = "Hom";
                            temp_alleles.push_back(a.alleles[k]);
                            fh_value = 0;
                        } else {
                            temp_alleles.push_back(a.alleles[k]);
                            het_hom_ind_temp = "Het";
                        }
                        // found = true;
                        // break;
                    }
                }

                // if (found) break;
            }

            if(temp_alleles.size() == 2){
                PairKey key(temp_alleles[0], temp_alleles[1]);
                auto* fh_map = (locus_name == "HLA_A") ? &fh_map_A
                              : (locus_name == "HLA_B") ? &fh_map_B
                              : (locus_name == "HLA_C") ? &fh_map_C
                                  : nullptr;

                    if (fh_map) {
                        auto it = fh_map->find(key);
                        if (it != fh_map->end()) {
                            fh_value = it->second;
                        }
                    }
                }
            std::string alleles_str;
            for (const auto& allele : temp_alleles) {
                alleles_str += allele + " ";
            }
            saver.push_back({s.sampleid, locus_name, het_hom_ind_temp, alleles_str, fh_value});
        }
    }

    return saver;
}

// int main(int arg, char* argv[]){
//     bool digits_choice = true;
//     auto test = map_het_samples(arg, argv, digits_choice);
//     for(auto &ttt :test){
//         std::cout << ttt.sampleid << '\t' << ttt.locus << "\t"  << ttt.het_hom << '\t' << ttt.alleles << '\t' << ttt.het_hom_fh << '\n';
//     }
//     return 0;
// }