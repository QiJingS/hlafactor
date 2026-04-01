#include "header.h"
// #include "input_files_parsers.cpp"
struct sample_variants_NK
{
	std::string sample_id;
	std::vector<std::string> var_name;
	std::vector<std::string> genotypes;
};

std::vector<sample_variants_NK> calculate_v_NK(int argc, char* argv[]) {
    auto samples = parse_file(argc, argv);
    std::vector<sample_variants_NK> def_NK;

    for (const auto& s : samples) {
        std::string sample_id = s.sampleid;
        std::vector<std::string> var_name;
        std::vector<std::string> genotypes;

        for (const auto& [locus_name, allele_calls] : s.loci) {
            for (const auto& call : allele_calls) {
                for (size_t k = 0; k < call.alleles.size(); ++k) {
                    int c_v = std::round(call.dosages[k]);

                    if (call.alleles[k].rfind("AA_B_-21", 0) == 0 && !call.alleles[k].empty() && call.alleles[k].back() == 'T') {
                        var_name.push_back("AA_B_-21");
                        genotypes.push_back(c_v == 1 ? "MT" : c_v == 2 ? "TT" : "MM");
                    }
                    else if (call.alleles[k].rfind("AA_B_82", 0) == 0 && !call.alleles[k].empty() && call.alleles[k].back() == 'L') {
                        var_name.push_back("AA_B_82");
                        genotypes.push_back(c_v == 1 ? "BW6/BW6" : c_v == 2 ? "BW4/BW4" : "BW4/BW6");
                    }
                    else if (call.alleles[k].rfind("AA_C_80", 0) == 0 && !call.alleles[k].empty() && call.alleles[k].back() == 'N') {
                        var_name.push_back("AA_C_80");
                        genotypes.push_back(c_v == 1 ? "C1/C2" : c_v == 2 ? "C1/C1" : "C2/C2");
                    }
                }
            }
        }

        def_NK.push_back({sample_id, var_name, genotypes});
    }

    return def_NK;
}
// int main(int argc, char* argv[]) {
//     auto test = calculate_v_NK(argc, argv);
// 	 for (const auto& rec : test) {
// 		std::cout << "Sample ID: " << rec.sample_id << "\n";
// 		for (size_t i = 0; i < rec.var_name.size(); ++i) {
// 			std::cout << "  Variant: " << rec.var_name[i] << ", Genotype: " << rec.genotypes[i] << "\n";
// 		}
// 	}
//     return 0;
// }