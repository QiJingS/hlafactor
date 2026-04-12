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
    def_NK.reserve(samples.size());

    auto decode_AA_B_minus21 = [](int c_v) -> std::string {
        switch (c_v) {
            case 0: return "MM";
            case 1: return "MT";
            case 2: return "TT";
            default: return "NA";
        }
    };

    auto decode_AA_B_82 = [](int c_v) -> std::string {
        switch (c_v) {
            case 0: return "BW6/BW6";
            case 1: return "BW4/BW6";
            case 2: return "BW4/BW4";
            default: return "NA";
        }
    };

    auto decode_AA_C_80 = [](int c_v) -> std::string {
        switch (c_v) {
            case 0: return "C2/C2";
            case 1: return "C1/C2";
            case 2: return "C1/C1";
            default: return "NA";
        }
    };

    for (const auto& s : samples) {
        std::string sample_id = s.sampleid;
        std::vector<std::string> var_name;
        std::vector<std::string> genotypes;

        for (const auto& [locus_name, allele_calls] : s.loci) {
            for (const auto& call : allele_calls) {
                for (size_t k = 0; k < call.alleles.size(); ++k) {
                    if (call.alleles[k].empty()) continue;

                    int c_v = static_cast<int>(std::round(call.dosages[k]));
                    const std::string& allele = call.alleles[k];

                    if (allele.rfind("AA_B_-21", 0) == 0 && allele.back() == 'T') {
                        var_name.push_back("AA_B_-21");
                        genotypes.push_back(decode_AA_B_minus21(c_v));
                    }
                    else if (allele.rfind("AA_B_82", 0) == 0 && allele.back() == 'L') {
                        var_name.push_back("AA_B_82");
                        genotypes.push_back(decode_AA_B_82(c_v));
                    }
                    else if (allele.rfind("AA_C_80", 0) == 0 && allele.back() == 'N') {
                        var_name.push_back("AA_C_80");
                        genotypes.push_back(decode_AA_C_80(c_v));
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