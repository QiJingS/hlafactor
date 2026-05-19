#include "header.h"
#include "tapasin.cpp"
#include "heterozygosity.cpp"
#include "expression.cpp"
#include "supertypes.cpp"
#include "input_files_parsers.cpp"
#include "HLA_class_I_heavy_chain.cpp"
#include "hlaqc.cpp"
#include "hla_amiaci.cpp"
#include "output_layer.cpp"
#include <filesystem>
namespace fs = std::filesystem;

bool isNumber(const std::string& s) {
    try {
        size_t pos;
        std::stod(s, &pos);
        return pos == s.size();
    } catch (...) {
        return false;
    }
}

std::string num_decro_stars(const std::string& msg) {
    return std::string(1 + msg.size(), '-');
}

std::string getArgValue(int argc, char* argv[], const std::string& flag) {
    for (int i = 1; i < argc - 1; ++i) {
        if (std::string(argv[i]) == flag) {
            return std::string(argv[i + 1]);
        }
    }
    return "";
}

bool pathExists(const std::string& path) {
    return !path.empty() && fs::exists(fs::u8path(path));
}

int main(int argc, char* argv[]) {
    
    std::string helpmsg = 
    "HLAfactor: a tool for HLA QC and analysis\n"
    "\n"
    "usage: hlafactor [options]\n"
    "  -h, --help        show this help message\n"
    "\n"
    "1. HLA QC:\n"
    "   hlafactor -i <input_file> -o <output_folder> -hlarg38|-hlarg19 [-af <min> [<max>]] [-impr2 <min> [<max>]]\n"
    "   -hlarg38        use hg38 (GRCh38) coordinates for QC\n"
    "   -hlarg19        use hg19 (GRCh19) coordinates for QC\n"
    "   -af              filter variants by allele frequency in the HLA region\n"
    "   -impr2           filter variants by imputation r2 in the HLA region\n"
    "\n"
    "2. HLA analysis:\n"
    "   hlafactor -i <input_file> -o <output_folder> [-tapasin | -het | -hetf | -exp | -sup | -lignk | -amac]\n"
    "   -tapasin         compute the tapasin dependence score\n"
    "   -het             compute general heterozygosity for HLA locus\n"
    "   -hetf            compute functional heterozygosity for HLA A/B/C locus\n"
    "   -exp             compute HLA A & C protein expression levels on the cell surface\n"
    "   -sup             compute HLA supertypes (A/B/class II)\n"
    "   -lignk           compute functional HLA markers linked to KIRs: HLA-B -21T/M, HLA-C C1/C2, HLA-B BW4/BW6\n"
    "   -amac            map HLA amino acid genotypes using IMGT alignments\n"
    "\n"
    "QC and analysis options cannot be combined. \n"
    "Supported input formats: .txt, .dosage, .vcf/.vcf.gz \n"
    "Output will be saved to the specified output folder, or current directory if not specified. \n"
    "Output files will be named based on the input file name with suffixes indicating the analysis type.\n";

    
    std::unordered_set<std::string> options;
    std::unordered_set<std::string> valid_options = {
        "-tapasin", "-hetf", "-i", "-o", "-het", "-exp", "-sup", "-lignk", "-af", "-amac", "-impr2", "-hlarg38", "-hlarg19"
    };

    if (argc < 2) {
        std::cerr << "!!ERROR: No options provided. Use -h or -help for usage information.\n";
        return 1;
    }

    if (std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "-help") == 0) {
        std::cout << helpmsg;
        return 0;
    }

    

    for (int i = 1; i < argc; ++i) {
        if (isNumber(argv[i])) continue;
        if (argv[i][0] == '-' && valid_options.find(argv[i]) == valid_options.end()) {
            std::string error_msg = "!!ERROR: Invalid option provided: " + std::string(argv[i]);
            std::cerr << num_decro_stars(error_msg) << '\n'
                      << error_msg << '\n'
                      << num_decro_stars(error_msg) << '\n'
                      << helpmsg;
            return 1;
        }
        options.insert(argv[i]);
    }
    if (options.find("-i") == options.end() || options.find("-o") == options.end()) {
        std::cerr << "!!ERROR: Missing required options: -i (input file) and -o (output folder) must both be provided.\n";
        return 1;
    }
    if (
        options.size() == 2 &&
        options.find("-i") != options.end() &&
        options.find("-o") != options.end()
    ) {
        std::cerr << "!!ERROR: At least one additional option must be provided besides -i and -o.\n" << helpmsg;
        return 1;
    }

    std::string input_path = getArgValue(argc, argv, "-i");
    std::string output_path = getArgValue(argc, argv, "-o");

    if (input_path.empty()) {
        std::cerr << "!!ERROR: -i must be followed by an input file path.\n";
        return 1;
    }
    if (output_path.empty()) {
        std::cerr << "!!ERROR: -o must be followed by an output folder path.\n";
        return 1;
    }

    if (!pathExists(input_path)) {
        std::cerr << "!!ERROR: Input file does not exist: " << input_path << '\n';
        return 1;
    }
    std::filesystem::path p(output_path);
    if (std::filesystem::exists(p) && std::filesystem::is_directory(p)) {
        
    }
    else {
        std::filesystem::path parent = p.parent_path();
        if (parent.empty()) {
            parent = ".";
        }

        if (!std::filesystem::exists(parent) || !std::filesystem::is_directory(parent)) {
            std::cerr << "!!ERROR: Output directory does not exist: "
                    << parent.string() << '\n';
            return 1;
        }
    }
    parse_file(argc, argv);
    bool do_qc = options.find("-af") != options.end() ||
                 options.find("-impr2") != options.end() ||
                 options.find("-hlarg38") != options.end() ||
                 options.find("-hlarg19") != options.end();

    if (do_qc) {
        if (options.find("-tapasin") != options.end() ||
            options.find("-het") != options.end() ||
            options.find("-exp") != options.end() ||
            options.find("-sup") != options.end() ||
            options.find("-lignk") != options.end()) {
            std::string error_msg =
                "!!ERROR: QC options (-af, -impr2, -hlarg38, -hlarg37, -amac) cannot be combined with analysis options (-tapasin, -het, -hetf, -exp, -sup, -lignk).";
            std::cerr << error_msg << '\n' << '\n';
            return 1;
        }
    }

    bool do_supertypes = options.find("-sup") != options.end();
    bool do_exp = options.find("-exp") != options.end();
    bool do_aa = options.find("-amac") != options.end();
    bool do_tap = options.find("-tapasin") != options.end();
    bool do_het = options.find("-het") != options.end();
    bool do_lignk = options.find("-lignk") != options.end();
    bool do_func_het = options.find("-hetf") != options.end();

    if (do_qc) {
        std::cout << "[INFO] Running QC pipeline...\n";
        function_hla_qc(argc, argv);
    }else if(do_aa) {
        std::cout << "\nThanks for using the HLAfactor :)\n";
        std::cout << "*********************************\n";
        std::cout << "[INFO] Program started.\n";

        std::cout << "[INFO] Running IMGT amino acid mapping...\n";
        auto test = map_imgt_amino_acids(argc, argv);
        auto output_info = get_output_spec(argc, argv, "amac");
        write_amino_acid_mapping_results(test, output_info);
    } else if (do_supertypes) {
        std::cout << "\nThanks for using the HLAfactor :)\n";
        std::cout << "*********************************\n";
        std::cout << "[INFO] Program started.\n";

        std::cout << "[INFO] Running supertypes mapping...\n";
        auto test = supertypes_mapping_data(argc, argv);
        auto output_info = get_output_spec(argc, argv, "supertypes");
        write_supertypes_app_results(test, output_info);
    } else if (do_exp) {

        std::cout << "\nThanks for using the HLAfactor :)\n";
        std::cout << "*********************************\n";
        std::cout << "[INFO] Program started.\n";

        std::cout << "[INFO] Running expression analysis...\n";
        auto summary = alleles_expression_A_C_summary(argc, argv);
        auto output_info = get_output_spec(argc, argv, "expression");
        write_expression_outputs(summary, output_info);
    } else if (do_tap) {
        std::cout << "\nThanks for using the HLAfactor :)\n";
        std::cout << "*********************************\n";
        std::cout << "[INFO] Program started.\n";

        std::cout << "[INFO] Running tapasin analysis...\n";
        auto summary = calculatoer_tapasin_A_B_C_globle(argc, argv);
        auto output_info = get_output_spec(argc, argv, "tapasin");
        write_taps_outputs(summary, output_info);
    } else if (do_het) {
        std::cout << "\nThanks for using the HLAfactor :)\n";
        std::cout << "*********************************\n";
        std::cout << "[INFO] Program started.\n";

        std::cout << "[INFO] Running HLA heterozygosity analysis...\n";
        auto test = map_het_samples(argc, argv, 4);
        auto output_info = get_output_spec(argc, argv, "het");
        write_hla_heterozygosity_results(test, output_info);
    } else if (do_lignk) {
        std::cout << "\nThanks for using the HLAfactor :)\n";
        std::cout << "*********************************\n";
        std::cout << "[INFO] Program started.\n";

        std::cout << "[INFO] Running functional HLA markers linked to KIRs analysis...\n";
        auto test = calculate_v_NK(argc, argv);
        auto output_info = get_output_spec(argc, argv, "lignk");
        write_HLA_variant_NK_results(test, output_info);
    } else if (do_func_het) {
        std::cout << "\nThanks for using the HLAfactor :)\n";
        std::cout << "*********************************\n";
        std::cout << "[INFO] Program started.\n";

        std::cout << "[INFO] Running functional heterozygosity analysis...\n";
        auto test = map_het_samples(argc, argv, 4);
        auto output_info = get_output_spec(argc, argv, "func_het");
        write_hla_functional_heterozygosity_results(test, output_info);
    } else {
        std::cout << helpmsg;
    }

    return 0;
}
