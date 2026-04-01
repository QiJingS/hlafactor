#include "header.h"
#include "tapasin.cpp"
#include "heterozygosity.cpp"
#include "expression.cpp"
#include "supertypes.cpp"
#include "input_files_parsers.cpp"
#include "HLA_class_I_heavy_chain.cpp"
#include "hlaqc.cpp"
#include "output_layer.cpp"
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
int main(int argc, char* argv[]) {
    
    parse_file(argc, argv);

    std::string helpmsg = "usage: hlafactor -tapasin | -het | -hetf | -exp | -sup | -lignk \n"
					 "  -tapasin : compute the tapasin depedence score\n"
					 "  -het : compute the general heterozygosity for HLA locus (default 2 resolution)\n"
                     "  -hetf : compute the functional heterozygosity for HLA A/B/C locus\n"
					 "  -exp : compute the HLA (A & C) expression levels on the cell surface \n"
					 "  -sup : compute the HLA supertypes (A/B/class II)\n"
                     "  -evet : included all the features above to the output results which included \n"
					 "  -lignk : suggest the important HLA variants resulted with KIRs\n";
    
    std::unordered_set<std::string> options;
    std::unordered_set<std::string> valid_options = {"-tapasin", "-hetf","-i" , "-o", "-het", "-exp", "-sup", "-lignk", "-af", "-impr2", "-hlarg38", "-hlarg37"};
    if(argc < 2) {
        std::string error_msg = "!!ERROR: No options provided. Use -h or --help for usage information.";
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg);
        return 1;
    }
    if (std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0) {
        std::cout << helpmsg;
        return 0;
    }
    for(int i = 1; i < argc; ++i) {
        if(isNumber(argv[i])) continue;
            if( argv[i][0] == '-' && valid_options.find(argv[i]) == valid_options.end()) {
                std::string error_msg =  "!!ERROR: Invalid option provided: " + std::string(argv[i]);
                std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << '\n' << helpmsg;
                return 1;
            }
        options.insert(argv[i]);
    }
    if(options.find("-i") == options.end() || options.find("-o") == options.end()) {

        std::string error_msg = "!!ERROR: Missing required options: -i (input file) and -o (output folder)";
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg);
        return 1;
    }

    bool do_qc = options.find("-af") != options.end() || options.find("-impr2") != options.end() || options.find("-hlarg38") != options.end() || options.find("-hlarg37") != options.end();

    if(do_qc) {
        if(options.find("-tapasin") != options.end() || options.find("-het") != options.end() || options.find("-exp") != options.end() || options.find("-sup") != options.end() || options.find("-lignk") != options.end()) {
             std::string error_msg = "!!ERROR: QC options (-af, -impr2, -hlarg38, -hlarg37) cannot be combined with analysis options (-tapasin, -het, -exp, -sup, -lignk).";
            std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg);
            return 1;
        }
    }

    bool do_supertypes = options.find("-sup") != options.end();
    bool do_exp = options.find("-exp") != options.end();
    bool do_tap = options.find("-tapasin") != options.end();
    bool do_het = options.find("-het") != options.end();
    bool do_lignk = options.find("-lignk") != options.end();
    bool do_func_het = options.find("-hetf") != options.end();

   if (do_qc) {
        function_hla_qc(argc, argv);   
    } else if (do_supertypes) {
        auto test = supertypes_mapping_data(argc, argv);
        auto output_info = get_output_spec(argc, argv, "supertypes");
        write_supertypes_app_results(test, output_info);
    } else if (do_exp) {
        auto summary = alleles_expression_A_C_summary(argc, argv);
        auto output_info = get_output_spec(argc, argv, "expression");
        write_expression_outputs(summary, output_info);
    } else if (do_tap) {
        auto summary = calculatoer_tapasin_A_B_C_globle(argc, argv);
        auto output_info = get_output_spec(argc, argv, "tapasin");
        write_taps_outputs(summary, output_info);
    } else if (do_het) {
        auto test = map_het_samples(argc, argv, 4);
        // auto output_info = get_output_spec(argc, argv, "het");
        // write_hla_heterozygosity_results(test, output_info);
    } else if (do_lignk) {
        auto test = calculate_v_NK(argc, argv);
        auto output_info = get_output_spec(argc, argv, "lignk");
        write_HLA_variant_NK_results(test, output_info);
    }else if (do_func_het)
    {
        auto test = map_het_samples(argc, argv, 4);
        auto output_info = get_output_spec(argc, argv, "func_het");
        write_hla_functional_heterozygosity_results(test, output_info);
    } else {
        std::cout << helpmsg;
    }
    return 0;
}