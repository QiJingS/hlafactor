#include "header.h"

namespace fs = std::filesystem;
bool in_HLA_region(int pos, int hg) {
    int pos1 = pos + 1;
    if (hg == 38) {
        return (pos1 >= 25726063 && pos1 <= 33400644);
    } else if (hg == 19) {
        return (pos1 >= 25726291 && pos1 <= 33368421);
    }
    return false;
}

bool is_chr6(const std::string& chr) {
    return (chr == "6" || chr == "chr6");
}

bool id_starts_with_HLA_or_AA(bcf1_t* rec) {
    bcf_unpack(rec, BCF_UN_STR);
    if (!rec->d.id) return false;

    std::string id = rec->d.id;

    return (id.rfind("HLA", 0) == 0) ||
           (id.rfind("AA_B_82", 0) == 0) ||
           (id.rfind("AA_B_-21", 0) == 0) ||
           (id.rfind("AA_C_80", 0) == 0);
}
static std::string resolve_output_path(const std::string& output) {
    fs::path p(output);
    if (fs::exists(p) && fs::is_directory(p)) {
        return (p / "filtered.vcf.gz").string();
    }
    if (!output.empty() && output.back() == '/') {
        return (p / "filtered.vcf.gz").string();
    }
    return output;
}
static std::string choose_hts_mode(const std::string& output) {
    if (output.size() >= 4 && output.rfind(".bcf") == output.size() - 4) {
        return "wb";
    }
    if (output.size() >= 3 && output.rfind(".gz") == output.size() - 3) {
        return "wz";
    }
    return "w";
}
void filter_HLA_region(int hg, const std::string& input, const std::string& output,
                       float af_min, float af_max,
                       bool use_af_min, bool use_af_max,
                       float r2_min, float r2_max,
                       bool use_r2_min, bool use_r2_max) {
    std::string outpath = resolve_output_path(output);
    std::string mode = choose_hts_mode(outpath);
    fs::path outp(outpath);
    if (outp.has_parent_path()) {
        std::error_code ec;
        fs::create_directories(outp.parent_path(), ec);
        if (ec) {
			std::string error_msg = "!!ERROR: cannot create output directory: " + outp.parent_path().string() + ". " + ec.message();
			std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
            return;
        }
    }
    htsFile* in = bcf_open(input.c_str(), "r");
    if (!in) {
		std::string error_msg = "!!ERROR: cannot open input file: " + input;
		std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        return;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(in);
    if (!hdr) {
		std::string error_msg = "!!ERROR: cannot read header from input file: " + input;
		std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        bcf_close(in);
        return;
    }
    int n_samples = bcf_hdr_nsamples(hdr);
    htsFile* out = bcf_open(outpath.c_str(), mode.c_str());
    if (!out) {
        std::string error_msg = "!!ERROR: cannot open output: " + outpath;
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        bcf_hdr_destroy(hdr);
        bcf_close(in);
        return;
    }
    if (bcf_hdr_write(out, hdr) != 0) {
        std::string error_msg = "!!ERROR: cannot write header\n";
        std::cerr << num_decro_stars(error_msg) << '\n' << error_msg << '\n' << num_decro_stars(error_msg) << "\n";
        bcf_hdr_destroy(hdr);
        bcf_close(in);
        bcf_close(out);
        return;
    }
    bcf1_t* rec = bcf_init();
    float* af_arr = NULL;
    int naf = 0;
    float* r2_arr = NULL;
    int nr2 = 0;

    int total_variants = 0;
    int kept_variants = 0;

    while (bcf_read(in, hdr, rec) == 0) {
        total_variants++;

        int rid = rec->rid;
        if (rid < 0) continue;

        const char* chrom = bcf_hdr_id2name(hdr, rid);
        int pos = rec->pos;

        if (!(chrom && is_chr6(chrom) && in_HLA_region(pos, hg))) continue;
        if (!id_starts_with_HLA_or_AA(rec)) continue;

        if (use_af_min || use_af_max) {
            if (bcf_get_info_float(hdr, rec, "AF", &af_arr, &naf) <= 0) continue;
            float af = af_arr[0];
            if (use_af_min && af < af_min) continue;
            if (use_af_max && af > af_max) continue;
        }

        if (use_r2_min || use_r2_max) {
            if (bcf_get_info_float(hdr, rec, "R2", &r2_arr, &nr2) <= 0) continue;
            float r2 = r2_arr[0];
            if (use_r2_min && r2 < r2_min) continue;
            if (use_r2_max && r2 > r2_max) continue;
        }

        if (bcf_write(out, hdr, rec) != 0) {
            std::cerr << "Error writing record\n";
            break;
        }
        kept_variants++;
    }

    free(af_arr);
    free(r2_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(in);
    bcf_close(out);

    std::cout << "\nThanks for using the HLAfactor :)\n";
    std::cout << "*********************************\n";
    std::cout << "Samples (individuals): " << n_samples << "\n";
    std::cout << "Total variants read:  " << total_variants << "\n";
    std::cout << "Variants kept:        " << kept_variants << "\n";
    std::cout << "Variants removed:     " << (total_variants - kept_variants) << "\n";
    std::cout << "*********************************\n";
    std::cout << "\nFiltering complete, please find the output file at: " << outpath << "\n";
}

void function_hla_qc(int argc, char* argv[]) {
    std::string input;
    std::string output;  
    bool has_output = false;
    int hg = 0;

    float af_min = 0, af_max = 0;
    bool use_af_min = false, use_af_max = false;

    float r2_min = 0, r2_max = 0;
    bool use_r2_min = false, use_r2_max = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-o" && i + 1 < argc) {
            output = argv[++i];
            has_output = true;
        }
        else if (arg == "-hlarg38") {
            hg = 38;
        }
        else if (arg == "-hlarg37") {
            hg = 19;
        }
        else if (arg == "-i" && i + 1 < argc) {
            input = argv[++i];
        }
        else if (arg == "-af" && i + 1 < argc) {
            af_min = std::stof(argv[++i]);
            use_af_min = true;

            if (i + 1 < argc && argv[i + 1][0] != '-') {
                af_max = std::stof(argv[++i]);
                use_af_max = true;
            }
        }
        else if (arg == "-impr2" && i + 1 < argc) {
            r2_min = std::stof(argv[++i]);
            use_r2_min = true;

            if (i + 1 < argc && argv[i + 1][0] != '-') {
                r2_max = std::stof(argv[++i]);
                use_r2_max = true;
            }
        }
    }

    if (input.empty()) {
        std::string error_msg = "!!Error: no input file (-i)";
        std::cerr << error_msg << '\n';
        return;
    }

    if (!has_output || output.empty()) {
        output = "filtered.vcf";
    }

    if (hg == 0) {
        std::string error_msg = "!!Error: specify -hlarg37 or -hlarg38";
        std::cerr <<  error_msg << '\n';
        return;
    }

    filter_HLA_region(hg, input, output,
                      af_min, af_max, use_af_min, use_af_max,
                      r2_min, r2_max, use_r2_min, use_r2_max);
}
// int main(int argc, char* argv[]) {
//     if (argc < 2) {
//         std::cout << "Usage:\n";
//         std::cout << "  ./hla_filter -hlarg38 -i input.vcf.gz -o out.vcf.gz\n";
//         std::cout << "Options:\n";
//         std::cout << "  -af x [y]      (AF >= x OR x <= AF <= y)\n";
//         std::cout << "  -impr2 x [y]   (R2 >= x OR x <= R2 <= y)\n";
//         return 0;
//     }
//
//     function_hla_qc(argc, argv);
//     return 0;
// }