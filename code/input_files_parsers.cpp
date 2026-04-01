#include "header.h"
std::vector<samplerecord> parse_txt_file(const std::string& filename) {
    std::ifstream txt_file(filename);
    if (!txt_file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;
    if (!std::getline(txt_file, line)) {
        throw std::runtime_error("Empty file or failed to read header: " + filename);
    }

    std::istringstream header_stream(line);
    std::string ID, second_token;
    header_stream >> ID; 

    std::vector<std::string> locus_names;
    std::vector<std::string> alleles_names;
    while (header_stream >> second_token) {
        alleles_names.push_back(second_token);
        std::string locus_name_i = second_token.substr(0, second_token.find('*'));
        locus_names.push_back(locus_name_i);
    }

    std::vector<samplerecord> result;
    while (std::getline(txt_file, line)) {
        if (line.empty()) continue;
        std::istringstream line_stream(line);
        std::string sample_id;
        line_stream >> sample_id;

        samplerecord rec;
        rec.sampleid = sample_id;

        int idx_locus = 0;
        std::string token;
        while (line_stream >> token && idx_locus < static_cast<int>(locus_names.size())) {
            double dosage;
            if (token == "NA" || token == "na" || token == "Na" || token == "nA") {
                dosage = std::numeric_limits<double>::quiet_NaN();
            } else {
                try {
                    dosage = std::stod(token);
                } catch (...) {
                    dosage = std::numeric_limits<double>::quiet_NaN();
                }
            }

            allelescall call;
            call.alleles.push_back(alleles_names[idx_locus]);
            call.dosages.push_back(dosage);
            rec.loci[locus_names[idx_locus]].push_back(std::move(call));

            ++idx_locus;
        }

        result.push_back(std::move(rec));
    }

    return result;
}

std::vector<samplerecord> parse_dosage_file(const std::string& filename) {
    std::ifstream dosage_file(filename);
    if (!dosage_file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;
    if (!std::getline(dosage_file, line)) {
        throw std::runtime_error("Empty file or failed to read header: " + filename);
    }

    std::istringstream header_stream(line);
    std::string chr, snpid, rsid, pos, ref, alt, sample_name;
    header_stream >> chr >> snpid >> rsid >> pos >> ref >> alt;

    std::vector<std::string> sample_names_vec;
    while (header_stream >> sample_name) {
        sample_names_vec.push_back(sample_name);
    }

    std::vector<samplerecord> samples;
    samples.reserve(sample_names_vec.size());
    for (const auto& s : sample_names_vec) {
        samples.push_back(samplerecord{s});
    }

    while (std::getline(dosage_file, line)) {
        if (line.empty()) continue;
        std::istringstream line_stream(line);
        std::string chromosome, SNPID, rsid2, position, alleleA, alleleB;
        line_stream >> chromosome >> SNPID >> rsid2 >> position >> alleleA >> alleleB;

        std::vector<double> dosages;
        double d;
        while (line_stream >> d) dosages.push_back(d);

        std::string locus = SNPID.substr(0, SNPID.find('*'));

        for (size_t i = 0; i < samples.size(); ++i) {
            allelescall call;
            call.alleles.push_back(SNPID);
            call.dosages.push_back(dosages[i]);
            samples[i].loci[locus].push_back(std::move(call));
        }
    }

    return samples;
}
void print_samples(const std::vector<samplerecord>& samples) {
    for (const auto& sample : samples) {
        std::cout << "Sample ID: " << sample.sampleid << '\n';
        for (const auto& [locus_name, allele_calls] : sample.loci) {
            for (const auto& call : allele_calls) {
                for (size_t k = 0; k < call.alleles.size(); ++k) {
                    double dosage = call.dosages[k];
                    std::cout << "  Locus: " << locus_name
                              << ", Allele: " << call.alleles[k]
                              << ", Dosage: " << dosage << '\n';
                }
            }
        }
    }
}



std::vector<samplerecord> parse_vcf_htslib(const std::string& filename) {
    htsFile* fp = bcf_open(filename.c_str(), "r");
    if (!fp) {
        throw std::runtime_error("Failed to open VCF: " + filename);
    }

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        bcf_close(fp);
        throw std::runtime_error("Failed to read VCF header");
    }

    int n_samples = bcf_hdr_nsamples(hdr);

    std::vector<samplerecord> samples;
    samples.reserve(n_samples);

    for (int i = 0; i < n_samples; ++i) {
        samples.emplace_back(samplerecord{hdr->samples[i]});
    }

    bcf1_t* rec = bcf_init();

    float* ds = nullptr;
    int nds = 0;

    int32_t* gt = nullptr;
    int ngt = 0;

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        std::string chrom = bcf_hdr_id2name(hdr, rec->rid);
        std::string pos = std::to_string(rec->pos + 1);
        std::string id = rec->d.id ? rec->d.id : ".";
        std::string locus = (id != "." ? id : chrom + ":" + pos);

        size_t star_pos = locus.find('*');
        if (star_pos != std::string::npos) {
            locus = locus.substr(0, star_pos);
        }

        int nret = bcf_get_format_float(hdr, rec, "DS", &ds, &nds);

        bool use_ds = (nret > 0);
        if (!use_ds) {
            nret = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
        }
        for (int i = 0; i < n_samples; ++i) {
            double dosage = 0.0;
            if (use_ds) {
                dosage = ds[i];
            } else if (nret > 0) {
                int32_t* ptr = &gt[i * 2];  
                int a1 = bcf_gt_allele(ptr[0]);
                int a2 = bcf_gt_allele(ptr[1]);
                dosage = (a1 == 1) + (a2 == 1);
            }
            allelescall call;
            call.alleles.push_back(id);
            call.dosages.push_back(dosage);
            samples[i].loci[locus].push_back(std::move(call));
        }
    }

    free(ds);
    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    return samples;
}

std::vector<samplerecord> parse_file(int argc, char* argv[]) {
    std::vector<samplerecord> samples;
    for (int i = 1; i < argc; ++i) {
        std::string input_s = argv[i];

        if (input_s == "-i") {
            std::string input_adress_names = argv[++i];
            std::filesystem::path filepath(input_adress_names);

            std::string ext = filepath.extension().string();
            std::cout << input_adress_names << std::endl;
            if (ext == ".txt") {
                samples = parse_txt_file(input_adress_names);
                std::cout << "Parsed txt samples: " << samples.size() << " rows\n";

            } else if (ext == ".dosage") {
                samples = parse_dosage_file(input_adress_names);
                // print_samples(samples);
                std::cout << "Parsed dosage samples: " << samples.size() << " columns\n";

            } else if (ext == ".vcf" || input_adress_names.ends_with(".vcf.gz")) {
                samples = parse_vcf_htslib(input_adress_names);
                // print_samples(samples);
                std::cout << "Parsed VCF samples: " << samples.size() << " columns\n";
            }

            break;
        }
    }
    return samples;
};
// int main(int argc, char* argv[]) {
//     auto test = parse_file(argc, argv);
//     print_samples(test);
//     return 0;
// }
// int main(){
//     // htsFile* fp = bcf_open("/Users/qijingshen/Desktop/Thesis and presentation/code/app_snps_HLA_4digits.vcf", "r");
//     auto test = parse_vcf_htslib("/Users/qijingshen/Desktop/Thesis and presentation/code/app_snps_HLA_4digits.vcf.gz");
//     print_samples(test);
//     // bcf_hdr_t* hdr = bcf_hdr_read(fp);
//     // bcf1_t* rec = bcf_init();
//     // while (bcf_read(fp, hdr, rec) == 0) {
//     //     bcf_unpack(rec, BCF_UN_STR);

//     //     std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t"
//     //               << rec->pos + 1 << "\t"
//     //               << rec->d.allele[0] << " -> "
//     //               << rec->d.allele[1] << "\n";
//     // }
//     // bcf_destroy(rec);
//     // bcf_hdr_destroy(hdr);
//     // bcf_close(fp);
//     return 0;
// }

