#include "header.h"

namespace {

bool is_missing_token(const std::string& token) {
    return token.empty() ||
           token == "."   ||
           token == "NA"  ||
           token == "na"  ||
           token == "Na"  ||
           token == "nA"  ||
           token == "NaN" ||
           token == "nan" ||
           token == "NAN";
}

double parse_dosage_token(const std::string& token) {
    if (is_missing_token(token)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    try {
        size_t pos = 0;
        double v = std::stod(token, &pos);

        if (pos != token.size() || !std::isfinite(v)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return v;
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

std::string locus_from_allele_name(const std::string& allele) {
    auto pos = allele.find('*');
    if (pos == std::string::npos) {
        return allele;
    }
    return allele.substr(0, pos);
}

bool starts_with_case_sensitive(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() &&
           std::equal(prefix.begin(), prefix.end(), s.begin());
}

bool is_hla_related_id(const std::string& id) {
    return starts_with_case_sensitive(id, "HLA") ||
           starts_with_case_sensitive(id, "AA");
}

} 

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

    std::vector<std::string> kept_locus_names;
    std::vector<std::string> kept_alleles_names;
    std::vector<bool> keep_flags;

    while (header_stream >> second_token) {
        bool keep = is_hla_related_id(second_token);
        keep_flags.push_back(keep);
        if (keep) {
            kept_alleles_names.push_back(second_token);
            kept_locus_names.push_back(locus_from_allele_name(second_token));
        }
    }

    std::vector<samplerecord> result;

    while (std::getline(txt_file, line)) {
        if (line.empty()) continue;

        std::istringstream line_stream(line);
        std::string sample_id;
        if (!(line_stream >> sample_id)) {
            continue;
        }

        samplerecord rec;
        rec.sampleid = sample_id;

        size_t idx_kept = 0;
        for (size_t idx_all = 0; idx_all < keep_flags.size(); ++idx_all) {
            std::string token;
            if (!(line_stream >> token)) {
                token = "NA";
            }

            if (!keep_flags[idx_all]) {
                continue;
            }

            double dosage = parse_dosage_token(token);

            allelescall call;
            call.alleles.push_back(kept_alleles_names[idx_kept]);
            call.dosages.push_back(dosage);
            rec.loci[kept_locus_names[idx_kept]].push_back(std::move(call));

            ++idx_kept;
        }

        result.push_back(std::move(rec));
    }

    return result;
}

// ===========================
// dosage parser
// ===========================
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
    if (!(header_stream >> chr >> snpid >> rsid >> pos >> ref >> alt)) {
        throw std::runtime_error("Malformed dosage header: " + filename);
    }

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
        if (!(line_stream >> chromosome >> SNPID >> rsid2 >> position >> alleleA >> alleleB)) {
            continue;
        }

        if (!is_hla_related_id(SNPID)) {
            continue;
        }

        std::string locus = locus_from_allele_name(SNPID);

        for (size_t i = 0; i < samples.size(); ++i) {
            std::string token;
            if (!(line_stream >> token)) {
                token = "NA";
            }

            double dosage = parse_dosage_token(token);

            allelescall call;
            call.alleles.push_back(SNPID);
            call.dosages.push_back(dosage);
            samples[i].loci[locus].push_back(std::move(call));
        }
    }

    return samples;
}

// ===========================
// VCF / BCF parser
// ===========================
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

        if (!is_hla_related_id(id)) {
            continue;
        }

        int nret_ds = bcf_get_format_float(hdr, rec, "DS", &ds, &nds);
        bool use_ds = (nret_ds > 0);

        int nret_gt = -1;
        if (!use_ds) {
            nret_gt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
        }

        for (int i = 0; i < n_samples; ++i) {
            double dosage = std::numeric_limits<double>::quiet_NaN();

            if (use_ds) {
                if (i < nret_ds) {
                    float v = ds[i];
                    if (v != bcf_float_missing &&
                        v != bcf_float_vector_end &&
                        std::isfinite(v)) {
                        dosage = static_cast<double>(v);
                    }
                }
            } else if (nret_gt > 0) {
                int base = i * 2;
                if (base + 1 < nret_gt) {
                    int32_t g1 = gt[base];
                    int32_t g2 = gt[base + 1];

                    if (!bcf_gt_is_missing(g1) && !bcf_gt_is_missing(g2)) {
                        int a1 = bcf_gt_allele(g1);
                        int a2 = bcf_gt_allele(g2);

                        if (a1 >= 0 && a2 >= 0) {
                            dosage = static_cast<double>((a1 == 1) + (a2 == 1));
                        }
                    }
                }
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

// ===========================
// unified parser entry
// ===========================
std::vector<samplerecord> parse_file(int argc, char* argv[]) {
    std::vector<samplerecord> samples;

    for (int i = 1; i < argc; ++i) {
        std::string input_s = argv[i];

        if (input_s == "-i" && i + 1 < argc) {
            std::string input_adress_names = argv[++i];
            std::filesystem::path filepath(input_adress_names);

            std::string ext = filepath.extension().string();
            std::cout << input_adress_names << std::endl;

            if (ext == ".txt") {
                samples = parse_txt_file(input_adress_names);
                std::cout << "Parsed txt samples: " << samples.size() << " rows\n";
            } else if (ext == ".dosage") {
                samples = parse_dosage_file(input_adress_names);
                std::cout << "Parsed dosage samples: " << samples.size() << " columns\n";
            } else if (ext == ".vcf" || input_adress_names.ends_with(".vcf.gz")) {
                samples = parse_vcf_htslib(input_adress_names);
                std::cout << "Parsed VCF samples: " << samples.size() << " columns\n";
            } else {
                throw std::runtime_error("Unsupported input format: " + input_adress_names);
            }

            break;
        }
    }

    return samples;
}