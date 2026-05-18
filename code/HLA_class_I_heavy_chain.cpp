#include "header.h"

struct sample_variants_NK
{
    std::string sample_id;
    std::vector<std::string> var_name;
    std::vector<std::string> genotypes;
};

namespace {

std::string normalize_nk_table_allele(const std::string& allele)
{
    if (allele.rfind("HLA_", 0) == 0) {
        return allele;
    }

    auto star_pos = allele.find('*');
    if (star_pos == std::string::npos) {
        return allele;
    }

    return "HLA_" + allele;
}

std::string resolve_nk_lookup_key(
    const std::string& allele,
    const std::unordered_map<std::string, std::pair<std::string, std::string>>& nk_lookup)
{
    auto exact = nk_lookup.find(allele);
    if (exact != nk_lookup.end()) {
        return exact->first;
    }

    size_t best_len = 0;
    std::string best_key;

    for (const auto& [key, value] : nk_lookup) {
        (void)value;
        if (allele.size() >= key.size() &&
            allele.compare(0, key.size(), key) == 0 &&
            key.size() > best_len) {
            best_len = key.size();
            best_key = key;
        }
    }

    return best_key;
}

}

std::vector<sample_variants_NK> calculate_v_NK(int argc, char* argv[])
{
    auto samples = parse_file(argc, argv);

    std::vector<sample_variants_NK> def_NK;
    def_NK.reserve(samples.size());

    std::unordered_map<std::string,
        std::pair<std::string, std::string>> nk_lookup;

    std::ifstream fin(resource_path("nk_table.txt"));

    if (!fin)
    {
        throw std::runtime_error(
            "Cannot open nk_table.txt");
    }

    std::string line;

    std::getline(fin, line);

    while (std::getline(fin, line))
    {
        if (line.empty()) continue;

        std::stringstream ss(line);

        std::string allele;
        std::string category;
        std::string ligand;

        std::getline(ss, allele, '\t');
        std::getline(ss, category, '\t');
        std::getline(ss, ligand, '\t');

        nk_lookup[normalize_nk_table_allele(allele)] = {category, ligand};
    }

    for (const auto& s : samples)
    {
        std::string sample_id = s.sampleid;

        std::vector<std::string> var_name;
        std::vector<std::string> genotypes;

        std::vector<std::string> c_ligands;
        std::vector<std::string> b_ligands;

        std::unordered_map<std::string, int> normalized_dosage;

        for (const auto& [locus_name, allele_calls] : s.loci)
        {
            for (const auto& call : allele_calls)
            {
                for (size_t k = 0; k < call.alleles.size(); ++k)
                {
                    if (call.alleles[k].empty()) continue;

                    int dosage =
                        static_cast<int>(
                            std::round(call.dosages[k]));

                    if (dosage == 0)
                    {
                        continue;
                    }

                    std::string allele = call.alleles[k];

                    auto current = normalized_dosage.find(allele);
                    if (current == normalized_dosage.end() || dosage > current->second)
                    {
                        normalized_dosage[allele] = dosage;
                    }
                }
            }
        }

        std::unordered_map<std::string, int> matched_dosage;

        for (const auto& [allele, dosage] : normalized_dosage)
        {
            std::string lookup_key = resolve_nk_lookup_key(allele, nk_lookup);
            if (lookup_key.empty())
            {
                continue;
            }

            auto current = matched_dosage.find(lookup_key);
            if (current == matched_dosage.end() || dosage > current->second)
            {
                matched_dosage[lookup_key] = dosage;
            }
        }

        for (const auto& [lookup_key, dosage] : matched_dosage)
        {
            auto it = nk_lookup.find(lookup_key);

            const std::string& category =
                it->second.first;

            const std::string& ligand =
                it->second.second;

            for (int d = 0; d < dosage; ++d)
            {
                if (category == "HLA_C")
                {
                    c_ligands.push_back(ligand);
                }
                else if (category == "HLA_B")
                {
                    b_ligands.push_back(ligand);
                }
            }
        }

        var_name.push_back("HLA_C_ligand");

        if (c_ligands.size() == 1)
        {
            genotypes.push_back(c_ligands[0]);
        }
        else if (c_ligands.size() == 2)
        {
            std::sort(
                c_ligands.begin(),
                c_ligands.end());

            genotypes.push_back(
                c_ligands[0] + "/" + c_ligands[1]);
        }
        else
        {
            genotypes.push_back("NA");
        }

        var_name.push_back("HLA_B_ligand");

        if (b_ligands.size() == 1)
        {
            genotypes.push_back(b_ligands[0]);
        }
        else if (b_ligands.size() == 2)
        {
            std::sort(
                b_ligands.begin(),
                b_ligands.end());

            genotypes.push_back(
                b_ligands[0] + "/" + b_ligands[1]);
        }
        else
        {
            genotypes.push_back("NA");
        }

        def_NK.push_back({
            sample_id,
            var_name,
            genotypes
        });
    }

    return def_NK;
}
