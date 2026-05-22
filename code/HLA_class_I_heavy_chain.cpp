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

std::string nk_category_from_allele(const std::string& allele)
{
    auto star_pos = allele.find('*');
    if (star_pos == std::string::npos) {
        return "";
    }

    std::string locus = allele.substr(0, star_pos);
    if (locus == "HLA_A" || locus == "HLA_B" || locus == "HLA_C") {
        return locus;
    }

    return "";
}

std::string normalize_nk_ligand(const std::string& ligand)
{
    std::string normalized = ligand;

    for (size_t pos = normalized.find("Bw");
         pos != std::string::npos;
         pos = normalized.find("Bw", pos + 2)) {
        normalized.replace(pos, 2, "BW");
    }

    return normalized;
}

bool allele_prefix_match(const std::string& lhs, const std::string& rhs)
{
    const std::string& shorter =
        lhs.size() <= rhs.size() ? lhs : rhs;
    const std::string& longer =
        lhs.size() <= rhs.size() ? rhs : lhs;

    return longer.size() > shorter.size() &&
           longer.compare(0, shorter.size(), shorter) == 0 &&
           longer[shorter.size()] == ':';
}

void add_ligand_call(
    std::vector<std::string>& ligands,
    const std::string& ligand,
    int dosage)
{
    for (int d = 0; d < dosage; ++d) {
        ligands.push_back(ligand);
    }
}

std::string format_ligands(std::vector<std::string> ligands)
{
    if (ligands.empty()) {
        return "NA";
    }

    std::sort(
        ligands.begin(),
        ligands.end(),
        [](const std::string& lhs, const std::string& rhs) {
            if (lhs == ".") return false;
            if (rhs == ".") return true;
            return lhs < rhs;
        });

    std::string out;
    for (const auto& ligand : ligands) {
        if (!out.empty()) {
            out += "/";
        }
        out += ligand;
    }

    return out;
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

        nk_lookup[normalize_nk_table_allele(allele)] = {
            category,
            normalize_nk_ligand(ligand)
        };
    }

    for (const auto& s : samples)
    {
        std::string sample_id = s.sampleid;

        std::vector<std::string> var_name;
        std::vector<std::string> genotypes;

        std::vector<std::string> a_ligands;
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
        std::unordered_map<std::string, std::pair<std::string, int>> unmatched_dosage;

        std::vector<std::pair<std::string, int>> allele_dosages(
            normalized_dosage.begin(),
            normalized_dosage.end());
        std::sort(
            allele_dosages.begin(),
            allele_dosages.end(),
            [](const auto& lhs, const auto& rhs) {
                if (lhs.first.size() != rhs.first.size()) {
                    return lhs.first.size() < rhs.first.size();
                }
                return lhs.first < rhs.first;
            });

        for (const auto& [allele, dosage] : allele_dosages)
        {
            std::string lookup_key = resolve_nk_lookup_key(allele, nk_lookup);
            if (lookup_key.empty())
            {
                std::string category = nk_category_from_allele(allele);
                if (category.empty()) {
                    continue;
                }

                std::string unmatched_key = allele;
                for (const auto& [existing_key, value] : unmatched_dosage) {
                    (void)value;
                    if (allele == existing_key ||
                        allele_prefix_match(allele, existing_key)) {
                        unmatched_key = existing_key;
                        break;
                    }
                }

                auto current = unmatched_dosage.find(unmatched_key);
                if (current == unmatched_dosage.end() || dosage > current->second.second)
                {
                    unmatched_dosage[unmatched_key] = {category, dosage};
                }

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

            if (category == "HLA_A")
            {
                add_ligand_call(a_ligands, ligand, dosage);
            }
            else if (category == "HLA_C")
            {
                add_ligand_call(c_ligands, ligand, dosage);
            }
            else if (category == "HLA_B")
            {
                add_ligand_call(b_ligands, ligand, dosage);
            }
        }

        for (const auto& [allele, value] : unmatched_dosage)
        {
            (void)allele;
            const std::string& category = value.first;
            int dosage = value.second;

            if (category == "HLA_A")
            {
                add_ligand_call(a_ligands, ".", dosage);
            }
            else if (category == "HLA_C")
            {
                add_ligand_call(c_ligands, ".", dosage);
            }
            else if (category == "HLA_B")
            {
                add_ligand_call(b_ligands, ".", dosage);
            }
        }

        var_name.push_back("HLA_A_ligand");
        genotypes.push_back(format_ligands(a_ligands));

        var_name.push_back("HLA_C_ligand");
        genotypes.push_back(format_ligands(c_ligands));

        var_name.push_back("HLA_B_ligand");
        genotypes.push_back(format_ligands(b_ligands));

        def_NK.push_back({
            sample_id,
            var_name,
            genotypes
        });
    }

    return def_NK;
}
