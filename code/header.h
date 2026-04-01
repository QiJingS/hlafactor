#ifndef INPUT_FILES_PARSERS_H
#define INPUT_FILES_PARSERS_H

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include <htslib/hts.h>
#include <htslib/vcf.h>

struct allelescall {
    std::vector<std::string> alleles;
    std::vector<double> dosages;
};

struct samplerecord {
    std::string sampleid;
    std::map<std::string, std::vector<allelescall>> loci;
};

std::vector<samplerecord> parse_txt_file(const std::string& filename);
std::vector<samplerecord> parse_dosage_file(const std::string& filename);
void error_message(const std::string& msg);
std::string num_decro_stars(const std::string& msg);
void print_samples(const std::vector<samplerecord>& samples);
std::vector<samplerecord> parse_file(int argc, char* argv[]);

#endif // INPUT_FILES_PARSERS_H