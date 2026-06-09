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

#if defined(__APPLE__)
  #include <mach-o/dyld.h>
#elif defined(__linux__)
  #include <unistd.h>
  #include <climits>
#endif

struct allelescall {
    std::vector<std::string> alleles;
    std::vector<double> dosages;
    std::vector<double> gt_dosages;
};

struct samplerecord {
    std::string sampleid;
    std::map<std::string, std::vector<allelescall>> loci;
};

inline std::filesystem::path get_executable_path() {
#if defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);

    std::string buf(size, '\0');
    if (_NSGetExecutablePath(buf.data(), &size) != 0) {
        throw std::runtime_error("Cannot locate executable path");
    }

    return std::filesystem::weakly_canonical(std::filesystem::path(buf.c_str()));

#elif defined(__linux__)
    char buf[PATH_MAX];         
    ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (len == -1) {
        throw std::runtime_error("Cannot locate executable path");
    }
    buf[len] = '\0';
    return std::filesystem::weakly_canonical(std::filesystem::path(buf));

#else
    throw std::runtime_error("Unsupported platform");
#endif
}

inline std::filesystem::path resource_path(const std::string& filename) {
    const auto exe_dir = get_executable_path().parent_path();
    const std::vector<std::filesystem::path> candidates = {
        (exe_dir / ".." / "data_input" / filename).lexically_normal(),
        (exe_dir / ".." / ".." / "data_input" / filename).lexically_normal(),
        (std::filesystem::current_path() / "data_input" / filename).lexically_normal(),
        (std::filesystem::current_path() / ".." / "data_input" / filename).lexically_normal()
    };

    for (const auto& candidate : candidates) {
        if (std::filesystem::exists(candidate)) {
            return candidate;
        }
    }

    return candidates.front();
}

std::vector<samplerecord> parse_txt_file(const std::string& filename);
std::vector<samplerecord> parse_dosage_file(const std::string& filename);
void error_message(const std::string& msg);
std::string num_decro_stars(const std::string& msg);
void print_samples(const std::vector<samplerecord>& samples);
std::vector<samplerecord> parse_file(int argc, char* argv[]);

#endif // INPUT_FILES_PARSERS_H
