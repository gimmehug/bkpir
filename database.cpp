#include "database.h"
#include <random>
#include <algorithm>
#include <stdexcept>
#include <cctype>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <iterator>

DatabaseGenerator::DatabaseGenerator(size_t N, 
                                     size_t num_keywords_1,
                                     size_t num_item,
                                     size_t MAX_KEYWORD_BYTELENGTH,
                                     size_t max_item_bytelength,
                                     size_t k)
    : N(N), num_keywords_1(num_keywords_1), num_item(num_item), MAX_KEYWORD_BYTELENGTH(MAX_KEYWORD_BYTELENGTH),
      max_item_bytelength(max_item_bytelength),
      k(k) {
    
    precompute_binom();
}

void DatabaseGenerator::generate_structure() {
    std::set<std::string> keyword_set;
    while (keyword_set.size() < num_keywords_1) {
        keyword_set.insert(generate_random_string(MAX_KEYWORD_BYTELENGTH));
    }
    keywords1.assign(keyword_set.begin(), keyword_set.end());
    
    std::vector<size_t> all_item_indices(num_item);
    std::iota(all_item_indices.begin(), all_item_indices.end(), 0);
    
    item_indices.resize(keywords1.size());
    
    for (size_t i = 0; i < keywords1.size(); ++i) {
        item_indices[i] = all_item_indices;
    }
    
    multi_hot_vectors.resize(keywords1.size(), std::vector<int>(num_item, 0));
    
    cw_keywords.clear();
    for (size_t i = 0; i < keywords1.size(); ++i) {
        uint64_t md5_int = convert_to_md5_int(keywords1[i]);
        
        std::string cw_key = generate_cw_code(md5_int, N, k);
        cw_keywords.push_back(cw_key);
        
        for (size_t j = 0; j < num_item; j++) {
            multi_hot_vectors[i][j] = 1;
        }
    }
}

void DatabaseGenerator::generate_content() {
    if (keywords1.empty()) {
        throw std::runtime_error("Database structure must be generated first");
    }
    
    global_items.clear();
    global_items.reserve(num_item);
    std::uniform_int_distribution<char> char_dis(32, 126);
    for (size_t i = 0; i < num_item; ++i) {
        std::string str(max_item_bytelength, '\0');
        for (size_t j = 0; j < max_item_bytelength; ++j) {
            str[j] = char_dis(gen);
        }
        global_items.push_back(str);
    }
    
    item_contents.resize(keywords1.size());
    for (size_t i = 0; i < keywords1.size(); ++i) {
        item_contents[i].clear();
        for (size_t index : item_indices[i]) {
            if (index < global_items.size()) {
                item_contents[i].push_back(global_items[index]);
            } else {
                throw std::out_of_range("Item index out of range");
            }
        }
    }
}
void DatabaseGenerator::precompute_binom() {
    binom.resize(N + 1, std::vector<long long>(k + 1, 0));
    for (size_t i = 0; i <= N; ++i) {
        for (size_t j = 0; j <= std::min(i, k); ++j) {
            if (j == 0 || j == i) {
                binom[i][j] = 1;
            } else {
                binom[i][j] = binom[i - 1][j - 1] + binom[i - 1][j];
            }
        }
    }
}

std::string DatabaseGenerator::generate_random_string(size_t length) const {
    static const char alphanum[] = 
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    
    std::string result;
    std::uniform_int_distribution<> dis(0, sizeof(alphanum) - 2);
    for (size_t i = 0; i < length; ++i) {
        result += alphanum[dis(gen)];
    }
    return result;
}

uint64_t DatabaseGenerator::convert_to_md5_int(const std::string& keyword) const {
    unsigned char digest[MD5_DIGEST_LENGTH];
    MD5((const unsigned char*)keyword.c_str(), keyword.size(), digest);
    
    uint64_t result = 0;
    for (int i = 0; i < 8; ++i) {
        result = (result << 8) | digest[i];
    }
    
    
    return result;
}

std::string DatabaseGenerator::generate_cw_code(uint64_t n, size_t m, size_t k_val) const {
    std::string y(m, '0');
    size_t k_temp = k_val;
    uint64_t n_temp = n;
    
    uint64_t total_combinations = binom[m][k_val];
    
    n_temp %= total_combinations;
    
    for (int m_prime = m - 1; m_prime >= 0; --m_prime) {
        if (k_temp == 0) break;
        
        if (static_cast<size_t>(m_prime) >= k_temp - 1) {
            uint64_t binom_val = binom[m_prime][k_temp];
            
            if (n_temp >= binom_val) {
                y[m_prime] = '1';
                n_temp -= binom_val;
                k_temp--;
            }
        }
    }
    
    return y;
}
