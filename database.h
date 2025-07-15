#ifndef DATABASE_GENERATOR_H
#define DATABASE_GENERATOR_H

#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <cstddef>
#include <random>
#include <openssl/md5.h>

class DatabaseGenerator {
public:
    DatabaseGenerator(size_t N, 
                      size_t num_keywords_1,
                      size_t num_item,
                      size_t MAX_KEYWORD_BYTELENGTH,
                      size_t max_item_bytelength,
                      size_t k);

    void generate_structure();
    
    void generate_content();

    const std::vector<std::string>& get_keywords1() const { return keywords1; }
    const std::vector<std::string>& get_cw_keywords() const { return cw_keywords; }
    const std::vector<std::vector<size_t>>& get_item_indices() const { return item_indices; }
    const std::vector<std::vector<int>>& get_multi_hot_vectors() const { return multi_hot_vectors; }
    const std::vector<std::vector<std::string>>& get_item_contents() const { return item_contents; }
    const std::vector<std::string>& get_global_items() const { return global_items; }

    size_t get_num_items() const { return num_item; }
    size_t get_max_item_bytelength() const { return max_item_bytelength; }

    void print_stats() const;
    void save_to_files(const std::string& base_path, int log_poly) const;

private:
    size_t N;
    size_t num_keywords_1;
    size_t num_item;
    size_t MAX_KEYWORD_BYTELENGTH;
    size_t max_item_bytelength;
    size_t k;
    size_t md5_length;
    
    std::vector<std::vector<long long>> binom;
    
    mutable std::mt19937 gen;
    
    std::vector<std::string> keywords1;
    std::vector<std::string> cw_keywords;
    std::vector<std::vector<size_t>> item_indices;
    std::vector<std::vector<int>> multi_hot_vectors;
    std::vector<std::vector<std::string>> item_contents;
    std::vector<std::string> global_items;

    void precompute_binom();
    std::string generate_random_string(size_t length) const;
    uint64_t convert_to_md5_int(const std::string& keyword) const;
    std::string generate_cw_code(uint64_t n, size_t m, size_t k_val) const;

};

#endif
