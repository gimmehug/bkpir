#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <unordered_map>
#include <seal/seal.h>

using namespace seal;

std::vector<std::vector<std::vector<uint64_t>>> transcoding(const std::vector<std::vector<std::string>>& vec, size_t maxBits);

std::vector<std::vector<uint64_t>> 
split_vector(const std::vector<int>& multiHotVector, int N);

std::vector<std::vector<std::string>> 
split_string_list(const std::vector<std::string>& string_list, int N);

std::vector<uint64_t> 
convert_to_int_vector(const std::string& str);

std::vector<uint64_t> 
generate_selection_vector(const std::vector<std::string>& input_keywords, 
                          std::vector<std::vector<uint64_t>>& previous_selections);

void 
reset_previous_selections(std::vector<std::vector<uint64_t>>& previous_selections);

std::vector<std::vector<uint64_t>> 
generate_vectors(size_t k, size_t n);

Ciphertext 
DAC_Product(const std::vector<Ciphertext>& ct_vector, 
            size_t left, 
            size_t right, 
            Evaluator& eval, 
            Decryptor& decrpt, 
            RelinKeys& relin_keys,
            int& mul_count, 
            int depth);

uint64_t factorial(size_t k);

int extended_gcd(int a, int b, int& x, int& y);

size_t multiplicative_inverse(int b, int t);

#endif // UTILS_H
