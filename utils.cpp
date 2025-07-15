#include "utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include <cctype>
#include <unistd.h>
#include <cstdint>
#include <stdexcept>

using namespace std;
using namespace seal;


std::vector<std::vector<std::vector<uint64_t>>> transcoding(
    const std::vector<std::vector<std::string>>& vec, size_t max_bits)
{
    const size_t bits_per_element = 16;
    const size_t elements_per_group = max_bits / bits_per_element;

    std::vector<std::vector<std::vector<uint64_t>>> grouped_int_vec_list(vec.size());

    for (int idx = 0; idx < vec.size(); ++idx) {
        const auto& inner_vec = vec[idx];
        std::vector<std::vector<uint64_t>> grouped_int_vectors;
        std::vector<uint64_t> current_int_vector;
        current_int_vector.reserve(elements_per_group);

        for (const std::string& str : inner_vec) {
            const char* cstr = str.data();
            size_t len = str.size();

            for (size_t i = 0; i < len; i += 2) {
                uint16_t combined = static_cast<uint8_t>(cstr[i]) |
                                   (i + 1 < len ? (static_cast<uint8_t>(cstr[i + 1]) << 8) : 0);
                current_int_vector.push_back(static_cast<uint64_t>(combined));

                if (current_int_vector.size() >= elements_per_group) {
                    grouped_int_vectors.push_back(std::move(current_int_vector));
                    current_int_vector.clear();
                    current_int_vector.reserve(elements_per_group);
                }
            }
        }

        if (!current_int_vector.empty()) {
            grouped_int_vectors.push_back(std::move(current_int_vector));
        }

        grouped_int_vec_list[idx] = std::move(grouped_int_vectors);
    }

    return grouped_int_vec_list;
}

std::vector<std::vector<uint64_t>> split_vector(const std::vector<int>& multiHotVector, int N) {
    int m = multiHotVector.size();
    int f = std::ceil(static_cast<double>(m) / N);
    std::vector<std::vector<uint64_t>> partitionedVectors(f, std::vector<uint64_t>(N, 0));

    for (int i = 0; i < m; ++i) {
        partitionedVectors[i / N][i % N] = static_cast<uint64_t>(multiHotVector[i]);
    }

    return partitionedVectors;
}

std::vector<vector<string>> split_string_list(const vector<string>& string_list, int N) {
    int f = ceil(static_cast<float>(string_list.size()) / N);
    vector<vector<string>> splitStrings(f, vector<string>(N, ""));
    for (int i = 0; i < string_list.size(); i++) {
        splitStrings[i / N][i % N] = string_list[i];
    }
    return splitStrings;
}


std::vector<uint64_t> convert_to_int_vector(const std::string& str) {
    std::vector<uint64_t> int_vector;
    for (char c : str) {
        int_vector.push_back(c - '0');
    }
    return int_vector;
}

std::vector<uint64_t> generate_selection_vector(const std::vector<std::string>& input_keywords, std::vector<std::vector<uint64_t>>& previous_selections) {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    std::vector<uint64_t> int_vector;
    std::vector<std::vector<uint64_t>> cw_vectors;
    for (const std::string& keyword : input_keywords) {
    int_vector = convert_to_int_vector(keyword);
    cw_vectors.push_back(int_vector);
    }
    std::vector<std::vector<uint64_t>> available_vectors;
    for (const auto& vector : cw_vectors) {
        if (std::find(previous_selections.begin(), previous_selections.end(), vector) == previous_selections.end()) {
            available_vectors.push_back(vector);
        }
    }
    if (available_vectors.empty()) {
        return {};
    }
    std::uniform_int_distribution<> dist(0, available_vectors.size() - 1);
    int index = dist(gen);
    std::vector<uint64_t> select_vector = available_vectors[index];
    previous_selections.push_back(select_vector);
    
    return select_vector;
}
void reset_previous_selections(std::vector<vector<uint64_t>>& previous_selections) {
    previous_selections.clear();
}

std::vector<std::vector<uint64_t>> generate_vectors(size_t k, size_t n) {
    std::vector<std::vector<uint64_t>> result;
    for (uint64_t i = 0; i < k; ++i) {
        std::vector<uint64_t> vec(n);
        fill(vec.begin(), vec.end(), i);
        result.push_back(vec);
    }
    return result;
}

seal::Ciphertext DAC_Product(const std::vector<Ciphertext>& ct_vector, size_t left, size_t right, Evaluator& eval, Decryptor& decrpt, RelinKeys& relin_keys,int& mul_count, int depth)  {
    Ciphertext result;
    if (left == right) {
        return ct_vector[left];
    } else {
        size_t mid = (left + right) / 2;
        Ciphertext left_product = DAC_Product(ct_vector, left, mid, eval, decrpt, relin_keys, mul_count, depth + 1);
        Ciphertext right_product = DAC_Product(ct_vector, mid + 1, right, eval, decrpt, relin_keys,mul_count, depth + 1);
        eval.multiply(left_product, right_product, result);
        eval.relinearize_inplace(result, relin_keys);
        mul_count ++;
        return result;
    }
}

uint64_t factorial(size_t k) {
    static const uint64_t table[] = {
        1ULL,                   // 0!
        1ULL,                   // 1!
        2ULL,                   // 2!
        6ULL,                   // 3!
        24ULL,                  // 4!
        120ULL,                 // 5!
        720ULL,                 // 6!
        5040ULL,                // 7!
        40320ULL,               // 8!
        362880ULL,              // 9!
        3628800ULL,             // 10!
        39916800ULL,            // 11!
        479001600ULL,           // 12!
        6227020800ULL,          // 13!
        87178291200ULL,         // 14!
        1307674368000ULL,       // 15!
        20922789888000ULL,      // 16!
        355687428096000ULL,     // 17!
        6402373705728000ULL,    // 18!
        121645100408832000ULL,  // 19!
        2432902008176640000ULL  // 20!
    };

    if (k > 20) {
        throw std::overflow_error("k exceeds maximum value (20) for uint64_t factorial");
    }
    return table[k];
}

int extended_gcd(int a, int b, int& x, int& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    } else {
        int x1, y1;
        int gcd = extended_gcd(b, a % b, x1, y1);
        x = y1;
        y = x1 - (a / b) * y1;
        return gcd;
    }
}
size_t multiplicative_inverse(int b, int t) {
    int x, y;
    int gcd = extended_gcd(b, t, x, y);
    if (gcd != 1) {
        return -1;
    } else {
        return (x % t + t) % t;
    }
}
