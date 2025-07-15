#include "utils.h"
#include "main.h"
#include "database.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <bitset>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <random>
#include <unistd.h>
#include <algorithm>

using namespace std;
using namespace seal;

void mmkpirl(size_t N, size_t num_n, size_t num_m, size_t max_keyword_bytelength, 
                   size_t max_item_bytelength, size_t k)
{
    print_example_banner("MMKPIR(L)");
    

    uint64_t frac_k = factorial(k);

    cout << "Using parameters:" << endl;
    cout << "   N (polynomial modulus): " << N << endl;
    cout << "   n (keywords): " << num_n << endl;
    cout << "   m (items): " << num_m << endl;
    cout << "   max_keyword_bytelength: " << max_keyword_bytelength << endl;
    cout << "   max_item_bytelength: " << max_item_bytelength << endl;
    cout << "   k (CW hamming weight): " << k << endl;

    DatabaseGenerator dbGenerator(
        N,
        num_n,
        num_m,
        max_keyword_bytelength,
        max_item_bytelength,
        k
    );

    dbGenerator.generate_structure();

    const auto& cw_keywords = dbGenerator.get_cw_keywords();
    const auto& item_indices = dbGenerator.get_item_indices();
    const auto& multi_hot_vectors = dbGenerator.get_multi_hot_vectors(); 

    dbGenerator.generate_content();
    const auto& item_contents = dbGenerator.get_item_contents();

    const int BITS_PER_CHAR = 8;
    std::vector<std::vector<int>> item_info;
    int max_total_bits = 0;
    int element_count = 0;
        
    for (const auto& vec : item_contents) {
        element_count = vec.size();
        int totalBits = 0;
        
        for (const std::string& item : vec) {
            totalBits += item.size() * BITS_PER_CHAR;
        }

        item_info.push_back({element_count, totalBits});

        max_total_bits = std::max(max_total_bits, totalBits);
    }
    
    cout << "   Maximum Valueset Size: " << max_total_bits/8/1024 << "KB" << std::endl;
    cout << "   Database Size: " << ((max_total_bits/8/1024)*item_info.size())/1024 << "MB" << std::endl;

    std::vector<uint64_t> cw_vector;
    std::vector<std::vector<uint64_t>> cw_vectors;
    for (const std::string& cw_key : cw_keywords) {
        cw_vector = convert_to_int_vector(cw_key);
        cw_vectors.push_back(cw_vector);
    }
    
    std::cout << "Generating selection vector..."<< std::endl;
    std::vector<vector<uint64_t>> previous_selections;
    vector<uint64_t> selection_vector1 = generate_selection_vector(cw_keywords, previous_selections);

    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = N;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(65537);
    SEALContext context(parms);

    print_parameters(context);
    cout << "Parameter validation (success): " << context.parameter_error_message() << endl;

    stringstream parms_stream;
    stringstream data_stream;
    stringstream sk_stream;

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    Encryptor encryptor(context, public_key);
    encryptor.set_secret_key(secret_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    BatchEncoder batch_encoder(context);
    GaloisKeys galois_keys; 
    keygen.create_galois_keys(galois_keys);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    size_t t_bits = parms.plain_modulus().bit_count();
    size_t t = parms.plain_modulus().value();
    size_t var_t = t_bits - 1;

    size_t slot_count = batch_encoder.slot_count();
    std::vector<uint64_t> zero_vec(slot_count,0UL);
    Plaintext plain_zero_vec;
    batch_encoder.encode(zero_vec, plain_zero_vec);

    std::vector<uint64_t> one_vec(slot_count,1UL);
    Plaintext plain_one_vec;
    batch_encoder.encode(one_vec, plain_one_vec);

    cout << "Generating query..." << endl;
    Plaintext plain_q1;
    Plaintext plain_y;
    std::vector<Plaintext> plain_ys;
    for (const auto& cw_vector : cw_vectors){
        batch_encoder.encode(cw_vector, plain_y);
        plain_ys.push_back(plain_y);
    }

    batch_encoder.encode(selection_vector1, plain_q1);

    Ciphertext tilde_q1;
    encryptor.encrypt_symmetric(plain_q1, tilde_q1);
    int b_init   = decryptor.invariant_noise_budget(tilde_q1); 
    cout << "   Fresh noise budget: " << b_init << " bits" << endl;
    auto size_sym_encrypted_selector = encryptor.encrypt_symmetric(plain_q1).save(data_stream);

    cout << "Computing tilde_u..." << endl;
    clock_t t_tu_start = clock();
    
    std::vector<Ciphertext> tilde_us(plain_ys.size());
    for (size_t i = 0; i < plain_ys.size(); ++i) {
        evaluator.multiply_plain(tilde_q1,plain_ys[i],tilde_us[i]); 
        evaluator.relinearize_inplace(tilde_us[i], relin_keys);
    }
    clock_t t_tu_end = clock();
    double t_tu = (static_cast<double>(t_tu_end) - t_tu_start) / CLOCKS_PER_SEC;

    Ciphertext temp;
    std::vector<Ciphertext> tilde_u_primes(tilde_us.size());
    int total_r_count = 0;
    int total_add_count = 0;

    clock_t t_tup_start = clock();
    for (size_t i = 0; i < tilde_us.size(); ++i) {
        temp = tilde_us[i];
        tilde_u_primes[i] = tilde_us[i];
        int m = slot_count;
        int r_count = 0;
        int add_count = 0;
        m = m/2;
        int l = 1;
        while (l < m){
        evaluator.rotate_rows_inplace(temp, l, galois_keys);
        r_count += 1;
        evaluator.add_inplace(tilde_u_primes[i],temp);         
        temp = tilde_u_primes[i];
        add_count += 1;
        l = l*2;
        }
        evaluator.rotate_columns_inplace(tilde_u_primes[i], galois_keys);
        r_count += 1;
        evaluator.add_inplace(tilde_u_primes[i],temp);
        add_count += 1;
        total_r_count += r_count;
        total_add_count += add_count;
    }
    clock_t t_tup_end = clock();
    double t_tup = (static_cast<double>(t_tup_end) - t_tup_start) / CLOCKS_PER_SEC;
    std::cout << "   Rotated " << total_r_count << " times and added " << total_add_count << " times" << std::endl;

    cout << "Computing tilde_s..." << endl;
    Plaintext plain_inv_i;
    std::vector<Plaintext> plain_inv_i_vectors;
    Ciphertext tilde_s;
    Ciphertext tilde_vp;
    std::vector<Ciphertext> tilde_vps;
    std::vector<Ciphertext> tilde_ss;
    std::vector<std::vector<uint64_t>> inv_i_vec = generate_vectors(k,slot_count);
    int mul_count = 0;
    for (const auto& inv_i : inv_i_vec){
        batch_encoder.encode(inv_i, plain_inv_i);
        plain_inv_i_vectors.push_back(plain_inv_i);
        }
    
    clock_t t_ts_start = clock();
    for (const auto& tilde_u_prime : tilde_u_primes) {
        Ciphertext temp;
        std::vector<Ciphertext> temp_vector;
        for (const auto& plain_inv_i : plain_inv_i_vectors){
            evaluator.sub_plain(tilde_u_prime,plain_inv_i,temp);  
            temp_vector.push_back(temp);
        }
        tilde_vp = DAC_Product(temp_vector, 0, temp_vector.size() - 1, evaluator, decryptor, relin_keys, mul_count, 1);
        evaluator.relinearize_inplace(tilde_vp, relin_keys);         
        tilde_vps.push_back(tilde_vp);
    }

    uint64_t inv_frac_k = multiplicative_inverse(frac_k,t);
    std::vector<uint64_t> inv_frac_k_vec(slot_count,inv_frac_k);
    Plaintext plain_inv_frac_k_vec;
    batch_encoder.encode(inv_frac_k_vec, plain_inv_frac_k_vec);

    for (auto& tilde_vp : tilde_vps){
        evaluator.multiply_plain(tilde_vp,plain_inv_frac_k_vec,tilde_s); 
        evaluator.relinearize_inplace(tilde_s, relin_keys);       
        tilde_ss.push_back(tilde_s);
    }
    for (auto& tilde_s : tilde_ss){
        evaluator.transform_to_ntt_inplace(tilde_s);
    }
    clock_t t_ts_end = clock();
    double t_ts = ((static_cast<double>(t_ts_end) - t_ts_start)) / CLOCKS_PER_SEC;

    cout << "Encoding database..." << endl;
    clock_t t_predb_start = clock(); 
    size_t total_element_bits = 0;
    for (int j=0; j<element_count; ++j ) {
        total_element_bits += max_item_bytelength* BITS_PER_CHAR;
    }

    size_t var_e = ceil(static_cast<float>(total_element_bits)/(N * var_t));
    std::cout << "   e = "<< var_e<< std::endl;

    auto payloads = transcoding(item_contents, var_t * N);
    vector<vector<Plaintext>> plain_payload_vectors_list(payloads.size());
    vector<vector<vector<uint64_t>>> payload_vector(payloads.size());

    for (size_t n = 0; n < payloads.size(); n++) {
        plain_payload_vectors_list[n].resize(var_e);
        payload_vector[n].resize(var_e);
        
        for (size_t e = 0; e < var_e; e++) {
            payload_vector[n][e].resize(N, 0);
        }
    }
    for (size_t n = 0; n < payloads.size(); n++) {
        const auto& payload_n = payloads[n];
        
        for (size_t e = 0; e < min(var_e, payload_n.size()); e++) {
            const auto& src = payload_n[e];
            auto& dst = payload_vector[n][e];
            
            if (src.size() > 0) {
                const size_t copy_size = min(src.size(), N);
                copy(src.begin(), src.begin() + copy_size, dst.begin());
            }
            
            batch_encoder.encode(dst, plain_payload_vectors_list[n][e]);
            evaluator.transform_to_ntt_inplace(plain_payload_vectors_list[n][e], context.first_parms_id());
        }
        
        for (size_t e = payload_n.size(); e < var_e; e++) {
            batch_encoder.encode(payload_vector[n][e], plain_payload_vectors_list[n][e]);
            evaluator.transform_to_ntt_inplace(plain_payload_vectors_list[n][e], context.first_parms_id());
        }
    }

    clock_t t_predb_end = clock();
    
    cout << "Computing inner product..." << endl;
    std::vector<std::vector<Ciphertext>> tilde_osml;

    clock_t t_ip_start = clock();

    if (plain_payload_vectors_list.size() == tilde_ss.size()){
        for (size_t i=0; i < plain_payload_vectors_list.size(); ++i){
            std::vector<Ciphertext> group_tilde_osml;
            for (size_t group_index = 0; group_index < plain_payload_vectors_list[i].size(); ++group_index){
                Ciphertext tilde_oml;
                evaluator.multiply_plain(tilde_ss[i],plain_payload_vectors_list[i][group_index],tilde_oml);         
                evaluator.relinearize_inplace(tilde_oml, relin_keys); 
                group_tilde_osml.push_back(tilde_oml);
            }
            tilde_osml.push_back(group_tilde_osml);
        }
    }
    else {
        throw std::invalid_argument("plain_payload_vectors_list must be of the same length with selection vectors.");
    } 

    size_t maxE = 0;
    for (const auto& vec : tilde_osml) {
        if (vec.size() > maxE) {
            maxE = vec.size();
        }
        
    }

    std::vector<Ciphertext> tilde_RM(maxE);
    for (Ciphertext& tilde_r : tilde_RM){
        encryptor.encrypt_symmetric(plain_zero_vec,tilde_r);
    }

    if (tilde_osml.size()==1){
        tilde_RM = tilde_osml[0];
        for (Ciphertext& tilde_r : tilde_RM){
            evaluator.transform_from_ntt_inplace(tilde_r);
        }
    }
    else{
        std::vector<Ciphertext> sum(maxE);   
        for (size_t x = 0; x < maxE; ++x) {
            encryptor.encrypt_symmetric(plain_zero_vec,sum[x]);
            evaluator.transform_to_ntt_inplace(sum[x]);
            for (const auto& vec : tilde_osml) {
                if (x < vec.size()) {
                    evaluator.add_inplace(sum[x],vec[x]);         
                    evaluator.relinearize_inplace(sum[x], relin_keys);                 
                }
            }

            evaluator.transform_from_ntt_inplace(sum[x]);
            
            tilde_RM[x] = sum[x];
        }
    }

    clock_t t_ip_end = clock();

    cout << "\n...Finished!" << endl;
    double t_ip = ((static_cast<double>(t_ip_end) - t_ip_start)) / CLOCKS_PER_SEC; 

    std::cout << "\nPreparation time: " << (static_cast<double>(t_predb_end) - t_predb_start) / CLOCKS_PER_SEC << " seconds" << std::endl; 
    std::cout << "Selection time: " << t_tu+t_tup+t_ts << " seconds" << std::endl; 
    std::cout << "Inner product time: " << t_ip << " seconds" << std::endl; 
    std::cout << "Total server time: " << t_tu+t_tup+t_ts+t_ip << " seconds" << std::endl; 

    cout << "Remaining noise budget: " << decryptor.invariant_noise_budget(tilde_RM[0]) << " bits" << endl;

    
    data_stream.str("");
    auto size_encrypted_answer = 0;

    for (size_t i = 0; i < tilde_RM.size(); i++) {
        while (context.get_context_data(tilde_RM[i].parms_id())->next_context_data()) {
            evaluator.mod_switch_to_next_inplace(tilde_RM[i]);
        }
        size_encrypted_answer += tilde_RM[i].save(data_stream);
    }

    std::cout << "Upload cost: " << size_sym_encrypted_selector/1024 <<" KB" << std::endl; 
    std::cout << "Download cost: "<< size_encrypted_answer / 1024 << " KB" << std::endl;

}
