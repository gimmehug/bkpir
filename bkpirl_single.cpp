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

void bkpirl_single(size_t N, size_t num_n, size_t num_m, size_t max_keyword_bytelength, 
                   size_t max_item_bytelength, size_t k)
{
    print_example_banner("BKPIR(L)-Single");
    

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

    cout << "Generating queries..." << endl;
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
    uint64_t inv_frac_k = multiplicative_inverse(frac_k,t);
    std::vector<uint64_t> inv_frac_k_vec(slot_count,inv_frac_k);
    Plaintext plain_inv_frac_k_vec;
    batch_encoder.encode(inv_frac_k_vec, plain_inv_frac_k_vec);

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

    for (auto& tilde_vp : tilde_vps){
        evaluator.multiply_plain(tilde_vp,plain_inv_frac_k_vec,tilde_s); 
        evaluator.relinearize_inplace(tilde_s, relin_keys);       
        tilde_ss.push_back(tilde_s);
    }
    clock_t t_ts_end = clock();
    double t_ts = ((static_cast<double>(t_ts_end) - t_ts_start)) / CLOCKS_PER_SEC;

    std::vector<std::vector<Plaintext>> plain_Vf;
    size_t len_mh_vecs = 0;
   
    for (const auto& vec : multi_hot_vectors) {  
        std::vector<std::vector<uint64_t>> re_multi_hot_vecs = split_vector(vec, N);
        len_mh_vecs = re_multi_hot_vecs.size();
        std::vector<Plaintext> plain_mh_vecs;
        for(auto& mh_vec : re_multi_hot_vecs){
            Plaintext plain_mh;
            batch_encoder.encode(mh_vec, plain_mh);
            plain_mh_vecs.push_back(plain_mh);
        }
        plain_Vf.push_back(plain_mh_vecs);
    }
    
    std::cout << "Multi-hot vector list size: " << len_mh_vecs << std::endl;

    cout << "Generating selection mask..." << endl;
    clock_t t_sm_start = clock();
    std::vector<std::vector<Ciphertext>> tilde_os;
    if (plain_Vf.size() == tilde_ss.size()){
        for (size_t i=0; i < plain_Vf.size(); i++){
            std::vector<Ciphertext> group_tilde_os;
            for (size_t group_index = 0; group_index < plain_Vf[i].size(); ++group_index){
                Ciphertext tilde_o;                
                evaluator.multiply_plain(tilde_ss[i],plain_Vf[i][group_index],tilde_o);         
                evaluator.relinearize_inplace(tilde_o, relin_keys); 
                group_tilde_os.push_back(tilde_o);
            }
            tilde_os.push_back(group_tilde_os);
        }
    }
    else {
        throw std::invalid_argument("plain_Vf must be of the same length with selection vectors.");
    } 

    std::vector<Ciphertext> tilde_RB(len_mh_vecs);

    for (Ciphertext& tilde_r : tilde_RB){
        encryptor.encrypt_symmetric(plain_zero_vec,tilde_r);
    }

    for (size_t i = 0; i < len_mh_vecs; i++) {        
        Ciphertext sum1;    
        encryptor.encrypt_symmetric(plain_zero_vec,sum1);
        for (auto& vec : tilde_os) {
            evaluator.add_inplace(sum1,vec[i]);         
            evaluator.relinearize_inplace(sum1, relin_keys);            
        }
        tilde_RB[i] = sum1;
    }

    clock_t t_sm_end = clock();

    int b_select = decryptor.invariant_noise_budget(tilde_RB[0]); 

    for (size_t i = 0; i < tilde_RB.size(); i++) {
        evaluator.transform_to_ntt_inplace(tilde_RB[i]);
    }

    double t_sm = (static_cast<double>(t_sm_end) - t_sm_start) / CLOCKS_PER_SEC;
    
    double t_select = t_sm + t_ts + t_tup + t_tu;

    Ciphertext tilde_rp1;
    const auto& global_items = dbGenerator.get_global_items();

    clock_t t_predb_start = clock();
    vector<vector<string>> splited_vl = split_string_list(global_items, N); 
    size_t largest_element_bits = max_item_bytelength * 8;
    size_t var_l = ceil(static_cast<float>(largest_element_bits) / var_t);
    size_t var_f = splited_vl.size();
    size_t var_N = splited_vl[0].size();
    
    
    cout << "Encoding database..." << endl;
    std::cout << "   l= "<< var_l<< " , f = "<<var_f<< std::endl;

    auto encoded_payloads = transcoding(splited_vl, var_t * N);
  
    vector<vector<Plaintext>> plain_PL(var_l, vector<Plaintext>(var_f));
    vector<vector<vector<uint64_t>>> PL(var_l, vector<vector<uint64_t>>(var_f, vector<uint64_t>(var_N, 0)));

    for (size_t i = 0; i < var_f; ++i) {
        if (i < encoded_payloads.size()) {
            const auto& encoded_group = encoded_payloads[i];
            
            for (size_t k = 0; k < var_l; ++k) {
                if (k < encoded_group.size()) {
                    const auto& payloads = encoded_group[k];
                    
                    for (size_t j = 0; j < var_N; ++j) {
                        if (j < payloads.size()) {
                            PL[k][i][j] = payloads[j];
                        } else {
                            PL[k][i][j] = 0;
                        }
                    }
                } else {
                    std::fill(PL[k][i].begin(), PL[k][i].end(), 0);
                }
                
                batch_encoder.encode(PL[k][i], plain_PL[k][i]);
                evaluator.transform_to_ntt_inplace(plain_PL[k][i], context.first_parms_id());
            }
        } else {
            for (size_t k = 0; k < var_l; ++k) {
                std::fill(PL[k][i].begin(), PL[k][i].end(), 0);
                batch_encoder.encode(PL[k][i], plain_PL[k][i]);
                evaluator.transform_to_ntt_inplace(plain_PL[k][i], context.first_parms_id());
            }
        }
    }

    clock_t t_predb_end = clock();
    double t_predb = (static_cast<double>(t_predb_end) - t_predb_start) / CLOCKS_PER_SEC;

    cout << "Computing inner product..." << endl;
    
    clock_t t_ip_start = clock();
    vector<vector<Ciphertext>> tilde_RBP(var_l,vector<Ciphertext>(var_f));    
    for (size_t k = 0; k < var_l; k++) {
        for (Ciphertext& rpl : tilde_RBP[k]){
            encryptor.encrypt_symmetric(plain_zero_vec,rpl);             
        }
    }
    for (size_t k = 0; k < var_l; k++) {
        for (size_t i = 0; i < var_f; i++) {
            evaluator.multiply_plain(tilde_RB[i],plain_PL[k][i],tilde_RBP[k][i]); 
            evaluator.relinearize_inplace(tilde_RBP[k][i], relin_keys);
            evaluator.transform_from_ntt_inplace(tilde_RBP[k][i]);
        }
    }
    clock_t t_ip_end = clock();
    
    int b_final  = decryptor.invariant_noise_budget(tilde_RBP[0][0]); 

    double t_ip =((static_cast<double>(t_ip_end) - t_ip_start)) / CLOCKS_PER_SEC;  
    
    double t_B = t_ip+t_select;

    cout << "\n...Finished!" << endl;

    std::cout << "\nPreparation time: " << t_predb << " seconds" << std::endl;
    std::cout << "Selection time: " << t_select << " seconds" << std::endl; 

    std::cout << "Inner product time: " << t_ip << " seconds" << std::endl;
    std::cout << "BKPIR(L) total server time: " << t_B << " seconds" << std::endl;

    int max_logic_depth = 0;
    bool valid_select = b_select > 0;
    bool valid_final  = b_final  > 0;

    if (valid_select && valid_final) {
        int select_cost = b_init - b_select;
        int inner_cost  = b_select  - b_final;
        int base_cost   = select_cost + inner_cost;
        std::cout << "Selection noise budget cost: "       << select_cost << "\n";
        std::cout << "Inner product noise budget cost: " << inner_cost << " bits\n";
        std::cout << "Basic budget (selection + inner product) cost: " << base_cost << " bits\n";
    } else {
        std::cout << "  Budget exhausted at some stage. Logic depth = 0.\n";
    }

    data_stream.str("");
    auto size_encrypted_answer = 0;
    for (size_t k = 0; k < var_l; k++) {
        for (size_t i = 0; i < var_f; i++) {
            while (context.get_context_data(tilde_RBP[k][i].parms_id())->next_context_data()) {
                evaluator.mod_switch_to_next_inplace(tilde_RBP[k][i]);
            }
            size_encrypted_answer += tilde_RBP[k][i].save(data_stream);
        }
    }
    
    std::cout << "Upload cost: " << size_sym_encrypted_selector/1024 <<" KB" << std::endl; 
    std::cout << "Download cost: "<< size_encrypted_answer / 1024 << " KB" << std::endl;

}
