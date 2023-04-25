#include <vector>
#include <iostream>
#include <stdint.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sstream>
#include <bit>

using namespace std;

vector<vector<int64_t>> binomials;

// Precalculates (n choose k) for 0 <= n <= N and 0 <= k <= n
void binomial_precalc(int64_t N){
    binomials.clear();
    binomials.push_back({1}); // 0 choose 0
    for(int64_t n = 1; n <= N; n++){
        vector<int64_t> row = {1}; // n choose 0
        for(int64_t i = 1; i <= n-1; i++){ // n choose 1..n-1
            row.push_back(binomials.back()[i-1] + binomials.back()[i]);
        }
        row.push_back(1); // n choose n
        binomials.push_back(row);
    }
}

int bits_needed(uint64_t x){
    return 64 - std::countl_zero(x);
}

// Returns the multinomial (n choose k_1, ... k_m) where K = (k_1, ... k_m)
// (32 choose 8,8,8,8) < 2^63
int64_t multinomial(int64_t n, const int64_t* K, int64_t m){
    if(n == 0) return 1;
    int64_t ans = 1;
    for(int64_t i = 0; i < m; i++){
        int64_t k = K[i];
        ans *= binomials[n][k];
        n -= k;
    }
    return ans;
}

// Takes a list counts of symbols
// Returns the string with the given rank from the lexicographically sorted
// list of strings that have the multiset of symbols with the given counts.
// The rank is zero-based.
vector<char> block_unrank(const int64_t* counts, int64_t n_counts, int64_t total_count, int64_t rank){
    int64_t sigma = n_counts;
    int64_t len = total_count;
    int64_t sum_so_far = 0;

    // Copy the counts into a local array for modification
    int64_t remaining_counts[4];
    for(int64_t i = 0; i < n_counts; i++)
        remaining_counts[i] = counts[i];

    vector<char> ans(len);
    for(int64_t i = 0; i < len; i++){
        int64_t remaining_total = len - i;
        for(char c = 0; c < sigma; c++){
            // Is c the next character?
            if(remaining_counts[c] == 0) continue;            
            remaining_counts[c]--; // Temporary
            int64_t howmany_from_c = multinomial(remaining_total-1, remaining_counts, sigma);
            remaining_counts[c]++; // Undo temporary
            if(rank < sum_so_far + howmany_from_c){
                // c is the next character
                ans[i] = c;
                remaining_counts[c]--; // Subtract permanently
                break;
            } else{
                sum_so_far += howmany_from_c;
            }
        }
    }
    return ans;
}

// Takes a list counts of symbols and a lexicographic rank.
// Decodes the block and counts the number of occurrences of `symbol` in the block indices [0..p)
int64_t decode_and_count_prefix(const int64_t* counts, int64_t n_counts, int64_t total_count, int64_t rank, int64_t symbol, int64_t p){
    int64_t sigma = n_counts;
    int64_t len = total_count;
    int64_t sum_so_far = 0;

    // Copy the counts into a local array for modification
    int64_t remaining_counts[4];
    for(int64_t i = 0; i < n_counts; i++)
        remaining_counts[i] = counts[i];

    int64_t ans = 0;
    for(int64_t i = 0; i < p; i++){
        int64_t remaining_total = len - i;
        for(char c = 0; c < sigma; c++){
            // Is c the next character?
            if(remaining_counts[c] == 0) continue;            
            remaining_counts[c]--; // Temporary
            int64_t howmany_from_c = multinomial(remaining_total-1, remaining_counts, sigma);
            remaining_counts[c]++; // Undo temporary
            if(rank < sum_so_far + howmany_from_c){
                // c is the next character
                ans += (c == symbol);
                remaining_counts[c]--; // Subtract permanently
                break;
            } else{
                sum_so_far += howmany_from_c;
            }
        }
    }
    return ans;
}

// Returns the 0-based lexicographic rank of the block among all blocks that have the same multiset of symbols
uint64_t block_rank(const vector<char>& block, const vector<int64_t>& counters){
    int64_t ans = 0;
    int64_t n = block.size();
    vector<int64_t> counters_copy = counters;
    for(char c : block){
        for(char d = 0; d < c; d++){
            if(counters_copy[d] > 0){
                counters_copy[d]--; // Subtract temporarily
                ans += multinomial(n-1, counters_copy.data(), counters_copy.size()); // remaining combinations
                counters_copy[d]++; // Add back temporarily
            }
        }
        n--;
        counters_copy[c]--;
    }
    return ans;
}

class Base4RRR{

public:

    typedef uint32_t A_int_t;
    typedef uint32_t B_int_t;

    static constexpr int64_t BLOCK_SIZE = 64; // bits
    static constexpr int64_t SUPERBLOCK_SIZE = 2048; // bits
    static constexpr int64_t SYMBOLS_IN_BLOCK = BLOCK_SIZE / 2; // Two bits per symbol
    static constexpr int64_t SYMBOLS_IN_SUPERBLOCK = SUPERBLOCK_SIZE / 2; // Two bits per symbol

    // Concatenated (sigma-1) counters per superblock, 32 bits per count.
    vector<A_int_t> A;

    // Pointers from superblock index to the start of the binary representation in D of the
    // lexicographic rank of the first block in the superblock.
    vector<B_int_t> B;

    // Block classes. Each class is a list of (sigma-1) 6-bit numbers.
    // The count for the last symbol can be inferred.
    sdsl::int_vector<6> C;

    // Concatenated variable-length lexicographic ranks of blocks within their classes
    sdsl::bit_vector D;

    // Total number of symbols in the dadta structure
    int64_t n_symbols;

    // Size of the alphabet. Must be 4 or 3
    int64_t sigma;

    Base4RRR(const vector<char>& original_seq) : n_symbols(original_seq.size()){

        sigma = *std::max_element(original_seq.begin(), original_seq.end()) + 1;

        assert(sigma == 3 || sigma == 4);
        assert(*std::min_element(original_seq.begin(), original_seq.end()) == 0);

        binomial_precalc(SYMBOLS_IN_BLOCK);

        // Pad the sequence to end at a superblock boundary
        vector<char> seq = original_seq;
        while(seq.size()%SYMBOLS_IN_SUPERBLOCK > 0)
            seq.push_back(0);

        int64_t n_blocks = seq.size() / SYMBOLS_IN_BLOCK;
        int64_t n_superblocks = seq.size() / SYMBOLS_IN_SUPERBLOCK;
        int64_t blocks_in_superblock = SYMBOLS_IN_SUPERBLOCK / SYMBOLS_IN_BLOCK;

        A.resize(n_superblocks*(sigma-1));
        B.resize(n_superblocks);
        C.resize((sigma-1)*n_blocks);
        // Dont' know the size of D yet
        
        vector<bool> D_bits;

        vector<int64_t> global_counters(sigma, 0);
        for(int64_t b = 0; b < n_blocks; b++){
            if(b % blocks_in_superblock == 0){
                // First block in superblock
                int64_t superblock_idx = b / blocks_in_superblock;
                for(int64_t i = 0; i < sigma-1; i++){
                    A[superblock_idx*(sigma-1) + i] = global_counters[i];
                }
                B[superblock_idx] = D_bits.size();
            }

            vector<int64_t> local_counters(sigma, 0); // Counters inside the block
            vector<char> block_symbols;
            for(int64_t i = 0; i < SYMBOLS_IN_BLOCK; i++){
                char c = seq[b*SYMBOLS_IN_BLOCK + i];
                block_symbols.push_back(c);
                global_counters[c]++;
                local_counters[c]++;
            }

            for(int64_t i = 0; i < sigma-1; i++){
                C[b*(sigma-1) + i] = local_counters[i];
                // The last counter for i = sigma is not stored
            }
            
            uint64_t class_offset = block_rank(block_symbols, local_counters);
            uint64_t n_bits = bits_needed(multinomial(SYMBOLS_IN_BLOCK, local_counters.data(), local_counters.size()));

            sdsl::bit_vector buf(64); // TODO do without memory allocation
            buf.set_int(0, class_offset, n_bits);

            // Extract the bits in the class offset and store to D_bits
            for(int64_t i = 0; i < n_bits; i++){
                D_bits.push_back(buf[i]);
            }
        }

        D.resize(D_bits.size());
        for(int64_t i = 0; i < D.size(); i++) D[i] = D_bits[i];
    }

    void debug_print_binary(uint64_t x) const{
        for(int64_t i = 0; i < 64; i++){
            cout << ((x >> (63-i)) & 1);
        }
        cout << endl;
    }

    // Rank of symbol in half-open interval [0..pos)
    int64_t rank(int64_t pos, char symbol) const{
        // TODO: does this work correctly for pos one past the end of the sequence?

        int64_t superblock_idx = pos / SYMBOLS_IN_SUPERBLOCK;
        int64_t block_idx = pos / SYMBOLS_IN_BLOCK;

        int64_t ans = 0;

        // Add the counts up to the superblock
        if(symbol == sigma-1){ // The count is not explicitly stored for the last symbol, but we can infer it
            int64_t count = superblock_idx * SYMBOLS_IN_SUPERBLOCK; // Count of all symbols up to but not including this superblock
            for(int64_t i = 0; i < sigma-1; i++) // Subtract the counts of all symbols except the last one
                count -= A[superblock_idx*(sigma-1) + i];
            ans += count;
        } else{ // The count is explicitly stored
            ans += A[superblock_idx*(sigma-1)+symbol];
        }

        // Add the counts in the block inside the superblock up to the block containing pos
        int64_t first_block = superblock_idx * SUPERBLOCK_SIZE / BLOCK_SIZE;
        int64_t D_offset = 0; // Sum of binary representations of block ranks
        int64_t counters[4]; // Assumes sigma <= 4
        for(int64_t b = first_block; b <= block_idx; b++){

            // Extract counters for the block
            int64_t sum_of_counters = 0;
            for(int64_t i = 0; i < sigma-1; i++){
                counters[i] = C[(sigma-1)*b + i];
                sum_of_counters += counters[i];
            }
            counters[sigma-1] = SYMBOLS_IN_BLOCK - sum_of_counters; // The last counter

            int64_t lex_rank_bit_size = bits_needed(multinomial(SYMBOLS_IN_BLOCK, counters, sigma));
            if(b <= block_idx - 1){ // Not the last one
                D_offset += lex_rank_bit_size;
                ans += counters[symbol];
            } else { // The last block
                // Decode the last block and count symbols
                int64_t lex_rank = D.get_int(B[superblock_idx] + D_offset, lex_rank_bit_size); // Decode the lex rank of the block
                ans += decode_and_count_prefix(counters, sigma, SYMBOLS_IN_BLOCK, lex_rank, symbol, pos % SYMBOLS_IN_BLOCK);
            }
        }

        return ans;
    }

    string size_report() const{
        
        // Bit sizes
        int64_t A_size = A.size() * sizeof(A_int_t) * 8;
        int64_t B_size = B.size() * sizeof(B_int_t) * 8;
        int64_t C_size = sdsl::size_in_bytes(C) * 8;
        int64_t D_size = sdsl::size_in_bytes(D) * 8;

        stringstream ss;
        ss << "A: " << (double) A_size / n_symbols << " bits per symbol" << endl;
        ss << "B: " << (double) B_size / n_symbols << " bits per symbol" << endl;
        ss << "C: " << (double) C_size / n_symbols << " bits per symbol" << endl;
        ss << "D: " << (double) D_size / n_symbols << " bits per symbol" << endl;
        ss << "Total: " << (double)(A_size + B_size + C_size + D_size) / n_symbols << " bits per symbols" << endl;
        ss << "Total: " << A_size + B_size + C_size + D_size << " bits"; // No endline
        return ss.str();
    }

};