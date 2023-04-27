#include <vector>
#include <iostream>
#include <stdint.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sstream>
#include <bit>

using namespace std;

// (n choose k) for 0 <= n,k <= 31
constexpr int32_t binomials[32][32] =
{
{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 3, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 4, 6, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 5, 10, 10, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 6, 15, 20, 15, 6, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 7, 21, 35, 35, 21, 7, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 8, 28, 56, 70, 56, 28, 8, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 9, 36, 84, 126, 126, 84, 36, 9, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 11, 55, 165, 330, 462, 462, 330, 165, 55, 11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 12, 66, 220, 495, 792, 924, 792, 495, 220, 66, 12, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 13, 78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78, 13, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 14, 91, 364, 1001, 2002, 3003, 3432, 3003, 2002, 1001, 364, 91, 14, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 15, 105, 455, 1365, 3003, 5005, 6435, 6435, 5005, 3003, 1365, 455, 105, 15, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 11440, 8008, 4368, 1820, 560, 120, 16, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 17, 136, 680, 2380, 6188, 12376, 19448, 24310, 24310, 19448, 12376, 6188, 2380, 680, 136, 17, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 18, 153, 816, 3060, 8568, 18564, 31824, 43758, 48620, 43758, 31824, 18564, 8568, 3060, 816, 153, 18, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 19, 171, 969, 3876, 11628, 27132, 50388, 75582, 92378, 92378, 75582, 50388, 27132, 11628, 3876, 969, 171, 19, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 20, 190, 1140, 4845, 15504, 38760, 77520, 125970, 167960, 184756, 167960, 125970, 77520, 38760, 15504, 4845, 1140, 190, 20, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 21, 210, 1330, 5985, 20349, 54264, 116280, 203490, 293930, 352716, 352716, 293930, 203490, 116280, 54264, 20349, 5985, 1330, 210, 21, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 22, 231, 1540, 7315, 26334, 74613, 170544, 319770, 497420, 646646, 705432, 646646, 497420, 319770, 170544, 74613, 26334, 7315, 1540, 231, 22, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 23, 253, 1771, 8855, 33649, 100947, 245157, 490314, 817190, 1144066, 1352078, 1352078, 1144066, 817190, 490314, 245157, 100947, 33649, 8855, 1771, 253, 23, 1, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 24, 276, 2024, 10626, 42504, 134596, 346104, 735471, 1307504, 1961256, 2496144, 2704156, 2496144, 1961256, 1307504, 735471, 346104, 134596, 42504, 10626, 2024, 276, 24, 1, 0, 0, 0, 0, 0, 0, 0},
{1, 25, 300, 2300, 12650, 53130, 177100, 480700, 1081575, 2042975, 3268760, 4457400, 5200300, 5200300, 4457400, 3268760, 2042975, 1081575, 480700, 177100, 53130, 12650, 2300, 300, 25, 1, 0, 0, 0, 0, 0, 0},
{1, 26, 325, 2600, 14950, 65780, 230230, 657800, 1562275, 3124550, 5311735, 7726160, 9657700, 10400600, 9657700, 7726160, 5311735, 3124550, 1562275, 657800, 230230, 65780, 14950, 2600, 325, 26, 1, 0, 0, 0, 0, 0},
{1, 27, 351, 2925, 17550, 80730, 296010, 888030, 2220075, 4686825, 8436285, 13037895, 17383860, 20058300, 20058300, 17383860, 13037895, 8436285, 4686825, 2220075, 888030, 296010, 80730, 17550, 2925, 351, 27, 1, 0, 0, 0, 0},
{1, 28, 378, 3276, 20475, 98280, 376740, 1184040, 3108105, 6906900, 13123110, 21474180, 30421755, 37442160, 40116600, 37442160, 30421755, 21474180, 13123110, 6906900, 3108105, 1184040, 376740, 98280, 20475, 3276, 378, 28, 1, 0, 0, 0},
{1, 29, 406, 3654, 23751, 118755, 475020, 1560780, 4292145, 10015005, 20030010, 34597290, 51895935, 67863915, 77558760, 77558760, 67863915, 51895935, 34597290, 20030010, 10015005, 4292145, 1560780, 475020, 118755, 23751, 3654, 406, 29, 1, 0, 0},
{1, 30, 435, 4060, 27405, 142506, 593775, 2035800, 5852925, 14307150, 30045015, 54627300, 86493225, 119759850, 145422675, 155117520, 145422675, 119759850, 86493225, 54627300, 30045015, 14307150, 5852925, 2035800, 593775, 142506, 27405, 4060, 435, 30, 1, 0},
{1, 31, 465, 4495, 31465, 169911, 736281, 2629575, 7888725, 20160075, 44352165, 84672315, 141120525, 206253075, 265182525, 300540195, 300540195, 265182525, 206253075, 141120525, 84672315, 44352165, 20160075, 7888725, 2629575, 736281, 169911, 31465, 4495, 465, 31, 1}
};

constexpr int bits_needed(uint64_t x){
    return 64 - __builtin_clzll(x); //std::countl_zero(x);
}

// Returns the multinomial (n choose k_1, ... k_m) where K = (k_1, ... k_m)
// (32 choose 8,8,8,8) < 2^63
constexpr int64_t multinomial(int64_t n, const int64_t* K, int64_t m){
    if(n == 0) return 1;
    int64_t ans = 1;
    for(int64_t i = 0; i < m; i++){
        int64_t k = K[i];
        ans *= binomials[n][k];
        n -= k;
    }
    return ans;
}

// First two counts are enough
constexpr int64_t trinomial(int64_t n, int64_t k1, int64_t k2){
    return (int64_t)binomials[n][k1] * binomials[n-k1][k2];
}

// First three counts are enough
constexpr int64_t quadrinomial(int64_t n, int64_t k1, int64_t k2, int64_t k3){
    return (int64_t)binomials[n][k1] * binomials[n-k1][k2] * binomials[n-k1-k2][k3];
}

// Takes a list counts of symbols and a lexicographic rank.
// Decodes the block and counts the number of occurrences of `symbol` in the block indices [0..p)
constexpr int64_t decode_and_count_prefix(const int64_t* counts, int64_t n_counts, int64_t total_count, int64_t rank, int64_t symbol, int64_t p){
    // 69% of the time is spent in this function by the current benchmark

    const int64_t sigma = n_counts;
    const int64_t len = total_count;
    int64_t sum_so_far = 0;

    // Copy the counts into a local array for modification
    int64_t remaining_counts[4] = {0,0,0,0};
    for(int64_t i = 0; i < n_counts; i++)
        remaining_counts[i] = counts[i];

    int64_t ans = 0;
    for(int64_t i = 0; i < p; i++){
        int64_t remaining_total = len - i;
        for(char c = 0; c < sigma; c++){
            // Is c the next character?
            if(remaining_counts[c] == 0) continue;            
            remaining_counts[c]--; // Temporary
            int64_t howmany_from_c = 0;
            if(sigma == 3) 
                howmany_from_c = trinomial(remaining_total-1, remaining_counts[0], remaining_counts[1]);
            else // sigma = 4
                howmany_from_c = quadrinomial(remaining_total-1, remaining_counts[0], remaining_counts[1], remaining_counts[2]);
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

// Takes a list counts of symbols and a lexicographic rank.
// Decodes the block and counts the number of occurrences of `symbol` in the block indices [0..p)
constexpr int64_t decode_and_count_pair_prefix(const int64_t* counts, int64_t n_counts, int64_t total_count, int64_t rank, int64_t symbol, int64_t p){
    // 69% of the time is spent in this function by the current benchmark

    const int64_t sigma = n_counts;
    const int64_t len = total_count;
    int64_t sum_so_far = 0;

    // Copy the counts into a local array for modification
    int64_t remaining_counts[4] = {0,0,0,0};
    for(int64_t i = 0; i < n_counts; i++)
        remaining_counts[i] = counts[i];

    int64_t ans = 0;
    for(int64_t i = 0; i < p; i++){
        int64_t remaining_total = len - i;
        for(char c = 0; c < sigma; c++){
            // Is c the next character?
            if(remaining_counts[c] == 0) continue;            
            remaining_counts[c]--; // Temporary
            int64_t howmany_from_c = 0;
            if(sigma == 3) 
                howmany_from_c = trinomial(remaining_total-1, remaining_counts[0], remaining_counts[1]);
            else // sigma = 4
                howmany_from_c = quadrinomial(remaining_total-1, remaining_counts[0], remaining_counts[1], remaining_counts[2]);
            remaining_counts[c]++; // Undo temporary
            if(rank < sum_so_far + howmany_from_c){
                // c is the next character
                ans += (c == symbol || c == sigma-1); // symbol OR LAST
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

class RRR_Generalization{

public:

    typedef uint32_t A_int_t;

    static constexpr int64_t BLOCK_SIZE = 62; // bits
    static constexpr int64_t SUPERBLOCK_SIZE = BLOCK_SIZE * 32; // bits
    static constexpr int64_t SYMBOLS_IN_BLOCK = BLOCK_SIZE / 2; // Two bits per symbol
    static constexpr int64_t SYMBOLS_IN_SUPERBLOCK = SUPERBLOCK_SIZE / 2; // Two bits per symbol

    // Concatenated tuples of sigma integers per superblock: sigma-1 counters of symbols
    // for the first sigma-1 symbols (last is inferred), and a pointer to the start of
    // the superblock in D.
    vector<A_int_t> A;

    // Block classes. Each class is a list of (sigma-1) 5-bit numbers.
    // The count for the last symbol can be inferred.
    sdsl::int_vector<5> C;

    // Concatenated variable-length lexicographic ranks of blocks within their classes
    sdsl::bit_vector D;

    // Total number of symbols in the dadta structure
    int64_t n_symbols;

    // Size of the alphabet. Must be 4 or 3
    int64_t sigma;

    RRR_Generalization(){}

    RRR_Generalization(const vector<char>& original_seq) : n_symbols(original_seq.size()){

        sigma = *std::max_element(original_seq.begin(), original_seq.end()) + 1;

        assert(sigma == 3 || sigma == 4);
        assert(*std::min_element(original_seq.begin(), original_seq.end()) == 0);

        // Pad the sequence to end at a superblock boundary
        vector<char> seq = original_seq;
        while(seq.size()%SYMBOLS_IN_SUPERBLOCK > 0)
            seq.push_back(0);

        int64_t n_blocks = seq.size() / SYMBOLS_IN_BLOCK;
        int64_t n_superblocks = seq.size() / SYMBOLS_IN_SUPERBLOCK;
        int64_t blocks_in_superblock = SYMBOLS_IN_SUPERBLOCK / SYMBOLS_IN_BLOCK;

        A.resize(n_superblocks*sigma);
        C.resize((sigma-1)*n_blocks);
        // Dont' know the size of D yet
        
        vector<bool> D_bits;

        sdsl::bit_vector buf(64); // Used for writing sdsl integers to a bit vector

        vector<int64_t> global_counters(sigma, 0);
        vector<int64_t> local_counters(sigma, 0); // Counters inside the block
        vector<char> block_symbols; // Sequence in block
        for(int64_t b = 0; b < n_blocks; b++){

            // Clear memory from previous iteration
            block_symbols.clear();
            for(int64_t i = 0; i < sigma; i++) local_counters[i] = 0;

            if(b % blocks_in_superblock == 0){
                // First block in superblock
                int64_t superblock_idx = b / blocks_in_superblock;
                for(int64_t i = 0; i < sigma-1; i++){
                    A[superblock_idx*sigma + i] = global_counters[i];
                }
                A[superblock_idx*sigma + (sigma-1)] = D_bits.size(); // "Pointer" to d
            }

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
                count -= A[superblock_idx*sigma + i];
            ans += count;
        } else{ // The count is explicitly stored
            ans += A[superblock_idx*sigma+symbol];
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
                int64_t lex_rank = D.get_int(A[superblock_idx*sigma + (sigma-1)] + D_offset, lex_rank_bit_size); // Decode the lex rank of the block
                ans += decode_and_count_prefix(counters, sigma, SYMBOLS_IN_BLOCK, lex_rank, symbol, pos % SYMBOLS_IN_BLOCK);
            }
        }

        return ans;
    }

    // Sums of ranks of symbol AND the last symbol in the alphabet in half-open interval [0..pos)
    int64_t rankpair(int64_t pos, char symbol) const{
        // TODO: does this work correctly for pos one past the end of the sequence?

        assert(symbol != sigma-1);

        int64_t superblock_idx = pos / SYMBOLS_IN_SUPERBLOCK;
        int64_t block_idx = pos / SYMBOLS_IN_BLOCK;

        int64_t ans = 0;

        // Add the counts up to the superblock        
        {
            // Add the count of the last symbol 
            int64_t count = superblock_idx * SYMBOLS_IN_SUPERBLOCK; // Count of all symbols up to but not including this superblock
            for(int64_t i = 0; i < sigma-1; i++) // Subtract the counts of all symbols except the last one
                count -= A[superblock_idx*sigma + i];
            ans += count;

            // Add the count for symbol
            ans += A[superblock_idx*sigma+symbol];
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
                ans += counters[sigma-1]; // Last one
            } else { // The last block
                // Decode the last block and count symbols
                int64_t lex_rank = D.get_int(A[superblock_idx*sigma + (sigma-1)] + D_offset, lex_rank_bit_size); // Decode the lex rank of the block
                ans += decode_and_count_pair_prefix(counters, sigma, SYMBOLS_IN_BLOCK, lex_rank, symbol, pos % SYMBOLS_IN_BLOCK);
            }
        }

        return ans;
    }

    string size_report() const{
        
        // Bit sizes
        int64_t A_size = A.size() * sizeof(A_int_t) * 8;
        int64_t C_size = sdsl::size_in_bytes(C) * 8;
        int64_t D_size = sdsl::size_in_bytes(D) * 8;

        stringstream ss;
        ss << "A: " << (double) A_size / n_symbols << " bits per symbol" << endl;
        ss << "C: " << (double) C_size / n_symbols << " bits per symbol" << endl;
        ss << "D: " << (double) D_size / n_symbols << " bits per symbol" << endl;
        ss << "Total: " << (double)(A_size + C_size + D_size) / n_symbols << " bits per symbols" << endl;
        ss << "Total: " << A_size + C_size + D_size << " bits"; // No endline
        return ss.str();
    }

    size_t size_in_bytes() const{
        int64_t A_size = A.size() * sizeof(A_int_t);
        int64_t C_size = sdsl::size_in_bytes(C);
        int64_t D_size = sdsl::size_in_bytes(D);
        return A_size + C_size + D_size;
    }

};
