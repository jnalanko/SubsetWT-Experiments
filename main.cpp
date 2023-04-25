#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "SBWT.hh"
#include "variants.hh"
#include "throwing_streams.hh"
#include "RRR_generalization.hh"

// A wrapper for the basic SDSL wavelet tree
class Basic_WT{
    public:

        sdsl::wt_blcd<sdsl::rrr_vector<>,
                      sdsl::rrr_vector<>::rank_1_type,
                      sdsl::rrr_vector<>::select_1_type,
                      sdsl::rrr_vector<>::select_0_type> wt;

        Basic_WT(const vector<char>& seq){
            // We have to translate 0..sigma-1 to ascii digits for the SDSL wavelet tree because
            // apparently we can't give a zero-byte because that will terminate the string?
            string ascii;
            for(char c : seq) ascii.push_back('0' + c);
            sdsl::construct_im(wt, ascii.c_str(), 1); // 1: file format is a sequence, not a serialized sdsl object
        }

        // Rank of symbol in half-open interval [0..pos)
        int64_t rank(int64_t pos, char symbol) const{
            return wt.rank(pos, '0' + symbol); // Translate to ascii for the query
        }
};

class SplitStructure{ // TODO Elena
    public:
        SplitStructure(const vector<char>& seq){
            // TODO
        } 

        // Rank of symbol in half-open interval [0..pos)
        int64_t rank(int64_t pos, char symbol) const{
            return 0; // TODO
        }
};

class BitMagic{ // TODO Simon

    public:
        BitMagic(const vector<char>& seq){
            // TODO
        } 

        // Rank of symbol in half-open interval [0..pos)
        int64_t rank(int64_t pos, char symbol) const{
            return 0; // TODO
        }
};

struct Query{
    int64_t pos;
    char symbol;
};

vector<Query> generate_queries(int64_t n_queries, int64_t sequence_length, int64_t alphabet_size){
    vector<Query> queries(n_queries);
    for(int64_t i = 0; i < n_queries; i++){
        queries[i] = {rand()%sequence_length, (char)(rand()%alphabet_size)};
    }
    return queries;
}

template<typename rank_support_t>
void benchmark(const rank_support_t& rank_support, const vector<Query>& queries){

    int64_t sum = 0; // Sum ranks to prevent the compiler optimizing everything away
    int64_t t0 = current_time_micros();
    for(const Query& Q : queries){
        sum += rank_support.rank(Q.pos, Q.symbol);
    }
    int64_t t1 = current_time_micros();

    cout << "Time: " << (double)(t1-t0) / queries.size() << " microseconds per query (sum " << sum << ")" << endl;
}

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Please give a plain-subsetwt sbwt file as input" << endl;
        return 1;
    }

    string sbwt_index_file = string(argv[1]);

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-subsetwt"){
        cerr << "Error: input is not a plain-subsetwt SBWT" << endl;
    }

    sbwt::plain_sswt_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain SubsetWT with " << sbwt.number_of_subsets() << " subsets" << endl;
    cerr << "Extracting sequences from sswt nodes" << endl;            
    auto sswt = sbwt.get_subset_rank_structure();

    vector<char> root, left_child, right_child;

    for(auto x : sswt.ACGT_wt) root.push_back(x - '0');
    for(auto x : sswt.AC_wt) left_child.push_back(x - '1'); // Symbol '0' does not exist in children
    for(auto x : sswt.GT_wt) right_child.push_back(x - '1'); // Symbol '0' does not exist in children

    cerr << "Generating queries" << endl;
    srand(1234);
    int64_t n_queries = 1e6;
    vector<Query> base4_queries_root = generate_queries(n_queries, root.size(), 4);
    vector<Query> base3_queries_left = generate_queries(n_queries, left_child.size(), 3);
    vector<Query> base3_queries_right = generate_queries(n_queries, right_child.size(), 3);

    cerr << "Building structures for the base 4 sequence" << endl;
    Basic_WT basic_wt(root);
    Base4RRR rrr_generalization(root);
    SplitStructure split_structure(root);
    BitMagic bit_magic(root);

    cout << rrr_generalization.size_report() << endl;
    cout << sdsl::size_in_bytes(basic_wt.wt) * 8 << endl;

    cerr << "Benchmarking basic WT" << endl;
    benchmark(basic_wt, base4_queries_root);
    cerr << "Benchmarking RRR generalization" << endl;
    benchmark(rrr_generalization, base4_queries_root);
    cerr << "Benchmarking Split structure" << endl;
    benchmark(split_structure, base4_queries_root);
    cerr << "Benchmarking Bit magic" << endl;
    benchmark(bit_magic, base4_queries_root);

    cerr << "Building structures for the base 3 sequence in the left child" << endl;
    Basic_WT basic_wt_left(left_child);
    Base4RRR rrr_generalization_left(left_child);
    SplitStructure split_structure_left(left_child);
    BitMagic bit_magic_left(left_child);

    cout << rrr_generalization_left.size_report() << endl;
    cout << sdsl::size_in_bytes(basic_wt_left.wt) * 8 << endl;    

    cerr << "Benchmarking basic WT" << endl;
    benchmark(basic_wt_left, base3_queries_left);
    cerr << "Benchmarking RRR generalization" << endl;
    benchmark(rrr_generalization_left, base3_queries_left);
    cerr << "Benchmarking Split structure" << endl;
    benchmark(split_structure_left, base3_queries_left);
    cerr << "Benchmarking Bit magic" << endl;
    benchmark(bit_magic_left, base3_queries_left);

    cerr << "Building structures for the base 3 sequence in the right child" << endl;
    Basic_WT basic_wt_right(right_child);
    Base4RRR rrr_generalization_right(right_child);
    SplitStructure split_structure_right(right_child);
    BitMagic bit_magic_right(right_child);

    cout << rrr_generalization_right.size_report() << endl;
    cout << sdsl::size_in_bytes(basic_wt_right.wt) * 8 << endl;    

    cerr << "Benchmarking basic WT" << endl;
    benchmark(basic_wt_right, base3_queries_right);
    cerr << "Benchmarking RRR generalization" << endl;
    benchmark(rrr_generalization_right, base3_queries_right);
    cerr << "Benchmarking Split structure" << endl;
    benchmark(split_structure_right, base3_queries_right);
    cerr << "Benchmarking Bit magic" << endl;
    benchmark(bit_magic_right, base3_queries_right);

}
