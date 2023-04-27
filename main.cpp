#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "SBWT.hh"
#include "variants.hh"
#include "throwing_streams.hh"
#include "RRR_generalization.hh"
#include "SplitStructure.hh"
#include "BitMagic.hh"
#include "SDSL_WT.hh"

struct Query{
    int64_t pos;
    char symbol;
};

vector<Query> generate_queries(int64_t n_queries, int64_t sequence_length, const vector<char>& query_alphabet){
    vector<Query> queries(n_queries);
    for(int64_t i = 0; i < n_queries; i++){
        queries[i] = {rand()%sequence_length, query_alphabet[rand() % query_alphabet.size()]};
    }
    return queries;
}

template<typename rank_support_t>
double benchmark_single_rank(const rank_support_t& rank_support, const vector<Query>& queries){

    int64_t sum = 0; // Sum ranks to prevent the compiler optimizing everything away
    int64_t t0 = current_time_micros();
    for(const Query& Q : queries){
        sum += rank_support.rank(Q.pos, Q.symbol);
    }
    int64_t t1 = current_time_micros();

    cerr << "Single rank final sum: " << sum << endl;
    cerr << "Single rank ns/query: " << (double)(t1-t0) / queries.size() * 1000 << endl;

    return (double)(t1-t0) / queries.size() * 1000; // Nanoseconds;
    
}

template<typename rank_support_t>
double benchmark_rankpair(const rank_support_t& rank_support, const vector<Query>& queries){

    int64_t sum = 0; // Sum ranks to prevent the compiler optimizing everything away
    int64_t t0 = current_time_micros();
    for(const Query& Q : queries){
        sum += rank_support.rankpair(Q.pos, Q.symbol);
    }
    int64_t t1 = current_time_micros();

    cerr << "Rank pair final sum: " << sum << endl;
    cerr << "Rank pair ns/query: " << (double)(t1-t0) / queries.size() * 1000 << endl;

    return (double)(t1-t0) / queries.size() * 1000; // Nanoseconds

}

template<typename rank_support_t>
void benchmark_structure(rank_support_t& rank_support, const string& name, const vector<Query>& single_queries, const vector<Query>& pair_queries, int64_t seq_size, ostream& csv_out){

    cerr << "Benchmarking " << name << endl;

    csv_out << name << ", ";
    csv_out << benchmark_single_rank(rank_support, single_queries) << ", ";
    csv_out << benchmark_rankpair(rank_support, pair_queries) << ", ";
    csv_out << rank_support.size_in_bytes()*8.0 / seq_size << endl;

    cerr << "bits/symbol: " << rank_support.size_in_bytes()*8.0 / seq_size << endl;

}

template<int64_t alphabet_size>
void benchmark_all(const vector<char>& seq, const vector<Query>& single_queries, const vector<Query>& pair_queries, const string& csv_outfile){

    typedef sdsl::wt_blcd<sdsl::bit_vector,
            sdsl::rank_support_v5<>,
            sdsl::select_support_scan<1>,
            sdsl::select_support_scan<0>> plain_wt_t; // No select sppport

    typedef sdsl::wt_blcd<sdsl::rrr_vector<>,
                          sdsl::rrr_vector<>::rank_1_type,
                          sdsl::rrr_vector<>::select_1_type,
                          sdsl::rrr_vector<>::select_0_type> rrr_wt_t; // No select support

    cerr << "Building structures" << endl;
    SDSL_WT<plain_wt_t> sdsl_plain_wt(seq);
    SDSL_WT<rrr_wt_t> sdsl_rrr_wt(seq);
    RRR_Generalization rrr_generalization(seq);
    SplitStructure<alphabet_size> split_structure(seq);
    BitMagic bit_magic(seq);

    ofstream out(csv_outfile);

    // Write csv header
    out << "structure, single_rank (ns/query), rank_pair (ns/query), space (bits/symbol)" << endl;

    benchmark_structure(sdsl_plain_wt, "plain-wt", single_queries, pair_queries, seq.size(), out);;
    benchmark_structure(sdsl_rrr_wt, "rrr-wt", single_queries, pair_queries, seq.size(), out);
    benchmark_structure(rrr_generalization, "rrr-generalization", single_queries, pair_queries,  seq.size(), out);
    benchmark_structure(split_structure, "split", single_queries, pair_queries, seq.size(), out);
    benchmark_structure(bit_magic, "bitmagic", single_queries, pair_queries, seq.size(), out);
}

// Returns root, left child, right child
tuple<vector<char>,vector<char>,vector<char>> get_sequences(const sbwt::plain_matrix_sbwt_t& sbwt){
    vector<char> root, left, right;
    int64_t n = sbwt.number_of_subsets();
    for(int64_t i = 0; i < n; i++){
        bool A = sbwt.get_subset_rank_structure().A_bits[i];
        bool C = sbwt.get_subset_rank_structure().C_bits[i];
        bool G = sbwt.get_subset_rank_structure().G_bits[i];
        bool T = sbwt.get_subset_rank_structure().T_bits[i];

        bool AC = A || C;
        bool GT = G || T;
        root.push_back(AC*2 + GT);
        if(AC) left.push_back(A*2 + C - 1); // -1: shift from [1,3] to [0,2]
        if(GT) right.push_back(G*2 + T - 1);  // -1: shift from [1,3] to [0,2]
    }
    return {root, left, right};
}

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Please give a plain-matrix sbwt file as input" << endl;
        return 1;
    }

    string sbwt_index_file = string(argv[1]);

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
    }

    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;
    cerr << "Building Base4 and Base3 sequences from SBWT" << endl;            

    auto sswt = sbwt.get_subset_rank_structure();

    vector<char> root, left_child, right_child;
    std::tie(root, left_child, right_child) = get_sequences(sbwt);

    cerr << "Generating queries" << endl;

    srand(1234);
    int64_t n_queries = 1e7;

    // The root has alphabet 0 = '00', 1 = '01', 2 = '10', 3 = '11'
    vector<Query> base4_single_queries_root = generate_queries(n_queries, root.size(), {0,1,2,3}); // '00', '01', '10', '11'
    vector<Query> base4_pair_queries_root = generate_queries(n_queries, root.size(), {1,2}); // '01', '10'

    // The children have alphabet 0 = '01', 1 = '10', 2 = '11'
    vector<Query> base3_single_queries_left = generate_queries(n_queries, left_child.size(), {0,1,2}); // '01', '10', '11'
    vector<Query> base3_pair_queries_left = generate_queries(n_queries, left_child.size(), {0,1}); // '01', '10'

    vector<Query> base3_single_queries_right = generate_queries(n_queries, right_child.size(), {0,1,2}); // '01', '10', '11'
    vector<Query> base3_pair_queries_right = generate_queries(n_queries, right_child.size(), {0,1}); // '01', '10'

    cerr << "Running queries for the root (base 4)" << endl;
    benchmark_all<4>(root, base4_single_queries_root, base4_pair_queries_root, "root_results.csv");

    cerr << "Running queries for the left child (base 3)" << endl;
    benchmark_all<3>(left_child, base3_single_queries_left, base3_pair_queries_left, "left_child_results.csv");

    cerr << "Running queries for the right child (base 3)" << endl;
    benchmark_all<3>(right_child, base3_single_queries_right, base3_pair_queries_right, "right_child_results.csv");

}
