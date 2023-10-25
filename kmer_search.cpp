#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "SBWT.hh"
#include "variants.hh"
#include "throwing_streams.hh"
#include "RRR_generalization.hh"
#include "SplitStructure.hh"
#include "BitMagic.hh"
#include "SDSL_WT.hh"
#include "NewSubsetWT.hh"
#include "SeqIO.hh"

template<typename variant_t>
variant_t build_variant(const sbwt::plain_matrix_sbwt_t& input){
    variant_t var(
        input.get_subset_rank_structure().A_bits,
        input.get_subset_rank_structure().C_bits,
        input.get_subset_rank_structure().G_bits,
        input.get_subset_rank_structure().T_bits,
        input.get_streaming_support(),
        input.get_k(),
        input.number_of_kmers(),
        input.get_precalc_k()
    );
    return var;
}

template<typename variant_t> 
void benchmark(const variant_t& index, const string& index_name, const string& query_file, ostream& csv_out){

    sbwt::SeqIO::Reader<> reader(query_file);

    int64_t sum = 0;

    int64_t total_time_micros = 0;
    int64_t n_kmers = 0;

    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        int64_t t0 = current_time_micros();
	for(int64_t i = 0; i < len-index.get_k()+1; i++){
            sum += index.search(reader.read_buf + i);
            n_kmers++;
        }
        int64_t t1 = current_time_micros();
        total_time_micros += t1 - t0; 
    }

    cerr << "Index kmers and subsets " << index.number_of_kmers() << " " << index.number_of_subsets() << endl;
    cerr << "Time: " << (double)(total_time_micros) / n_kmers << " us/kmer" << endl;
    cerr << "Sum: " << sum << endl;

    int64_t space_bits = index.get_subset_rank_structure().size_in_bytes() * 8;
    space_bits += index.get_precalc().size() * sizeof(pair<int64_t, int64_t>) * 8; // Precalc
    csv_out << index_name << ", " << (double)(total_time_micros) / n_kmers << ", " << (double)space_bits / index.number_of_kmers() << endl;

    cerr << "Total k-mers queried: " << n_kmers << endl;
}

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Please give a plain-matrix sbwt file and a fastq query file as input" << endl;
        return 1;
    }

    string sbwt_index_file = string(argv[1]);
    string query_file = string(argv[2]);

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
    }

    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;

    typedef sdsl::wt_blcd<sdsl::bit_vector,
            sdsl::rank_support_v5<>,
            sdsl::select_support_scan<1>,
            sdsl::select_support_scan<0>> plain_wt_t; // No select sppport

    typedef sdsl::wt_blcd<sdsl::rrr_vector<>,
                          sdsl::rrr_vector<>::rank_1_type,
                          sdsl::rrr_vector<>::select_1_type,
                          sdsl::rrr_vector<>::select_0_type> rrr_wt_t;

    auto rrr_generalization = build_variant<sbwt::SBWT<NewSubsetWT<RRR_Generalization, RRR_Generalization>>>(sbwt);
    auto bitmagic = build_variant<sbwt::SBWT<NewSubsetWT<BitMagic, BitMagic>>>(sbwt);
    auto split = build_variant<sbwt::SBWT<NewSubsetWT<SplitStructure<4>, SplitStructure<3>>>>(sbwt);
    auto basic_wt = build_variant<sbwt::SBWT<NewSubsetWT<SDSL_WT<plain_wt_t>, SDSL_WT<plain_wt_t>>>>(sbwt);
    auto rrr_wt = build_variant<sbwt::SBWT<NewSubsetWT<SDSL_WT<rrr_wt_t>, SDSL_WT<rrr_wt_t>>>>(sbwt);

    ofstream results("kmer_search_results.csv");
    benchmark(sbwt, "plain-matrix", query_file, results);
    benchmark(rrr_generalization, "rrr-generalization", query_file, results);
    benchmark(basic_wt, "basic-wt", query_file, results);
    benchmark(rrr_wt, "rrr-wt", query_file, results);
    benchmark(bitmagic, "bitmagic", query_file, results);
    benchmark(split, "split", query_file, results);

}
