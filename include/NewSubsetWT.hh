#pragma once 


template<typename base4_rank_t, typename base3_rank_t>
class NewSubsetWT{

public:

    static const int64_t ROOT_LEFT = 2; // '10'
    static const int64_t ROOT_RIGHT = 1; // '01'
    static const int64_t CHILD_LEFT = 1; // '10'
    static const int64_t CHILD_RIGHT = 0; // '01'

    base4_rank_t root; 
    base3_rank_t left_child;
    base3_rank_t right_child;

    NewSubsetWT(){}

    NewSubsetWT(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits){
        vector<char> root_seq, left_seq, right_seq;
        int64_t n = A_bits.size();
        for(int64_t i = 0; i < n; i++){
            bool A = A_bits[i];
            bool C = C_bits[i];
            bool G = G_bits[i];
            bool T = T_bits[i];

            bool AC = A || C;
            bool GT = G || T;
            root_seq.push_back(AC*2 + GT);
            if(AC) left_seq.push_back(A*2 + C - 1); // -1: shift from [1,3] to [0,2]
            if(GT) right_seq.push_back(G*2 + T - 1);  // -1: shift from [1,3] to [0,2]
        }

        root = base4_rank_t(root_seq);
        left_child = base3_rank_t(left_seq);
        right_child = base3_rank_t(right_seq);

    }

    // Count of character c in subsets up to pos, not including pos
    int64_t rank(int64_t pos, char c) const{
        assert(c == 'A' || c == 'C' || c == 'G' || c == 'T');
        char root_sym = (c == 'A' || c == 'C') ? ROOT_LEFT : ROOT_RIGHT; // G and T go right
        char child_sym = (c == 'A' || c == 'G') ? CHILD_LEFT : CHILD_RIGHT; // C and T go right

        int64_t x = root.rankpair(pos, root_sym);

        if(root_sym == ROOT_LEFT) return left_child.rankpair(x, child_sym);
        else return right_child.rankpair(x, child_sym);
    }

    int64_t size_in_bytes() const{
        return root.size_in_bytes() + left_child.size_in_bytes() + right_child.size_in_bytes();
    }

    int64_t serialize(ostream& os) const{
        return 0; // TODO
    }

    void load(istream& is){
        return; // TODO
    }


};
