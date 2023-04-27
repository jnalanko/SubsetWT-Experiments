#pragma once

#include <vector>
#include <iostream>
#include "globals.hh"
#include "SBWT.hh"
using namespace std;

class BitMagic{

public:

    BitMagic(){
        _bits = new uint64_t[0];
    }

    BitMagic(const vector<char>& seq){
        _n = seq.size();
        uint64_t symsPerWord = 32; // seems to be unnecessary
        uint64_t nblocks = _n/_b + 1;
        //_b is the number of symbols per block
        _N = nblocks*(2*_b + 128);//length in bits of data structure; +128 is for 4*32-bit ints per block of _b 2-bit symbols
        _bits = new uint64_t[_N/64];
        uint32_t *psums = new uint32_t[4];
        psums[0] = psums[1] = psums[2] = psums[3] = 0;
        _sigma = 0;
        for(uint64_t i=0;i<_n;i++){
            if(!(psums[seq[i]])) {_sigma++;}
            psums[seq[i]]++;
        }
        //we now know the alphabet (_sigma) which is assumed to be either 3 or 4
        psums[0] = psums[1] = psums[2] = psums[3] = 0;
        uint64_t bi = 0;
        for(uint64_t i=0;i<_n;){
            if(i%_b == 0){
                //we are at the start of a block - set its prefix sums
                ((uint32_t *)(_bits+bi))[0] = psums[0];
                ((uint32_t *)(_bits+bi))[1] = psums[1];
                ((uint32_t *)(_bits+bi))[2] = psums[2];
                ((uint32_t *)(_bits+bi))[3] = psums[3];
                bi+=2; //move past the words containing the prefix sums for this block
            }
            uint64_t j=0;
            uint64_t w=0;
            while(j<32 && (i+j)<_n){
                uint8_t sym = seq[i+j];
                if(_sigma == 3) sym++; //for ternary sequences, remap the alphabet from {0,1,2} to {1,2,3}
                psums[sym]++;
                w = w | (((uint64_t)sym) << (2*j));
                j++;
            }
            _bits[bi] = w;
            bi++;
            i+=j;
        }
        delete [] psums;
    }

    BitMagic(const BitMagic& other){
        assert(&other != this); // What on earth are you trying to do?
        operator=(other);
    }

    BitMagic& operator =(const BitMagic& other) {
        if(&other != this){
            this->_b = other._b; //number of symbols per block
            this->_logb = other._logb;
            this->_n = other._n;
            this->_N = other._N;
            this->_sigma = other._sigma;

            this-> _bits = new uint64_t[_N/64];
            std::copy(other._bits, other._bits + _N/64, _bits);
            //this->_MASK = other._MASK; //??
            std::copy(other._MASK, other._MASK + 4, _MASK);
            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

    ~BitMagic(){
        if (_bits != nullptr){
            delete [] _bits;
            _bits = nullptr;
            //delete [] _MASK;
        }
    }


    size_t size_in_bytes() const{
        size_t sz = 0;
        sz += (sizeof(uint64_t)*(_N/64)); //_bits
        sz += (sizeof(uint64_t)*4); //_b, _logb, _n, _N
        sz += (sizeof(uint64_t*)); //_bits's pointer
        return sz;
    }

    //Useful during debugging
    void print64bitword(uint64_t w) const{
        for(uint64_t i=0; i<64; i++){
            cout << ((w>>(63-i))&1);
        }
        cout << '\n';
    }

    inline void count1110(uint64_t w, uint64_t &countpA, uint64_t &countpB) const{
        uint64_t w2 = w & 0xAAAAAAAAAAAAAAAA; //popcnt(w2) == # of 11s + # of 10s
        uint64_t w3 = (w2 >> 1);
        w3 = w3 & w;
        countpA = __builtin_popcountll(w3);
        countpB = __builtin_popcountll(w2) - countpA;
    }


    inline void count1110Partial(uint64_t w, uint64_t &countpA, uint64_t &countpB, uint32_t shift) const{
        uint64_t w2 = w & 0xAAAAAAAAAAAAAAAA; //popcnt(w2) == # of 11s + # of 10s
        uint64_t w3 = (w2 >> 1);
        w3 = w3 & w;
        countpA = __builtin_popcountll((w3 << shift)>>(shift-2));
        countpB = __builtin_popcountll((w2 << shift)>>(shift-2)) - countpA;
    }

    //return the sum of the 11s and (10s or 01s, depending on sym)
    inline int64_t rankpairBase4Word(uint64_t w, uint8_t sym) const{
        w = w & _MASK[sym-1];
        return __builtin_popcountll(w);
    }
    inline int64_t rankpairBase4WordPartial(uint64_t w, uint8_t sym, uint32_t shift) const{
        w = w & _MASK[sym-1];
        return __builtin_popcountll((w << shift)>>(shift-2));
    }

    //TODO: currently only correct for Base-4 sequences
    int64_t rankpair(int64_t pos, char symbol) const{
        uint64_t sym = (uint64_t)symbol;
        if(_sigma == 3){
            sym++;
        }
        //return (rank(pos,symbol) + rank(pos,3));
        uint64_t blockstart = (pos>>_logb)*(2*_b + 128)/64;
        //blockstart is the word offset of the start of the block containing position i
        //2*_b is the number of bits from symbols in a block, because there are _b items per block and 2 bits per symbol
        //128 = 4*32 is the number of bits needed for the preblock ranks for each symbol
        uint64_t preBlockRank = ((uint32_t *)(_bits+blockstart))[sym]; //retrieve the appropriate preblock rank
        preBlockRank += ((uint32_t *)(_bits+blockstart))[3]; //retrieve the appropriate preblock rank
        uint64_t *blockwords = _bits + blockstart + 2; //move past the prefix sums
        uint64_t blocki = (pos&(_b-1))/32; //index of word in this block containing the query position

        uint64_t wholeWordRank = 0, leftOverRank = 0;

        for(uint64_t i=0; i<blocki; i++){
            uint64_t w = blockwords[i];
            wholeWordRank += rankpairBase4Word(w,sym);
        }
        if(pos%32){ //possibly inspect part of the next word
            uint64_t w = blockwords[blocki];
            //compute an appropriate shift
            uint32_t shift = 64 - 2*(pos%32); //pos%32 is never 0 inside this if statement
            leftOverRank = rankpairBase4WordPartial(w,sym,shift);
        }
        return preBlockRank + wholeWordRank + leftOverRank;
    }

    // Rank of symbol in half-open interval [0..pos)
    int64_t rank(int64_t pos, char symbol) const{
        uint64_t sym = (uint64_t)symbol;
        uint64_t blockstart = (pos>>_logb)*(2*_b + 128)/64;
        //blockstart is the word offset of the start of the block containing position i
        //2*_b is the number of bits from symbols in a block, because there are _b items per block and 2 bits per symbol
        //128 = 4*32 is the number of bits needed for the preblock ranks for each symbol
        uint64_t preBlockRank = ((uint32_t *)(_bits+blockstart))[sym]; //retrieve the appropriate preblock rank
        uint64_t *blockwords = _bits + blockstart + 2; //move past the prefix sums
        uint64_t blocki = (pos&(_b-1))/32; //index of word in this block containing the query position

        uint64_t countpA = 0, countpB = 0, wholeWordRank = 0, leftOverRank = 0;
        if(sym < 2){
            for(uint64_t i=0; i<blocki; i++){
                uint64_t w = blockwords[i];
                w = ~w;
                count1110(w,countpA,countpB);
                uint64_t count = countpB;
                bool x = (sym == 0 || sym == 3);
                if(x) count = countpA;
                wholeWordRank += count;
            }
            if(pos%32){ //possibly inspect part of the next word
                uint64_t w = blockwords[blocki];
                //compute an appropriate shift
                uint32_t shift = 64 - 2*(pos%32); //pos%32 is never 0 inside this if statement
                w = ~w;
                count1110Partial(w,countpA,countpB,shift);
                uint64_t count = countpB;
                bool x = (sym == 0 || sym == 3);
                if(x) count = countpA;
                leftOverRank = count;
            }
        }
        else{ //sym >= 2
            for(uint64_t i=0; i<blocki; i++){
                uint64_t w = blockwords[i];
                count1110(w,countpA,countpB);
                uint64_t count = countpB;
                bool x = (sym == 0 || sym == 3);
                if(x) count = countpA;
                wholeWordRank += count;
            }
            if(pos%32){ //possibly inspect part of the next word
                uint64_t w = blockwords[blocki];
                //compute an appropriate shift
                uint32_t shift = 64 - 2*(pos%32); //pos%32 is never 0 inside this if statement
                count1110Partial(w,countpA,countpB,shift);
                uint64_t count = countpB;
                bool x = (sym == 0 || sym == 3);
                if(x) count = countpA;
                leftOverRank = count;
            }
        }

        return preBlockRank + wholeWordRank + leftOverRank;
    }

    uint64_t _b = 512; //number of symbols per block
    uint64_t _logb = 9;
    uint64_t _n;
    uint64_t _N;
    uint8_t _sigma;
    uint64_t *_bits;
    uint64_t _MASK[4] = {0x5555555555555555,0xAAAAAAAAAAAAAAAA,0,0};
};
