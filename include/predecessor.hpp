#ifndef _PRED_H_
#define _PRED_H_

#include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <chrono>
#include <ctime>
#include <ratio>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <random>

using namespace std;
using namespace std::chrono;

class Predecessor {

public:
    Predecessor() {
    }

    Predecessor(const vector<uint64_t>& data) {
        _r = data.size();
        _gaps = new uint16_t[_r]; // allocates an array _gaps of _r ints, whose values are initially undefined
        _gaps[0] = data[0];
        for(uint64_t i=1;i<_r;i++){
            uint64_t g = data[i] - data[i-1];
            if(g > 65535){ // bigger than 16 bit integers
                cerr << "Uh oh: big gap\n";
            }
            _gaps[i] = (uint16_t)g;
        }

        uint32_t maxGap = 0;
        uint32_t minGap = 1<<31;
        _n = 0;
        for (size_t j = 0; j < _r; j++){
            if(maxGap < _gaps[j]){
                maxGap = _gaps[j];
            }
            if(minGap > _gaps[j]){
                minGap = _gaps[j];
            }
            _n += _gaps[j];
        }
        //build predecessor data structure for phrase starting positions
        _numpblocks = (_r / _pbs);
        if (_r % _pbs) {
            _numpblocks++;
        }
        _partialSums = new uint64_t[_numpblocks + 1];
        _partialSums[_numpblocks] = UINT64_MAX; //_n+1;
        uint64_t sum = 0;
        for (uint64_t i = 0; i < _r; i++) {
            if (i % _pbs == 0) {
                _partialSums[i / _pbs] = sum;
            }
            sum += _gaps[i];
        }

        _msbn = __builtin_clzll(_n); // Returns the number of leading 0-bits in _n, starting at the most significant bit position. If _n is 0, the result is undefined.
        _jtShift = (64 - _msbn) - _logjtSize;
        _jumpTable = new uint32_t[_jtSize];

        for(uint32_t j=0; j<_jtSize; j++){
            _jumpTable[j] = UINT32_MAX; //was useful for debugging
        }

        _jumpTable[0] = 0;
        uint32_t jtp = 0;
        uint32_t prevKey = 0;
        uint32_t run = 1, maxRun = 0;
        for(uint32_t i=1; i<_numpblocks; i++){
            uint32_t key = getJumpTableKey(_partialSums[i]);
            if(key != prevKey){
                if(run > maxRun) maxRun = run;
                run = 0;
                while(jtp < key){
                    _jumpTable[jtp++] = _jumpTable[prevKey];
                }
                _jumpTable[key] = i;
            }
            //jtp++;
            prevKey = key;
            run++;
        }
        while(jtp < _jtSize){
            _jumpTable[jtp++] = _jumpTable[prevKey];
        }
//       cerr << "maxRun (in jt): " << maxRun << '\n';
    }

    Predecessor(const Predecessor& other){
        assert(&other != this); // What on earth are you trying to do?
        operator=(other);
    }

    Predecessor& operator =(const Predecessor& other) {
        if(&other != this){
            this->_n = other._n;   //length of BWT and text
            this->_r = other._r;    //number of runs
            this->_pbs = other._pbs;
            this->_numpblocks = other._numpblocks;
            this-> _partialSums = new uint64_t[_numpblocks + 1];
            std::copy(other._partialSums, other._partialSums + _numpblocks + 1, _partialSums);
            this-> _gaps = new uint16_t[_r];
            std::copy(other._gaps, other._gaps + _r, _gaps);
            this->_msbn = other._msbn;
            this->_jtShift = other._jtShift;
            this->_jtSize = other._jtSize;
            this-> _jumpTable = new uint32_t[_jtSize];
            std::copy(other._jumpTable, other._jumpTable + _jtSize, _jumpTable);
            this->_logjtSize = other._logjtSize;
            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

    ~Predecessor() {
       delete [] _partialSums;
        delete [] _gaps;
        delete [] _jumpTable;
       _partialSums = nullptr;
        _gaps = nullptr;
        _jumpTable = nullptr;
    }

    inline uint32_t getJumpTableKey(uint64_t index){
        return (index >> _jtShift);
    }

    //outputs a pair <p,bool>
    // p is the index of the predecessor in the set
    // bool is 1 if the value at index p is equal to key, else 0
    pair<int64_t,bool> inline getPred(int64_t key){
        if(key < _partialSums[0]){
            return {-1,0};
        }
        uint64_t *x = std::lower_bound(_partialSums, _partialSums + _numpblocks , key); // Returns an iterator pointing to the first element in the range [first, last) that does not satisfy element < value
        if (x - _partialSums) { //x > _partialSum start
            x--; //x now pointing at predecessor in _partialSums
        }
        uint64_t sum = *x;
        uint64_t p = _pbs * (x - _partialSums);
        while (p < _r) { // sum < key &&
            uint64_t len = _gaps[p];
            if (sum + len > key) break;
            p++;
            sum += len;
        }
        return {p-1, (sum==key)};
    }

    pair<int64_t,bool> inline getPredWithJumpTable(uint64_t pkey){
        if(pkey < _partialSums[0]){
            return {-1,0};
        }
        if (pkey > _partialSums[_numpblocks]){
            pkey = _partialSums[_numpblocks];
        }
        uint32_t jkey = getJumpTableKey(pkey);
        uint64_t *x = _partialSums + _jumpTable[jkey];
        while(*x < pkey){
            x++;}
        if (x - _partialSums) {
            x--; //x now pointing at predecessor in _partialSums
        }
        uint64_t sum = *x;
        int64_t p = _pbs * (x - _partialSums); //should surely be replaced by a shift (?)
        while (p < _r) {
            uint64_t len = _gaps[p];
            if (sum + len > pkey) { break; }
            p++;
            sum += len;
        }
        return {p-1, (sum == pkey)};
    }

    size_t getn(){
        return _n;
    }

    uint64_t sizeInBytes(){
        uint64_t sz = 0;
        sz += _r*sizeof(uint16_t);
        sz += _numpblocks*sizeof(uint64_t);
        return sz;
    }

private:

    inline uint64_t getGap(uint64_t index) {
        return (_gaps[index]);
    }

    size_t _n = 0;   //length of BWT and text
    size_t _r = 0;   //number of runs

    uint32_t _pbs = 64; //predecessor block size
    //uint32_t _pbs = 8; //predecessor block size
    uint32_t _numpblocks; //number of blocks in predecessor structure
    uint64_t *_partialSums; //every _pbs^th
    uint16_t *_gaps;

    uint32_t _msbn;
    uint32_t _jtShift;
    uint32_t *_jumpTable;
    //uint32_t _jtSize = 65536;
    //uint32_t _logjtSize = 16;
    uint32_t _jtSize = 131072;
    uint32_t _logjtSize = 17;
    //uint32_t _jtSize = 262144;
    //uint32_t _logjtSize = 18;

//    double _totalPSscanning = 0;
};

#endif
