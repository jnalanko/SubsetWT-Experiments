
#include "globals.hh"
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

int64_t current_time_micros(){
    return (std::chrono::duration_cast< microseconds >(high_resolution_clock::now().time_since_epoch())).count();
}

// Returns the number of bytes written
int64_t serialize_string(const string& S, ostream& out){
    int64_t size = S.size();
    out.write((char*)&size, sizeof(size));
    out.write(S.data(), size);
    return sizeof(size) + size;
}

string load_string(istream& in){
    int64_t size;
    in.read((char*)&size, sizeof(size));
    string S(size, '\0');
    in.read((char*)&S[0], size); // The C++ standard guarantees that std::string is stored contiguously in memory
    return S;
}