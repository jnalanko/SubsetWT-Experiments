microbenchmark:
	g++ main.cpp src/SeqIO.cpp src/globals.cpp ./sdsl-lite/build/lib/libsdsl.a -std=c++20 -I sdsl-lite/include/ -O3 -I include -I ./sdsl-lite/build/external/libdivsufsort/include/ -g -o benchmark -Wno-deprecated-declarations -march=native -DNDEBUG

kmer_search_benchmark:
	g++ kmer_search.cpp src/SeqIO.cpp src/globals.cpp ./sdsl-lite/build/lib/libsdsl.a -std=c++20 -I sdsl-lite/include/ -O3 -I include -I ./sdsl-lite/build/external/libdivsufsort/include/ -g -o kmer_search -Wno-deprecated-declarations -march=native -DNDEBUG

