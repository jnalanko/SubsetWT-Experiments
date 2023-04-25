all:
	g++ main.cpp ./sdsl-lite/build/lib/libsdsl.a -std=c++20 -I sdsl-lite/include/ -I include -I ./sdsl-lite/build/external/libdivsufsort/include/ -g -o benchmark -Wno-deprecated-declarations -O3

Base4RRR:
	g++ Base4RRR_main.cpp ./sdsl-lite/build/lib/libsdsl.a -std=c++20 -I sdsl-lite/include/ -I include -I ./sdsl-lite/build/external/libdivsufsort/include/ -g -o Base4RRR -Wno-deprecated-declarations
