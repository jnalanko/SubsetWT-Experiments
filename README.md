# Building

```
git submodule init
git submodule update
cd sdsl-lite
cd build
cmake .. -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc)
make
cd ../..
make
```

This creates an executable called `benchmark`.

# Running

The code takes in a plain-subsetwt SBWT file. There is one small example containing 3 E. coli genomes at `index.sswt`. To benchmark on it, run:

```
./benchmark index.sswt
```