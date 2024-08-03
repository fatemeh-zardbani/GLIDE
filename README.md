# GLIDE
Adaptive Spatial Indexing in Dynamic Environments 

This is an implementation of the proposed update module for spatial adaptive index, GLIDE, implemented on top of the AIR index.

## Dependencies

- `gcc`

## Usage

The index is defined and implemented in the `GLIDE_*.h` files. The runner files, `general_testing_*.cpp` includes the header, creates an instance of the tree, and tests its use. This runner file performs a workload with a mix of insertions and queries. every line in the workload file is regarding an action, which is determined by the first value in the line: 0 for queries, 1 for insertions, and -1 for deletions.
 This can be tuned using many flags at compile time, and input values at runtime. An example of the compilation command is as follows:  

```
g++ -D data_size6M -D versionGSR -D fo16 -D maxt64 -D SSL8 -D dh0 -D mt128 -mcmodel=large general_testing_glide_deletion.cpp -O3 -o glide
```

Then the running code will be:

```
./glide data.txt 100000 query.txt times.txt 100 2000 
```

## Compile time configurations

The flags that you can use at compile time to set parameters or perform other tasks, are as follows:

- `-D datasizeX`: from the values defined in the beginning of the header file, which are the values used in the experiments of the study, this determines the size of the pre-loaded data.
- `-D versoionX`: from the choices defined in the beginning of the runner file, you can choose which header file to include. This determines which version of the method you want to choose. The choices are: `GS`, `GSM`, `GSQ`, `GSS`, `CQ`, `CQM`, `CQripple`, `CS`, `CSM`, `GQ`, `GQM`.
- `-D foX`: from the values defined in the beginning of the header file, which are the values used in the experiments of the study, this parameter sets the fan-out value for the tree. This is the maximum number of children an internal node can have.
- `-D maxtX`: from the values defined in the beginning of the header file, which are the values used in the experiments of the study, this paarameter sets the maximum cracking threshold. This determines whether a leaf is regular or irregular and consequently whether it should be cracked or not.
- `-D SSLX`: from the values defined in the beginning of the header file, which are the values used in the experiments of the study, this parameter determines the maximum number of spares each node(internal and leaves) can hold.
- `-D dhX`: from the values defined in the beginning of the header file, which are the values used in the experiments of the study, this parameter determines the number of default holes a leaf will be given upon being "sling"ed.
- `-D mt128`: from the values defined in the beginning of the header file, which are the values used in the experiments of the study, this paramter determines whether a leaf should be cracked before being "sling"ed. In the study it is set as 128, which is double 64 the crackign threshold.
- `-mcmodel=large`: since we use compile-time allocated memory, when using the code for larger datasets, we need more space the standard. Using this flag allows the compilation to go through.
- `-O3`: compiler optimization.

## Runtime configurations

The parameters given at runtime are as follows:

- data file address: there should be enough data for the re-loading and all the insertions in this file. The data should be stored as follows for a 2D example: $minx miny maxx maxy$ 
- size of query workload: count of how many range queries 
- query file address
- time file address
- $\iota$: the frequency at which insertions will be made
- insertion count: how many insertions to make every $\iota$ query 
