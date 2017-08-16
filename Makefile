simname := surface_decoder
so_name := square_octagonal_mapping
dec_name:= square_octagonal_decoder

blossom := ./blossom/blossom5

std_mat := ./src/standard_matching.cpp
rn_gen  := ./src/rng-omp.cpp
oct_map := ./src/square_octagonal.cpp
sq_dec  := ./src/sq_oct_matching.cpp


CXX := g++
CXXFLAGS := -std=c++0x -I./include/ -fopenmp -lgsl -lgslcblas

srcfiles := $(std_mat) $(rn_gen)
sqoctfiles := $(oct_map) $(rn_gen)
objects := $(blossom) $(simname) $(so_name) $(dec_name)

decoder: $(objects)

$(blossom):
	make -C blossom

$(simname): $(srcfiles)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(simname) $(srcfiles) $(LDLIBS)

clean:
	make -C blossom clean
	rm -f $(objects)

sqoct: $(blossom)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(so_name) $(sqoctfiles) $(LDLIBS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(dec_name) $(rn_gen) $(sq_dec) $(LDLIBS)

