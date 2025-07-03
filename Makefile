CXX = g++-15
CXXFLAGS = -O3 -std=c++17 -fopenmp

SRC_DIR = src
NUM_LIB = $(SRC_DIR)/numlib/*.cpp

TARGETS = apwd3-c1.x apwd3-c3.x apwd3-c4.x apwd3-cD.x apwd3-cE.x apwd3-benchmark.x

all: $(TARGETS)

apwd3-c1.x: $(NUM_LIB) $(SRC_DIR)/aPWD3_part_c1.cpp $(SRC_DIR)/main_part_c1.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

apwd3-c3.x: $(NUM_LIB) $(SRC_DIR)/aPWD3_part_c3.cpp $(SRC_DIR)/main_part_c3.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

apwd3-c4.x: $(NUM_LIB) $(SRC_DIR)/aPWD3_part_c4.cpp $(SRC_DIR)/main_part_c4.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

apwd3-cD.x: $(NUM_LIB) $(SRC_DIR)/aPWD3_part_cD.cpp $(SRC_DIR)/main_part_cD.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

apwd3-cE.x: $(NUM_LIB) $(SRC_DIR)/aPWD3_part_cE.cpp $(SRC_DIR)/main_part_cE.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

apwd3-benchmark.x: $(NUM_LIB) $(SRC_DIR)/aPWD3_part_c1.cpp $(SRC_DIR)/aPWD3_part_c3.cpp $(SRC_DIR)/aPWD3_part_c4.cpp $(SRC_DIR)/main_benchmark.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(TARGETS)

