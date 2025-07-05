CXX = g++-15
CXXFLAGS = -O3 -std=c++17 -fopenmp

SRC_DIR := src
BUILD_DIR := build

NUM_LIB_SRCS := $(wildcard $(SRC_DIR)/numlib/*.cpp)
NUM_LIB_OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(NUM_LIB_SRCS))

TARGETS := apwd3-c1.x apwd3-c3.x apwd3-c4.x apwd3-cD.x apwd3-cE.x apwd3-benchmark.x
TARGETS_IN_BUILD := $(addprefix $(BUILD_DIR)/, $(TARGETS))

.PHONY: all clean

all: $(TARGETS_IN_BUILD)

$(BUILD_DIR)/apwd3-c1.x: $(NUM_LIB_OBJS) $(BUILD_DIR)/aPWD3_part_c1.o $(BUILD_DIR)/main_part_c1.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/apwd3-c3.x: $(NUM_LIB_OBJS) $(BUILD_DIR)/aPWD3_part_c3.o $(BUILD_DIR)/main_part_c3.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/apwd3-c4.x: $(NUM_LIB_OBJS) $(BUILD_DIR)/aPWD3_part_c4.o $(BUILD_DIR)/main_part_c4.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/apwd3-cD.x: $(NUM_LIB_OBJS) $(BUILD_DIR)/aPWD3_part_cD.o $(BUILD_DIR)/main_part_cD.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/apwd3-cE.x: $(NUM_LIB_OBJS) $(BUILD_DIR)/aPWD3_part_cE.o $(BUILD_DIR)/main_part_cE.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/apwd3-benchmark.x: $(NUM_LIB_OBJS) $(BUILD_DIR)/aPWD3_part_c1.o $(BUILD_DIR)/aPWD3_part_c3.o $(BUILD_DIR)/aPWD3_part_c4.o $(BUILD_DIR)/main_benchmark.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@echo "Cleaning up the build directory..."
	@rm -rf $(BUILD_DIR)
	@echo "Done."