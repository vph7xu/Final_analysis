#############################################
# Final‑analysis Makefile                   #
#                                           #
# • Builds the executable `raw_asymmetry`   #
# • Requires CERN ROOT and g++ ≥ C++17      #
#############################################

# ---- ROOT flags ----
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

# ---- Compiler & flags ----
CXXFLAGS  := -std=c++17 -O2 -Iinclude $(ROOTCFLAGS)
LDFLAGS   := $(ROOTLIBS)
CXX       := g++

# ---- Source layout ----
SRC_DIR   := src
SRCS      := $(wildcard $(SRC_DIR)/*.cpp)
OBJS      := $(SRCS:.cpp=.o)
MAIN      := main.cpp
EXEC      := raw_asymmetry

# ---- Default target ----
all: $(EXEC)

# Link step
$(EXEC): $(OBJS) $(MAIN)
	$(CXX) $(CXXFLAGS) $(MAIN) $(OBJS) -o $@ $(LDFLAGS)

# Compile each object file
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---- Convenience targets ----
#   make debug   → recompiles with -g -O0
#   make clean   → removes objects & executable

debug: CXXFLAGS += -g -O0
debug: clean all

clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean debug

