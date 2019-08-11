
# -----------------------------------------------------------------
#   Makefile for GCTA 
#   
#   Supported platforms: Unix / Linux
# ---------------------------------------------------------------------

# Directory of the target
OUTPUT = build/gcta64

# Source directoru
SOURCE = src

# Compiler
CXX = g++

# EIGEN library
EIGEN_PATH = ../eigen

# Intel MKL library
MKL_PATH = /opt/intel/mkl

# Compiler flags
CXXFLAGS = -w -m64 -g -static -fopenmp -o3 -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG -I $(MKL_PATH)/include 
LIB += -static -lz -Wl,--start-group  $(MKL_PATH)/lib/intel64/libmkl_intel_lp64.a $(MKL_PATH)/lib/intel64/libmkl_gnu_thread.a $(MKL_PATH)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -lcppunit

#EXCLUDE = UnusedCode.cpp read_imput.cpp gbat.cpp
#SRC := $(filter-out $(SOURCE)/$(EXCLUDE), $(wildcard $(SOURCE)/*.cpp))
SRC := $(wildcard $(SOURCE)/*.cpp)
OBJ := $(SRC:.cpp=.o)

HEADER := $(wildcard $(SOURCE)/*.h)

all: $(OUTPUT)

$(OUTPUT): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LIB) $(HEADER)

clean: 
	rm -f $(SOURCE)/*.o
	rm -f $(OUTPUT)
