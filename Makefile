# NOTE: needs boost, tclap, STXXL, KMC, and sdsl

CXX=g++
CPP_FLAGS=-pipe -m64 -std=c++14  -W -Wall -Wextra -Wpointer-arith \
					-Wunused -Wwrite-strings #-pedantic-errors -Wcast-qual #-Werror

CLANG_WARNINGS=-Wbool-conversions -Wshift-overflow -Wliteral-conversion # CLANG ONLY
BOOST_PATH=../3rd_party_inst# #/s/chopin/h/proj/soma/boost_1_54

COSMO_PATH=../
DEP_PATH=$(COSMO_PATH)/3rd_party_inst
KMC_PATH=$(COSMO_PATH)/3rd_party_src/KMC
#RANKSELECT_PATH=./3rd_party_inst/include/rankselect
INC=-isystem $(DEP_PATH)/include -isystem $(BOOST_PATH)/include -isystem ./Sux -isystem $(COSMO_PATH) -isystem include
LIB=-L$(DEP_PATH)/lib -L./ -L$(BOOST_PATH)/lib
BOOST_FLAGS= -lboost_system -lboost_filesystem
DEP_FLAGS=$(INC) $(LIB) $(BOOST_FLAGS) -isystem $(KMC_PATH) -lsdsl -fopenmp #openmp is needed for logging
DEBUG_FLAGS=-pg -g -gstabs
NDEBUG_FLAGS=-DNDEBUG
OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing
NOPT_FLAGS=-O0

# Using Semantic Versioning: http://semver.org/
VERSION=0.7.0
BANNER='Copyright Alex Bowe (c) 2016'
CPP_FLAGS+=-DVERSION=\"$(VERSION)\" -DBANNER=\"$(BANNER)\"

k?=32
# TODO: quantize k
CPP_FLAGS+=-DK_LEN=$(k)

colors?=6000
CPP_FLAGS+=-DNUM_COLS=$(colors)

ifeq ($(asm),1)
CPP_FLAGS+=-S -fverbose-asm
endif

ifeq ($(optimise),0)
CPP_FLAGS+=$(NOPT_FLAGS)
else
CPP_FLAGS+=$(OPT_FLAGS)
endif

ifeq ($(debug),1)
CPP_FLAGS+=$(DEBUG_FLAGS)
else
CPP_FLAGS+=$(NDEBUG_FLAGS)
endif

ifeq ($(verbose),1)
CPP_FLAGS+=-DVERBOSE
endif

KMC_OBJS=$(KMC_PATH)/kmc_api/kmc_file.o $(KMC_PATH)/kmc_api/kmer_api.o $(KMC_PATH)/kmc_api/mmer.o
RANKSELECT_OBJS=$(RANKSELECT_PATH)/bitmap.o
BUILD_REQS=$(COSMO_PATH)/lut.hpp $(COSMO_PATH)/debug.hpp $(COSMO_PATH)/utility.hpp $(COSMO_PATH)/io.hpp $(COSMO_PATH)/sort.hpp $(COSMO_PATH)/kmer.hpp $(COSMO_PATH)/dummies.hpp $(COSMO_PATH)/debruijn_graph.hpp $(COSMO_PATH)/debruijn_graph_shifted.hpp $(COSMO_PATH)/pack-color.hpp $(COSMO_PATH)/cosmo-color-pd.hpp $(COSMO_PATH)
COLOR_REQS=$(COSMO_PATH)/colored_debruijn_graph.hpp $(COSMO_PATH)/io.hpp $(COSMO_PATH)/debug.hpp
RANK9SEL_SRC=./Sux/*.cpp
RB_VEC_SRC=rb-vec.cpp rb-vec.hpp bit_array.c $(RANK9SEL_SRC)
RB_QUERY_SRC=rb-query.cpp $(RB_VEC_SRC)
RB_FS_SRC=rb-filesystem.cpp

BINARIES=rb-pack-color rb-find-bubble rb-validate

default: all

# Stores last compiled options to know whether we need to recompile when comp variables change
.PHONY: force
compiler_flags: force
	@echo $(CPP_FLAGS) | cmp -s - $@ || echo $(CPP_FLAGS) > $@

rb-pack-color: rb-pack-color.cpp $(RB_VEC_SRC) $(BUILD_REQS) compiler_flags $(RB_FS_SRC)
	$(CXX) $(CPP_FLAGS) rb-pack-color.cpp $(RB_VEC_SRC) $(RB_FS_SRC) -o $@ $(KMC_OBJS) $(DEP_FLAGS)

rb-find-bubble: rb-find-bubble.cpp $(RB_QUERY_SRC) $(BUILD_REQS) compiler_flags $(RB_FS_SRC)
	$(CXX) $(CPP_FLAGS) $(RB_QUERY_REQS) -o $@ rb-find-bubble.cpp $(RB_QUERY_SRC) $(RB_FS_SRC) $(KMC_OBJS) $(DEP_FLAGS)

rb-validate: rb-validate.cpp $(RB_QUERY_SRC) $(BUILD_REQS) compiler_flags $(RB_FS_SRC)
	$(CXX) $(CPP_FLAGS) $(RB_QUERY_REQS) -o $@ rb-validate.cpp $(RB_QUERY_SRC) $(RB_FS_SRC) $(KMC_OBJS) $(DEP_FLAGS)



all: $(BINARIES)

clean:
		rm -rf $(BINARIES) *.o *.dSYM



