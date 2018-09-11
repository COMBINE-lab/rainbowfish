cmake_minimum_required(VERSION 3.0)

project(rainbowfish)

set(VSTRING 0.7.0)
set(BSTRING 'Copyright-people-2017')
set(NCOLORS 6000)

set(KMER_LEN 32)

message(STATUS "vstring = ${VSTRING}")
set(RBF_CPP_FLAGS "-pipe -m64 -std=c++14 -W -Wall -Wextra -Wpointer-arith -Wunused -Wwrite-strings -openmp")
set(RBF_CPP_FLAGS "${RBF_CPP_FLAGS} -DVERSION=\"${VSTRING}\" -DBANNER=\"${BSTRING}\" -DNUM_COLS=${NCOLORS} -DK_LEN=${KMER_LEN}")
set(OPT_FLAGS "-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing")
set(DEBUG_FLAGS "-pg -g -gstabs")

set(RBF_CPP_FLAGS "${RBF_CPP_FLAGS} ${OPT_FLAGS}")

set(CMAKE_CXX_FLAGS ${RBF_CPP_FLAGS})

#CLANG_WARNINGS=-Wbool-conversions -Wshift-overflow -Wliteral-conversion # CLANG ONLY
set(BOOST_PATH "${PROJECT_SOURCE_DIR}/../3rd_party_inst")
set(COSMO_PATH "${PROJECT_SOURCE_DIR}/../")
set(DEP_PATH "${COSMO_PATH}/3rd_party_inst")
set(KMC_PATH "${COSMO_PATH}/3rd_party_src/KMC")

include_directories(${KMC_PATH}
		    ${DEP_PATH}/include 
		    ${BOOST_PATH}/include 
		    ${COSMO_PATH} 
		    ${PROJECT_SOURCE_DIR}/include
		    )
	
link_directories(${DEP_PATH}/lib 
       		 ${PROJECT_SOURCE_DIR}
		 ${BOOST_PATH}/lib)

	 #file(GLOB SUX_SRC ${PROJECT_SOURCE_DIR}/Sux/*.cpp)
	 #message(STATUS ${SUX_SRC})
	 #add_library(sux STATIC ${SUX_SRC})
add_library(rainbowfish STATIC rb-vec.cpp rb-filesystem.cpp rb-query.cpp xxhash.c)

add_executable(rb-pack-color rb-pack-color.cpp)
target_link_libraries(rb-pack-color rainbowfish boost_system boost_filesystem sdsl)

add_executable(rb-find-bubble rb-find-bubble.cpp)
target_link_libraries(rb-find-bubble rainbowfish boost_system boost_filesystem sdsl)

add_executable(rb-validate rb-validate.cpp)
target_link_libraries(rb-validate rainbowfish boost_system boost_filesystem sdsl)

#add_executable(test-select test_select.cpp)
#target_link_libraries(test-select rainbowfish sux sdsl)


#RB_VEC_SRC=rb-vec.cpp rb-vec.hpp bit_array.c $(RANK9SEL_SRC)
#
#target_link_libraries(
#
#rb-pack-color: rb-pack-color.cpp $(RB_VEC_SRC) $(BUILD_REQS) compiler_flags $(RB_FS_SRC)
#	$(CXX) $(CPP_FLAGS) rb-pack-color.cpp $(RB_VEC_SRC) $(RB_FS_SRC) -o $@ $(KMC_OBJS) $(DEP_FLAGS)
#
#rb-find-bubble: rb-find-bubble.cpp $(RB_QUERY_SRC) $(BUILD_REQS) compiler_flags $(RB_FS_SRC)
#	$(CXX) $(CPP_FLAGS) $(RB_QUERY_REQS) -o $@ rb-find-bubble.cpp $(RB_QUERY_SRC) $(RB_FS_SRC) $(KMC_OBJS) $(DEP_FLAGS)
#
#rb-validate: rb-validate.cpp $(RB_QUERY_SRC) $(BUILD_REQS) compiler_flags $(RB_FS_SRC)
#	$(CXX) $(CPP_FLAGS) $(RB_QUERY_REQS) -o $@ rb-validate.cpp $(RB_QUERY_SRC) $(RB_FS_SRC) $(KMC_OBJS) $(DEP_FLAGS)
#
#
#
#
#BOOST_FLAGS= -lboost_system -lboost_filesystem
#DEP_FLAGS=$(INC) $(LIB) $(BOOST_FLAGS) -isystem $(KMC_PATH) -lsdsl -fopenmp #openmp is needed for logging
#DEBUG_FLAGS=-pg -g -gstabs
#NDEBUG_FLAGS=-DNDEBUG
#OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing
#NOPT_FLAGS=-O0
#
## Using Semantic Versioning: http://semver.org/
#VERSION=0.7.0
#BANNER='Copyright Alex Bowe (c) 2016'
#CPP_FLAGS+=-DVERSION=\"$(VERSION)\" -DBANNER=\"$(BANNER)\"
#
#k?=32
## TODO: quantize k
#CPP_FLAGS+=-DK_LEN=$(k)
#
#colors?=6000
#CPP_FLAGS+=-DNUM_COLS=$(colors)
#
#ifeq ($(asm),1)
#CPP_FLAGS+=-S -fverbose-asm
#endif
#
#ifeq ($(optimise),0)
#CPP_FLAGS+=$(NOPT_FLAGS)
#else
#CPP_FLAGS+=$(OPT_FLAGS)
#endif
#
#ifeq ($(debug),1)
#CPP_FLAGS+=$(DEBUG_FLAGS)
#else
#CPP_FLAGS+=$(NDEBUG_FLAGS)
#endif
#
#ifeq ($(verbose),1)
#CPP_FLAGS+=-DVERBOSE
#endif
#
#KMC_OBJS=$(KMC_PATH)/kmc_api/kmc_file.o $(KMC_PATH)/kmc_api/kmer_api.o $(KMC_PATH)/kmc_api/mmer.o
#RANKSELECT_OBJS=$(RANKSELECT_PATH)/bitmap.o
#BUILD_REQS=$(COSMO_PATH)/lut.hpp $(COSMO_PATH)/debug.hpp $(COSMO_PATH)/utility.hpp $(COSMO_PATH)/io.hpp $(COSMO_PATH)/sort.hpp $(COSMO_PATH)/kmer.hpp $(COSMO_PATH)/dummies.hpp $(COSMO_PATH)/debruijn_graph.hpp $(COSMO_PATH)/debruijn_graph_shifted.hpp $(COSMO_PATH)/pack-color.hpp $(COSMO_PATH)/cosmo-color-pd.hpp $(COSMO_PATH)
#COLOR_REQS=$(COSMO_PATH)/colored_debruijn_graph.hpp $(COSMO_PATH)/io.hpp $(COSMO_PATH)/debug.hpp
#RANK9SEL_SRC=./Sux/*.cpp
#RB_VEC_SRC=rb-vec.cpp rb-vec.hpp bit_array.c $(RANK9SEL_SRC)
#RB_QUERY_SRC=rb-query.cpp $(RB_VEC_SRC)
#RB_FS_SRC=rb-filesystem.cpp
#
#BINARIES=rb-pack-color rb-find-bubble rb-validate
#
#
