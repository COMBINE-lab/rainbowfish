#ifndef _RBVECS_
#define _RBVECS_

#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <math.h>

//#include "bitmap.h"
//#include "shared.h"
//#include "rank9sel.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/select_support.hpp>

#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

#define VEC_EXT ".bitvec"
#define COMPRESSEDVEC_EXT ".rrr"

class RBVec {
	private:
			bool hasSelect_;
			sdsl::bit_vector bitvec_;
      sdsl::bit_vector::select_1_type selbitvec_;
			size_t bvsize_;
	public:
			RBVec(){}
  //~RBVec();
			RBVec(std::string fileName, bool hasSelect);
			RBVec(uint64_t bitSize);
			RBVec& operator=(const RBVec& other);
			bool operator[](uint64_t idx) const;
			void set(uint64_t idx);
			uint64_t select(uint64_t rnk);
			uint64_t getInt(uint64_t offset, uint64_t bitLen);
			bool setInt(uint64_t offset, uint64_t num, uint8_t bitlen);
			bool serialize(std::string fileName, uint64_t finalSize);
};

class RBVecCompressed {
	using compressed_vec = sdsl::rrr_vector<63>;
	private:
			bool hasSelect_;
			size_t bvsize_;
			sdsl::bit_vector bitvec_;
			compressed_vec rrr_bitvec_;
			decltype(rrr_bitvec_)::select_1_type selbitvec_;
	public:
			RBVecCompressed(){}
			RBVecCompressed(std::string fileName, bool hasSelect);
			RBVecCompressed(uint64_t bitSize);
			RBVecCompressed& operator=(const RBVecCompressed& other);
			bool operator[](uint64_t idx) const;
			void set(uint64_t idx);
			uint64_t select(uint64_t rnk);
			uint64_t getInt(uint64_t offset, uint64_t bitLen);
			bool setInt(uint64_t offset, uint64_t num, uint8_t bitlen);
			bool serialize(std::string fileName, uint64_t finalSize);
};

#endif
