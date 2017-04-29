#ifndef _RBVECS_
#define _RBVECS_

#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <math.h>

//#include "bitmap.h"
//#include "shared.h"
#include "rank9sel.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/select_support.hpp>

#define VEC_EXT ".bitvec"
#define COMPRESSEDVEC_EXT ".rrr"

class RBVec {
	private:
			bool hasSelect;
			sdsl::bit_vector bitvec_;
			//sdsl::select_support_mcl<1,1>* selbitvec_{nullptr};
			rank9sel* selbitvec_{nullptr};
			//sdsl::select_support* selbitvec_{nullptr};
			//boost::dynamic_bitset<> bitvec_;
			boost::dynamic_bitset<> readBitset(std::string infileName);
			size_t bvsize_;
			//BitmapPoppy* bitselvec{nullptr};
			//BitmapPoppy* readRSBitset(std::string infileName);
			rank9sel* readRSBitset(std::string infileName);

	public:
			RBVec(std::string fileName, bool hasSelect);
			RBVec(uint64_t bitSize);
			bool operator[](uint64_t idx) const;
			void set(uint64_t idx);
			uint64_t select(uint64_t rnk);
			uint64_t getInt(uint64_t offset, uint64_t bitLen);
			bool setInt(uint64_t offset, uint64_t num, uint8_t bitlen);
			bool serialize(std::string fileName);
			//uint64_t getNextOne(uint64_t offset);
};

class RBVecCompressed {
	using compressed_vec = sdsl::rrr_vector<15, sdsl::int_vector<>, 8>;
	private:
			bool hasSelect;
			size_t bvsize_;
			sdsl::bit_vector bitvec_;
			compressed_vec rrr_bitvec_;
			decltype(rrr_bitvec_)::select_1_type selbitvec_;
			//sdsl::rrr_vector<63, sdsl::int_vector<>, 8> readRSbitset(std::string fileName);
	public:
			RBVecCompressed(std::string fileName, bool hasSelect);
			RBVecCompressed(uint64_t bitSize);
			bool operator[](uint64_t idx) const;
			void set(uint64_t idx);
			uint64_t select(uint64_t rnk);
			uint64_t getInt(uint64_t offset, uint64_t bitLen);
			bool setInt(uint64_t offset, uint64_t num, uint8_t bitlen);
			bool serialize(std::string fileName);
			//uint64_t getNextOne(uint64_t offset);

};

#endif
