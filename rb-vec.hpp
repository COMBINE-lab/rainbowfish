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

class RBVec {
	private:
			bool hasSelect;
			boost::dynamic_bitset<> bitvec;
			boost::dynamic_bitset<> readBitset(std::string infileName);
			//BitmapPoppy* bitselvec{nullptr};
			//BitmapPoppy* readRSBitset(std::string infileName);
			rank9sel* bitselvec{nullptr};
			rank9sel* readRSBitset(std::string infileName);

	public:
			RBVec(std::string fileName, bool hasSelect);
			bool operator[](uint64_t idx) const;
			uint64_t select(uint64_t rnk);
			uint64_t getInt(uint64_t offset, uint64_t bitLen);
			//uint64_t getNextOne(uint64_t offset);
};

class RBVecCompressed {
	private:
			bool hasSelect;
			sdsl::rrr_vector<63> bitvec_;
			sdsl::rrr_vector<63>::select_1_type selbitvec_;
	public:
			RBVecCompressed(std::string fileName, bool hasSelect);
			bool operator[](uint64_t idx) const;
			uint64_t select(uint64_t rnk);
			uint64_t getInt(uint64_t offset, uint64_t bitLen);
			//uint64_t getNextOne(uint64_t offset);

};

#endif
