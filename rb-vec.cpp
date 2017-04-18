#include "rb-vec.hpp"
#include "bit_array.h"

RBVec::RBVec(std::string fileName, bool hasSelect) {
		this->hasSelect = hasSelect;
		if (!hasSelect) bitvec_ = readBitset(fileName);
		else selbitvec_ = readRSBitset(fileName);
}

RBVec::RBVec(uint64_t bitSize) {
		boost::dynamic_bitset<> tmp_bitvec(bitSize);
		bitvec_ = tmp_bitvec;
}

//TODO implement it for rank9sel
bool RBVec::operator[](uint64_t idx) const {
	return bitvec_[idx];
}

//TODO will get so slow if I set if/else for bitvec vs bitselvec
void RBVec::set(uint64_t idx) {
	bitvec_.set();
}

// will just return a valid number if hasSelect is true
uint64_t RBVec::select(uint64_t rnk) {
	//no +1 for rank9sel but +1 for rankselect
	return selbitvec_->select(rnk);
}

//todo implement it for BitPoppy
uint64_t RBVec::getInt(uint64_t offset, uint64_t bitLen) {
	uint64_t shifter = 0;
	uint64_t inted = 0;
	//assumption: labels are put in A from least significant digit to most (e.g. label = 4 is put in A in the order 001)
	for (uint64_t ctr = 0; ctr < bitLen; ctr++) {
		inted |= (bitvec_[offset+ctr] << shifter++);
	}
	return inted;
}

bool RBVec::serialize(std::string fileName) {
	std::ofstream out(fileName+VEC_EXT, std::ios::out | std::ios::binary);
	int bitctr{7};
	unsigned char byte{0};
	size_t bs_size = bitvec_.size();
	out.write(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
	for (size_t i = 0; i < bs_size; ++i) {
		byte |= (bitvec_[i] << bitctr);
		if (bitctr == 0) {
			out.write(reinterpret_cast<char*>(&byte), sizeof(byte));
			byte = 0;
			bitctr = 8;
		}
		bitctr--;	
	}
	if (bitctr > 0) {
		out.write(reinterpret_cast<char*>(&byte), sizeof(byte));
	}
	return true;
}

boost::dynamic_bitset<> RBVec::readBitset(std::string infileName){
  std::ifstream infile(infileName+VEC_EXT, std::ios::binary);
  size_t bs_size;
  infile.read(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
  boost::dynamic_bitset<> bs(bs_size);
  size_t bytes_len = (size_t)ceil(bs_size/8.0);
  std::cout << infileName << " -- # of bytes: " << bytes_len <<std::endl;
  unsigned char bytes[bytes_len] = {0};
  infile.read(reinterpret_cast<char*>(bytes), bytes_len);
  size_t ctr = 0;
  for (size_t bytectr=0; bytectr < bytes_len; ++bytectr) {
    for (int i = 7; i >= 0; --i) {
      bs[ctr] = (bytes[bytectr] >> i) & 1;
      ctr++;
      if(ctr == bs_size) {bytectr = bytes_len;break;}
    }
  }
  infile.close(); 
  return bs;
}

rank9sel* RBVec::readRSBitset(std::string infileName) {
	std::ifstream infile(infileName+VEC_EXT, std::ios::in|std::ios::binary);
  	size_t bs_size;
	infile.read(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
	size_t byteCnt = ceil(bs_size/8.0);
    std::cout << infileName << " -- # of bytes: " << byteCnt <<std::endl;
    BIT_ARRAY* ba = bit_array_create(bs_size);
	for (size_t i = 0; i < byteCnt; ++i) {
		char b;
		infile.read(&b, sizeof(b));
		for (size_t j = 0; j < 8; ++j) {
			if ((b >> (7-j)) & 1) bit_array_set(ba, i * 8 + j); 
		}
	}
	infile.close();
	rank9sel* bitmap = new rank9sel(ba->words, ba->num_of_bits);
	return bitmap;
}



/* Compressed bitvector */

// loads a compressed rrr_vector from file
RBVecCompressed::RBVecCompressed(std::string fileName, bool hasSelect) {
	sdsl::load_from_file(rrr_bitvec_, fileName+COMPRESSEDVEC_EXT);
	if (hasSelect) {
		selbitvec_ = sdsl:: rrr_vector<63>::select_1_type(&rrr_bitvec_);
	}
}

// creates an uncompressed bitvector (to be converted to compressed version later)
RBVecCompressed::RBVecCompressed(uint64_t bitSize) {
	sdsl::bit_vector tmp_bitvec(bitSize);
	bitvec_ = tmp_bitvec;
}

// gets the rrr_vector value at position idx
bool RBVecCompressed::operator[](uint64_t idx) const {
	return (bool)rrr_bitvec_[idx];
}

//TODO ask Rob about the inconsistnecy of set and []
// sets bitvector value at position idx (cannot change rrr_vector values, so we have to assume, set is requested for uncompressed bitvector)
void RBVecCompressed::set(uint64_t idx) {
		bitvec_[idx] = 1;
}

//Assume that it is just called when hasSelect is set
uint64_t RBVecCompressed::select(uint64_t rnk) {
	//rank starts at 1 in sdsl::rrr_vector::select_1_type so we +1 the input
	rnk++;
	return selbitvec_(rnk);
}

// gets int value from rrr_vector
uint64_t RBVecCompressed::getInt(uint64_t offset, uint64_t bitLen) {
	uint64_t shifter = 0;
	uint64_t inted = 0;
	//assumption: labels are put in A from least significant digit to most (e.g. label = 4 is put in A in the order 001)
	for (uint64_t ctr = 0; ctr < bitLen; ctr++) {
		inted |= (rrr_bitvec_[offset+ctr] << shifter++);
	}
	return inted;
}

// compresses bitvector into rrr_vector and serializes it
bool RBVecCompressed::serialize(std::string fileName) {
	sdsl::rrr_vector<63> rrr_bitvec(bitvec_);
	return sdsl::store_to_file(rrr_bitvec, fileName+COMPRESSEDVEC_EXT);
}

sdsl::rrr_vector<63> RBVecCompressed::readRSbitset(std::string fileName) {
  std::ifstream infile(fileName+COMPRESSEDVEC_EXT, std::ios::binary);
  size_t bs_size;
  infile.read(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
  boost::dynamic_bitset<> bs(bs_size);
  size_t bytes_len = (size_t)ceil(bs_size/8.0);
  std::cout << fileName << " -- # of bytes: " << bytes_len <<std::endl;
  sdsl::bit_vector bitvec(bs_size);
  unsigned char bytes[bytes_len] = {0};
  infile.read(reinterpret_cast<char*>(bytes), bytes_len);
  size_t ctr = 0;
  for (size_t bytectr=0; bytectr < bytes_len; ++bytectr) {
    for (int i = 7; i >= 0; --i) {
      bitvec[ctr] = (bytes[bytectr] >> i) & 1;
      ctr++;
      if(ctr == bs_size) {bytectr = bytes_len;break;}
    }
  }
  infile.close(); 
  return bitvec;
}
