#include "rb-vec.hpp"
#include "bit_array.h"

RBVec::RBVec(std::string fileName, bool hasSelect) {
		this->hasSelect = hasSelect;
		if (!hasSelect) bitvec = readBitset(fileName);
		else bitselvec = readRSBitset(fileName);
}

//todo implement it for BitPoppy
bool RBVec::operator[](uint64_t idx) const {
	return bitvec[idx];
}

// will just return a valid number if hasSelect is true
uint64_t RBVec::select(uint64_t rnk) {
	//no +1 for rank9sel but +1 for rankselect
	return bitselvec->select(rnk);
}

//todo implement it for BitPoppy
uint64_t RBVec::getInt(uint64_t offset, uint64_t bitLen) {
	uint64_t shifter = 0;
	uint64_t inted = 0;
	//assumption: labels are put in A from least significant digit to most (e.g. label = 4 is put in A in the order 001)
	for (uint64_t ctr = 0; ctr < bitLen; ctr++) {
		inted |= (bitvec[offset+ctr] << shifter++);
	}
	return inted;
}

boost::dynamic_bitset<> RBVec::readBitset(std::string infileName){
  std::ifstream infile(infileName, std::ios::binary);
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

rank9sel* RBVec::readRSBitset(std::string infileName){
	std::ifstream infile(infileName, std::ios::in|std::ios::binary);
  	size_t bs_size;
	infile.read(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
	size_t byteCnt = ceil(bs_size/8.0);
    std::cout << infileName << " -- # of bytes: " << byteCnt <<std::endl;
    BIT_ARRAY* ba = bit_array_create(bs_size);
	// read byte by byte from file, go over bits of each byte and set bit_array if the bit is 1
	for (size_t i = 0; i < byteCnt; ++i) {
			char b;
			infile.read(&b, sizeof(b));
			for (size_t j = 0; j < 8; ++j) {
					if ((b >> (7-j)) & 1) bit_array_set(ba, i * 8 + j); 
			}
	}
	infile.close();
	//char s[65];
	//s[64] = '\0';
    //bit_array_to_substr(ba, 0, 64, s, '1', '0', 1);
	//std::cout << "[" << s << "]\n";
	rank9sel* bitmap = new rank9sel(ba->words, ba->num_of_bits);
	return bitmap;
}



/* Compressed bitvector */

/*RBVecCompressed::RBVecCompressed(std::string fileName, bool hasSelect) {
  std::ifstream infile(fileName, std::ios::binary);
  size_t bs_size;
  infile.read(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
  boost::dynamic_bitset<> bs(bs_size);
  size_t bytes_len = (size_t)ceil(bs_size/8.0);
  std::cout << fileName << " -- # of bytes: " << bytes_len <<std::endl;
  sdsl::bit_vector* bitvec = new sdsl::bit_vector(bs_size, 0);
  unsigned char bytes[bytes_len] = {0};
  infile.read(reinterpret_cast<char*>(bytes), bytes_len);
  size_t ctr = 0;
  for (size_t bytectr=0; bytectr < bytes_len; ++bytectr) {
    for (int i = 7; i >= 0; --i) {
      (*bitvec)[ctr] = (bytes[bytectr] >> i) & 1;
      ctr++;
      if(ctr == bs_size) {bytectr = bytes_len;break;}
    }
  }
  infile.close(); 
  bitvec_ = sdsl::rrr_vector<63>(*bitvec);
}
*/
