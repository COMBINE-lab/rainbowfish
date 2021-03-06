#include "rb-vec.hpp"

RBVec::RBVec(std::string fileName, bool hasSelect) {
		hasSelect_ = hasSelect;
		sdsl::load_from_file(bitvec_, fileName+VEC_EXT);
		if (hasSelect) {
			//selbitvec_ = new rank9sel(bitvec_.data(), bitvec_.size());
			selbitvec_ = sdsl::bit_vector::select_1_type(&bitvec_);
    }
		bvsize_ = bitvec_.size();
	  	std::cerr << fileName+VEC_EXT << " -- size " << bitvec_.size() <<std::endl;
}

/*
RBVec::~RBVec() {
		if (selbitvec_) { delete selbitvec_; }
}
*/

RBVec::RBVec(uint64_t bitSize) {
		hasSelect_ = false;
		bitvec_ = sdsl::bit_vector(bitSize);
		std::cerr << " -- size " << bitvec_.size() << std::endl;
}

RBVec& RBVec::operator=(const RBVec& other) {
		hasSelect_ = other.hasSelect_;
		bitvec_ = other.bitvec_;
		if (hasSelect_) { 
      selbitvec_ = sdsl::bit_vector::select_1_type(&bitvec_);
		}
		bvsize_ = other.bvsize_;		
		std::cerr << "Ran RBVec operator=\n";
		return *this;
}

//TODO implement it for rank9sel
bool RBVec::operator[](uint64_t idx) const {
	return bitvec_[idx];
}

//TODO will get so slow if I set if/else for bitvec vs bitselvec
void RBVec::set(uint64_t idx) {
	if (bitvec_.size() <= idx) bitvec_.resize((uint64_t)(bitvec_.size()*1.5));
	bitvec_[idx] = 1;
}

// will just return a valid number if hasSelect is true
uint64_t RBVec::select(uint64_t rnk) {
	//no +1 for rank9sel but +1 for rankselect
	//return (rnk < 28273951) ? selbitvec_->select(rnk+1) : bitvec_.size() - 1;
	//return selbitvec_->select(rnk); 
  return selbitvec_(rnk+1);
}

//todo implement it for BitPoppy
uint64_t RBVec::getInt(uint64_t offset, uint64_t bitLen) {
	return bitvec_.get_int(offset, static_cast<uint8_t>(bvsize_-offset>bitLen?bitLen:bvsize_-offset));
}

bool RBVec::setInt(uint64_t offset, uint64_t num, uint8_t bitlen) {
	if (bitvec_.size() <= offset+bitlen-1) bitvec_.resize((uint64_t)(bitvec_.size()*1.5));
	bitvec_.set_int(offset, num, bitlen);
	return true;
}

bool RBVec::serialize(std::string fileName, uint64_t finalSize) {
	fileName = fileName + VEC_EXT;
	if (finalSize != bitvec_.size()) bitvec_.resize(finalSize);
	std::cerr << fileName << " -- size " << bitvec_.size() <<std::endl;
	return sdsl::store_to_file(bitvec_, fileName);
}

/* Compressed bitvector */

// loads a compressed rrr_vector from file
RBVecCompressed::RBVecCompressed(std::string fileName, bool hasSelect) {
	hasSelect_ = hasSelect;
	fileName = fileName+COMPRESSEDVEC_EXT;
	sdsl::load_from_file(rrr_bitvec_, fileName);
	if (hasSelect_) {
	  selbitvec_ = decltype(rrr_bitvec_)::select_1_type(&rrr_bitvec_);
	}
	bvsize_ = rrr_bitvec_.size();
  	std::cerr << fileName << " -- size " << rrr_bitvec_.size() <<std::endl;
}

// creates an uncompressed bitvector (to be converted to compressed version later)
RBVecCompressed::RBVecCompressed(uint64_t bitSize) {
	hasSelect_ = false;
	sdsl::bit_vector tmp_bitvec(bitSize);
	bitvec_ = tmp_bitvec;
	std::cerr << " -- size " << bitvec_.size() << std::endl;
}

RBVecCompressed& RBVecCompressed::operator=(const RBVecCompressed& other) {
		hasSelect_ = other.hasSelect_;
		bitvec_ = other.bitvec_;
		rrr_bitvec_ = other.rrr_bitvec_;
		if (hasSelect_) { 
				selbitvec_ = decltype(rrr_bitvec_)::select_1_type(&rrr_bitvec_);
		}
		bvsize_ = other.bvsize_;		
		std::cerr << "Ran RBVecCompressed operator=\n";
		return *this;
}

// gets the rrr_vector value at position idx
bool RBVecCompressed::operator[](uint64_t idx) const {
	return (bool)rrr_bitvec_[idx];
}

//TODO ask Rob about the inconsistnecy of set and []
// sets bitvector value at position idx (cannot change rrr_vector values, so we have to assume, set is requested for uncompressed bitvector)
void RBVecCompressed::set(uint64_t idx) {
	if (bitvec_.size() <= idx) bitvec_.resize((uint64_t)(bitvec_.size()*1.5));
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
	return rrr_bitvec_.get_int(offset, bvsize_-offset>bitLen?bitLen:bvsize_-offset);
}

bool RBVecCompressed::setInt(uint64_t offset, uint64_t num, uint8_t bitlen) {
	if (bitvec_.size() <= offset+bitlen-1) bitvec_.resize((uint64_t)(bitvec_.size()*1.5));
	bitvec_.set_int(offset, num, bitlen);
	return true;
}

// compresses bitvector into rrr_vector and serializes it
bool RBVecCompressed::serialize(std::string fileName, uint64_t finalSize) {
	fileName = fileName+COMPRESSEDVEC_EXT;
	if (finalSize != bitvec_.size()) bitvec_.resize(finalSize);
	compressed_vec rrr_bitvec(bitvec_);
  	std::cerr << fileName << " -- size " << rrr_bitvec.size() <<std::endl;
	return sdsl::store_to_file(rrr_bitvec, fileName);
}

/*
sdsl::rrr_vector<63, sdsl::int_vector<>, 8> RBVecCompressed::readRSbitset(std::string fileName) {
  std::ifstream infile(fileName+COMPRESSEDVEC_EXT, std::ios::binary);
  size_t bs_size;
  infile.read(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
  boost::dynamic_bitset<> bs(bs_size);
  size_t bytes_len = (size_t)ceil(bs_size/8.0);
  std::cerr << fileName << " -- size " << bs_size <<std::endl;
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
*/
