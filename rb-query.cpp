#include "rb-query.hpp"

ColorDetector::ColorDetector(std::string Afile, std::string bfile, std::string eqfile, size_t colorCnt) {
	A = readBitset(Afile);
	b = readRSBitset(bfile);	
	eqT = readBitset(eqfile);
	this->colorCnt = colorCnt;
}

bool ColorDetector::contains(unsigned int color, uint64_t edge) {
	edge++;
	uint64_t start = b->select(edge);
	uint64_t end = b->select(edge+1);//todo: what if the edge is the last one?
	uint64_t colorIdx = 0;
	size_t shifter = 0;
	//assumption: labels are put in A from least significant digit to most (e.g. label = 4 is put in A in the order 001)
	for (uint64_t ctr = start; ctr < end; ctr++) {
		colorIdx |= (A[ctr] << shifter++);
	}
	//assumption: color i comes in position i in the original bv_color data structure
	return eqT[colorCnt*colorIdx + color];
}

size_t ColorDetector::getColorCnt() {
	return colorCnt;
}

boost::dynamic_bitset<> ColorDetector::readBitset(std::string infileName){
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

BitmapPoppy* ColorDetector::readRSBitset(std::string infileName){
	std::ifstream infile(infileName, std::ios::in|std::ios::binary);
  	size_t bs_size;
	infile.read(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
	size_t byteCnt = ceil(bs_size/8.0);
    std::cout << infileName << " -- # of bytes: " << byteCnt <<std::endl;
	uint64_t* bitsInUint = new uint64_t[(int)ceil(bs_size/64.0)];
	//infile.read(reinterpret_cast<char*>(bitsInUint), byteCnt);
	//todo: resolve the case when byteCnt = 0 in a clean way
	// we have at least one byte in the file
	unsigned char* bytes;
	size_t byteCtr = 0;
	size_t uintCtr = 0;
	do {
		// either read 8 bytes (uint64_t) or if it is the end of the file, the remaining bytes.
		size_t readBytes = byteCnt-byteCtr > 8? 8:byteCnt-byteCtr; 		
		bytes = new unsigned char[readBytes];
		infile.read(reinterpret_cast<char*>(bytes), readBytes);
		unsigned char convertedBytes[8];
		// flip byte positions in the desired way!!!
		for (size_t i = 0; i < readBytes; i++) convertedBytes[7-i] = bytes[i];
		// create a uint64_t out of array of 8 bytes and put it in the result array
		bitsInUint[uintCtr++] = *((uint64_t*)convertedBytes);
 		// move the counter to the first unread byte in file
		byteCtr += readBytes;		
	} while(byteCtr < byteCnt);
	infile.close();
	BitmapPoppy* bitmap = new BitmapPoppy(bitsInUint, bs_size);
	return bitmap;
}

bool writeBitset(boost::dynamic_bitset<>& bs, std::string outfileName) {
  std::ofstream out(outfileName, std::ios::out|std::ios::binary);
  int bitctr{8};
  unsigned char byte{0};
  size_t bs_size = bs.size();
  out.write(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
  for (size_t i = 0; i < bs.size(); ++i) {
    bitctr--;
    byte |= (bs[i] << bitctr);
    if (bitctr == 0) {
       out.write(reinterpret_cast<char*>(&byte), sizeof(byte));
       byte = 0;
       bitctr = 8;
    }    
  }
  if (bitctr < 8) {
    out.write(reinterpret_cast<char*>(&byte), sizeof(byte));
  }
  return true;
}

/*int main(int, char*[]) {
	ColorDetector* cd = new ColorDetector("A.bitvec", "B.bitvec", "eqTable.bitvec", 6);
	std::ofstream f{"ours.res"};
	size_t num_colors = 6;
	for (int edge = 0; edge < 9254208; edge++) {
		short our_color_mask = 0;
		for (size_t c = 0; c < num_colors; c++) {
			our_color_mask |= (cd->contains(c, edge) << c);
		}
		f<<our_color_mask<<std::endl;
	}
	f.close();
	return EXIT_SUCCESS;
	
}*/

