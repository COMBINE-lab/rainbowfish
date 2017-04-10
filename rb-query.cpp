#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <math.h>
//#include <cstdint>

#include "bitmap.h"
#include "shared.h"

bool writeBitset(boost::dynamic_bitset<>& bs, std::string outfileName) {
  std::ofstream out(outfileName, std::ios::out|std::ios::binary);
  int bitctr{8};
  unsigned char byte{0};
  size_t bs_size = bs.size();
  out.write(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
  //std::cout<<"size:"<<bs_size<<"\n";
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

boost::dynamic_bitset<> readBitset(std::string infileName){
  std::ifstream infile(infileName, std::ios::binary);
  size_t bs_size;
  infile.read(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
  //std::cout<<"size:"<<bs_size<<"\n";
  boost::dynamic_bitset<> bs(bs_size);
  size_t bytes_len = (size_t)ceil(bs_size/8.0);
  std::cout << bytes_len <<std::endl;
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
std::pair<uint64*, uint64> readRSBitset(std::string infileName){
	std::ifstream infile(infileName, std::ios::in|std::ios::binary);
  	size_t bs_size;
	infile.read(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
	size_t byteCnt = ceil(bs_size/8.0);
	uint64* bitsInUint = new uint64[(int)ceil(bs_size/64.0)];
	//infile.read(reinterpret_cast<char*>(bitsInUint), byteCnt);
	//todo: resolve the case when byteCnt = 0 in a clean way
	// we have at least one byte in the file
	unsigned char* bytes;
	size_t byteCtr = 0;
	size_t uintCtr = 0;
	do {
		size_t readBytes = byteCnt-byteCtr > 8? 8:byteCnt-byteCtr; // either read 8 bytes (uint64) or if it is the end of the file, the remaining bytes.
		bytes = new unsigned char[readBytes];
		infile.read(reinterpret_cast<char*>(bytes), readBytes);
		unsigned char convertedBytes[8];
		for (size_t i = 0; i < readBytes; i++) convertedBytes[7-i] = bytes[i]; // flip byte positions in the desired way!!!
		bitsInUint[uintCtr++] = *((uint64*)convertedBytes); // create a uint64 out of array of 8 bytes and put it in the result array
		byteCtr += readBytes; // move the counter to the first unread byte in file		
	} while(byteCtr < byteCnt);
	infile.close();
	return std::make_pair(bitsInUint, *(reinterpret_cast<uint64*>(&bs_size)));
}
int main(int, char*[]) {
  boost::dynamic_bitset<> x(120); // all 0's by default
  x[8] = 1;
  x[15] = 1;
  x[119] = 1;
  x[115] = 1;
  std::cout << "Original bitset: " << x << "\n";
  if (writeBitset(x, "test.bs")) std::cout<<"Successfully written.\n";
  
  // test readBitset to boost::dynamic_bitset
  boost::dynamic_bitset<> b = readBitset("test.bs");
  std::cout<<"After serializing & deserializing: "<<b<<"\n";
  
  // test readBitset to uint64*
  std::pair<uint64*, uint64> bits = readRSBitset("test.bs");
  BitmapPoppy* bitmap = new BitmapPoppy(bits.first, bits.second);
  std::cout<<bitmap->select(1)<<"\n";
  std::cout<<bitmap->select(2)<<"\n";
  std::cout<<bitmap->select(3)<<"\n";
  std::cout<<bitmap->select(4)<<"\n";
  return EXIT_SUCCESS;
}

