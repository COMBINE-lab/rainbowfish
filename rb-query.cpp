#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <math.h>

bool writeBitset(boost::dynamic_bitset<>& bs, std::string outfileName) {
  std::ofstream out(outfileName, std::ios::out|std::ios::binary);
  int bitctr{7};
  unsigned char byte{0};
  size_t bs_size = bs.size();
  out.write(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
  for (size_t i = 0; i < bs.size(); ++i) {
    byte |= (bs[i] << bitctr);
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

boost::dynamic_bitset<> readBitset(std::string infileName){
  std::ifstream infile(infileName, std::ios::binary);
  size_t bs_size;
  infile.read(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
  boost::dynamic_bitset<> bs(bs_size);
  size_t bytes_len = (size_t)ceil(bs_size/8.0);
  unsigned char bytes[bytes_len] = {0};
  infile.read(reinterpret_cast<char*>(bytes), bs_size);
  size_t ctr = 0;
  for (size_t bytectr=0; bytectr < bytes_len; ++bytectr) { 
    for (size_t i = 7; i >= 0; --i) {
      bs[ctr] = (bytes[bytectr] >> i) & 1;
      ctr++;
      if(ctr == bs_size) {bytectr = bytes_len;break;}
    }
  }
  infile.close();
  return bs;
}
/*void readBitset(std::string infileName, uint64* bits){
	std::ifstream infile(infileName, std::ios::in|std::ios::binary);
	infile.seekg (0,infile.end);
	long size = infile.tellg();
	infile.seekg (0);
	bits = new unint64[size/64];
	infile.read(bits, sizeof(bits));
	infile.close();
}*/
int main(int, char*[]) {
  boost::dynamic_bitset<> x(10); // all 0's by default
  x[0] = 1;
  x[1] = 1;
  x[4] = 1;
  x[6] = 1;
  std::cout << "Original bitset: " << x << "\n";
  if (writeBitset(x, "test.bs")) std::cout<<"Successfully written.\n";
  boost::dynamic_bitset<> b = readBitset("test.bs");
  std::cout<<"After serializing & deserializing: "<<b<<"\n";
  return EXIT_SUCCESS;
}

