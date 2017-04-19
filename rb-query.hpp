#ifndef _RB_QUERY
#define _RB_QUERY

#include "rb-vec.hpp"

template <class T>
class ColorDetector {
  private:
	T A;
	T b;
	T eqT;
	size_t colorCnt;
	uint64_t prevEdge_;
	uint64_t prevColor_;
	bool prevContains_(unsigned int color);
  public:
	ColorDetector(std::string dir, size_t colorCnt); // todo: we should keep count of colors or k-mers in top of the file
	bool contains(unsigned int color, uint64_t edge);
	size_t getColorCnt();
  //private:
	//boost::dynamic_bitset<> readBitset(std::string filename);
	//BitmapPoppy* readRSBitset(std::string filename);	
};

#endif
