#ifndef _RB_QUERY
#define _RB_QUERY

#include "rb-vec.hpp"

template <class T1, class T2, class T3>
class ColorDetector {
  private:
	T1 A;
	T2 b;
	T3 eqT;
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
