#ifndef _RB_QUERY
#define _RB_QUERY

#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <math.h>

#include "bitmap.h"
#include "shared.h"


class ColorDetector {
	boost::dynamic_bitset<> A;
	BitmapPoppy* b;
	boost::dynamic_bitset<> eqT;
	size_t colorCnt;
  public:
	ColorDetector(std::string Afile, std::string bfile, std::string eqfile, size_t colorCnt); // todo: we should keep count of colors or k-mers in top of the file
	bool contains(unsigned int color, uint64_t edge);
	size_t getColorCnt();
  private:
	boost::dynamic_bitset<> readBitset(std::string filename);
	BitmapPoppy* readRSBitset(std::string filename);	
};

#endif
