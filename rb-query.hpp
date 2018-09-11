#ifndef _RB_QUERY
#define _RB_QUERY

#include "rb-vec.hpp"

#include <sys/timeb.h>



template <class T1, class T2, class T3>
class ColorDetector {
  private:
	T1 A;
	T2 b;
	T3 eqT;
	uint64_t select_t;
	uint64_t getInt_t;
	uint64_t getColor_t;
	
	uint64_t colorCnt_;
	uint64_t prevEdge_;
	uint64_t prevColor_;
	uint64_t* prevColorVal_;
	bool isDynamicLblLength_;
	uint64_t lblFixedLength_;
	bool prevContains_(unsigned int color);
  public:
	ColorDetector(std::string dir, uint64_t colorCnt, bool isLblDynamicLength, uint64_t lblFixedLength); 
	bool contains(unsigned int color, uint64_t edge);
	uint64_t getColorCnt();
	void printStatistics();
int getMilliCountt(){
  timeb tb;
  ftime(&tb);
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}


int getMilliSpant(int nTimeStart){
  int nSpan = getMilliCountt() - nTimeStart;
  if(nSpan < 0)
    nSpan += 0x100000 * 1000;
  return nSpan;
}


  //private:
	//boost::dynamic_bitset<> readBitset(std::string filename);
	//BitmapPoppy* readRSBitset(std::string filename);	
};

#endif
