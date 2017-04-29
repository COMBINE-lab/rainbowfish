
#include "rb-query.hpp"

template <class T1, class T2, class T3>
ColorDetector<T1, T2, T3>::ColorDetector(std::string dir, uint64_t colorCnt) :
		A(dir + "/lbl", false), b(dir + "/rnk", true), eqT(dir + "/eqTable", false) {
	colorCnt_ = colorCnt;
	prevEdge_=std::numeric_limits<uint64_t>::max();
	select_t = 0;
	getInt_t = 0;
	getColor_t = 0;
}

template <class T1, class T2, class T3>
bool ColorDetector<T1, T2, T3>::prevContains_(unsigned int color) {
	//uint64_t st = getMilliCountt();
	bool res = eqT[colorCnt_*prevColor_+color];
	//getColor_t += getMilliSpant(st);
	return res;
}

inline int64_t bitscanforward(uint64_t val)
{
	if (val == 0) {
		return -1;
	} else {
		asm("bsf %[val], %[val]"
				: [val] "+r" (val)
				:
				: "cc");
		return val;
	}
}

template <class T1, class T2, class T3>
bool ColorDetector<T1, T2, T3>::contains(unsigned int color, uint64_t edge) {
	if (edge == prevEdge_) { return prevContains_(color); }
	//uint64_t st = getMilliCountt();
	uint64_t start = b.select(edge);

	uint64_t start_index = start / 64;
	uint64_t start_offset = start % 64;
	// get the machine word corresponding to the start index and mask it using
	// the start_offset
	uint64_t masked_word = (b.get_word(start) << (start_offset+1));
	// if the masked_word is 0 that means the next bit set is in the next word
	if (!masked_word) {
		start_index++;
		masked_word = (b.get_word(start+64));
	}
	uint64_t end = 0;
	int64_t end_offset = bitscanforward(masked_word);
	if (end_offset != -1)
		end = start_index * 64 + end_offset;
	else	// start is the last bit set in the bit vector.
		end = start_index * 64 + 63;
		
	//uint64_t end = b.select(edge+1);//todo: what if the edge is the last one?

	//select_t += getMilliSpant(st);
	//st = getMilliCountt();
	uint64_t colorIdx = A.getInt(start, end-start);
	//getInt_t += getMilliSpant(st);
	prevEdge_ = edge;
	prevColor_ = colorIdx;
	//st = getMilliCountt();
	bool res = eqT[colorCnt_*colorIdx+color];
	//getColor_t += getMilliSpant(st);
	return res;
}

template <class T1, class T2, class T3>
size_t ColorDetector<T1, T2, T3>::getColorCnt() {
	return colorCnt_;
}

template <class T1, class T2, class T3>
void ColorDetector<T1, T2, T3>::printStatistics() {
	std::cout<<"------------Some STATISTICS-------------\n";
	std::cout<<"		select timing: "<<select_t<<" ms";
	std::cout<<"		getColor timing: "<<getColor_t<<" ms";
	std::cout<<"		getInt timing: "<<getInt_t<<" ms"<<std::endl;
}

template class ColorDetector<RBVec, RBVec, RBVec>;
template class ColorDetector<RBVecCompressed, RBVecCompressed, RBVecCompressed>;
template class ColorDetector<RBVec, RBVecCompressed, RBVecCompressed>;
template class ColorDetector<RBVec, RBVec, RBVecCompressed>;

