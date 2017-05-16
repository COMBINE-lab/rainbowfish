
#include "rb-query.hpp"

template <class T1, class T2, class T3>
ColorDetector<T1, T2, T3>::ColorDetector(std::string dir, uint64_t colorCnt, bool isLblDynamicLength, uint64_t lblFixedLength) {
	A = T1(dir + "/lbl", false);
	isDynamicLblLength_ = isLblDynamicLength;
	lblFixedLength_ = lblFixedLength;
	if (isDynamicLblLength_) { 
		b = T2(dir + "/rnk", true); 
	}
	eqT = T3(dir + "/eqTable", false);
	colorCnt_ = colorCnt;
	prevEdge_=std::numeric_limits<uint64_t>::max();
	select_t = 0;
	getInt_t = 0;
	getColor_t = 0;
	prevColorVal_ = new uint64_t[ (colorCnt+63)/64 ];
}

template <class T1, class T2, class T3>
bool ColorDetector<T1, T2, T3>::prevContains_(unsigned int color) {
	//uint64_t st = getMilliCountt();
	bool res = eqT[colorCnt_*prevColor_+color];
	//uint64_t mask = 1;
	//bool res1 = prevColorVal_[color / 64] & (mask<<(color % 64));
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
		//we need the distance between two 1s, so I add 1 to the index which starts from 0
		return val + 1;
	}
}

template <class T1, class T2, class T3>
bool ColorDetector<T1, T2, T3>::contains(unsigned int color, uint64_t edge) {
	if (edge == prevEdge_) { return prevContains_(color); }
	uint64_t colorIdx = 0;
	if (isDynamicLblLength_) {
			uint64_t start = b.select(edge);
			//if (start > 83381100)
			//std::cerr<<edge<<" : "<<start<<"\n";				
			uint64_t next = b.getInt(start, 64);
			if (!(next & 0x0000000000000001)) {
				std::cerr << "lowest bit should always be set!\n";
			}
			//next = (next & 0x7FFFFFFFFFFFFFFF);
			next = (next & 0xFFFFFFFFFFFFFFFE);
			//std::cerr << "after\n";
			auto nzero = __builtin_ctzll(next);
			uint64_t end = start + nzero;
			//uint64_t end = b.select(edge+1);
			/*uint64_t end_check = b.select(edge+1);//todo: what if the edge is the last one?
			if (end != end_check) {
				std::cerr << "start = " << start << ", end with __builtin_clz = " << end << ", but with select = " << end_check << "\n";
			}*/

			/** Alternative bitscan forward **/
			/*
			//uint64_t end = b.select(edge+1);//todo: what if the edge is the last one?
			uint64_t next_word = b.getInt(start+1, 64);
			uint64_t len = bitscanforward(next_word); 
			*/
			/** Alternative bitscan forward **/

			colorIdx = A.getInt(start, end-start);
			colorIdx = colorIdx + ((1<<(end-start))-2);
	}
	else {
			colorIdx = A.getInt(edge*lblFixedLength_, lblFixedLength_);
	}
	//getInt_t += getMilliSpant(st);
	prevEdge_ = edge;
	prevColor_ = colorIdx;
	//st = getMilliCountt();
	bool res = eqT[colorCnt_*colorIdx+color];
	/*for (uint64_t c=0, cntr=0; c<colorCnt_; c+=64,cntr++) {
		prevColorVal_[cntr] = eqT.getInt(colorCnt_*colorIdx+c, std::min((int)(colorCnt_-c), 64));
	}*/
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

