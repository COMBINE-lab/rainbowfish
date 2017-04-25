
#include "rb-query.hpp"

template <class T1, class T2, class T3>
ColorDetector<T1, T2, T3>::ColorDetector(std::string dir, size_t colorCnt) :
		eqT(dir + "/eqTable", false), A(dir + "/lbl", false), b(dir + "/rnk", true) {
	this->colorCnt = colorCnt;
	prevEdge_=std::numeric_limits<uint64_t>::max();
}

template <class T1, class T2, class T3>
bool ColorDetector<T1, T2, T3>::prevContains_(unsigned int color) {
	return eqT[colorCnt*prevColor_+ color];
}

template <class T1, class T2, class T3>
bool ColorDetector<T1, T2, T3>::contains(unsigned int color, uint64_t edge) {
	if (edge == prevEdge_) { return prevContains_(color); }
	uint64_t start = b.select(edge);
	uint64_t end = b.select(edge+1);//todo: what if the edge is the last one?
	uint64_t colorIdx = A.getInt(start, end-start);
	prevEdge_ = edge;
	prevColor_ = colorIdx;
	return eqT[colorCnt*colorIdx + color];
}

template <class T1, class T2, class T3>
size_t ColorDetector<T1, T2, T3>::getColorCnt() {
	return colorCnt;
}

template class ColorDetector<RBVec, RBVec, RBVec>;
template class ColorDetector<RBVecCompressed, RBVecCompressed, RBVecCompressed>;
template class ColorDetector<RBVec, RBVecCompressed, RBVecCompressed>;
template class ColorDetector<RBVec, RBVec, RBVecCompressed>;

/*int main(int, char*[]) {
	ColorDetector<RBVecCompressed> cd("A.bitvec", "B.bitvec", "eqTable.bitvec", 6);
	std::ofstream f{"ours.res"};
	size_t num_colors = 6;
	uint64_t edge_cnt = 9256522;
	uint64_t edge_cnt = 12;
	for (uint64_t edge = 0; edge < edge_cnt; edge++) {
		short our_color_mask = 0;
		std::cout<<edge<<":";
		for (size_t c = 0; c < num_colors; c++) {
			std::cout<<cd.contains(c, edge);
			our_color_mask |= (cd.contains(c, edge) << c);
		}
		std::cout<< " --> "<<our_color_mask<<"\n";
		f<<our_color_mask<<std::endl;
	}
	f.close();
	return EXIT_SUCCESS;
}*/
