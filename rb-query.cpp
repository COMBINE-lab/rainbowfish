
#include "rb-query.hpp"

template <class T>
ColorDetector<T>::ColorDetector(std::string dir, size_t colorCnt) :
		A(dir + "/lbl", false), b(dir + "/rnk", true), eqT(dir + "/eqTable", false) {
	this->colorCnt = colorCnt;
}

template <class T>
bool ColorDetector<T>::prevContains_(unsigned int color) {
	return eqT[colorCnt*prevColor_+ color];
}

template <class T>
bool ColorDetector<T>::contains(unsigned int color, uint64_t edge) {
	if (edge == prevEdge_) { return prevContains_(color); }
	uint64_t start = b.select(edge);
	uint64_t end = b.select(edge+1);//todo: what if the edge is the last one?
	uint64_t colorIdx = A.getInt(start, end-start);
	prevEdge_ = edge;
	prevColor_ = colorIdx;	
	return eqT[colorCnt*colorIdx + color];
}

template <class T>
size_t ColorDetector<T>::getColorCnt() {
	return colorCnt;
}

template class ColorDetector<RBVec>;
template class ColorDetector<RBVecCompressed>;

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
