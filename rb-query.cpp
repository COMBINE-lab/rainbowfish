
#include "rb-query.hpp"

template <class T>
ColorDetector<T>::ColorDetector(std::string Afile, std::string bfile, std::string eqfile, size_t colorCnt) :
		A(Afile, false), b(bfile, true), eqT(eqfile, false) {
	this->colorCnt = colorCnt;
}

template <class T>
bool ColorDetector<T>::contains(unsigned int color, uint64_t edge) {
	//edge++;
	uint64_t start = b.select(edge);
	uint64_t end = b.select(edge+1);//todo: what if the edge is the last one?
	//if (end - start > 2)
//		std::cout<<edge<<":s"<<start<<":e"<<end<<"\n";
	uint64_t colorIdx = A.getInt(start, end-start);
	//assumption: labels are put in A from least significant digit to most (e.g. label = 4 is put in A in the order 001)
	//for (uint64_t ctr = start; ctr < end; ctr++) {
		//colorIdx |= (A[ctr] << shifter++);
	//}
	//assumption: color i comes in position i in the original bv_color data structure
	//std::cout<<"B select:\n";
	//for (int i=0;i<10;i++) std::cout<<b.select(i);
	//std::cout<<"\n";
	//std::cout<<"Eq:\n";
	//for (int i =0;i<20;i++) std::cout<<eqT[i];
	//std::cout<<"\n";
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
	//uint64_t edge_cnt = 9256522;
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
