
#include <iostream>
#include <fstream>
#include <utility>
#include <ctime>

// TCLAP
#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <cstdio>

#include <cstdlib>

#include <libgen.h>
#include <sparsepp/spp.h>
#include <boost/dynamic_bitset.hpp>
using spp::sparse_hash_map;
// Custom Headers
//#include "uint128_t.hpp"
//#include "debug.h"
#include "kmer.hpp"

//using namespace std;
//using namespace sdsl;

#include <cstdlib>
#include <sys/timeb.h>
#include "rb-pack-color.hpp"
#include <bitset>
#include "rb-vec.hpp"
int getMilliCount()
{
    timeb tb;
    ftime(&tb);
    int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
    return nCount;
}


int getMilliSpan(int nTimeStart)
{
    int nSpan = getMilliCount() - nTimeStart;
    if(nSpan < 0)
        nSpan += 0x100000 * 1000;
    return nSpan;
}

std::string file_extension = ".<extension>";


void parse_arguments(int argc, char **argv, parameters_t & params)
{
    TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
    TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
                                                             "Input file. Currently only supports DSK's binary format (for k<=64).", true, "", "input_file", cmd);
    TCLAP::UnlabeledValueArg<std::string> num_colors_arg("num_colors",
                                                         "Number of colors", true, "", "num colors", cmd);
    TCLAP::UnlabeledValueArg<std::string> res_dir_arg("dir",
                                                             "Result directory. Should have created the directory first.", true, "", "res_dir", cmd);
	//TCLAP::UnlabeledValueArg<std::string> compress_arg("compress", "say compress if you want the data saved in a compressed format.", false, "", "compress", cmd);
    cmd.parse( argc, argv );
    params.input_filename  = input_filename_arg.getValue();
    params.num_colors  = atoi(num_colors_arg.getValue().c_str());
	params.res_dir = res_dir_arg.getValue();
}

void deserialize_color_bv(std::ifstream &colorfile, color_bv &value)
{
    colorfile.read((char *)&value, sizeof(color_bv));
}

bool serialize_info(uint64_t num_colors, uint64_t num_edges, std::string res_dir) {
	std::ofstream outfile(res_dir+"/info", std::ios::binary);
	outfile.write(reinterpret_cast<char*>(&num_colors), sizeof(num_colors));
	outfile.write(reinterpret_cast<char*>(&num_edges), sizeof(num_edges));
	outfile.close();
}

template <class T1, class T2, class T3>
class ColorPacker {
	public:
			T1 lblvec;
			T2 rnkvec;
			T3 eqTvec;
	public:
			ColorPacker(uint64_t eqBitSize, uint64_t lblBitSize) :
				eqTvec(eqBitSize), lblvec(lblBitSize), rnkvec(lblBitSize) {
				//TODO fillout info file
				// 	color cnt
				// 	kmer cnt
				// 	...
			}

			size_t insertColorLabel(unsigned num, uint64_t pos) {
				// most significant bit of number goes down to the end of the bitset
				do {
					if (num&1) {
						lblvec.set(pos);
					}
					num>>=1;  
					pos++;
				} while(num);
				return pos;
			}
			
			bool storeAll(std::string dir) {
				return eqTvec.serialize(dir + "/eqTable") && lblvec.serialize(dir + "/lbl") && rnkvec.serialize(dir + "/rnk");
			}
};

int main(int argc, char * argv[])
{
	bool sort = true;
	bool compress = true;
    std::cerr << "pack-color compiled with supported colors=" << NUM_COLS << std::endl;
    std::cerr <<"Starting" << std::endl;
    parameters_t params;
    parse_arguments(argc, argv, params);

    const char * file_name = params.input_filename.c_str();
    std::cerr << "file name: " << file_name << std::endl;

	const char * res_dir = params.res_dir.c_str();//"bitvectors";
	// Open File
    std::ifstream colorfile(file_name, std::ios::in|std::ios::binary);

    colorfile.seekg(0, colorfile.end);
    size_t end = colorfile.tellg();
    colorfile.seekg(0, colorfile.beg);
    std::cerr << "file size: " << end << std::endl;
    std::cerr << "sizeof(color_bv): " << sizeof(color_bv) << std::endl;
    size_t num_color = params.num_colors;
    size_t num_edges = end / sizeof(color_bv);

	int startTime = getMilliCount();
	int checkPointTime = getMilliCount();
	//FIRST ROUND going over all edges
	// Read file and fill out equivalence classes
    std::cerr << "edges: " << num_edges << " colors: " << num_color << " Total: " << num_edges * num_color << std::endl;
    sparse_hash_map<color_bv, int> eqCls;
    bool first = true;
	for (size_t i=0; i < num_edges; i++) {
        if (i % 100000000 == 0) {
            std::cerr << "deserializing edge " << i << std::endl;
        }
        color_bv value;
        deserialize_color_bv(colorfile, value);
		//std::cerr<<"value.size() : "<< value.size() << "\n";
		//std::cerr<<i << " : ";
		//for (int j = 0; j < value.size(); j++ ) 
		//		if (value[j] && first) {std::cerr<<i<<","<<j<<"\n";first = false;}
		//if (i == 100)
		//for (int j =0; j < value.size(); j++) std::cerr<<value[j];
		//std::cerr<<" ";
		//for (int j =5000; j < 6000; j++) std::cerr<<value[j];
		//std::cerr<<"\n";
		if (eqCls.find(value)==eqCls.end())
	    	eqCls[value] = 1;
		else
			eqCls[value] = eqCls[value]+1;
    }
    std::cerr << getMilliSpan(checkPointTime) << " ms : Succinct builder object allocated" << std::endl;
	checkPointTime = getMilliCount();
	// Put data in hashmap to vector for further probable sorting!!
    std::vector<std::pair<color_bv, int>> eqClsVec;
    eqClsVec.reserve(eqCls.size());
    for (const auto& c : eqCls) { eqClsVec.push_back(c); }	
	
	// sort the hashmap
	if (sort) {
		auto cmp = [](std::pair<color_bv, int> const & a, std::pair<color_bv, int> const & b) 
		{	 
    	 	return a.second > b.second;
		};
		std::sort(eqClsVec.begin(), eqClsVec.end(), cmp);
	}
	// replacing labels instead of k-mer counts as the hash map values
    int lbl = 0;
    size_t totalBits = 0;
	uint64_t total_edges = 0;
	std::cerr<<"Starting from hereeeee\n";
    for (const auto& c : eqClsVec) {
	 	std::cout <<"cls"<<lbl<< ", "<< c.second<<": ";
		total_edges += c.second;
		for (int k=0;k<num_color;k++) if (c.first[k] == true) std::cout<<k<<" ";
		std::cout<<"\n";
		totalBits += (lbl==0?c.second:ceil(log2(lbl+1))*c.second);
		eqCls[c.first] = lbl++;
	}
	std::cerr << getMilliSpan(checkPointTime) << " ms : (Sorting eq vector and ) assigning a label to each eq class." << std::endl;
	std::cerr << " Total edge vs total edge: "<<total_edges<<" vs "<< num_edges<<"\n";
    size_t vecBits = totalBits;
    totalBits *= 2;
    totalBits += num_color * eqCls.size();
	std::cerr << "total bits: " << totalBits << " or " << totalBits/(8*pow(1024,2)) << " MB\n";
	
	ColorPacker<RBVecCompressed, RBVecCompressed, RBVecCompressed> * cp = new ColorPacker<RBVecCompressed, RBVecCompressed, RBVecCompressed>(eqCls.size()*num_color, vecBits);
	//ColorPacker<RBVec, RBVec, RBVec> * cp = new ColorPacker<RBVec, RBVec, RBVec>(eqCls.size()*num_color, vecBits);
	//ColorPacker<RBVec, RBVecCompressed, RBVecCompressed> * cp = new ColorPacker<RBVec, RBVecCompressed, RBVecCompressed>(eqCls.size()*num_color, vecBits);

//	if (compress) cp = new ColorPacker<RBVecCompressed>(res_dir);
//	else cp = new ColorPacker<RBVec>(res_dir);

	checkPointTime = getMilliCount();
	// pack eqTable in bitvector
	uint64_t i = 0;
	for (const auto& c : eqClsVec) {
		for (size_t j = 0; j < num_color; ++j) {
			if (c.first[j]) (cp->eqTvec).set(i);
			i++;
		}
	}
	std::cerr << getMilliSpan(checkPointTime) << " ms : Packing eq. table into bitvector." << std::endl;

	// SECOND ROUND going over all edges
	checkPointTime = getMilliCount();
	int packStartTime = getMilliCount();
	//startTime = getMilliCount();
	// create label & rank vectors
   	colorfile.seekg(0, colorfile.beg);
	uint64_t curPos = 0;
	for (size_t i=0; i < num_edges; i++) {
		if (i % 100000000/*1000000*/ == 0) {
				std::cerr<<getMilliSpan(checkPointTime) << " ms : "<<i<<std::endl;
				checkPointTime = getMilliCount();
		}
		color_bv value;
		deserialize_color_bv(colorfile, value);
		if (i <= 2000)
			std::cout<<"e"<<i <<":" << eqCls[value] << " ";
	//	std::cout<<value<<"\n";
		//if (i < 100) {
			//std::cerr<<"pos"<<curPos<<", eq"<<eqCls[value]<<": ";
			//for (int j =0; j < num_color/*value.size()*/; j++) std::cerr<<value[j];
			//std::cerr<<" ";
		//}
		(cp->rnkvec).set(curPos);
		curPos = cp->insertColorLabel(eqCls[value], curPos);
		// if we want to set the end, here we should say b.set(curPos-1);
		// TODO if for last edge we need an extra 1 we should have 1 extra bit at the end of the vector and set it here
	}	
	std::cerr << "\n" << getMilliSpan(packStartTime) << " ms : Packing label & rank into bitvector." << std::endl;
	
	checkPointTime = getMilliCount();
	cp->storeAll(res_dir);
	std::cerr << getMilliSpan(checkPointTime) << " ms : Storing all three bitvectors." << std::endl << std::endl;
	serialize_info(num_color, num_edges, res_dir);
	std::cerr << getMilliSpan(startTime)/1000.0 << " s : Total Time." << std::endl;
}
