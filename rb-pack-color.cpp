
#include <iostream>
#include <fstream>
#include <utility>
#include <ctime>

// TCLAP
#include "tclap/CmdLine.h"

#include "cereal/archives/json.hpp"

#include "rb-filesystem.hpp"

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
	TCLAP::UnlabeledValueArg<std::string> pass_arg("pass", "1pass or 2pass", false, "", "pass", cmd);
    cmd.parse( argc, argv );
    params.input_filename  = input_filename_arg.getValue();
    params.num_colors  = atoi(num_colors_arg.getValue().c_str());
	params.res_dir = res_dir_arg.getValue();
	params.pass = pass_arg.getValue();
}

void deserialize_color_bv(std::ifstream &colorfile, color_bv &value)
{
    colorfile.read((char *)&value, sizeof(color_bv));
}

bool serialize_info(uint64_t num_colors, 
		    uint64_t num_edges, 
			uint64_t num_eqCls,
		    std::string label_type,
		    std::string select_type,
		    std::string eqtable_type,
		    std::string res_dir) {

	std::string jsonFileName = res_dir + "/info.json";
	std::ofstream jsonFile(jsonFileName);
	{	
	cereal::JSONOutputArchive archive(jsonFile);
	archive(cereal::make_nvp("label_type", label_type)); 
	archive(cereal::make_nvp("select_type", select_type));
	archive(cereal::make_nvp("eqtable_type", eqtable_type));
	archive(cereal::make_nvp("num_colors", num_colors));
	archive(cereal::make_nvp("num_edges", num_edges));
	archive(cereal::make_nvp("num_eqCls", num_eqCls));
	}
	jsonFile.close();
	return true;
}

template <class T1, class T2, class T3>
class ColorPacker {
	public:
			T1 lblvec;
			T2 rnkvec;
			T3 eqTvec;
	public:
			ColorPacker(uint64_t eqBitSize, uint64_t lblBitSize) :
				lblvec(lblBitSize), rnkvec(lblBitSize+1), eqTvec(eqBitSize) {
				//TODO fillout info file
				// 	color cnt
				// 	kmer cnt
				// 	...
			}

			size_t insertColorLabel(uint64_t num, uint64_t pos) {
				// most significant bit of number goes down to the end of the bitset
				uint8_t nbits = static_cast<uint8_t>(num==0?1:ceil(log2(num+1)));
				lblvec.setInt(pos, num, nbits);
				return pos + nbits;
			}
			
			bool storeAll(std::string dir, uint64_t bitvecSize, uint64_t eqClsSize) {
				return eqTvec.serialize(dir + "/eqTable", eqClsSize) && lblvec.serialize(dir + "/lbl", bitvecSize) && rnkvec.serialize(dir + "/rnk", bitvecSize+1);
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
	if (params.pass == "1pass") sort = false; //Anything else means apply sorting!! HeHe!!
    if (!rainbowfish::fs::FileExists(params.res_dir.c_str())) {
	rainbowfish::fs::MakeDir(params.res_dir.c_str());
    } 
    const char * file_name = params.input_filename.c_str();
    std::cerr << "file name: " << file_name << std::endl;

	const char * res_dir = params.res_dir.c_str();//"bitvectors";
	// Open File
    std::ifstream colorfile(file_name, std::ios::in|std::ios::binary);

    colorfile.seekg(0, colorfile.end);
    size_t end = colorfile.tellg();
    std::cerr << "file size: " << end << std::endl;
    std::cerr << "sizeof(color_bv): " << sizeof(color_bv) << std::endl;
    size_t num_color = params.num_colors;
    size_t num_edges = end / sizeof(color_bv);

	int startTime = getMilliCount();
	int checkPointTime = getMilliCount();
	int allocationTime = getMilliCount();
	//FIRST ROUND going over all edges
	// Read file and fill out equivalence classes
    std::cerr << "edges: " << num_edges << " colors: " << num_color << " Total: " << num_edges * num_color << std::endl;
    sparse_hash_map<color_bv, int> eqCls;
    std::vector<std::pair<color_bv, int>> eqClsVec;
    //bool first = true;
//	ColorPacker<RBVecCompressed, RBVecCompressed, RBVecCompressed> * cp;
	ColorPacker<RBVec, RBVec, RBVec> * cp;
	uint64_t curPos = 0;
	if (sort) {
    		colorfile.seekg(0, colorfile.beg);
			for (size_t i=0; i < num_edges; i++) {
				if (i % 100000000 == 0) {
					std::cerr<<getMilliSpan(checkPointTime) << " ms : "<<i<<std::endl;
					checkPointTime = getMilliCount();
				}
				color_bv value;
				deserialize_color_bv(colorfile, value);
				if (eqCls.find(value)==eqCls.end())
					eqCls[value] = 1;
				else
					eqCls[value] = eqCls[value]+1;
			}
			std::cerr << getMilliSpan(allocationTime) << " ms : Succinct builder object allocated" << std::endl;
			checkPointTime = getMilliCount();
			// Put data in hashmap to vector for further probable sorting!!
			eqClsVec.reserve(eqCls.size());
			for (const auto& c : eqCls) { eqClsVec.push_back(c); }	
			
			// sort the hashmap
			auto cmp = [](std::pair<color_bv, int> const & a, std::pair<color_bv, int> const & b) 
			{	 
				return a.second > b.second;
			};
			std::sort(eqClsVec.begin(), eqClsVec.end(), cmp);
			
			// replacing labels instead of k-mer counts as the hash map values
			int lbl = 0;
			size_t totalBits = 0;
			uint64_t total_edges = 0;
			for (const auto& c : eqClsVec) {
				std::cout <<"cls"<<lbl<< ", "<< c.second<<": ";
				total_edges += c.second;
				for (uint64_t k=0;k<num_color;k++) if (c.first[k] == true) std::cout<<k<<" ";
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
			
//			cp = new ColorPacker<RBVecCompressed, RBVecCompressed, RBVecCompressed>(eqCls.size()*num_color, vecBits);
			cp = new ColorPacker<RBVec, RBVec, RBVec>(eqCls.size()*num_color, vecBits);


			// SECOND ROUND going over all edges
			checkPointTime = getMilliCount();
			int packStartTime = getMilliCount();
			// create label & rank vectors
			colorfile.seekg(0, colorfile.beg);
			for (size_t i=0; i < num_edges; i++) {
				if (i % 100000000/*1000000*/ == 0) {
						std::cerr<<getMilliSpan(checkPointTime) << " ms : "<<i<<std::endl;
						checkPointTime = getMilliCount();
				}
				color_bv value;
				deserialize_color_bv(colorfile, value);
				(cp->rnkvec).set(curPos);
				curPos = cp->insertColorLabel(eqCls[value], curPos);
				// if we want to set the end, here we should say b.set(curPos-1);
				// TODO if for last edge we need an extra 1 we should have 1 extra bit at the end of the vector and set it here
			}
			(cp->rnkvec).set(curPos);
			std::cerr << "\n" << getMilliSpan(packStartTime) << " ms : Packing label & rank into bitvector." << std::endl;
	}
	else {
		checkPointTime = getMilliCount();
		int packStartTime = getMilliCount();
		// create label & rank vectors
		colorfile.seekg(0, colorfile.beg);
//		cp = new ColorPacker<RBVecCompressed, RBVecCompressed, RBVecCompressed>(num_color*num_color, 2*num_edges*num_color);
		cp = new ColorPacker<RBVec, RBVec, RBVec>(num_color*num_color, 2*num_edges*num_color);

		for (size_t i=0; i < num_edges; i++) {
				if (i % 100000000/*1000000*/ == 0) {
						std::cerr<<getMilliSpan(checkPointTime) << " ms : "<<i<<std::endl;
						checkPointTime = getMilliCount();
				}
				color_bv value;
				deserialize_color_bv(colorfile, value);
				if (eqCls.find(value)==eqCls.end()) {
					eqCls[value] = eqClsVec.size();
					eqClsVec.push_back(std::make_pair(value, eqCls[value]));
				}
				(cp->rnkvec).set(curPos);
				curPos = cp->insertColorLabel(eqCls[value], curPos);
				// if we want to set the end, here we should say b.set(curPos-1);
				// TODO if for last edge we need an extra 1 we should have 1 extra bit at the end of the vector and set it here
			}
			(cp->rnkvec).set(curPos);
			std::cerr << "\n" << getMilliSpan(packStartTime) << " ms : Packing label & rank into bitvector." << std::endl;	
	}
	// pack eqTable in bitvector
	checkPointTime = getMilliCount();
	uint64_t i = 0;
	for (const auto& c : eqClsVec) {
		for (size_t j = 0; j < num_color; ++j) {
			if (c.first[j]) (cp->eqTvec).set(i);
			i++;
		}
	}
	std::cerr << getMilliSpan(checkPointTime) << " ms : Packing eq. table into bitvector." << std::endl;


	checkPointTime = getMilliCount();
	cp->storeAll(res_dir, curPos, eqCls.size()*num_color);
	std::cerr << getMilliSpan(checkPointTime) << " ms : Storing all three bitvectors." << std::endl << std::endl;
	serialize_info(num_color, num_edges, eqCls.size(), "uncompressed", "compressed", "compressed", res_dir);
	std::cerr << getMilliSpan(startTime)/1000.0 << " s : Total Time." << std::endl;
}
