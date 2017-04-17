
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
#include "pack-color.hpp"
#include <bitset>
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
	//TCLAP::UnlabeledValueArg<std::string> compress_arg("compress", "say compress if you want the data saved in a compressed format.", false, "", "compress", cmd);
    cmd.parse( argc, argv );
    params.input_filename  = input_filename_arg.getValue();
    params.num_colors  = atoi(num_colors_arg.getValue().c_str());
		
}

void deserialize_color_bv(std::ifstream &colorfile, color_bv &value)
{
    colorfile.read((char *)&value, sizeof(color_bv));
}

/*template <class T>
class ColorPacker {
	private:
			T eqTvec;
			T lblvec;
			T rnkvec;
	public:
			ColorPacker(std::string dir) :
				eqTvec(dir + "/eqTable"), lblvec(dir + "/lbl"), rnkvec(dir + "/rnk") {
				//TODO fillout info file
				// 	color cnt
				// 	kmer cnt
				// 	...
			}

			size_t insertColorLabel(T& bs, unsigned num, size_t pos) {
				// most significant bit of number goes down to the end of the bitset
				do {
					if (num&1) {
						bs.set(pos);
					}
					num>>=1;  
					pos++;
				} while(num);
				return pos;
			}

}*/
size_t insertColorLabel(boost::dynamic_bitset<>& bs, unsigned num, size_t pos) {
	// most significant bit of number goes down to the end of the bitset
    do {
    	if (num&1) {
	      bs.set(pos);
    	}
	num>>=1;  
    	pos++;
    } while(num);
    return pos;
}
size_t insertColorLabel(sdsl::bit_vector* bs, unsigned num, size_t pos) {
	// most significant bit of number goes down to the end of the bitset
    do {
    	if (num&1) {
	      (*bs)[pos] = 1;
    	}
	num>>=1;  
    	pos++;
    } while(num);
    return pos;
}


bool writeBitset(boost::dynamic_bitset<>& bs, std::ofstream& out) {
  int bitctr{7};
  unsigned char byte{0};
  size_t bs_size = bs.size();
  std::cout<<bs_size<<"\n";
  out.write(reinterpret_cast<char*>(&bs_size), sizeof(bs_size));
  for (size_t i = 0; i < bs.size(); ++i) {
    byte |= (bs[i] << bitctr);
    if (bitctr == 0) { 
       out.write(reinterpret_cast<char*>(&byte), sizeof(byte));
       byte = 0;
       bitctr = 8; 
    }
    bitctr--;
  }
  if (bitctr > 0) {
    out.write(reinterpret_cast<char*>(&byte), sizeof(byte));
  }
  return true;
}

int main(int argc, char * argv[])
{
	bool sort = true;
	bool compress = true;
	std::string res_dir = "bitvectors";
    //std::cerr << "pack-color compiled with supported colors=" << NUM_COLS << std::endl;
    std::cerr <<"Starting" << std::endl;
    parameters_t params;
    parse_arguments(argc, argv, params);

    const char * file_name = params.input_filename.c_str();
    std::cerr << "file name: " << file_name << std::endl;

	// Open File
    std::ifstream colorfile(file_name, std::ios::in|std::ios::binary);

    colorfile.seekg(0, colorfile.end);
    size_t end = colorfile.tellg();
    colorfile.seekg(0, colorfile.beg);
    std::cerr << "file size: " << end << std::endl;
    std::cerr << "sizeof(color_bv): " << sizeof(color_bv) << std::endl;
    size_t num_color = params.num_colors;
    size_t num_edges = end / sizeof(color_bv);

	//FIRST ROUND
	// Read file and fill out equivalence classes
    std::cerr << "edges: " << num_edges << " colors: " << num_color << " Total: " << num_edges * num_color << std::endl;
    sparse_hash_map<color_bv, int> eqCls;
    for (size_t i=0; i < num_edges; i++) {
        if (i % 100000000 == 0) {
            std::cerr << "deserializing edge " << i << std::endl;
        }
        color_bv value;
        deserialize_color_bv(colorfile, value);
		if (eqCls.find(value)==eqCls.end())
	    	eqCls[value] = 1;
		else
			eqCls[value] = eqCls[value]+1;
    }
    std::cerr << "Succinct builder object allocated." << std::endl;
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
    for (const auto& c : eqClsVec) {
	 	//std::cerr <<"eq cls "<<lbl++<< ": "<<c.first<<" : "<<c.second << "\n";
		totalBits += (lbl==0?c.second:ceil(log2(lbl+1))*c.second);
		//std::cerr << lbl << " : " << c.second << std::endl;
		eqCls[c.first] = lbl++;
	}
    size_t vecBits = totalBits;
    totalBits *= 2;
    totalBits += num_color * eqCls.size();
	std::cerr << "total bits: " << totalBits << " or " << totalBits/(8*pow(1024,2)) << " MB\n";

	if (compress) {
		// pack eqTable in bitvector
		sdsl::bit_vector* eqTable_nc = new sdsl::bit_vector(eqCls.size()*num_color);
		std::cerr << "eqTable allocated.\n";
		sdsl::bit_vector lbl_nc(vecBits);
		std::cerr << "lbl allocated.\n";
		sdsl::bit_vector rnk_nc(vecBits);
		std::cerr << "rnk allocated.\n";
		size_t i = 0;
		for (const auto& c : eqClsVec) {
			for (size_t j = 0; j < num_color; ++j) {
				(*eqTable_nc)[i++] = c.first[j];
			}
		}
		// compress eqTable
		sdsl::rrr_vector<63> eqTablevec(*eqTable_nc);
		std::cerr << "eqTable size :\n " 
				<< "		Not Compressed : " << sdsl::size_in_mega_bytes(*eqTable_nc) << " MB\n"
				<< "		Compressed     : " << sdsl::size_in_mega_bytes(eqTablevec) << " MB\n";
		delete eqTable_nc;
		std::cerr << " eqTable deleted.\n";
		// creating bitvectors A and b
		// A contains eq. class label for each k-mer in bits
		// b is set in the start position of each k-mer in A
		int sysTime = getMilliCount();
    	colorfile.seekg(0, colorfile.beg);
		size_t curPos = 0;
		for (size_t i=0; i < num_edges; i++) {
			//if (i % 1000000 == 0) std::cerr<<i<<std::endl;
			color_bv value;
			deserialize_color_bv(colorfile, value);
			rnk_nc[curPos] = 1;
			curPos = insertColorLabel(&lbl_nc, eqCls[value], curPos);
			// if we want to set the end, here we should say b.set(curPos-1);
			// TODO if for last edge we need an extra 1 we should have 1 extra bit at the end of the vector and set it here
		}
		// compress lbl & rnk
		// lbl
		std::cerr<<" Finished packing lbl and rnk.\n";
		sdsl::rrr_vector<63> lblvec(lbl_nc);
		std::cerr << "lbl size :\n " 
				<< "		Not Compressed : " << sdsl::size_in_mega_bytes(lbl_nc) << " MB\n"
				<< "		Compressed     : " << sdsl::size_in_mega_bytes(lblvec) << " MB\n";
		//delete lbl_nc;
		// rnk
		sdsl::rrr_vector<63> rnkvec(rnk_nc);
		std::cerr << "rnk size :\n " 
				<< "		Not Compressed : " << sdsl::size_in_mega_bytes(rnk_nc) << " MB\n"
				<< "		Compressed     : " << sdsl::size_in_mega_bytes(rnkvec) << " MB\n";
		//delete rnk_nc;
		
		// Serialize eqTable, lbl, and rnk bit_vectors
		sdsl::store_to_file(eqTablevec, res_dir + "/" + "eqTable.rrr");
		sdsl::store_to_file(lblvec, res_dir + "/" + "lbl.rrr");
		sdsl::store_to_file(rnkvec, res_dir + "/" + "rnk.rrr");
	}
	else {
		// put eqTable in bitvector
		boost::dynamic_bitset<> eqTable(eqCls.size()*num_color);
		size_t i = 0;
		for (const auto& c : eqClsVec) {
			for (size_t j = 0; j < num_color; ++j) {
				eqTable[i] = c.first[j];
				i++;
			}
		}

		// creating bitvectors A and b
		// A contains eq. class label for each k-mer in bits
		// b is set in the start position of each k-mer in A
		int sysTime = getMilliCount();
    	colorfile.seekg(0, colorfile.beg);
		boost::dynamic_bitset<> A(vecBits);
		boost::dynamic_bitset<> rnk(vecBits);
		size_t curPos = 0;
		for (size_t i=0; i < num_edges; i++) {
			color_bv value;
			deserialize_color_bv(colorfile, value);
			rnk.set(curPos);
			curPos = insertColorLabel(A, eqCls[value], curPos);
			// if we want to set the end, here we should say b.set(curPos-1);
			// TODO if for last edge we need an extra 1 we should have 1 extra bit at the end of the vector and set it here
		}
	
		std::cerr << "\nA, rnk & eq. cls BVs creation time : " << getMilliSpan(sysTime) << " ms\n";
    	std::cout<<"A.bitvec:";
	    std::ofstream Aout("A.bitvec", std::ios::out|std::ios::binary);
    	if (!writeBitset(A, Aout)) {
       		std::cerr << "Oh noes; couldn't write A!\n";
	    }
    	Aout.close();
	    std::cout<<"B.bitvec:";
    	std::ofstream Bout("B.bitvec", std::ios::out|std::ios::binary);
	    if (!writeBitset(rnk, Bout)) {
    	   std::cerr << "Oh noes; couldn't write B!\n";
    	}
	    Bout.close();
    	std::cout<<"eqTable.bitvec:";
		std::ofstream eqTableout("eqTable.bitvec", std::ios::out|std::ios::binary);
    	if (!writeBitset(eqTable, eqTableout)) {
	       std::cerr << "Oh noes; couldn't write eqTable!\n";
    	}
	    eqTableout.close();
	}	
	//std::cerr << "**********************    A    ***********************\n" << A;
	//std::cerr << "**********************   rank  ***********************\n" << rnk;
}
