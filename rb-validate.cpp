#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debruijn_graph_shifted.hpp"
#include "algorithm.hpp"
#include "rb-query.hpp"

#include <sys/timeb.h>
#include <random>

#include "cereal/archives/json.hpp"

struct parameters_t {
  std::string input_filename = "";
  std::string color_filename = "";
  std::string res_dir = "";
  std::string bvs_type = "";
  std::string validation_type = "";
};


int getMilliCount(){
  timeb tb;
  ftime(&tb);
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}


int getMilliSpan(int nTimeStart){
  int nSpan = getMilliCount() - nTimeStart;
  if(nSpan < 0)
    nSpan += 0x100000 * 1000;
  return nSpan;
}

void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Rainbowfish validate", ' ', "0.1.0");
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", ".dbg file.", true, "", "graph_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> color_filename_arg("color", ".rrr file.", true, "", "color_file", cmd);
    TCLAP::UnlabeledValueArg<std::string> res_dir_arg("dir", "Result directory. Should have created the directory first.", true, "", "res_dir", cmd);
    TCLAP::UnlabeledValueArg<std::string> bitvectors_type_arg("bv_type","format is like ccc, uuu, ucc, .. c=compressed, u=uncompressed. order = label, rank, eqTable", true, "", "bv_type", cmd);	
    TCLAP::UnlabeledValueArg<std::string> validation_type_arg("validation_type","Validation Type: Accepted values:compare, query, random-query, cosmo-query", true, "", "validation_type", cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.color_filename  = color_filename_arg.getValue();
  params.res_dir		 = res_dir_arg.getValue();
  params.bvs_type = bitvectors_type_arg.getValue();
  params.validation_type = validation_type_arg.getValue();
}

void deserialize_info(uint64_t& num_colors, uint64_t& num_edges, std::string res_dir) {

	std::string jsonFileName = res_dir + "/info.json";
	std::ifstream jsonFile(jsonFileName);
	{	
		cereal::JSONInputArchive archive(jsonFile);
		archive(cereal::make_nvp("num_colors", num_colors));
		archive(cereal::make_nvp("num_edges", num_edges));
	}
	jsonFile.close();
}

class MainBase {
	public:
			MainBase(){}
			virtual void run(parameters_t& p){ std::cout<<p.bvs_type<<"\n";}
};

template <class T1, class T2, class T3>
class MainTemplatized : public MainBase {
	public:
			MainTemplatized(){}
			void run(parameters_t& p) {
  					cerr << typeid(T1).name() << " " << typeid(T2).name() << " " << typeid(T3).name() << endl;
					cerr << "pack-color compiled with supported colors=" << NUM_COLS << std::endl;
					//ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
					//Can add this to save a couple seconds off traversal - not really worth it.
					cerr << "loading dbg" << std::endl;
					debruijn_graph_shifted<> dbg;
					load_from_file(dbg, p.input_filename);
					//input.close();
					cerr << "loading colors" << std::endl;
					sd_vector<> colors;
					load_from_file(colors, p.color_filename);
					uint64_t num_colors = 0;
					uint64_t num_edges = 0;
					deserialize_info(num_colors, num_edges, p.res_dir);
					cerr << "k             : " << dbg.k << endl;
					cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
					cerr << "num_edges()   : " << dbg.num_edges() << " or " << num_edges <<  endl;
					cerr << "colors        : " << colors.size() / dbg.size() << " or " << num_colors << endl; 
					cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
					cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;

				  std::string res_dir = p.res_dir;
				  ColorDetector<T1, T2, T3> cd(res_dir, num_colors);
				  uint64_t checkPointTime = getMilliCount();
				
				  if (p.validation_type == "compare") {
				  	uint64_t startTime = getMilliCount();
				  	uint64_t rbsum = 0;
				  	uint64_t cosmosum = 0;
				  	bool allTheSame = true;
					for (uint64_t edge = 0; edge < num_edges; edge++) {
					   bool first = true;
					   for (size_t c = 0; c < num_colors; c++) {			 			   

							 short rb = cd.contains(c, edge);
							 short cosmo = colors[edge*num_colors+c];
							 rbsum += rb;
							 cosmosum += cosmo;
							 if (rb != cosmo) {
									 allTheSame = false;
									 if (first) {
											 std::cout << "rbVSvari e" <<edge << "--> ";
											 first = false;					
									 }
									 std::cout << "c"<<c << ":" << cd.contains(c, edge) << "," << colors[edge*num_colors+c] << " ";
							 }
					   }
					   if (!first) std::cout << "\n";
					   if (edge % 10000000 == 0) {
							   std::cerr << getMilliSpan(checkPointTime) << " ms : " << edge << " out of " << num_edges << " edges were compared.\n";
							   checkPointTime = getMilliCount();
					   }
				  }
				  std::cerr << "rbsum: " << rbsum << " cosmosum: "<<cosmosum;
				  std::cerr << "\n\n" << getMilliSpan(startTime) << " ms : Time for total of " << num_edges * num_colors << " comparisons.\n";  
				  if (allTheSame) std::cerr<<" HURRAAAAY! Validation Test Passed.\n";
				  else std::cerr<<"NOT GOOD! Validation Test Failed.\n";
				  }		
				  if (p.validation_type == "query" || p.validation_type == "random-query") {
						  bool temp;
						  uint64_t totalsetbit = 0;
						  uint64_t rainbow_st = getMilliCount();
						  checkPointTime = getMilliCount();
						  if (p.validation_type == "random-query") {
								std::random_device rd;
								std::mt19937_64 gen(rd());
								std::uniform_int_distribution<unsigned long long> dis;
								uint64_t *edgeIdx = new uint64_t[num_edges];
								uint64_t *colorIdx = new uint64_t[num_colors];
						  		for (uint64_t edge = 0; edge < num_edges; edge++) {
									edgeIdx[edge] = dis(gen)%num_edges;
							   }
									for (size_t c = 0; c < num_colors; c++) {
											colorIdx[c] = rand() % num_colors;
									}


							   for (uint64_t edge = 0; edge < num_edges; edge++) {
							   for (size_t c = 0; c < num_colors; c++) {
									   temp = cd.contains(colorIdx[c], edgeIdx[edge]);
									   totalsetbit += temp;
							   }
							   
							   if (edge % 10000000 == 0) {
									 std::cerr << "rb " << getMilliSpan(checkPointTime)/1000 << " s : " << edge << " of " << num_edges << "\n";
									 checkPointTime = getMilliCount();
							   }
									 //  std::cerr << "rb " << getMilliSpan(checkPointTime)/1000 << " s : " <<  num_edges << " for color " <<c< "\n";
									  // checkPointTime = getMilliCount();
							  // }
						  	  }
						  }
						  for (uint64_t edge = 0; edge < num_edges; edge++) {
							   for (size_t c = 0; c < num_colors; c++) {
									   temp = cd.contains(c, edge);
									   totalsetbit += temp;
							   }
							   if (edge % 10000000 == 0) {
									   std::cerr << "rb " << getMilliSpan(checkPointTime)/1000 << " s : " << edge << " of " << num_edges << "\n";
									   //cd.printStatistics();
									   checkPointTime = getMilliCount();
							   }
						  }
						  std::cerr<<"\n\n\n";
						  rainbow_st = getMilliSpan(rainbow_st);
						  std::cerr << "\n\n Total of " << num_edges * num_colors << " comparisons with "<<totalsetbit<<" set bits:\n";
						  std::cerr << "		" << rainbow_st << " ms : RAINBOWFISH\n";
						  cd.printStatistics();
				 }
		if (p.validation_type == "cosmo-query") {			
		  uint64_t vari_st = getMilliCount();
				  uint64_t totalsetbit = 0;
				  bool temp;
		  checkPointTime = getMilliCount();
		  for (uint64_t edge = 0; edge < dbg.num_edges(); edge++) {
			   for (size_t c = 0; c < 1000; c++) {
					   temp = colors[edge*num_colors+c];			   
							   totalsetbit += temp;
			   }
		//	   if (edge % 100000 == 0) {
		//			   std::cerr << "v " << getMilliSpan(checkPointTime)/1000 << " s : "<< edge << " of " << num_edges <<"\n";
		//			   checkPointTime = getMilliCount();
		//	   }
		  }
		  std::cerr<<"\n\n\n";
		  vari_st = getMilliSpan(vari_st);
		  std::cerr << "\n\n Total of " << num_edges * num_colors << " comparisons with "<<totalsetbit<<" set bits:\n";
		  //std::cerr << "		" << rainbow_st << " ms : RAINBOWFISH\n";
		  std::cerr << "		" << vari_st << " ms : VARI\n";
		}
	}			
};

template class MainTemplatized<RBVec, RBVec, RBVec>;
template class MainTemplatized<RBVecCompressed, RBVecCompressed, RBVecCompressed>;
template class MainTemplatized<RBVec, RBVecCompressed, RBVecCompressed>;
template class MainTemplatized<RBVec, RBVec, RBVecCompressed>;


int main(int argc, char* argv[]) {
	parameters_t p;
	parse_arguments(argc, argv, p);
	MainBase* m{nullptr};		
	if (p.bvs_type == "ccc")
				m = new MainTemplatized<RBVecCompressed, RBVecCompressed, RBVecCompressed>();
	else if (p.bvs_type == "uuu")
				m = new  MainTemplatized<RBVec, RBVec, RBVec>();
	else if (p.bvs_type == "ucc")
				m = new  MainTemplatized<RBVec, RBVecCompressed, RBVecCompressed>();
	else if (p.bvs_type == "uuc")
				m = new  MainTemplatized<RBVec, RBVec, RBVecCompressed>();
	if (m)
		m->run(p);
	else std::cout<<"Initialization failed.\n";

}
