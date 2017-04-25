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

struct parameters_t {
  std::string input_filename = "";
  std::string color_filename = "";
  std::string res_dir = "";
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
  TCLAP::CmdLine cmd(BANNER, ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", ".dbg file.", true, "", "graph_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> color_filename_arg("color", ".rrr file.", true, "", "color_file", cmd);
    TCLAP::UnlabeledValueArg<std::string> res_dir_arg("dir",
					                                                             "Result directory. Should have created the directory first.", true, "", "res_dir", cmd);
	
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.color_filename  = color_filename_arg.getValue();
  params.res_dir		 = res_dir_arg.getValue();
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);
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
  uint64_t num_colors = 10;//00;

  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  cerr << "colors        : " << colors.size() / dbg.size() << endl; 
  cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;

  uint64_t startTime = getMilliCount();
  uint64_t checkPointTime = getMilliCount();
  uint64_t num_edges = dbg.num_edges();
  std::string res_dir = p.res_dir;
  ColorDetector<RBVecCompressed, RBVecCompressed, RBVecCompressed> cd(res_dir, num_colors);
  //ColorDetector<RBVec, RBVec, RBVec> cd(res_dir, num_colors);
  //ColorDetector<RBVec, RBVecCompressed, RBVecCompressed> cd(res_dir, num_colors);

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
			 //std::cout << cd.contains(c, edge) << "," << colors[edge*num_colors+c] << " ";
//			 if (cd.contains(c, edge) != colors[edge * num_colors + c]) {
			 if (rb != cosmo) {
					 allTheSame = false;
					 if (first) {
							 std::cout << "rbVSvari e" <<edge << "--> ";
							 first = false;					
					 }
					 std::cout << "c"<<c << ":" << cd.contains(c, edge) << "," << colors[edge*num_colors+c] << " ";
			 }
			 //if (edge < 100) std::cerr<<cd.contains(c, edge);
       }
	   //if (edge < 100) std::cerr<<" ";
	   if (!first) std::cout << "\n";
	   if (edge % 1000000000/*1000000*/ == 0) {
			   std::cerr << getMilliSpan(checkPointTime) << " ms : " << edge << " out of " << num_edges << " edges were compared.\n";
			   checkPointTime = getMilliCount();
	   }
  }
  std::cerr << "rbsum: " << rbsum << " cosmosum: "<<cosmosum;
  std::cerr << "\n\n" << getMilliSpan(startTime) << " ms : Time for total of " << num_edges * num_colors << " comparisons.\n";  
  if (allTheSame) std::cerr<<" HURRAAAAY! Validation Test Passed.\n";
  else std::cerr<<"NOT GOOD! Validation Test Failed.\n";

  /*
  bool temp;
  uint64_t rainbow_st = getMilliCount();
  checkPointTime = getMilliCount();
  for (uint64_t edge = 0; edge < num_edges; edge++) {
       for (size_t c = 0; c < num_colors; c++) {
			   temp = cd.contains(c, edge);
	   }
//	   if (edge % 100000 == 0) {
//			   std::cerr << "rb " << getMilliSpan(checkPointTime)/1000 << " s : " << edge << " of " << num_edges << "\n";
//			   checkPointTime = getMilliCount();
//	   }
  }
  std::cerr<<"\n\n\n";
  rainbow_st = getMilliSpan(rainbow_st);

  uint64_t vari_st = getMilliCount();
  checkPointTime = getMilliCount();
  for (uint64_t edge = 0; edge < num_edges; edge++) {
       for (size_t c = 0; c < num_colors; c++) {
			   temp = colors[edge*num_colors+c];			   
	   }
//	   if (edge % 100000 == 0) {
//			   std::cerr << "v " << getMilliSpan(checkPointTime)/1000 << " s : "<< edge << " of " << num_edges <<"\n";
//			   checkPointTime = getMilliCount();
//	   }
  }
  std::cerr<<"\n\n\n";
  vari_st = getMilliSpan(vari_st);
  std::cerr << "\n\n Total of " << num_edges * num_colors << " comparisons:\n";
  std::cerr << "		" << rainbow_st << " ms : RAINBOWFISH\n";
  std::cerr << "		" << vari_st << " ms : VARI\n";
*/ 
}
