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
#include "rb-find-bubble.hpp"

//using namespace std;
//using namespace sdsl;

#include <sys/timeb.h>

#include "cereal/archives/json.hpp"
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
  TCLAP::CmdLine cmd("Rainbowfish", ' ', "0.1.0");
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", ".dbg file.", true, "", "graph_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> color_filename_arg("color", ".rrr file.", true, "", "color_file", cmd);
  string color_mask1 = "color_mask1";
  TCLAP::ValueArg<std::string> color_mask1_arg("a", "color_mask1",
	    "Color mask 1, color1 [" + color_mask1 + "]", false, "", color_mask1, cmd);
  string color_mask2 = "color_mask2";
  TCLAP::ValueArg<std::string> color_mask2_arg("b", "color_mask2",
	    "Color mask 2, color2 [" + color_mask2 + "]", false, "", color_mask2, cmd);
    TCLAP::UnlabeledValueArg<std::string> res_dir_arg("dir", "Result directory. Should have created the directory first.", true, "", "res_dir", cmd);	
    TCLAP::UnlabeledValueArg<std::string> bitvectors_type_arg("bv_type","format is like ccc, uuu, ucc, .. c=compressed, u=uncompressed. order = label, rank, eqTable", true, "", "bv_type", cmd);

  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.color_filename  = color_filename_arg.getValue();
  params.color_mask1     = color_mask1_arg.getValue();
  params.color_mask2     = color_mask2_arg.getValue();
  params.res_dir		 = res_dir_arg.getValue();
  params.bvs_type         = bitvectors_type_arg.getValue();
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
static char base[] = {'?','A','C','G','T'};


void test_symmetry(debruijn_graph_shifted<> dbg) {
  for (unsigned long x = 0; x<dbg.sigma+1;x++) {
    ssize_t in = dbg.incoming(43, x);
    if (in == -1)
      continue;
    for (unsigned long y = 0; y<dbg.sigma+1;y++) {
      ssize_t out = dbg.outgoing(in, y);
      if (out == -1)
	continue;
      cout << "Incoming " << in <<  ":" << out <<"\n";
    }
  }
}



void dump_nodes(debruijn_graph_shifted<> dbg, uint64_t * colors) {
  for (size_t i = 0; i < dbg.num_nodes(); i++) {
    cout << i << ":" << dbg.node_label(i) << colors[dbg._node_to_edge(i)] << "\n";
  }
}


void dump_edges(debruijn_graph_shifted<> dbg, uint64_t * colors) {
  for (size_t i = 0; i < dbg.size(); i++) {
    cout << i << "e:" << dbg.edge_label(i) << colors[i] << "\n";
  }
}

const char *const starts[] = {"GCCATACTGCGTCATGTCGCCCTGACGCGC","GCAGGTTCGAATCCTGCACGACCCACCAAT","GCTTAACCTCACAACCCGAAGATGTTTCTT","AAAACCCGCCGAAGCGGGTTTTTACGTAAA","AATCCTGCACGACCCACCAGTTTTAACATC","AGAGTTCCCCGCGCCAGCGGGGATAAACCG","GAATACGTGCGCAACAACCGTCTTCCGGAG"};

void print_color(color_bv& color)
{
    std::string colstr = color.to_string();
    for (unsigned int first1 = 0; first1 < colstr.size() ; first1++) {
        if (colstr[first1] == '1') {
            std::string outstring = colstr.substr(first1, colstr.size());
            cout << outstring;
            return;
        }
    }
    cout << "0";

}

template <class T1, class T2, class T3>
void find_bubbles(const debruijn_graph_shifted<> &dbg, ColorDetector<T1, T2, T3>& colors, /*sd_vector<> colors2,*/ color_bv color_mask1, color_bv color_mask2)
{
	uint64_t qt = getMilliCount();
	uint64_t tqt = 0;
    int t = getMilliCount();
    uint64_t num_colors = colors.getColorCnt();

    sdsl::bit_vector visited = sdsl::bit_vector(dbg.num_nodes(), 0);
	std::cerr << "Starting to look for bubbles for "<<num_colors<< " colors.\n";
    std::vector<std::string> branch_labels(2);

    // for each candidate start nodein the graph
    for (size_t start_node = 0; start_node < dbg.num_nodes(); start_node++) {
		if (start_node % 10000000 == 0) std::cerr<< start_node << " out of " << dbg.num_nodes() <<"\n";
        // if its out degree is two and we haven't encountered it already, start processing it like it's the start of a bubble
        if (!visited[start_node] && dbg.outdegree(start_node) == 2) { //FIXME: why do we only care about outdegree == 2?
            // initialize bubble tracking variables

 

            color_bv branch_color[2];
            size_t end_nodes[2]; // AKA right flank start. place to store end of branch node
            // start of a bubble handling
            int branch_num = -1;
            for (unsigned long x = 1; x < dbg.sigma + 1; x++) { // iterate through the DNA alphabet looking for *the* two outgoing edges from node i
                // follow each strand or supernode
                ssize_t edge = dbg.outgoing_edge(start_node, x);
                if (edge == -1)
                    continue;
                branch_num++;
                branch_labels[branch_num].clear();
                
                branch_labels[branch_num] += base[x];
                // build color mask
                color_bv color_mask = 0;
				//color_bv color_mask2 = 0;
				//if (edge >= 9254208) std::cout<<edge<<":";
				//std::cout<<edge<<":";
				qt = getMilliCount();
                for (uint64_t c = 0; c < num_colors; c++) {
						bool ours = colors.contains(c, edge);					
						color_mask.set(c, ours);
						//if (ours != colors2[edge*num_colors+c]) 
						//		std::cout<<" e"<<edge<<",c"<<c<<":"<<ours<<","<<colors2[edge*num_colors+c]<<"\n";

                     	//color_mask2 |= color_bv(colors2[edge * num_colors + c]) << c;
				}

				//if (color_mask != color_mask2) std::cout <<"\n dif e"<<edge<<":"<<color_mask.count()<<","<<color_mask2.count()<<"\n";
				tqt += getMilliSpan(qt);				
                branch_color[branch_num] = color_mask;

                // walk along edges until we encounter 
                ssize_t pos = dbg._edge_to_node(edge);
                while (dbg.indegree(pos) == 1 && dbg.outdegree(pos) == 1) {
                    visited[pos] = 1;
                    ssize_t next_edge = 0;
                    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) { // iterate through the alphabet
                        next_edge = dbg.outgoing_edge(pos, x2);
                        if (next_edge != -1) {
                            branch_labels[branch_num] += base[x2];
                            break;
                        }
                    }
                    pos = dbg._edge_to_node(next_edge);
                }

                // if we stopped walking along the bubble on a node where indegree > 1, then record this new node
                end_nodes[branch_num] =  (dbg.indegree(pos) > 1) ? pos : 0;
            }
                                    
            // check same end node
            if ((end_nodes[0] && end_nodes[0] == end_nodes[1]) ) {
                // check color:
                if (((color_mask1 & branch_color[0]).any() && (~color_mask1 & branch_color[0]).none() &&
                     (color_mask2 & branch_color[1]).any() && (~color_mask2 & branch_color[1]).none()) || 
                    ((color_mask1 & branch_color[1]).any() && (~color_mask1 & branch_color[1]).none() &&
                     (color_mask2 & branch_color[0]).any() && (~color_mask2 & branch_color[0]).none())) {
                    cout << "\nStart flank: " << dbg.node_label(start_node) << " c: ";
                    print_color ( branch_color[0]);
                    cout << ":";
                    print_color( branch_color[1]);
                    cout << "\n";
                    cout << "Branch: " << branch_labels[0] << "\n";
                    cout << "Branch: " << branch_labels[1] << "\n";
                    cout << "End flank: " << dbg.node_label(end_nodes[0]) << "\n";

                }
            }
        }
    }
	std::cerr << "Find bubbles time: " << getMilliSpan(t) << " ms , query time: " << tqt << " ms" << std::endl;
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
  cerr << "pack-color compiled with supported colors=" << NUM_COLS << endl;
  //ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
   //Can add this to save a couple seconds off traversal - not really worth it.
  cerr << "loading dbg" << std::endl;
  debruijn_graph_shifted<> dbg;
  load_from_file(dbg, p.input_filename);
  cerr << "loading colors" << std::endl;
  uint64_t num_colors = 0;
  uint64_t num_edges = 0;
  deserialize_info(num_colors, num_edges, p.res_dir);
  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << " or " << num_edges << endl;
  cerr << "num colors    : " << num_colors << endl;
  cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;

  //dump_nodes(dbg, colors);
  //dump_edges(dbg, colors);
  color_bv mask1 = (p.color_mask1.length() > 0) ? atoi(p.color_mask1.c_str()) : -1;
  color_bv mask2 = (p.color_mask2.length() > 0) ? atoi(p.color_mask2.c_str()) : -1;

  uint64_t startTime = getMilliCount();
  //uint64_t checkPointTime = getMilliCount();
  std::string res_dir = p.res_dir;
  ColorDetector<T1, T2, T3> cd(res_dir, num_colors);
  //sd_vector<> colors;
  //load_from_file(colors, p.color_filename);
  find_bubbles(dbg, cd,/*colors,*/ mask1, mask2);
  cerr<<" FIND BUBBLE Total time (including loading bitvectors but not dbg) : " << getMilliSpan(startTime) << " ms\n";
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
