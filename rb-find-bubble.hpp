#ifndef RBBUBBKE_COLOR_H
#define RBBUBBLE_COLOR_H

struct parameters_t {
  std::string input_filename = "";
  std::string color_filename = "";
  std::string output_prefix = "";
  std::string color_mask1 = "";
  std::string color_mask2 = "";
  std::string res_dir = "";
};

int getMilliCount();
int getMilliSpan(int nTimeStart);
void parse_arguments(int argc, char **argv, parameters_t & params);
void test_symmetry(debruijn_graph_shifted<> dbg);
void dump_nodes(debruijn_graph_shifted<> dbg, uint64_t * colors);
void dump_edges(debruijn_graph_shifted<> dbg, uint64_t * colors);
template <class T1, class T2, class T3>
void find_bubbles(debruijn_graph<> dbg, ColorDetector<T1, T2, T3> &colors, uint64_t color_mask1, uint64_t color_mask2);













#endif
