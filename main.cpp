#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "graphs.h"

using namespace std;

int main(int argc, char** argv){
  Graph* G = new Graph("Graph.in");
  /*
  G->DFS();
  G->print_vertices();
  if(G->bellman_ford()){
    std::cout << "There are no negative-weight cycles in G." << std::endl;
  }
  else{
    std::cout << "There IS a negative-weight cycle in G." << std::endl;
  }

  delete G;
  */

  if(argc < 3){
    return 1;
  }

  //MultiColors** MC = multi_color_graph_with_simple_latin_square_seed(G, atoi(argv[1]), atoi(argv[2]));
  MultiColors** MC = multi_color_graph_with_mod_sum_latin_square_seed(G, atoi(argv[1]));
  //print_latin_squares(G, MC, (long long unsigned)atoi(argv[1]));
  //MultiColors** MC = multi_color_graph_with_xor_latin_square_seed(G, atoi(argv[1]));
  return 0;
}