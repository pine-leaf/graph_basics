#ifndef __GRAPHS_H__
#define __GRAPHS_H__
#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <stdlib.h>
#include <string>

class ColorTupleHash{
public:
  ColorTupleHash(uint64_t array_size);
  ~ColorTupleHash();
  void insert_color_tuple(uint64_t*, uint64_t);
  uint64_t array_size;
  class ColorTupleNode* get_color_tuple(uint64_t*, uint64_t);
  void remove_color_tuple(uint64_t*, uint64_t);
  class ColorTupleNode** array;
};

class ColorTupleNode{
public:
  ColorTupleNode(uint64_t*, uint64_t);
  ~ColorTupleNode();
  uint64_t* tuple;
  uint64_t length;
  class ColorTupleNode* child;
  class ColorTupleNode* parent;
  uint64_t key;
  void print_me();
};

class MultiColors{
public:
  MultiColors(uint64_t, uint64_t, uint64_t, class Vertex*, class ColorTupleHash*);
  ~MultiColors();
  uint64_t index;
  uint64_t length;
  uint64_t chrom_number;
  uint64_t* color_array;
  uint64_t* color_possibilities_count_array;
  uint64_t** color_possibilities_array;
  class ColorTupleHash* color_tuple_hash;
  class Vertex* vertex;
  class MultiColors** MC;
  int try_next_possibility();
  void zero_me();
  void reset_me();
  void re_restrict_me();
  void print_me();
  void operator=(const class MultiColors&);
private:
  int increment_color_array();
};

class VertexFilter{
public:
  VertexFilter(uint64_t);
  void insert_vertex(class Vertex*);
  void remove_vertex(class Vertex*);
  class VertexNode* is_vertex_in_filter(class Vertex*);
  uint64_t get_size();
private:
  uint64_t array_size;
  uint64_t a;
  uint64_t b;
  uint64_t size;
  class VertexNode** array;
};

class Edge{
friend class Graph;
public:
  Edge(class Vertex*, class Vertex*, double long);
  class Vertex* get_head_vertex();
  class Vertex* get_tail_vertex();
  double long get_weight();
  uint64_t get_color();
  void set_color(uint64_t);
private:
  class Vertex* head;
  class Vertex* tail;
  double long weight;
  uint64_t color; // Let's allow for edge colorings as well!
};

class EdgeNode{
friend class Vertex;
friend class Edge;
EdgeNode(class Edge*);
public:
private:
  class Edge* edge;
  class EdgeNode* child;
};

class VertexNode{
friend class Vertex;
friend class VertexFilter;
public:
VertexNode(class Vertex*);
private:
class Vertex* vertex;
class VertexNode* child;
class VertexNode* parent;
};

class Vertex{
friend class Edge;
friend class Graph;
friend class MultiColors;
public:
  Vertex(uint64_t);
  void print_me();
  void set_color(uint64_t);
  class Vertex* get_next_out();
  void reset_next_out();
  class Vertex* get_next_in();
  void reset_next_in();
  uint64_t get_graph_index();
  uint64_t get_color();
  uint64_t get_in_degree();
  uint64_t get_out_degree();
  void block_color_for_neighbors();
  void un_block_color_for_neighbors();
private:
  uint64_t graph_index;
  uint64_t color;
  uint64_t in_degree;
  uint64_t out_degree;
  uint64_t start_time;
  uint64_t end_time;
  class EdgeNode* next_out;
  class EdgeNode* next_in;
  class EdgeNode* out_edges;
  class EdgeNode* in_edges;
  double long distance;
  class Vertex* previous;
};

class Graph{
public:
  //Graph();
  Graph(const std::string&);
  ~Graph();
  void print_matrix();
  void print_vertices();
  void print_distances();
  void DFS();
  int bellman_ford();
  void vertex_color_graph(uint64_t);
  uint64_t get_n_vertices();
  uint64_t get_graph_index();
  class Vertex* get_ith_vertex(uint64_t i);
private:
  uint64_t n_vertices;
  uint64_t n_edges;
  class Vertex** vertices;
  class Edge** edges;
};

class Set{

};

class SetGraph{

};

#include "graphs.cpp"
#endif