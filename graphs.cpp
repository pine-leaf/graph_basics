#ifndef __GRAPHS_CPP__
#define __GRAPHS_CPP__
#include "graphs.h"
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <omp.h>
#include <math.h>

class Vertex* find_next_child_to_color(class Vertex* V, class MultiColors** MC);

void print_latin_squares(class Graph* G, class MultiColors** MC, uint64_t dimension);

uint64_t key_for_color_tuple(uint64_t* array, uint64_t length){
  uint64_t output = 0ULL;
  uint64_t i = 0ULL;
  uint64_t R = 42197ULL;
  for(i; i < length; i++){
    output ^= array[i];
    output += array[i]*(R + i*2);
  }
  return output;
}

uint64_t index_for_color_array(uint64_t* array, uint64_t length, uint64_t chrom_number){
  uint64_t output = 0ULL;
  uint64_t multiplier = 1ULL;
  for(uint64_t i = 0; i < length; i++){
    output += multiplier*(array[0] - 1ULL);
    multiplier *= chrom_number;
  }
  return output;
}

uint64_t b_power_p(uint64_t b, uint64_t p){
  uint64_t output = 1ULL;
  uint64_t tmp = b;
  uint64_t exponent = 1ULL;
  for(int i = 0; i < 64; i++){
    if(p&exponent){
      output *= tmp;
    }
    tmp *= tmp;
    exponent <<= 1;
  }
  return output;
}

uint64_t gcd(uint64_t a, uint64_t b){
  if(a < b){
    return gcd(b,a);
  }
  uint64_t r = a%b;
  if(!r){
    return b;
  }
  return gcd(b, r);
}

Graph::Graph(const std::string& fileName){
  std::string word;
  std::string line;
  std::ifstream file;
  file.open(fileName);
  const char* c_string;
  if(!file.is_open()){
    std::cout << "Coudln't open file: " << fileName << std::endl;
    return;
  }
  uint64_t* edges_array = NULL;
  uint64_t* degree_array = NULL;
  long double* weights_array = NULL;
  while(std::getline(file, line)){
    std::stringstream iss_line(line);
    if(!(iss_line >> word)){continue;}
    c_string = word.c_str();
    if(!strcmp(c_string,"V:")){
      iss_line >> word;
      c_string = word.c_str();
      this->n_vertices = strtoull(c_string, NULL, 10);
      degree_array = new uint64_t[this->n_vertices+1];
      degree_array[0] = 0;
      this->vertices = new Vertex*[this->n_vertices];
      }
    else if(!strcmp(c_string,"E:")){
      iss_line >> word;
      c_string = word.c_str();
      this->n_edges = strtoull(c_string, NULL, 10);
      this->edges = new Edge*[this->n_edges];
      edges_array = new uint64_t[this->n_edges];
      weights_array = (long double*) malloc(sizeof(long double)*this->n_edges);
    }
    else if(!strcmp(c_string, "Edges:")){
      for(uint64_t i = 0ULL; i < this->n_edges; i++){
        if(!(iss_line >> word)){
          break;
        }
        //std::stringstream iss_line(line);
        c_string = word.c_str();
        edges_array[i] = strtoull(c_string, NULL, 10);
        //printf("%s %u ", c_string, edges_array[i]);
      }
    }
    else if(!strcmp(c_string, "Weights:")){
      for(uint64_t i = 0ULL; i < this->n_edges; i++){
        if(!(iss_line >> word)){
          break;
        }
        //std::stringstream iss_line(line);
        c_string = word.c_str();
        //std::cout << strtold(c_string, NULL) << " ";
        weights_array[i] = (long double) strtold(c_string, NULL);
        //weights_array[i] = 10;
        //printf("%s **%f** ", c_string, (float) weights_array[i]);
      }
    }
    else if(!strcmp(c_string, "Degrees:")){
      for(uint64_t i = 0ULL; i < this->n_vertices; i++){
        if(!(iss_line >> word)){
          break;
        }
        //std::stringstream iss_line(line);
        c_string = word.c_str();
        degree_array[i+1] = strtoull(c_string, NULL, 10) + degree_array[i];
      }
    }
  }
  uint64_t h_index;
  uint64_t t_index;
  class Vertex* H;
  class Vertex* T;
  //std::cout << "V: " << this->n_vertices << " E: " << this->n_edges << std::endl;
  class Vertex** vertices_ref = this->vertices;
  int tid;
  #pragma omp parallel for shared(vertices_ref)
  for(uint64_t i = 0ULL; i < this->n_vertices; i++){
    vertices_ref[i] = new Vertex(i);
  }
  #pragma end parallel

  uint64_t e = 0ULL;
  uint64_t j;
  //#pragma omp parallel for shared(vertices_ref, weights_array)
  for(uint64_t i = 0ULL; i < this->n_vertices; i++){
    for(j = degree_array[i]; j < degree_array[i+1]; j++){
      h_index = i;
      t_index = edges_array[j];
      T = vertices_ref[t_index];
      H = vertices_ref[h_index];
      this->edges[e] = new Edge(T,H,weights_array[e]);
      e++;
    }
  }
  //printf("\n");
  delete degree_array;
  delete edges_array;
  free(weights_array);
}

Graph::~Graph(){
  uint64_t i;
  for(i = 0ULL; i < this->n_vertices; i++){
    delete this->vertices[i];
  }
  for(i = 0ULL; i < this->n_edges; i++){
    delete this->edges[i];
  }
  delete[] this->vertices;
  delete[] this->edges;
}

Vertex::Vertex(uint64_t graph_index){
  this->graph_index = graph_index;
  this->color = 0ULL;
  this->in_degree = 0ULL;
  this->out_degree = 0ULL;
  this->start_time = 0ULL;
  this->end_time = 0ULL;
  this->next_out = NULL;
  this->next_in = NULL;
  this->out_edges = NULL;
  this->in_edges = NULL;
  this->distance = 0.0f;
  this->previous = NULL;
}

void Vertex::print_me(){
  std::cout << this->graph_index << ": Color: " << this->color << " in_degree: " << this->in_degree << " out_degree: " << this->out_degree << " start_time: " << this->start_time << " end_time: " << this->end_time << " distance: " << this->distance << " previous: ";
    if(this->previous){
      std::cout << this->previous->graph_index << std::endl;
    }else{
      std::cout << "NULL" << std::endl;
    }
}

Edge::Edge(class Vertex* H, class Vertex* T, double long weight){
  // This function is NOT safe! I.e. it will create a duplicate edge.
  H->in_degree++;
  T->out_degree++;
  this->head = H;
  this->tail = T;
  this->weight = weight;
  class EdgeNode* H_node = new EdgeNode(this);
  class EdgeNode* T_node = new EdgeNode(this);
  class EdgeNode* tmp = H->in_edges;
  H->in_edges = H_node;
  H->next_in = H_node;
  H_node->child = tmp;
  tmp = T->out_edges;
  T->out_edges = T_node;
  T->next_out = T_node;
  T_node->child = tmp;
  this->color = 0ULL;
}

void Edge::set_color(uint64_t c){
  this->color = c;
}

uint64_t Edge::get_color(){
  return this->color;
}

EdgeNode::EdgeNode(class Edge* E){
  this->child = NULL;
  this->edge = E;
}

void Vertex::set_color(uint64_t c){
  this->color = c;
}

class Vertex* Vertex::get_next_out(){
  class EdgeNode* N = this->next_out;
  class Vertex* output = NULL;
  if(N){
    output = N->edge->get_head_vertex();
    this->next_out = N->child;
  }
  else{
    this->next_out = this->out_edges;
  }
  return output;
}

void Vertex::reset_next_out(){
  this->next_out = this->out_edges;
}

class Vertex* Vertex::get_next_in(){
  class EdgeNode* N = this->next_in;
  class Vertex* output = NULL;
  if(N){
    output = N->edge->get_tail_vertex();
    this->next_in = N->child;
  }
  else{
    this->next_in = this->in_edges;
  }
  return output;
}

void Vertex::reset_next_in(){
  this->next_in = this->in_edges;
}

uint64_t Vertex::get_graph_index(){
  return this->graph_index;
}

uint64_t Vertex::get_color(){
  return this->color;
}

uint64_t Vertex::get_out_degree(){
  return this->out_degree;
}

uint64_t Vertex::get_in_degree(){
  return this->in_degree;
}

class Vertex* Edge::get_head_vertex(){
  return this->head;
}

class Vertex* Edge::get_tail_vertex(){
  return this->tail;
}

long double Edge::get_weight(){
  return this->weight;
}

void Graph::DFS(){ // can't really be done in parallel
  uint64_t i;
  uint64_t nv = this->n_vertices;
  uint64_t* color_array = new uint64_t[nv];
  class Vertex* V;
  class Vertex* tmp;
  #pragma omp parallel for shared(color_array) private(V)
  for(i = 0ULL; i < nv; i++){
    V = this->vertices[i];
    color_array[i] = V->color;
    V->color = 0ULL;
    V->start_time = 0ULL;
    V->end_time = 0ULL;
    V->previous = NULL;
  }
  uint64_t index = 0ULL;
  uint64_t time = 0ULL;
  uint64_t finished = 0ULL;
  V = NULL;
  tmp = NULL;
  while(finished < nv){
    time++;
    if(!V){
      V = this->vertices[index];
      index++;
      continue;
    }
    switch(V->color){
      case 0ULL:
        V->color = 1ULL;
        V->start_time = time;
        continue;
      case 1ULL:
        tmp = V->get_next_out();
        if(tmp){
          if(!tmp->color){
            tmp->previous = V;
            V = tmp;
          }
          continue;
        }
        else{
          V->color = 2ULL;
          V->end_time = time;
          V = V->previous;
          finished++;
          continue;
        }
      case 2ULL:
        V = V->previous;
    }
  }
}

int Graph::bellman_ford(){
  uint64_t nv = this->n_vertices;
  uint64_t ne = this->n_edges;
  if(!nv || !ne){
    return 1;
  }
  uint64_t* color_array = new uint64_t[nv];
  class Edge** edge_array = this->edges;
  class Vertex** vertex_array = this->vertices;
  //class Vertex** vertex_heap = new Vertex*[nv];
  uint64_t i;
  #pragma omp parallel for shared(color_array)
  for(i = 0ULL; i < nv; i++){
    //vertex_heap[i] = NULL;
    color_array[i] = this->vertices[i]->color;
    this->vertices[i]->color = 0ULL;
    this->vertices[i]->distance = 0.0f;
  }
  #pragma end parallel

  this->vertices[0]->color = 1ULL;

  uint64_t e = 0ULL;
  uint64_t h_index;
  class Edge* E;
  class Vertex* T;
  class Vertex* H;
  int tid;
  for(i = 0ULL; i < nv-1; i++){
    #pragma omp parallel for shared(edge_array, vertex_array) private(E, tid, T, H, h_index)
    for(e = 0ULL; e < ne; e++){
      //tid = omp_get_thread_num();
      E = edge_array[e];
      T = E->tail;
      H = E->head;
      if(!T->color){
        continue;
      }
      if(!H->color){
        h_index = H->graph_index;
        vertex_array[h_index]->distance = T->distance + E->weight;
        H->previous = T;
        H->color = 1ULL;
      }
      else if(H->distance > T->distance + E->weight){
        h_index = H->graph_index;
        vertex_array[h_index]->distance = T->distance + E->weight;
        H->previous = T;
      }
    }
  }

  int output = 1;
  #pragma omp parallel for shared(vertex_array, output) private(E)
  for(e = 0ULL; e < ne; e++){
    E = edge_array[e];
    H = E->head;
    if(!H->color){
      continue;
    }
    T = E->tail;
    if(H->distance > T->distance + E->weight){
      output = 0;
    }
  }

  #pragma omp parallel for shared(vertex_array, color_array) private(T)
  for(i = 0ULL; i < nv; i++){
    T = vertex_array[i];
    T->color = color_array[i];
  }

  return output;
}

void vertex_heap_distance_min_heapify_up(class Vertex** vertex_heap, uint64_t start, uint64_t end){
  if(!vertex_heap){
    return;
  }

  uint64_t L = 2*start + 1;
  uint64_t R = L + 1;
  uint64_t swap_index = start;

  if(L < end){
    if(vertex_heap[L] < vertex_heap[swap_index]){
      swap_index = L;
    }
  }
  if(R < end){
    if(vertex_heap[R] < vertex_heap[swap_index]){
      swap_index = R;
    }
  }

  if(swap_index != start){
    class Vertex* tmp = vertex_heap[start];
    vertex_heap[start] = vertex_heap[swap_index];
    vertex_heap[swap_index] = tmp;
    vertex_heap_distance_min_heapify_up(vertex_heap, swap_index, end);
  }
}

void Graph::print_distances(){
  uint64_t i = 0ULL;
  class Vertex* V;
  for(;i < this->n_vertices; i++){
    V = this->vertices[i];
    printf("%i -> %f\n", i, (float)V->distance);
  }
}

void Graph::print_vertices(){
  uint64_t i = 0ULL;
  for(; i < this->n_vertices; i++){
    this->vertices[i]->print_me();
  }
}

void Graph::vertex_color_graph(uint64_t n){
  // Finish this later, I guess...
}

class Vertex* find_next_vertex(class Graph* G, uint64_t length, uint64_t** color_grid){
  uint64_t max_set = 0ULL;
  uint64_t tmp_set;
  uint64_t max_degrees = 0ULL;
  uint64_t i,j;
  uint64_t output = 0ULL;
  uint64_t nv = G->get_n_vertices();
  if(!G->get_n_vertices() || !length){
    return NULL;
  }
  class Vertex* V;
  for(i = 0ULL; i < nv; i++){
    V = G->get_ith_vertex(i);
    tmp_set = 0ULL;
    for(j = 0; j < length; j++){
      if(color_grid[i][j]){
        tmp_set++;
      }
    }
    if(tmp_set > max_set){
      max_set = tmp_set;
      output = i;
    }
    else if(tmp_set == max_set){
      if(V->get_out_degree() + V->get_in_degree() > max_degrees){
        max_degrees = V->get_out_degree() + V->get_in_degree();
        output = i;
      }
    }
  }
  return G->get_ith_vertex(output);
}

void print_previous_array(class Vertex** previous_array, uint64_t length){
  uint64_t i = 0ULL;
  for(; i < length; i++){
    std::cout << i << "-> ";
    if(previous_array[i]){
      std::cout << previous_array[i]->get_graph_index() << " ";
    }
    else{
      std::cout << "NULL ";
    }
  }
  std::cout << std::endl;
}

class MultiColors** multi_color_graph_with_multi_colors(class Graph* G, uint64_t chrom_number, uint64_t length, class MultiColors** multi_colors){
  uint64_t nv = G->get_n_vertices();
  uint64_t i,j,k;

  if(!multi_colors){
    return NULL;
  }
  class MultiColors** MC = new MultiColors*[nv];
  for(i = 0; i < nv; i++){
    MC[i] = new MultiColors(0,0,0,NULL,NULL);
    *MC[i] = *multi_colors[i];
    MC[i]->MC = MC;
  }
  class Vertex* V;
  class Vertex* previous_array[nv];
  for(i = 0ULL; i < nv; i++){
    previous_array[i] = NULL;
  }

  // Do a bunch of stuff here with DFS-type things...
  // Initialize V:
  V = G->get_ith_vertex(0ULL);

  class Vertex* C = NULL;
  class Vertex* P = NULL;
  class Vertex* Latest = NULL;
  class MultiColors* M = NULL;
  uint64_t index = 0ULL;
  uint64_t dimension;
  uint64_t colored = 0ULL;
  uint64_t max_colored = 0ULL;
  uint64_t time = 0ULL;
  int back_tracking = 0;

  for(i = 0ULL; i < nv; i++){
    if(MC[i]->color_array[0]){
      colored++;
    }
  }
  Latest = NULL;
  dimension = colored;

  while(colored < nv){
    time++;
    if(!V){
      back_tracking = 0;
      V = G->get_ith_vertex(index);
      M = MC[index];
      if(M->color_array[0]){
        V = NULL;
      }
      index++;
      if(index >= nv){
        break;
      }
      continue;
    }
    if(colored > max_colored){
      max_colored = colored;
      std::cout << max_colored << ", " << time << std::endl;
    }
    // We know that V exists here!
    M = MC[V->get_graph_index()];

    if(back_tracking){
      if(M->color_array[0]){
        colored--;
      }
      M->reset_me();
      if(M->try_next_possibility()){
        // we got a coloring that worked, so let's go with it!
        colored++;
        if(colored == nv){
          std::cout << "Time: " << time << ", V: " << V->get_graph_index() << std::endl; 
          print_latin_squares(G, MC, dimension);
          std::cout << "---------------------" << std::endl;
          continue;
        }
        back_tracking = 0; // we got a coloring that works, so don't back-track now
        C = find_next_child_to_color(V, MC);
        if(C){
          // C exists, so let's go there and try to color it!
          previous_array[C->get_graph_index()] = V;
          Latest = V;
          V = C;
          continue;
        }
        else{
          // There are no more children, so let's just go back and look for more vertices
          P = previous_array[V->get_graph_index()];
          V = P;
          continue;
        }
      }
      else{
        // Still can't find a good coloring... must back-track even further!
        *M = *multi_colors[V->get_graph_index()];
        V->reset_next_out();
        P = previous_array[V->get_graph_index()];
        if(Latest){
          Latest = previous_array[Latest->get_graph_index()];
        }
        else{
          //std::cout << "Latest is NULL. V: " << V->get_graph_index() << std::endl;
          V = NULL;
          continue;
        }
        if(P){
          P->reset_next_out();
        }
        previous_array[V->get_graph_index()] = NULL;
        V = P;
        continue;
      }
    }
    else{
      // We know we're going forward!
      if(!M->color_array[0]){
        // M hasn't been colored, so let's try to color it!
        M->reset_me();
        if(M->try_next_possibility()){
          // We colored M!!
          colored++;
          if(colored == nv){
            std::cout << "Time: " << time << ", V: " << V->get_graph_index() << std::endl; 
            print_latin_squares(G, MC, dimension);
            std::cout << "---------------------" << std::endl;
            back_tracking = 1;
            continue;
          }
          C = find_next_child_to_color(V, MC);
          if(C){
            // C exists, so let's go forward with C
            previous_array[C->get_graph_index()] = V;
            Latest = V;
            V = C;
            continue;
          }
          else{
            // V has no more uncolored children, so let's go back and look for some
            P = previous_array[V->get_graph_index()];
            V = P;
            continue;
          }
        }
        else{
          // Can't color M, so we need to back-track!!
          *M = *multi_colors[V->get_graph_index()];
          V->reset_next_out();
          back_tracking = 1;
          P = previous_array[V->get_graph_index()];
          if(Latest){
            Latest = previous_array[Latest->get_graph_index()];
          }
          else{
            V = NULL;
            continue;
          }
          if(P){
            P->reset_next_out();
          }
          previous_array[V->get_graph_index()] = NULL;
          V = P;
          continue;
        }
      }
      else{
        // V HAS been colored, but we're okay, so let's just look for the next child;
        C = find_next_child_to_color(V, MC);
        if(C){
          // C exists, so let's go forward with C
          previous_array[C->get_graph_index()] = Latest;
          //Latest = V;
          V = C;
          continue;
        }
        else{
          // V has no more uncolored children, so let's go back and look for some
          P = previous_array[V->get_graph_index()];
          V = P;
          continue;
        }
      }
    }
  }
  std::cout << "Time: " << time << std::endl;
  std::cout << "Colored: " << colored << std::endl;

  for(i = 0ULL; i < nv; i++){
    MC[i]->color_tuple_hash = NULL;
  }

  return MC;
}

class MultiColors** multi_color_graph(class Graph* G, uint64_t chrom_number, uint64_t length){
  uint64_t nv = G->get_n_vertices();
  uint64_t i,j,k;
  class ColorTupleHash* H = new ColorTupleHash(nv);
  class MultiColors** MC = new MultiColors*[nv];
  class Vertex* V;
  class Vertex* previous_array[nv];
  for(i = 0ULL; i < nv; i++){
    V = G->get_ith_vertex(i);
    MC[i] = new MultiColors(length, chrom_number, i, V, H);
    MC[i]->MC = MC;
    previous_array[i] = NULL;
  }

  return multi_color_graph_with_multi_colors(G, chrom_number, length, MC);

}

class MultiColors** simple_latin_square_seed(class Graph* G, uint64_t chrom_number, uint64_t length){
  if(!G || !chrom_number || !length){
    return NULL;
  }
  uint64_t nv = G->get_n_vertices();
  uint64_t i,j,k;
  class ColorTupleHash* H = new ColorTupleHash(nv);
  class MultiColors** MC = new MultiColors*[nv];
  class Vertex* V;
  //class Vertex* previous_array[nv];
  for(i = 0ULL; i < nv; i++){
    V = G->get_ith_vertex(i);
    MC[i] = new MultiColors(length, chrom_number, i, V, H);
    MC[i]->MC = MC;
    //previous_array[i] = NULL;
  }
  for(i = 0ULL; i < chrom_number; i++){
    class MultiColors* M = MC[i];
    for(j = 0; j < length; j++){
      M->color_array[j] = i + 1ULL;
      M->color_possibilities_count_array[j] = 1ULL;
      M->color_possibilities_array[j][0] = 0ULL;
      for(k = 1ULL; k <= chrom_number; k++){
        if(k != i + 1ULL){
          M->color_possibilities_array[j][k] = 0ULL;
        }
      }
      H->insert_color_tuple(M->color_array, M->length);
    }
  }
  return MC;
}

class MultiColors** mod_sum_latin_square_seed(class Graph* G, uint64_t chrom_number){
  uint64_t length = 2ULL;
  if(!G || !chrom_number){
    return NULL;
  }
  uint64_t nv = G->get_n_vertices();
  uint64_t i,j,k,s,l,m,n,c;
  class ColorTupleHash* H = new ColorTupleHash(nv);
  class MultiColors** MC = new MultiColors*[nv];
  class Vertex* V;
  class MultiColors* M;

  for(i = 0ULL; i < nv; i++){
    V = G->get_ith_vertex(i);
    MC[i] = new MultiColors(length, chrom_number, i, V, H);
    MC[i]->MC = MC;
  }

  for(i = 0ULL; i < chrom_number; i++){
    M = MC[i];
    for(j = 0; j < length; j++){
      M->color_array[j] = i + 1ULL;
      M->color_possibilities_count_array[j] = 1ULL;
      M->color_possibilities_array[j][0] = 0ULL;
      for(k = 1ULL; k <= chrom_number; k++){
        if(k != i + 1ULL){
          M->color_possibilities_array[j][k] = 0ULL;
        }
      }
      H->insert_color_tuple(M->color_array, M->length);
    }
  }
  for(i = chrom_number; i < nv; i++){
    M = MC[i];
    m = i/chrom_number;
    n = i%chrom_number;
    s = (m + n)%chrom_number + 1ULL;
    for(c = 1ULL; c <= chrom_number; c++){
      if(c == s){
        M->color_possibilities_array[0][c] = 1ULL;
      }
      else{
        M->color_possibilities_array[0][c] = 0ULL;
      }
    }
    M->color_possibilities_count_array[0] = 1ULL;
  }
  return MC;
}

class MultiColors** xor_latin_square_seed(class Graph* G, uint64_t chrom_number){
  uint64_t length = 2ULL;
  if(!G || !chrom_number){
    return NULL;
  }
  uint64_t nv = G->get_n_vertices();
  uint64_t i,j,k,s,l,m,n,c;
  class ColorTupleHash* H = new ColorTupleHash(nv);
  class MultiColors** MC = new MultiColors*[nv];
  class Vertex* V;
  class MultiColors* M;

  for(i = 0ULL; i < nv; i++){
    V = G->get_ith_vertex(i);
    MC[i] = new MultiColors(length, chrom_number, i, V, H);
    MC[i]->MC = MC;
  }

  for(i = 0ULL; i < chrom_number; i++){
    M = MC[i];
    for(j = 0; j < length; j++){
      M->color_array[j] = i + 1ULL;
      M->color_possibilities_count_array[j] = 1ULL;
      M->color_possibilities_array[j][0] = 0ULL;
      for(k = 1ULL; k <= chrom_number; k++){
        if(k != i + 1ULL){
          M->color_possibilities_array[j][k] = 0ULL;
        }
      }
      H->insert_color_tuple(M->color_array, M->length);
    }
  }
  for(i = chrom_number; i < nv; i++){
    M = MC[i];
    m = i/chrom_number;
    n = i%chrom_number;
    s = (m ^ n)%chrom_number + 1ULL;
    for(c = 1ULL; c <= chrom_number; c++){
      if(c == s){
        M->color_possibilities_array[0][c] = 1ULL;
      }
      else{
        M->color_possibilities_array[0][c] = 0ULL;
      }
    }
    M->color_possibilities_count_array[0] = 1ULL;
  }
  return MC;
}

class MultiColors** multi_color_graph_with_simple_latin_square_seed(class Graph* G, uint64_t chrom_number, uint64_t length){
  if(!G){
    return NULL;
  }
  class MultiColors** MC = simple_latin_square_seed(G, chrom_number, length);
  return multi_color_graph_with_multi_colors(G, chrom_number, length, MC);
}

class MultiColors** multi_color_graph_with_mod_sum_latin_square_seed(class Graph* G, uint64_t chrom_number){
  if(!G){
    return NULL;
  }
  class MultiColors** MC = mod_sum_latin_square_seed(G, chrom_number);
  return multi_color_graph_with_multi_colors(G, chrom_number, 2ULL, MC);
}

class MultiColors** multi_color_graph_with_xor_latin_square_seed(class Graph* G, uint64_t chrom_number){
  if(!G){
    return NULL;
  }
  class MultiColors** MC = xor_latin_square_seed(G, chrom_number);
  return multi_color_graph_with_multi_colors(G, chrom_number, 2ULL, MC);
}

VertexNode::VertexNode(class Vertex* V){
  this->vertex = V;
  this->child = NULL;
  this->parent = NULL;
}

VertexFilter::VertexFilter(uint64_t array_size){
  this->array_size = array_size;
  this->array = new VertexNode*[array_size];
  #pragma omp parallel for shared(array)
  for(uint64_t i = 0ULL; i < array_size; i++){
    this->array[i] = NULL;
  }
  uint64_t a, b;
  a = array_size/5;
  b = a;
  while(gcd(array_size, a) != 1ULL){
    a++;
  }
  this->a = a;
  this->b = b;
}

void VertexFilter::insert_vertex(class Vertex* V){
  class VertexNode* N = new VertexNode(V);
  uint64_t index = V->get_graph_index();
  index %= this->array_size;
  index *= this->a;
  index %= this->array_size;
  index += this->b;
  index %= this->array_size;
  class VertexNode* tmpNode = this->array[index];
  this->array[index] = N;
  N->child = tmpNode;
  if(N->child){
    N->child->parent = N;
  }
  this->size++;
}

class VertexNode* VertexFilter::is_vertex_in_filter(class Vertex* V){
  uint64_t index = V->get_graph_index();
  index %= this->array_size;
  index *= this->a;
  index %= this->array_size;
  index += this->b;
  index %= this->array_size;
  class VertexNode* N = this->array[index];
  if(!N){
    return NULL;
  }
  while(N){
    if(N->vertex == V){
      return N;
    }
    N = N->child;
  }
  return NULL;
}

void VertexFilter::remove_vertex(class Vertex* V){
  class VertexNode* N = is_vertex_in_filter(V);
  if(!N){
    return;
  }
  class VertexNode* P = N->parent;
  class VertexNode* C = N->child;
  if(P){
    P->child = C;
  }
  if(C){
    C->parent = P;
  }
  delete N;
}

class Vertex* Graph::get_ith_vertex(uint64_t i){
  if(i >= this->n_vertices){
    return NULL;
  }
  return this->vertices[i];
}

uint64_t Graph::get_n_vertices(){
  return this->n_vertices;
}

ColorTupleNode::ColorTupleNode(uint64_t* array, uint64_t length){
  this->parent = NULL;
  this->child = NULL;
  this->key = key_for_color_tuple(array, length);
  this->length = length;
  this->tuple = new uint64_t[length];
  for(uint64_t i = 0ULL; i < length; i++){
    this->tuple[i] = array[i];
  }
}

ColorTupleNode::~ColorTupleNode(){
  delete this->tuple;
}

ColorTupleHash::ColorTupleHash(uint64_t array_size){
  this->array_size = array_size;
  this->array = new ColorTupleNode*[array_size];
  for(uint64_t i = 0ULL; i < array_size; i++){
    this->array[i] = NULL;
  }
}

ColorTupleHash::~ColorTupleHash(){
  delete this->array;
}

ColorTupleNode* ColorTupleHash::get_color_tuple(uint64_t* tuple, uint64_t length){
  if(length <= 1ULL){
    return NULL;
  }
  uint64_t key = key_for_color_tuple(tuple, length);
  uint64_t i;
  uint64_t index = key%(this->array_size);
  int equal = 1;
  class ColorTupleNode* N = this->array[index];
  class ColorTupleNode* tmp;
  while(N){
    tmp = N;
    N = N->child;
    equal = 1;
    for(i = 0ULL; i < length; i++){
      if(tmp->tuple[i] != tuple[i]){
        equal = 0;
        break;
      }
    }
    if(!equal){
      continue;
    }
    N = tmp;
    break;
  }
  return N;
}

void ColorTupleHash::insert_color_tuple(uint64_t* tuple, uint64_t length){
  if(length <= 1ULL){
    return;
  }
  class ColorTupleNode* N = new ColorTupleNode(tuple, length);
  uint64_t index = (N->key)%(this->array_size);
  class ColorTupleNode* H = this->array[index];
  this->array[index] = N;
  N->child = H;
  if(H){
    H->parent = N;
  }
}

void ColorTupleHash::remove_color_tuple(uint64_t* tuple, uint64_t length){
  uint64_t i;
  class ColorTupleNode* N = this->get_color_tuple(tuple, length);
  if(!N){
    return;
  }
  uint64_t index = (N->key)%(this->array_size);
  if(N == this->array[index]){
    if(!N->child){
      this->array[index] = NULL;
      delete N;
      return;
    }
    else{
      this->array[index] = N->child;
      N->child->parent = NULL;
      delete N;
      return;
    }
  }
  class ColorTupleNode* P = N->parent;
  class ColorTupleNode* C = N->child;
  if(P){
    P->child = C;
  }
  if(C){
    C->parent = P;
  }
  delete N;
}

MultiColors::MultiColors(uint64_t length, uint64_t chrom_number, uint64_t index, class Vertex* V, class ColorTupleHash* H){
  this->vertex = V;
  this->color_tuple_hash = H;
  this->color_array = new uint64_t[length];
  this->color_possibilities_count_array = new uint64_t[length];
  this->index = index;
  this->length = length;
  this->chrom_number = chrom_number;
  this->MC = NULL;
  uint64_t i,j;
  for(i = 0ULL; i < length; i++){
    this->color_array[i] = 0ULL;
    this->color_possibilities_count_array[i] = chrom_number;
  }
  this->color_possibilities_array = new uint64_t*[length];
  for(i = 0; i < length; i++){
    this->color_possibilities_array[i] = new uint64_t[chrom_number+1ULL];
    this->color_possibilities_array[i][0] = 0ULL;
    for(j = 1ULL; j < chrom_number + 1ULL; j++){
      this->color_possibilities_array[i][j] = 1ULL;
    }
  }
}

MultiColors::~MultiColors(){
  uint64_t i,j;
  if(this->color_possibilities_array){
    for(i = 0ULL; i <= this->length; i++){
      if(this->color_possibilities_array[i]){
        delete this->color_possibilities_array[i];
      }
    }
    delete this->color_possibilities_array;
  }
  if(this->color_array){
    delete this->color_array;
  }
  if(this->color_possibilities_count_array){
    delete this->color_possibilities_count_array;
  }
}

int MultiColors::increment_color_array(){
  uint64_t i,j;
  if(!this->length){
    return 1;
  }
  // this loop will find the first possible color_array
  if(!this->color_array[0]){
    for(i = 0ULL; i < this->length; i++){
      for(j = 1ULL; j <= this->chrom_number; j++){
        if(this->color_possibilities_array[i][j]){
          this->color_array[i] = j;
          break;
        }
      }
      if(this->color_array[i] == 0ULL){
        // Making it here means that there are NO legal colorings.
        for(j = 0ULL; j < this->length; j++){
          this->color_array[j] = 0ULL;
        }
        return 0;
      }
    }
    return 1;
  }
  // add 1
  uint64_t high, low;
  for(i = 0ULL; i < this->length; i++){
    high = this->color_array[i];
    for(j = high + 1ULL; j <= this->chrom_number; j++){
      if(this->color_possibilities_array[i][j]){
        high = j;
        break;
      }
    }
    if(this->color_array[i] != high){
      this->color_array[i] = high;
      return 1;
    }
    low = this->color_array[i];
    for(j = 1ULL; j <= this->color_array[i]; j++){
      if(this->color_possibilities_array[i][j]){
        low = j;
        break;
      }
    }
    this->color_array[i] = low;
  }
  return 0;
}

/*
  try_next_possibility will search for the next "multi-coloring"
  that is consistent with all of the given vertex's neighbors
  multi-colorings. If it finds one, it returns 1, else it returns 0
*/
int MultiColors::try_next_possibility(){
  uint64_t* colors = this->color_array;
  uint64_t i,j;
  class ColorTupleHash* H = this->color_tuple_hash;
  H->remove_color_tuple(this->color_array, this->length);
  if(!this->increment_color_array())
  {
    return 0;
  }
  class ColorTupleNode* N = H->get_color_tuple(this->color_array, this->length);
  while(N){
    if(!this->increment_color_array()){
      return 0;
    }
    N = H->get_color_tuple(this->color_array, this->length);
  }
  H->insert_color_tuple(this->color_array, this->length);
  return 1;
}

void MultiColors::re_restrict_me(){
  if(!this->MC){
    return;
  }
  class Vertex* V = this->vertex;
  class EdgeNode* EN = V->next_out;
  class Vertex* neighbor = V->get_next_out();
  class MultiColors* M = NULL;
  uint64_t neighbor_index;
  uint64_t i,j;
  for(i = 0ULL; i < this->length; i++){
    this->color_possibilities_count_array[i] = this->chrom_number;
    for(j = 1ULL; j <= this->chrom_number; j++){
      this->color_possibilities_array[i][j] = 1ULL;
    }
  }
  while(V->next_out != EN){
    if(neighbor){
      neighbor_index = neighbor->graph_index;
      M = this->MC[neighbor_index];
      for(i = 0ULL; i < this->length; i++){
        if(this->color_possibilities_array[i][M->color_array[i]]){
          this->color_possibilities_array[i][M->color_array[i]] = 0ULL;
          this->color_possibilities_count_array[i]--;
        }
      }
    }
    neighbor = V->get_next_out();
  }
}

void MultiColors::reset_me(){
  if(!this->MC){
    return;
  }
  class Vertex* V = this->vertex;
  class EdgeNode* EN = V->next_out;
  class Vertex* neighbor = V->get_next_out();
  class MultiColors* M = NULL;
  uint64_t neighbor_index;
  uint64_t neighbor_color;
  uint64_t i;
  while(V->next_out != EN){ // This way we leave V unchanged when we exit.
    if(neighbor){
      neighbor_index = neighbor->graph_index;
      M = this->MC[neighbor_index];
      for(i = 0ULL; i < this->length; i++){
        neighbor_color = M->color_array[i];
        if(!neighbor_color){
          i = this->length;
          continue;
        }
        if(this->color_possibilities_array[i][neighbor_color]){
          this->color_possibilities_count_array[i]--;
          this->color_possibilities_array[i][neighbor_color] = 0ULL;
        }
      }
    }
    neighbor = V->get_next_out();
  }
}

void MultiColors::zero_me(){
  uint64_t i, j;
  class ColorTupleHash* H = this->color_tuple_hash;
  for(i = 0ULL; i < this->length; i++){
    this->color_array[i] = 0ULL;
    this->color_possibilities_count_array[i] = this->chrom_number;
    this->color_possibilities_array[i][0] = 0ULL;
    for(j = 1ULL; j <= this->chrom_number; j++){
      this->color_possibilities_array[i][j] = 1ULL;
    }
  }
  class Vertex* V = this->vertex;
  class EdgeNode* EN = V->next_out;
  class Vertex* neighbor = V->get_next_out();
  class MultiColors* M = NULL;
  uint64_t neighbor_index;
  while(V->next_out != EN){ // This way we leave V unchanged when we exit.
    if(neighbor){
      neighbor_index = neighbor->graph_index;
      M = this->MC[neighbor_index];
      M->re_restrict_me();
    }
    neighbor = V->get_next_out();
  }
}

uint64_t score_multi_color(class MultiColors* M){
  uint64_t output = 0;
  uint64_t i;
  class Vertex* V = M->vertex;
  if(!V){
    std::cout << "V is NULL." << std::endl;
    return 0ULL;
  }

  output += V->get_out_degree() + V->get_in_degree();

  for(i = 0ULL; i < M->length; i++){
    switch(M->color_possibilities_count_array[i]){
      case 1ULL:
        output += 100ULL;
        break;
      case 2ULL:
        output += 30ULL;
        break;
      case 3ULL:
        output += 10ULL;
        break;
      case 4ULL:
        output += 3ULL;
        break;
      case 5ULL:
        output += 1ULL;
        break;
    }
  }

  return output;
}

class Vertex* find_next_child_to_color(class Vertex* V, class MultiColors** MC){
  if(!V || !MC){
    return NULL;
  }
 
  V->reset_next_out();
  class Vertex* C = V->get_next_out();
  class Vertex* B = NULL;
  if(!C){
    return NULL;
  }
  uint64_t C_index = C->get_graph_index();
  uint64_t max_score = 0ULL;
  uint64_t tmp_score;
  class MultiColors* M = MC[C_index];
  while(C){
    C_index = C->get_graph_index();
    M = MC[C_index];
    if(!M->color_array[0]){
      tmp_score = score_multi_color(M);
      if(tmp_score >= max_score){
        max_score = tmp_score;
        B = C;
      }
    }
    C = V->get_next_out();
  }
  return B;
}

void MultiColors::print_me(){
  if(!this->length){
    return;
  }
  uint64_t i,j;
  std::cout << "color_array: " << std::endl;
  for(i = 0ULL; i < this->length; i++){
    std::cout << this->color_array[i] << " ";
  }
  std::cout << "color_possibilities_array:" << std::endl;
  for(i = 0;i < length; i++){
    for(j = 0; j <= chrom_number; j++){
      std::cout << this->color_possibilities_array[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void ColorTupleNode::print_me(){
  uint64_t i;
  std::cout << "Node: ";
  for(i = 0ULL; i < this->length; i++){
    std::cout << this->tuple[i] << " ";
  }
  std::cout << std::endl;
}

void print_latin_squares(class Graph* G, class MultiColors** MC, uint64_t dimension){
  uint64_t i, j, s;
  if(!G || !MC || !dimension){
    return;
  }
  if(!MC[0]){
    return;
  }
  if(G->get_n_vertices() != dimension*dimension){
    std::cout << "Cannot print Latin Squares. Dimension is wrong." << std::endl;
    return;
  }
  class MultiColors* M;
  uint64_t length = MC[0]->length;
  for(s = 0ULL; s < length; s++){
    printf("Latin Square %llu:\n", s);
    for(i = 0ULL; i < dimension; i++){
      for(j = 0ULL; j < dimension; j++){
        M = MC[dimension*i + j];
        printf("%3llu", M->color_array[s]);
      }
      printf("\n");
    }
  }
}

void MultiColors::operator=(const class MultiColors& M){
  if(&M == this){
    return;
  }
  uint64_t i,j;
  if((this->chrom_number != M.chrom_number) || (this->length != M.length)){
    if(this->color_possibilities_array){
      for(i = 0ULL; i < this->length; i++){
        if(this->color_possibilities_array[i]){
          delete this->color_possibilities_array[i];
        }
      }
      delete this->color_possibilities_array;
    }
    if(this->color_array){
      delete this->color_array;
    }
    if(this->color_possibilities_count_array){
      delete this->color_possibilities_count_array;
    }
    this->color_possibilities_array = new uint64_t*[M.length];
    this->color_array = new uint64_t[M.length];
    this->color_possibilities_count_array = new uint64_t[M.length];
    for(i = 0ULL; i < M.length; i++){
      this->color_possibilities_array[i] = new uint64_t[M.chrom_number+1ULL];
    }
    this->length = M.length;
    this->chrom_number = M.chrom_number;
  }
  for(i = 0ULL; i < this->length; i++){
    this->color_possibilities_array[i][0] = 0ULL;
    for(j = 1ULL; j <= M.chrom_number; j++){
      this->color_possibilities_array[i][j] = M.color_possibilities_array[i][j];
    }
    this->color_array[i] = M.color_array[i];
    this->color_possibilities_count_array[i] = M.color_possibilities_count_array[i];
  }
  this->color_tuple_hash = M.color_tuple_hash;
  this->vertex = M.vertex;
  //this->MC = M.MC;
  this->index = M.index;
}

#endif