#ifndef _GRAPH_H
#define _GRAPH_H

typedef unsigned int uint;

#include <vector>

#define CHUNK_SIZE (8*(sizeof(int))) // chunk size
class graph
{
  private:
  std::vector<std::vector<int> >  adj;
  uint mask[CHUNK_SIZE];
  std::vector<uint> mapka;
  protected:
  void reorder_custom(const std::vector<uint>& order); // order[i] = new pos of i in [0;n)
  public:
  uint nr_nodes;
  graph(uint n);
  ~graph();
  void add_edge(uint i, uint j);
  inline bool is_edge(uint i, uint j) {return adj[i][j];} 
  inline uint weight(uint i) { return 1; }
  // note: reordering might take a bit of time
  // do before RDS
  void reorder_degree(); // degree order from large to small
  void reorder_weight(); // weight from large to small
  void reorder_random(); // permute at random
  void reorder_2nb(); // order based on the size of 2-neigborhood, from large to small
  void reorder_color(uint s); // See S. Trukhanov et al.
  void reorder_none(); // don't reorder anything
  void reorder_rev(); // revert the order of vertices (usually used to change from small to large)
  void restore_order(std::vector<uint>& v);

  //debug
  void print_mapka() const;
};

graph* from_dimacs(const char* fname); // don't forget to delete

#endif
