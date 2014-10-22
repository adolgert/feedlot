#include "boost/graph/graphml.hpp"
#include "pen.hpp"


// This kind of adjacency list uses integers for the vertex id, which
// is special to those using vecS and vecS.
using PenContactGraph=boost::adjacency_list<boost::vecS,
    boost::vecS,boost::undirectedS>;
/*! Feedlots look like suburbs. There are long blocks with streets
 *  between. Each block has row_cnt pens, sitting back-to-back.
 *  Returns block_cnt*row_cnt*2 pens.
 */
PenContactGraph BlockStructure(int block_cnt, int row_cnt) {
  PenContactGraph g(block_cnt*row_cnt*2);
  for (int bidx=0; bidx<2*block_cnt; ++bidx) {
    int base=bidx*row_cnt;
    for (int right_idx=0; right_idx<row_cnt-1; ++right_idx) {
      add_edge(base+right_idx, base+right_idx+1, g);
    }
  }
  return g;
}

bool AdjacentPens(int i, int j, const PenContactGraph& g) {
  using AdjIter=boost::graph_traits<PenContactGraph>::adjacency_iterator;
  AdjIter start, end;
  assert(i<num_vertices(g));
  std::tie(start, end)=adjacent_vertices(i, g);
  for ( ; start!=end; ++start) {
    if (*start==j) {
      return true;
    }
  }
  return false;
}

int64_t pen_of(int64_t individual, int64_t per_pen) {
  return individual/per_pen;
}
