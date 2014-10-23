#ifndef _PEN_HPP_
#define _PEN_HPP_ 1
#include <cstdint>
#include <cstddef>
#include "boost/graph/adjacency_list.hpp"

// This kind of adjacency list uses integers for the vertex id, which
// is special to those using vecS and vecS.
using PenContactGraph=boost::adjacency_list<boost::vecS,
    boost::vecS,boost::undirectedS>;

/*! Feedlots look like suburbs. There are long blocks with streets
 *  between. Each block has row_cnt pens, sitting back-to-back.
 *  Returns block_cnt*row_cnt*2 pens.
 */
PenContactGraph BlockStructure(int block_cnt, int row_cnt);
PenContactGraph DisconnectedPens(int pen_cnt);
bool AdjacentPens(int i, int j, const PenContactGraph& g);
int64_t pen_of(int64_t individual, int64_t per_pen);

#endif
