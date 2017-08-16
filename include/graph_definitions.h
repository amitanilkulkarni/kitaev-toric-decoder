#ifndef GRAPH_DEFINITIONS_H
#define GRAPH_DEFINITIONS_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;

typedef adjacency_list<vecS,vecS,undirectedS,surface_node,surface_edge> surface_code_t;

typedef graph_traits<surface_code_t>::edge_descriptor qubit;
    
typedef graph_traits<surface_code_t>::vertex_descriptor Node;

typedef graph_traits<surface_code_t>::edge_iterator qubit_i;

//typedef graph_traits<surface_code_t>::out_edge_iterator out_ei;

typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_index1_t, int>, 
property<edge_weight_t, int> > distance_graph_t;

typedef graph_traits<distance_graph_t>::edge_descriptor dist_edge;
    
typedef graph_traits<distance_graph_t>::vertex_descriptor dist_node;
    
typedef graph_traits<distance_graph_t>::edge_iterator dist_edge_i;
    
std::pair<dist_edge_i,dist_edge_i> dist_edge_p;

#endif
