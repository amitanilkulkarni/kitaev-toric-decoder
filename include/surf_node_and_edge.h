#ifndef NODE_AND_EDGE_H
#define NODE_AND_EDGE_H

#include <string>

struct surface_edge{
    int z_error;
    int qubit_no;
    std::vector<std::string> assoc_stabs_surf;
};

struct surface_node{
    int syndrome_surf;
    std::string surf_vert_name;
};


#endif

