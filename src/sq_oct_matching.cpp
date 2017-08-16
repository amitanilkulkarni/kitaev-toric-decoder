
/* ------------------------------------------------------------------


   INPUT : <Path to map.txt> <mapID (1 or 2)>

   OUTPUT: /tmp/error_path<mapID.txt
   

*/ //-----------------------------------------------------------------


#include <iostream>                  // for std::cout
#include <iomanip>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <fstream>                   // for file i/o
#include <cstdlib>                   // for system command
#include <omp.h>		     // for multithreading
#include <chrono>
#include <stdexcept>
#include <sys/time.h>
#include "math.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include "surf_node_and_edge.h"
#include "graph_definitions.h"

#include "rng-omp.h"

//--------------------------------------------------------------------
// DEFAULT PARAMETERS
//--------------------------------------------------------------------

#define OUTEDGES_PER_NODE 2	// Two for square, three for triangular
#define FACE_SIZE 4		// Four for square, three for triangular

using namespace boost;


//--------------------------------------------------------------------
// Function to compute distances on a periodic toric code lattice
//--------------------------------------------------------------------

int toric_lattice_dist(int node1_curr,int node2_curr, int toric_width) {
  
  int node1_x=node1_curr%toric_width;
  int node1_y=node1_curr/toric_width;
  int node2_x=node2_curr%toric_width;
  int node2_y=node2_curr/toric_width;

  //Have to consider periodicity for errors!
  //Distance along X-axis
  int x_dist_direct=abs(node1_x-node2_x);
  int x_distance=std::min(x_dist_direct,toric_width-x_dist_direct);

  //Distance along Y-axis
  int y_dist_direct=abs(node1_y-node2_y);
  int y_distance=std::min(y_dist_direct,toric_width-y_dist_direct);

  return (x_distance+y_distance);
}


//--------------------------------------------------------------------
// Function to create a toric code graph with errors
//--------------------------------------------------------------------

void toric_code_creator(int num_nodes, int code_dimension, surface_code_t& surf_code, int num_threads) {

  
  //Taking periodic boundaries into consideration. Modulo counter
  int periodic_index=0;

  //Two edges from every vertex. Forward and downwards.
  qubit temp_qubit_1,temp_qubit_2;

  //Adding all the edges to the graph
  int node_count; //Node counter

#pragma omp parallel for num_threads(num_threads)
  for(node_count=0; node_count<num_nodes; node_count++) {
    //Adding forward edge
    if((node_count+1)%code_dimension==0) {
      periodic_index=node_count+1-code_dimension;
      temp_qubit_1 = add_edge(node_count,periodic_index,surf_code).first;
      surf_code[temp_qubit_1].qubit_no=node_count;
    }
    else {
      temp_qubit_1 = add_edge(node_count,node_count+1,surf_code).first;
      surf_code[temp_qubit_1].qubit_no=node_count;
    }

    //Adding downward edge
    if((node_count+1+code_dimension)>num_nodes) {
      periodic_index=node_count+code_dimension-num_nodes;
      temp_qubit_2 = add_edge(node_count,periodic_index,surf_code).first;
      surf_code[temp_qubit_2].qubit_no=node_count+num_nodes;
    }
    else {
      temp_qubit_2 = add_edge(node_count,node_count+code_dimension,surf_code).first;
      surf_code[temp_qubit_2].qubit_no=node_count+num_nodes;
    }
  }
}



//--------------------------------------------------------------------
// Function to create the distance graph
//--------------------------------------------------------------------

void distance_graph_creator(int code_dim_dist, std::list<int>& syndrome_indices_list, int num_syndromes_dist, distance_graph_t& dist_graph, property_map<distance_graph_t, vertex_index1_t>::type old_index, int ID) {

  
  int dist_graph_count1=0;

  // Mapping old indices to new ones
  for (std::list<int>::iterator it = syndrome_indices_list.begin(); it != syndrome_indices_list.end(); it++) {
    put(old_index,dist_graph_count1,*it);
    dist_graph_count1=dist_graph_count1+1;
  }

  int dist_graph_count2=0; //to access other vertices to add edges

  //Adding edges with appropriate weights (lattice distances)
  for (dist_graph_count1=0;dist_graph_count1<num_syndromes_dist;dist_graph_count1++) {

    dist_graph_count2=dist_graph_count1+1;//Starting from the next dist_vertex
    int curr_node1_index=old_index[dist_graph_count1];

    //file_p<<curr_node1_index<<"\n";

    while(dist_graph_count2<num_syndromes_dist) {
      int curr_node2_index=old_index[dist_graph_count2];

      //Calling function to compute lattice distance (check the function)
      int curr_weight=toric_lattice_dist(curr_node1_index,curr_node2_index,code_dim_dist);
      add_edge(dist_graph_count1,dist_graph_count2,curr_weight,dist_graph);
      dist_graph_count2=dist_graph_count2+1;
    }
  }
}





//--------------------------------------------------------------------
// MAIN FUNCTION
//--------------------------------------------------------------------
int main(int argc, char *argv[]) {

  int code_dim;
  int num_threads=1;

  
  // Reading the map
  std::ifstream fp;
  fp.open(argv[1]);
  fp >> code_dim;
  int mapID = atoi(argv[2]);
  
  //std::cout << code_dim << "\n";
    
  
  // Fixed value variables --------------------------------------------------
  int no_of_qubits= OUTEDGES_PER_NODE*code_dim*code_dim;
  int no_of_edges=no_of_qubits;
  int no_of_nodes=no_of_qubits/OUTEDGES_PER_NODE;
  int face_size=FACE_SIZE;
  // -------------------------------------------------------------------------

  //std::cout << std::setprecision(4);

  // Creating and populating surface code graphs
  surface_code_t surf_code(no_of_nodes);
    
  //Creating toric code on its own
  toric_code_creator(no_of_nodes,code_dim,surf_code, num_threads);  
    
  graph_traits<surface_code_t>::out_edge_iterator out_i, out_end;


  
  qubit temp_qubit_error;
  int error_weight=0;
  int ID=omp_get_thread_num();

     
  //For size of distance graph
  int no_of_syndromes=0;
  //To hold the list of indices from the original graph having a non-trivial syndrome
  std::list<int> syndrome_indices_list;
  //Edge iterator for each vertex
  //typename graph_traits<surface_code_t>::out_edge_iterator out_i, out_end;


  int x=0, y=0, node_count=0;
  bool temp_syn;
  while ( y < code_dim ) {

    fp >> temp_syn;
    node_count = x + y*code_dim;
    if(temp_syn) {
      surf_code[node_count].syndrome_surf=-1;
      syndrome_indices_list.push_back (node_count);
      no_of_syndromes++;
    }

    x++;
    
    fp >> temp_syn;
    node_count = x + y*code_dim;
    if(temp_syn) {
      surf_code[node_count].syndrome_surf=-1;
      syndrome_indices_list.push_back (node_count);
      no_of_syndromes++;
    }

    y++; x--;
    
    fp >> temp_syn;
    node_count = x + y*code_dim;
    if(temp_syn) {
      surf_code[node_count].syndrome_surf=-1;
      syndrome_indices_list.push_back (node_count);
      no_of_syndromes++;
    }

    x++;

    fp >> temp_syn;
    node_count = x + y*code_dim;
    if(temp_syn) {
      surf_code[node_count].syndrome_surf=-1;
      syndrome_indices_list.push_back (node_count);
      no_of_syndromes++;
    }

    x++; y--;

    if(x == code_dim) {
      x=0;
      y++; y++;
    }

  }

  
  if(no_of_syndromes==0) {
    std::cout<<"No syndromes\n";
  }

  //--------------------------------------------------------------------
  // Creating a new complete graph from the vertices with syndromes
  //--------------------------------------------------------------------
 
  distance_graph_t dist_graph(no_of_syndromes);

  property_map<distance_graph_t, vertex_index1_t>::type old_index = get(vertex_index1, dist_graph);
  property_map<distance_graph_t, edge_weight_t>::type weight = get(edge_weight, dist_graph);
  graph_traits<distance_graph_t>::out_edge_iterator out_dist_i, out_dist_end;

  distance_graph_creator(code_dim, syndrome_indices_list, no_of_syndromes, dist_graph, old_index, ID);

  //--------------------------------------------------------------------
  //Printing out to file for Blossom 5
  //--------------------------------------------------------------------

  char str[32];
  sprintf(str, "/tmp/trial_graph_%d.txt",ID);
  std::ofstream myfile;

  #pragma omp critical
  {
  myfile.open(str);

  myfile<<no_of_syndromes<<" "<<num_edges(dist_graph)<<"\n";
  dist_edge output_edge;
  for(dist_edge_p=edges(dist_graph); dist_edge_p.first != dist_edge_p.second; ++dist_edge_p.first) {
    output_edge = *dist_edge_p.first;
    dist_node src_out = source(output_edge, dist_graph), targ_out = target(output_edge, dist_graph);
    myfile << src_out<< " " << targ_out << " "<<weight[output_edge]<<"\n";
  }

  myfile.close();
  }
  //blossom5_input_creator(no_of_syndromes, dist_graph, weight, ID);

  //-------------------------------------------------------------------
  // Invoking blossom5 from within the program
  //-------------------------------------------------------------------
  char command[150];      
  sprintf (command, "./blossom/blossom5 -e %s%d.txt -w %s%d.txt > /dev/null", "/tmp/trial_graph_", ID, "/tmp/output_", ID);
  //std::cout <<ID<<" "<< command<<"\n";
  system(command);


  //--------------------------------------------------------------------
  // Reading in matching from file and storing it
  //--------------------------------------------------------------------

  // char str[32];
  sprintf(str, "/tmp/output_%d.txt",ID);
      
  std::ifstream infile(str);
  int in_1, in_2;
  int matching_weight, syndromes_computed;
  int no_of_lines=0;
  std::vector<int> matched_indices;

  while(infile>>in_1>>in_2) {
    if(no_of_lines==0) {
      matching_weight=in_1;
      syndromes_computed=2*in_2;
      try{
	if(syndromes_computed!=no_of_syndromes) {
	  throw std::invalid_argument( "Incorrect matching!" );
	}
      }
      catch (const std::invalid_argument& e) {
	std::cout<<"Incorrect matching!"<<std::endl;
	//return -1;
	exit(-1);
      }
    }

    else if(no_of_lines!=0) {
      matched_indices.push_back(old_index[in_1]);
      matched_indices.push_back(old_index[in_2]);
    }
    no_of_lines=no_of_lines+1;
  }

  infile.close();

  //--------------------------------------------------------------------
  // Create geodesics & writing error to a file
  //--------------------------------------------------------------------
  int indices_count;
  int g1, g2;
  int g1_x,g2_x,g1_y,g2_y;
  qubit temp_qubit_geo;
  Node targ_geo;
  int geo_count, geo_count_alt, geo_count_x, geo_count_y;

  sprintf(str, "/tmp/error_path%d.txt",mapID);
  std::ofstream errfile(str);

  
  for(indices_count=0;indices_count<no_of_syndromes;indices_count+=2) {

    geo_count=0;
    geo_count_x=0;
    geo_count_y=0;

    // Pair of matched indices
    g1=matched_indices[indices_count];
    g2=matched_indices[indices_count+1];

    // X and Y coordinates of matched indices
    g1_x=g1%code_dim;
    g1_y=g1/code_dim;
    g2_x=g2%code_dim;
    g2_y=g2/code_dim;

    // Absolute X and Y distances
    int x_min=std::min(g1_x,g2_x);
    int x_max=std::max(g1_x,g2_x);
    int x_dist=x_max-x_min;
    int y_min=std::min(g1_y,g2_y);
    int y_max=std::max(g1_y,g2_y);
    int y_dist=y_max-y_min;


    // ----------------------------------------------------------------------
    // Adding geodesic along x direction, taking care of crossing geodesics
    // ----------------------------------------------------------------------

    /*
      MOVING PURELY IN X DIRECTION

      First, we move only along the x-direction, considering only
      geo_count_x as the criterion.
      Three cases are possible.
    */


    // 1- When x distance is zero
    //   a- If the shortest y distance is direct (unwrapped),
    //      note the point with the smaller Y and move on.
    //   b- If the smallest distance is wraped,
    //      note the point with the larger Y and move on.
    if(x_dist==0) {
      if(y_dist<=code_dim-y_dist) {
	if(y_min==g1_y) {geo_count=g1;}
	else if(y_min==g2_y) {geo_count=g2;}
      }
      else if(y_dist>code_dim-y_dist) {
	if(y_max==g1_y) {geo_count=g1;}
	else if(y_max==g2_y) {geo_count=g2;}
      }
    }


    // 2- When the shortest X-distance is unwrapped
    //
    //    Start from the point with the smaller X,
    //    Move in the direction of the other,
    //    Toggle Z-error along the path.
    //    Ultimately, x-distance becomes zero.
    //
    else if(x_dist<=code_dim-x_dist) {

      geo_count_x=x_min;
      if(x_min==g1_x) {geo_count=g1;}
      else if(x_min==g2_x) {geo_count=g2;}

      while(geo_count_x!=x_max) {
	if((geo_count+1)%code_dim==0) {
	  for (tie(out_i, out_end) = out_edges(geo_count, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);
	    if(int (targ_geo)==geo_count+1-code_dim) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }
	  }
	  geo_count=geo_count+1-code_dim;
	  geo_count_x=(geo_count_x+1)%code_dim;
	}

	else if((geo_count+1)%code_dim!=0) {
	  for (tie(out_i, out_end) = out_edges(geo_count, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);
	    if(int (targ_geo)==geo_count+1) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }
	  }
	  geo_count=geo_count+1;
	  geo_count_x=(geo_count_x+1)%code_dim;
	  //std::cout<<geo_count_x<<","<<x_max<<std::endl;
	}


      }
    }


    // 3- When the shortest X-distance is wrapped
    //
    //   Similar procedure as case 2,
    //   Just start with the point with max-X,
    //   And wrap around to reach the point with min-X.
    //
    else if(x_dist>code_dim-x_dist) {

      //std::cout<<indices_count<<",x2"<<std::endl;
      geo_count_x=x_max;

      if(x_max==g1_x) {geo_count=g1;}
      else if(x_max==g2_x) {geo_count=g2;}

      while(geo_count_x!=x_min) {

	if((geo_count+1)%code_dim==0) {
	  for (tie(out_i, out_end) = out_edges(geo_count, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);
	    if(int (targ_geo)==geo_count+1-code_dim) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }
	  }
	  geo_count=geo_count+1-code_dim;
	  geo_count_x=(geo_count_x+1)%code_dim;
	}

	else if((geo_count+1)%code_dim!=0) {
	  for (tie(out_i, out_end) = out_edges(geo_count, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);
	    if(int (targ_geo)==geo_count+1) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }
	  }

	  geo_count=geo_count+1;
	  geo_count_x=(geo_count_x+1)%code_dim;
	}
      }
    }

    /*
      MOVING PURELY IN Y-DIRECTION

      Having moved in the x-direction, we are now at some vertex
      in the lattice.

      Depending on whether we are at the Y-max point or Y-min
      point, two cases arise.

      Depending on whether the shortest distance is wrapped or
      direct, two cases arise again.

      If the shortest path to the destination is an increasing
      path from this vertex, we move from it to the destination
      and use the periodicity.

      If the shortest path is an increasing path from the
      destination vertex to the current one, we move there and
      switch the roles of current vertex and destination vertex.
    */

    geo_count_y=geo_count/code_dim;

    //
    // 1- If at y minimum and shortest is from min to max
    //
    if(geo_count_y==y_min && y_dist<=code_dim-y_dist) {
      //std::cout<<indices_count<<",y1"<<std::endl;
      while(geo_count_y!=y_max) {
	if(geo_count/code_dim==(code_dim-1)) {
	  for (tie(out_i, out_end) = out_edges(geo_count, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);
	    if(int (targ_geo)==geo_count-code_dim*(code_dim-1)) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }
	  }
	  geo_count=geo_count-code_dim*(code_dim-1);
	  geo_count_y=(geo_count_y+1)%code_dim;
	}
	else if(geo_count/code_dim!=(code_dim-1)) {
	  for (tie(out_i, out_end) = out_edges(geo_count, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);
	    if(int (targ_geo)==geo_count+code_dim) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }
	  }
	  geo_count=geo_count+code_dim;
	  geo_count_y=(geo_count_y+1)%code_dim;
	}
      }
    }

    //
    // 2- If at y maximum and shortest is from min to max
    //
    else if(geo_count_y==y_max && y_dist<=code_dim-y_dist) {
      // std::cout<<indices_count<<",y2"<<std::endl;
      geo_count_y=y_min;
      if(y_min==g1_y) {geo_count_alt=g1;}
      else if(y_min==g2_y) {geo_count_alt=g2;}

      while(geo_count_y!=y_max) {

	if(geo_count_alt/code_dim==(code_dim-1)) {
	  for (tie(out_i, out_end) = out_edges(geo_count_alt, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);
	    if(int (targ_geo)==geo_count_alt-code_dim*(code_dim-1)) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }
	  }
	  geo_count_alt=geo_count_alt-code_dim*(code_dim-1);
	  geo_count_y=(geo_count_y+1)%code_dim;
	}

	else if(geo_count_alt/code_dim!=(code_dim-1)) {

	  for (tie(out_i, out_end) = out_edges(geo_count_alt, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);

	    if(int (targ_geo)==geo_count_alt+code_dim) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }

	  }

	  geo_count_alt=geo_count_alt+code_dim;
	  geo_count_y=(geo_count_y+1)%code_dim;
	}
      }
    }

    //
    // 3- If at y maximum and shortest is from max to min
    //
    else if(geo_count_y==y_max && y_dist>code_dim-y_dist) {

      // std::cout<<indices_count<<",y3"<<std::endl;
      while(geo_count_y!=y_min) {

	if(geo_count/code_dim==(code_dim-1)) {

	  for (tie(out_i, out_end) = out_edges(geo_count, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);

	    if(int (targ_geo)==geo_count-code_dim*(code_dim-1)) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }

	  }

	  geo_count=geo_count-code_dim*(code_dim-1);
	  geo_count_y=(geo_count_y+1)%code_dim;
	}

	else if(geo_count/code_dim!=(code_dim-1)) {

	  for (tie(out_i, out_end) = out_edges(geo_count, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);

	    if(int (targ_geo==geo_count+code_dim)) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }
	  }
	  geo_count=geo_count+code_dim;
	  geo_count_y=(geo_count_y+1)%code_dim;
	}
      }
    }

    //
    // 4- If at y minimum and shortest is from max to min
    //
    else if(geo_count_y==y_min && y_dist>code_dim-y_dist) {
      // std::cout<<indices_count<<",y4"<<std::endl;
      geo_count_y=y_max;
      if(y_max==g1_y) {geo_count_alt=g1;}
      else if(y_max==g2_y) {geo_count_alt=g2;}
      while(geo_count_y!=y_min) {
	if(geo_count_alt/code_dim==(code_dim-1)) {
	  for (tie(out_i, out_end) = out_edges(geo_count_alt, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);
	    if(int (targ_geo)==geo_count_alt-code_dim*(code_dim-1)) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }
	  }
	  geo_count_alt=geo_count_alt-code_dim*(code_dim-1);
	  geo_count_y=(geo_count_y+1)%code_dim;
	}
	else if(geo_count_alt/code_dim!=(code_dim-1)) {
	  for (tie(out_i, out_end) = out_edges(geo_count_alt, surf_code); out_i != out_end; ++out_i) {
	    temp_qubit_geo = *out_i;
	    targ_geo = target(temp_qubit_geo, surf_code);
	    if(int (targ_geo)==geo_count_alt+code_dim) {
	      int temp_geo_error=0;
	      errfile<<surf_code[temp_qubit_geo].qubit_no<<" ";
	    }
	  }
	  geo_count_alt=geo_count_alt+code_dim;
	  geo_count_y=(geo_count_y+1)%code_dim;
	}
      }
    }

  }
  errfile<<"\n";
  errfile.close();
  fp.close();
}
