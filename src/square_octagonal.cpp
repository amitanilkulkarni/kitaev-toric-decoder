/* ---------------------------------------------------------------------

               DECODER FOR SQUARE OCTAGONAL COLOR CODE
	       ***************************************

   INPUT:  Error probability
   
   OUTPUT: Error rate for given error probability

   KNOWN BUGS:
   - Z error syndrome calculation is faulty ('2m' instead of '1').

--------------------------------------------------------------------- */


#include <iostream>
#include <iomanip>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <fstream>                   // for file i/o
#include <cstdlib>                   // for system command
#include <omp.h>		     // for multithreading
#include <chrono>
#include <stdexcept>
#include <sys/time.h>
#include <string>
#include "math.h"
#include <gsl/gsl_rng.h>             // GNU Scientific Lib's RNG
#include "rng-omp.h"


#define MIN_TRIALS 1
#define MIN_FAILS  1




// Defining the node class
class ColorNode {
public:
  bool z_error;
  bool x_error;
};


// Each codeblock is a tile of square octagonal lattice with 16 nodes
class CodeBlock {
public:
  int x, y;
  ColorNode node[16];
};


double RandomNum;
RNG rng(1);	// Syntax: rng(num_of_threads)





// -----------------------------------------------------------------------
//
//                             MAIN FUNCTION
//
// -----------------------------------------------------------------------

int main(int argc, char *argv[]) {


  // PARAMETERS ----------------------------------------------------------
  int qubits = 64; 	       // Must be of the form 16N^2 (default: N=2)
  int trials = MIN_TRIALS;     // Minimum number of trials
  int errors = MIN_FAILS;      // Minimum number of decoder failures
  int num_of_faces=qubits/16;
  int blockDim = sqrt(num_of_faces);
  int mPos = 4;	    // Position of '2m' when '1' is point 1. Must be even.
  bool verbose = 1;
  float errorProb = atof(argv[1]);
  // ---------------------------------------------------------------------

  
  // Checking if sufficient arguments are provided
  if(argc!=2) {
    std::cout << "Usage: ./square_octagonal_mapping <error probability>\n";
    return 1;
  }




  
  // ---------------------------------------------------------------------
  // Reading the LUTs for error lifting
  // ---------------------------------------------------------------------
  char filename1[32];
  // Reading the LUT-1 for lifting X errors (Xmap1.txt)
  sprintf(filename1, "LUT/d%d_pos%d_Xmap1.txt", 2*blockDim, mPos);
  std::ifstream lut1(filename1);
  std::vector<int> Xmap1rev[8*num_of_faces];
  
  int edgeNum, subGraph, vNum;
  char separator;
  
  // Storing the reverse map in a vector
  for(int i=0; i<8*num_of_faces; i++) {

    lut1 >> edgeNum;
    
    while(1) {

      // Identifying the action based on separator
      lut1 >> separator;
      if(separator=='-') {
	lut1 >> subGraph;
      }
      else if (separator == ',') {
	lut1 >> vNum; vNum--;
	Xmap1rev[edgeNum].push_back(vNum+16*subGraph);
      }
      else if(separator==';')
	break;
    }
    
  }

  lut1.close();


  // Reading the LUT-2 for lifting X errors (Xmap2.txt)
  sprintf(filename1, "LUT/d%d_pos%d_Xmap2.txt", 2*blockDim, mPos);
  std::ifstream lut2(filename1);
  std::vector<int> Xmap2rev[8*num_of_faces];
  // Storing the reverse map in a vector
  for(int i=0; i<8*num_of_faces; i++) {

    lut2 >> edgeNum;
    
    while(1) {

      // Identifying the action based on separator
      lut2 >> separator;
      if(separator=='-') {
	lut2 >> subGraph;
      }
      else if (separator == ',') {
	lut2 >> vNum; vNum--;
	Xmap2rev[edgeNum].push_back(vNum+16*subGraph);
      }
      else if(separator==';')
	break;
    }
    
  }

  lut2.close();
  




  
  CodeBlock graph[num_of_faces];
  int fail_count=0;
  
  // -----------------------------------------------------------------------------------------------
  // BEGIN TRIALS ----------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------------

  for (int trial_no=1; trial_no<=20000; trial_no++) {

    double error_rate;
    int Xerror_weight, Zerror_weight;
    int syn_check;

    Xerror_weight=0;
    Zerror_weight=0;
    
    if(verbose) std::cout<<"Begin trial: "<<trial_no<<"\n";
  
    // Tiling the 16-qubit blocks
    for (int i=0; i<num_of_faces; i++) {

      graph[i].x = i%blockDim;
      graph[i].y = i/blockDim;

      for (int j=0; j<16; j++) {

	// Introducing random X-error
	RandomNum = rng();
	if(RandomNum < errorProb) {
	  graph[i].node[j].x_error = 1; Xerror_weight++;
	  if(verbose) std::cout << "X at Subgraph " << i << " vertex " << j+1 << "\n";
	}
	else {
	  graph[i].node[j].x_error = 0;
	}

	// Introducing RandomNum Z-error (Will worry about it later)
	RandomNum = rng();
	if(RandomNum < errorProb) {
	  graph[i].node[j].z_error = 1; Zerror_weight++;
	  //if(verbose) std::cout << "Z at Subgraph " << i << " vertex " << j+1 << "\n";
	}
	else {
	  graph[i].node[j].z_error = 0;
	}
      
      }
    }

    if(Xerror_weight==0) {
      if(verbose) std::cout<<"Zero error weight graph; skipped.\n";
      continue;
    }
  

  

    // ---------------------------------------------------------------------
    // Calculating X syndromes
    // ---------------------------------------------------------------------
    bool XoctSyndrome[4*num_of_faces]; // For octahedral faces
    bool  XsqSyndrome[4*num_of_faces]; // For square faces

    // 'i' represents block number; remember 0-base for arrays!
    for (int i=0; i<num_of_faces; i++) {

      int j;

      // Syndrome of octahedral face 1 (0,0)
      XoctSyndrome[4*i+0] = graph[i].node[11].x_error ^ graph[i].node[12].x_error;

      j = (graph[i].x-1+blockDim)%blockDim + graph[i].y*blockDim;
      XoctSyndrome[4*i+0] ^= graph[j].node[13].x_error ^ graph[j].node[14].x_error;

      j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + graph[i].x;
      XoctSyndrome[4*i+0] ^= graph[j].node[9].x_error ^ graph[j].node[10].x_error;

      j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + (graph[i].x-1+blockDim)%blockDim;
      XoctSyndrome[4*i+0] ^= graph[j].node[8].x_error ^ graph[j].node[15].x_error;


      // Syndrome of octahedral face 2 (1,0)
      XoctSyndrome[4*i+1] = graph[i].node[11].x_error ^ graph[i].node[0].x_error
	^ graph[i].node[1].x_error ^ graph[i].node[14].x_error;

      j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + graph[i].x;
      XoctSyndrome[4*i+1] ^= graph[j].node[10].x_error ^ graph[j].node[5].x_error
	^ graph[j].node[4].x_error ^ graph[j].node[15].x_error;

      // Syndrome of octahedral face 3 (0,1)
      XoctSyndrome[4*i+2] = graph[i].node[12].x_error ^ graph[i].node[7].x_error
	^ graph[i].node[6].x_error ^ graph[i].node[9].x_error;

      j = ((graph[i].x-1+blockDim)%blockDim) + graph[i].y*blockDim;
      XoctSyndrome[4*i+2] ^= graph[j].node[13].x_error ^ graph[j].node[8].x_error
	^ graph[j].node[3].x_error ^ graph[j].node[2].x_error;

      // Syndrome of octahedral face 4 (1,1)
      XoctSyndrome[4*i+3] = graph[i].node[0].x_error ^ graph[i].node[1].x_error
	^ graph[i].node[2].x_error ^ graph[i].node[3].x_error^ graph[i].node[4].x_error
	^ graph[i].node[5].x_error^ graph[i].node[6].x_error^ graph[i].node[7].x_error;
    }

  
    for (int i=0; i<num_of_faces; i++){

      // Syndrome of square face 1 (0,0)
      XsqSyndrome[4*i+0] = graph[i].node[0].x_error ^ graph[i].node[7].x_error
	^ graph[i].node[11].x_error ^ graph[i].node[12].x_error;

      // Syndrome of square face 2 (0,1)
      XsqSyndrome[4*i+1] = graph[i].node[1].x_error ^ graph[i].node[2].x_error
	^ graph[i].node[13].x_error ^ graph[i].node[14].x_error;

      // Syndrome of square face 3 (1,0)
      XsqSyndrome[4*i+2] = graph[i].node[5].x_error ^ graph[i].node[6].x_error
	^ graph[i].node[9].x_error ^ graph[i].node[10].x_error;

      // Syndrome of square face 4 (1,1)
      XsqSyndrome[4*i+3] = graph[i].node[3].x_error ^ graph[i].node[4].x_error
	^ graph[i].node[8].x_error ^ graph[i].node[15].x_error;

    }


    // ---------------------------------------------------------------------
    // Calculating Z syndromes
    // ---------------------------------------------------------------------
    bool ZoctSyndrome[4*num_of_faces]; // For octahedral faces
    bool  ZsqSyndrome[4*num_of_faces]; // For square faces

    // 'i' represents block number; remember 0-base for arrays!
    for (int i=0; i<num_of_faces; i++) {

      int j;

      // Syndrome of octahedral face 1 (0,0)
      ZoctSyndrome[4*i+0] = graph[i].node[11].z_error ^ graph[i].node[12].z_error;

      j = (graph[i].x-1+blockDim)%blockDim + graph[i].y*blockDim;
      ZoctSyndrome[4*i+0] ^= graph[j].node[13].z_error ^ graph[j].node[14].z_error;

      j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + graph[i].x;
      ZoctSyndrome[4*i+0] ^= graph[j].node[9].z_error ^ graph[j].node[10].z_error;

      j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + (graph[i].x-1+blockDim)%blockDim;
      ZoctSyndrome[4*i+0] ^= graph[j].node[8].z_error ^ graph[j].node[15].z_error;


      // Syndrome of octahedral face 2 (0,1)
      ZoctSyndrome[4*i+1] = graph[i].node[11].z_error ^ graph[i].node[0].z_error
	^ graph[i].node[1].z_error ^ graph[i].node[14].z_error;

      j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + graph[i].x;
      ZoctSyndrome[4*i+1] ^= graph[j].node[10].z_error ^ graph[j].node[5].z_error
	^ graph[j].node[4].z_error ^ graph[j].node[15].z_error;

      // Syndrome of octahedral face 3 (1,0)
      ZoctSyndrome[4*i+2] = graph[i].node[12].z_error ^ graph[i].node[7].z_error
	^ graph[i].node[6].z_error ^ graph[i].node[9].z_error;

      j = ((graph[i].x-1+blockDim)%blockDim) + graph[i].y*blockDim;
      ZoctSyndrome[4*i+2] ^= graph[j].node[13].z_error ^ graph[j].node[8].z_error
	^ graph[j].node[3].z_error ^ graph[j].node[2].z_error;

      // Syndrome of octahedral face 4 (1,1)
      ZoctSyndrome[4*i+3] = graph[i].node[0].z_error ^ graph[i].node[1].z_error
	^ graph[i].node[2].z_error ^ graph[i].node[3].z_error^ graph[i].node[4].z_error
	^ graph[i].node[5].z_error^ graph[i].node[6].z_error^ graph[i].node[7].z_error;
    }

  
    for (int i=0; i<num_of_faces; i++){

      // Syndrome of square face 1 (0,0)
      ZsqSyndrome[4*i+0] = graph[i].node[0].z_error ^ graph[i].node[7].z_error
	^ graph[i].node[11].z_error ^ graph[i].node[12].z_error;

      // Syndrome of square face 2 (0,1)
      ZsqSyndrome[4*i+1] = graph[i].node[1].z_error ^ graph[i].node[2].z_error
	^ graph[i].node[13].z_error ^ graph[i].node[14].z_error;

      // Syndrome of square face 3 (1,0)
      ZsqSyndrome[4*i+2] = graph[i].node[5].z_error ^ graph[i].node[6].z_error
	^ graph[i].node[9].z_error ^ graph[i].node[10].z_error;

      // Syndrome of square face 4 (1,1)
      ZsqSyndrome[4*i+3] = graph[i].node[3].z_error ^ graph[i].node[4].z_error
	^ graph[i].node[8].z_error ^ graph[i].node[15].z_error;

    }



    // ---------------------------------------------------------------------
    // Taking care of the indirectly-mapped syndromes of square tiles
    // ---------------------------------------------------------------------
  
    switch( mPos ) {
    
    case 2:
      for (int i=0; i<num_of_faces; i++) {
	XsqSyndrome[4*i+1] ^= XoctSyndrome[4*i+3];
	XsqSyndrome[4*i+2] ^= XoctSyndrome[4*i+0];
      }

      for (int i=0; i<num_of_faces; i++) {
	ZsqSyndrome[4*i+1] ^= ZoctSyndrome[4*i+0];
	ZsqSyndrome[4*i+2] ^= ZoctSyndrome[4*i+3];
      }
      break;

    case 4:
      for (int i=0; i<num_of_faces; i++) {
	XsqSyndrome[4*i+0] ^= XoctSyndrome[4*i+0];
	XsqSyndrome[4*i+3] ^= XoctSyndrome[4*i+3];
      }

      for (int i=0; i<num_of_faces; i++) {
	ZsqSyndrome[4*i+0] ^= ZoctSyndrome[4*i+3];
	ZsqSyndrome[4*i+3] ^= ZoctSyndrome[4*i+0];
      }
      break;
    
    case 6:
      for (int i=0; i<num_of_faces; i++) {
	XsqSyndrome[4*i+1] ^= XoctSyndrome[4*i+0];
	XsqSyndrome[4*i+2] ^= XoctSyndrome[4*i+3];
      }

      for (int i=0; i<num_of_faces; i++) {
	ZsqSyndrome[4*i+1] ^= ZoctSyndrome[4*i+3];
	ZsqSyndrome[4*i+2] ^= ZoctSyndrome[4*i+0];
      }
      break;

    case 8:
      for (int i=0; i<num_of_faces; i++) {
	XsqSyndrome[4*i+0] ^= XoctSyndrome[4*i+3];
	XsqSyndrome[4*i+3] ^= XoctSyndrome[4*i+0];
      }

      for (int i=0; i<num_of_faces; i++) {
	ZsqSyndrome[4*i+0] ^= ZoctSyndrome[4*i+0];
	ZsqSyndrome[4*i+3] ^= ZoctSyndrome[4*i+3];
      }
      break;

    default:  std::cout<< "Error: mPos must be an even number.\n"; return 1;
    
    }


    // ---------------------------------------------------------------------
    // Printing out the mapping to a file
    // ---------------------------------------------------------------------

    // Map 1 for X errors
    std::ofstream fp0;
    sprintf(filename1, "/tmp/Xmap1.txt");
    fp0.open(filename1);
    fp0 << 2*blockDim << "\n";			// Dimension of sq lattice

  
    // Map 2 for X errors
    std::ofstream fp1;
    char filename2[32];
    sprintf(filename2, "/tmp/Xmap2.txt");
    fp1.open(filename2);
    fp1 << 2*blockDim << "\n";			// Dimension of sq lattice

  
    // Map 1 for Z errors
    std::ofstream fp2;
    char filename3[32];
    sprintf(filename3, "/tmp/Zmap1.txt");
    fp2.open(filename3);
    fp2 << 2*blockDim << "\n";			// Dimension of sq lattice
  
    // Map 2 for Z errors
    std::ofstream fp3;
    char filename4[32];
    sprintf(filename4, "/tmp/Zmap2.txt");
    fp3.open(filename4);
    fp3 << 2*blockDim << "\n";			// Dimension of sq lattice
  

    syn_check=0;
    // Printing X syndromes of octagonal faces to X map 1
    for (int i=0; i<num_of_faces; i++) {		// Vertex Syndromes
      for (int j=0; j<4; j++) {
	fp0 << XoctSyndrome[4*i+j] << " ";	// Note the order
	if(XoctSyndrome[4*i+j]) syn_check++;
      }
    }
    fp0 << "\n";  

    bool skip_map1=0;
    if(syn_check==0) {
      if(verbose) std::cout << "Skipped matching on map 1.\n";
      skip_map1=1;
    }

    syn_check=0;
    // Printing X syndromes of square faces to X map 2
    for (int i=0; i<num_of_faces; i++) {		// Vertex Syndromes
      for (int j=0; j<4; j++) {
	fp1 << XsqSyndrome[4*i+j] << " ";		// Note the order
	if(XsqSyndrome[4*i+j]) syn_check++;
      }
    }
    fp1 << "\n";

    bool skip_map2=0;
    if(syn_check==0) {
      if(verbose) std::cout << "Skipped matching on map 2.\n";
      skip_map2=1;
    }

    /*    
    // Printing Z syndromes of octagonal faces to Z map 1
    for (int i=0; i<num_of_faces; i++) {		// Vertex Syndromes
      for (int j=0; j<4; j++) {
	fp2 << ZoctSyndrome[4*i+j] << " ";	// Note the order
      }
    }
    fp2 << "\n";


    // Printing Z syndromes of octagonal faces to Z map 2
    for (int i=0; i<num_of_faces; i++) {		// Vertex Syndromes
      for (int j=0; j<4; j++) {
	fp3 << ZsqSyndrome[4*i+j] << " ";		// Note the order
      }
    }
    fp3 << "\n";
    */

    fp0.close();
    fp1.close();
    /*fp2.close();
      fp3.close();*/



  
    // ---------------------------------------------------------------------
    // Running perfect matching on Xmap1.txt
    // ---------------------------------------------------------------------  

    char command[150];
    sprintf(command, "./square_octagonal_decoder /tmp/Xmap1.txt 1");

    if(!skip_map1) { // Taking care of the case where Xmap1 has no syndromes
      system(command);
      // Reading the error path
      sprintf(filename1, "/tmp/error_path1.txt");
      std::ifstream errorPath1(filename1);

      while(errorPath1 >> edgeNum) {
	for(auto It=Xmap1rev[edgeNum].begin(); It!=Xmap1rev[edgeNum].end(); It++) {
	  graph[(*It)/16].node[(*It)%16].x_error = !graph[(*It)/16].node[(*It)%16].x_error;
	}
	if(verbose) std::cout<<edgeNum<<" ";
      }
      if(verbose) std::cout<<"\n";
      errorPath1.close();
    }

  
    // ---------------------------------------------------------------------
    // Running perfect matching on Xmap2.txt
    // ---------------------------------------------------------------------  

    //char command[150];
    sprintf(command, "./square_octagonal_decoder /tmp/Xmap2.txt 2");

    if(!skip_map2){ // Taking care of the case where Xmap2 has no syndromes
      system(command);
      // Reading the error path ----------------------------------------------
      sprintf(filename1, "/tmp/error_path2.txt");
      std::ifstream errorPath2(filename1);

      while(errorPath2 >> edgeNum) {
	for(auto It=Xmap2rev[edgeNum].begin(); It!=Xmap2rev[edgeNum].end(); It++) {
	  graph[(*It)/16].node[(*It)%16].x_error = !graph[(*It)/16].node[(*It)%16].x_error;
	}
	if(verbose) std::cout<<edgeNum<<" ";
      }
      if(verbose) std::cout<<"\n";
      errorPath2.close();
    }
    
    if(verbose) std::cout<<"Errors present after correction:\n";

    // Printing final path
    for (int i=0; i<num_of_faces*16; i++) {
      if( graph[i/16].node[i%16].x_error ) {
	if(verbose) std::cout << "X at subgraph " << i/16 << " vertex " << i%16+1 << "\n";
      }
    }

  
  

    // ---------------------------------------------------------------------
    // Calculating the syndromes after error correction
    // ---------------------------------------------------------------------

  

    // Calculating X syndromes ---------------------------------------------
  
    bool PostXoctSyndrome[4*num_of_faces]; // For octahedral faces
    bool  PostXsqSyndrome[4*num_of_faces]; // For square faces

    // 'i' represents block number; remember 0-base for arrays!
    for (int i=0; i<num_of_faces; i++) {

      int j;

      // Syndrome of octahedral face 1 (0,0)
      PostXoctSyndrome[4*i+0] = graph[i].node[11].x_error ^ graph[i].node[12].x_error;

      j = (graph[i].x-1+blockDim)%blockDim + graph[i].y*blockDim;
      PostXoctSyndrome[4*i+0] ^= graph[j].node[13].x_error ^ graph[j].node[14].x_error;

      j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + graph[i].x;
      PostXoctSyndrome[4*i+0] ^= graph[j].node[9].x_error ^ graph[j].node[10].x_error;

      j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + (graph[i].x-1+blockDim)%blockDim;
      PostXoctSyndrome[4*i+0] ^= graph[j].node[8].x_error ^ graph[j].node[15].x_error;


      // Syndrome of octahedral face 2 (1,0)
      PostXoctSyndrome[4*i+1] = graph[i].node[11].x_error ^ graph[i].node[0].x_error
	^ graph[i].node[1].x_error ^ graph[i].node[14].x_error;

      j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + graph[i].x;
      PostXoctSyndrome[4*i+1] ^= graph[j].node[10].x_error ^ graph[j].node[5].x_error
	^ graph[j].node[4].x_error ^ graph[j].node[15].x_error;

      // Syndrome of octahedral face 3 (0,1)
      PostXoctSyndrome[4*i+2] = graph[i].node[12].x_error ^ graph[i].node[7].x_error
	^ graph[i].node[6].x_error ^ graph[i].node[9].x_error;

      j = ((graph[i].x-1+blockDim)%blockDim) + graph[i].y*blockDim;
      PostXoctSyndrome[4*i+2] ^= graph[j].node[13].x_error ^ graph[j].node[8].x_error
	^ graph[j].node[3].x_error ^ graph[j].node[2].x_error;

      // Syndrome of octahedral face 4 (1,1)
      PostXoctSyndrome[4*i+3] = graph[i].node[0].x_error ^ graph[i].node[1].x_error
	^ graph[i].node[2].x_error ^ graph[i].node[3].x_error^ graph[i].node[4].x_error
	^ graph[i].node[5].x_error^ graph[i].node[6].x_error^ graph[i].node[7].x_error;
    }

  
    for (int i=0; i<num_of_faces; i++){

      // Syndrome of square face 1 (0,0)
      PostXsqSyndrome[4*i+0] = graph[i].node[0].x_error ^ graph[i].node[7].x_error
	^ graph[i].node[11].x_error ^ graph[i].node[12].x_error;

      // Syndrome of square face 2 (0,1)
      PostXsqSyndrome[4*i+1] = graph[i].node[1].x_error ^ graph[i].node[2].x_error
	^ graph[i].node[13].x_error ^ graph[i].node[14].x_error;

      // Syndrome of square face 3 (1,0)
      PostXsqSyndrome[4*i+2] = graph[i].node[5].x_error ^ graph[i].node[6].x_error
	^ graph[i].node[9].x_error ^ graph[i].node[10].x_error;

      // Syndrome of square face 4 (1,1)
      PostXsqSyndrome[4*i+3] = graph[i].node[3].x_error ^ graph[i].node[4].x_error
	^ graph[i].node[8].x_error ^ graph[i].node[15].x_error;

    }

    /*

    // Calculating Z syndromes ---------------------------------------------

    bool PostZoctSyndrome[4*num_of_faces]; // For octahedral faces
    bool  PostZsqSyndrome[4*num_of_faces]; // For square faces

    // 'i' represents block number; remember 0-base for arrays!
    for (int i=0; i<num_of_faces; i++) {

    int j;

    // Syndrome of octahedral face 1 (0,0)
    PostZoctSyndrome[4*i+0] = graph[i].node[11].z_error ^ graph[i].node[12].z_error;

    j = (graph[i].x-1+blockDim)%blockDim + graph[i].y*blockDim;
    PostZoctSyndrome[4*i+0] ^= graph[j].node[13].z_error ^ graph[j].node[14].z_error;

    j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + graph[i].x;
    PostZoctSyndrome[4*i+0] ^= graph[j].node[9].z_error ^ graph[j].node[10].z_error;

    j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + (graph[i].x-1+blockDim)%blockDim;
    PostZoctSyndrome[4*i+0] ^= graph[j].node[8].z_error ^ graph[j].node[15].z_error;


    // Syndrome of octahedral face 2 (0,1)
    PostZoctSyndrome[4*i+1] = graph[i].node[11].z_error ^ graph[i].node[0].z_error
    ^ graph[i].node[1].z_error ^ graph[i].node[14].z_error;

    j = ((graph[i].y-1+blockDim)%blockDim)*blockDim + graph[i].x;
    PostZoctSyndrome[4*i+1] ^= graph[j].node[10].z_error ^ graph[j].node[5].z_error
    ^ graph[j].node[4].z_error ^ graph[j].node[15].z_error;

    // Syndrome of octahedral face 3 (1,0)
    PostZoctSyndrome[4*i+2] = graph[i].node[12].z_error ^ graph[i].node[7].z_error
    ^ graph[i].node[6].z_error ^ graph[i].node[9].z_error;

    j = ((graph[i].x-1+blockDim)%blockDim) + graph[i].y*blockDim;
    PostZoctSyndrome[4*i+2] ^= graph[j].node[13].z_error ^ graph[j].node[8].z_error
    ^ graph[j].node[3].z_error ^ graph[j].node[2].z_error;

    // Syndrome of octahedral face 4 (1,1)
    PostZoctSyndrome[4*i+3] = graph[i].node[0].z_error ^ graph[i].node[1].z_error
    ^ graph[i].node[2].z_error ^ graph[i].node[3].z_error^ graph[i].node[4].z_error
    ^ graph[i].node[5].z_error^ graph[i].node[6].z_error^ graph[i].node[7].z_error;
    }

  
    for (int i=0; i<num_of_faces; i++){

    // Syndrome of square face 1 (0,0)
    PostZsqSyndrome[4*i+0] = graph[i].node[0].z_error ^ graph[i].node[7].z_error
    ^ graph[i].node[11].z_error ^ graph[i].node[12].z_error;

    // Syndrome of square face 2 (0,1)
    PostZsqSyndrome[4*i+1] = graph[i].node[1].z_error ^ graph[i].node[2].z_error
    ^ graph[i].node[13].z_error ^ graph[i].node[14].z_error;

    // Syndrome of square face 3 (1,0)
    PostZsqSyndrome[4*i+2] = graph[i].node[5].z_error ^ graph[i].node[6].z_error
    ^ graph[i].node[9].z_error ^ graph[i].node[10].z_error;

    // Syndrome of square face 4 (1,1)
    PostZsqSyndrome[4*i+3] = graph[i].node[3].z_error ^ graph[i].node[4].z_error
    ^ graph[i].node[8].z_error ^ graph[i].node[15].z_error;

    }
    */


    // ---------------------------------------------------------------------
    // Taking care of the indirectly-mapped syndromes of square tiles
    // ---------------------------------------------------------------------
  
    switch( mPos ) {
    
    case 2:
      for (int i=0; i<num_of_faces; i++) {
	PostXsqSyndrome[4*i+1] ^= PostXoctSyndrome[4*i+3];
	PostXsqSyndrome[4*i+2] ^= PostXoctSyndrome[4*i+0];
      }
      /*
	for (int i=0; i<num_of_faces; i++) {
	PostZsqSyndrome[4*i+1] ^= PostZoctSyndrome[4*i+0];
	PostZsqSyndrome[4*i+2] ^= PostZoctSyndrome[4*i+3];
	}*/
      break;
    
    case 4:
      for (int i=0; i<num_of_faces; i++) {
	PostXsqSyndrome[4*i+0] ^= PostXoctSyndrome[4*i+0];
	PostXsqSyndrome[4*i+3] ^= PostXoctSyndrome[4*i+3];
      }
      /*
	for (int i=0; i<num_of_faces; i++) {
	PostZsqSyndrome[4*i+0] ^= PostZoctSyndrome[4*i+3];
	PostZsqSyndrome[4*i+3] ^= PostZoctSyndrome[4*i+0];
	}*/
      break;
    
    case 6:
      for (int i=0; i<num_of_faces; i++) {
	PostXsqSyndrome[4*i+1] ^= PostXoctSyndrome[4*i+0];
	PostXsqSyndrome[4*i+2] ^= PostXoctSyndrome[4*i+3];
      }
      /*
	for (int i=0; i<num_of_faces; i++) {
	PostZsqSyndrome[4*i+1] ^= PostZoctSyndrome[4*i+3];
	PostZsqSyndrome[4*i+2] ^= PostZoctSyndrome[4*i+0];
	}*/
      break;
    
    case 8:
      for (int i=0; i<num_of_faces; i++) {
	PostXsqSyndrome[4*i+0] ^= PostXoctSyndrome[4*i+3];
	PostXsqSyndrome[4*i+3] ^= PostXoctSyndrome[4*i+0];
      }
      /*
	for (int i=0; i<num_of_faces; i++) {
	PostZsqSyndrome[4*i+0] ^= PostZoctSyndrome[4*i+0];
	PostZsqSyndrome[4*i+3] ^= PostZoctSyndrome[4*i+3];
	}*/
      break;
    
    default: if(verbose) std::cout<< "Error: mPos must be an even number.\n"; return 1;
    
    }

  
    // Recalculating Syndromes. Should be Zero. If not, matching failed.
    int synSum = 0;
    for (int i=0; i<num_of_faces; i++) {
      for (int j=0; j<4; j++) {
	synSum += PostXoctSyndrome[4*i+j];
      }
    }
    for (int i=0; i<num_of_faces; i++) {
      for (int j=0; j<4; j++) {
	synSum+= PostXsqSyndrome[4*i+j];
      }
    }

    if(synSum) {
      /*if(verbose)*/ std::cout<<"Syndromes nonzero after error correction. Matching failed!\n";
      return 1;
    }


    // ---------------------------------------------------------------------
    // Logical operators: Checking for non-trivial cycles
    // ---------------------------------------------------------------------

    bool cycleCheck=0;

    // First cycle (Vertical)
    for (int i=0; i<num_of_faces; i+=blockDim) {
      cycleCheck ^= graph[i].node[11].x_error^graph[i].node[7].x_error
                   ^graph[i].node[6].x_error^graph[i].node[10].x_error;
    }

    // Second cycle (Vertical)
    for (int i=0; i<num_of_faces; i+=blockDim) {
      cycleCheck ^= graph[i].node[11].x_error^graph[i].node[0].x_error
                   ^graph[i].node[6].x_error^graph[i].node[9].x_error;
    }

    // Third cycle (Horizontal)
    for (int i=0; i<blockDim; i++) {
      cycleCheck ^= graph[i].node[12].x_error^graph[i].node[0].x_error
                   ^graph[i].node[1].x_error^graph[i].node[13].x_error;
    }
    
    // Fourth cycle (Horizontal)
    for (int i=0; i<blockDim; i++) {
      cycleCheck ^= graph[i].node[12].x_error^graph[i].node[7].x_error
                   ^graph[i].node[1].x_error^graph[i].node[14].x_error;
    }
    
    // ---------------------------------------------------------------------
    // GATHERING TRIAL RESULTS & CHECKING TERMINATE CONDITION
    // ---------------------------------------------------------------------
    if (cycleCheck) { fail_count++; }

    error_rate = 1.0 * fail_count / trial_no;    

    if ( (trial_no > trials) && (fail_count > errors) ) {
      std::cout << "Final error rate for p = " << errorProb << " is " << error_rate << " .\n";
      break;
    }

    else if ( trial_no == 20000 ) {
      std::cout << "Final error rate for p = " << errorProb << " is " << error_rate << " .\n";
      break;
    }

    if( trial_no/1500*1500 == trial_no ) {
      std::cout << trial_no << "th trial: Errors = "<< fail_count << "\n";
    }

  }
  
  // -----------------------------------------------------------------------------------------------
  // END TRIALS ------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------------
}

