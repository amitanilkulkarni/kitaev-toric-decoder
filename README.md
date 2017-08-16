
----------------------------------------------------------------------
     CSS SURFACE CODES : DECODER SIMULATION FOR SQUARE LATTICE
----------------------------------------------------------------------

		AMIT A KULKARNI (EE11B069)
		DEPARTMENT OF ELECTRICAL ENGG
		IIT MADRAS, CHENNAI
		
----------------------------------------------------------------------

1) REQUIRED PACKAGES AND SYSTEM

   To compile the code, the following packages must be installed:

   gcc		Version should support C++11	(tested on v5.3)
   openmp	LLVM OpenMP Runtime Library	(tested on v3.7)
   boost	For boost graph libraries	(tested on v1.60)
   gsl		GNU Scientific Library		(tested on v2.1)

   BLOSSOM-V	Check README_FIRST.txt in ./blossom/

2) COMPILATION

    For toric decoder simulation     :  make decoder
    For square octagonal decoder     :  make sqoct
    To clean the compiled files      :  make clean

    Before building the code, it is recommended to clean the present
    object files using 'make clean'.

   
3) USAGE
 
   A) Toric code decoder simulation:

      ./surface_decoder [L] [P_START] [P_END] [P_STEP] [NUM_THREADS]
  
      Example:
          ./surface_decoder 8 0.05 0.18 0.02 12

      Where,

          L = Dimension of the square lattice (LxL)
          P_START = Starting point of error probability
          P_END = End point of error probability
          P_STEP = Step between values of error probability
          NUM_THREADS = Number of process threads. Optimal is $(nproc).

   B) Square octagonal decoder simulation:

      ./square_octagonal
