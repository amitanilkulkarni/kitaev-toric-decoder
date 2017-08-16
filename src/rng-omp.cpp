/* rng-omp.cpp */

#include "rng-omp.h"
#include <time.h>
#include <gsl/gsl_rng.h>             // GNU Scientific Lib's RNG

/* constructor */
RNG::RNG(unsigned int threads): _engines()
{
  for(unsigned int i = 0; i < threads; ++i) {
    const gsl_rng_type * T;	// Refer: http://is.gd/c4H81S
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, time(NULL)+i+1);  
    _engines.push_back(r);
  }
}

/* get uniform distributed random number in (0,1) */
auto RNG::operator()() -> double
{
	return gsl_rng_uniform( _engines[ omp_get_thread_num() ] );
}

/* End of File */
