/* rng-omp.h */

#ifndef RNG_H
#define RNG_H

#include <iostream>
#include <random>
#include <vector>

#include "omp.h"

#include <gsl/gsl_rng.h>             // GNU Scientific Lib's RNG

class RNG
{
	public:
		RNG(unsigned int);

		auto operator()()	-> double;

		RNG(const RNG&)							= delete;
		RNG& operator=(const RNG&)	= delete;
		RNG(RNG&&)									= delete;
		RNG& operator=(RNG&&)				= delete;

	private:
		std::vector<gsl_rng*> _engines;

};
#endif

/* End Of File */
