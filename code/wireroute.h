/**
 * Parallel VLSI Wire Routing via OpenMP
 * Jack Kasbeer (jkasbeer), Qifang Cai (qcai)
 */

#ifndef __WIREOPT_H__
#define __WIREOPT_H__

#include <omp.h>

typedef struct
{
  /* Define the data structure for wire here */
	int numBends;  // 0, 1, or 2
	int bends[4];  // bend 1, bend 2 (or no empty)
	int bounds[4]; // start point, end point ([x y x y])
} wire_t;

typedef int cost_t;

const char *get_option_string(const char *option_name, const char *default_value);
int get_option_int(const char *option_name, int default_value);
float get_option_float(const char *option_name, float default_value);

#endif
