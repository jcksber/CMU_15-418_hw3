/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1: Jack Kasbeer (jkasbeer), Name 2: Charlie Cai (qcai)
 */

#ifndef __WIREOPT_H__
#define __WIREOPT_H__

#include <omp.h>

typedef struct
{
  /* Define the data structure for wire here */
} wire_t;

typedef int cost_t;

const char *get_option_string(const char *option_name, const char *default_value);
int get_option_int(const char *option_name, int default_value);
float get_option_float(const char *option_name, float default_value);

#endif
