/**
 * Parallel VLSI Wire Routing via OpenMP
 * Jack Kasbeer (jkasbeer), Qifang Cai (qcai)
 */

#ifndef __WIREOPT_H__
#define __WIREOPT_H__

#include <omp.h>

/* Path struct
 * Defined by start & end points and 
 */
typedef struct 
{
	int numBends;   // 0, 1, or 2
	int bends[4];   // bend 1, bend 2 (or no empty)
	int bounds[4];  // start point, end point ([x y x y]) CONSTANT VALUES
} path_t;

/* Wire struct
 * Defined as two paths
 */
typedef struct
{
  // Current path 'definition'
	path_t *currentPath;
	// Previous path 'definition'
	path_t *prevPath;
} wire_t;

/* Cost_t definition 
 * Just an integer, but with a lock for cost array writes per cell
 */
typedef struct 
{
	int lock;
	int val;
} cost_t;

const char *get_option_string(const char *option_name, const char *default_value);
int get_option_int(const char *option_name, int default_value);
float get_option_float(const char *option_name, float default_value);
/* Our helper functions */
wire_t initWire();

#endif
