/**
 * Parallel VLSI Wire Routing via OpenMP
 * Jack Kasbeer (jkasbeer), Qifang Cai (qcai)
 */

#ifndef __WIREOPT_H__
#define __WIREOPT_H__

#include <omp.h>

/* Path struct - describe a single wire
 * Defined by start & end points and
 */
typedef struct
{
	int numBends;   // 0, 1, or 2
	int bends[4];   // bend 1, bend 2 (or no empty)
	int bounds[4];  // start point, end point ([x y x y]) CONSTANT VALUES
} path_t;

/* wire_t *
 * Wire struct - define a single wire as two path's
 */
typedef struct
{
  // Current path 'definition'
	path_t *currentPath;
	// Previous path 'definition'
	path_t *prevPath;
} wire_t;

/* cost_cell_t *
 * Just an integer, but with a lock for cost array writes per cell
 */
typedef struct
{
	omp_lock_t lock;
	int val;
} cost_cell_t;

/* cost_t *
 * the struct defines the board;
 * contains both the previous record and the current board
 */
typedef struct
{
  int prevMax;
  int prevAggrTotal;
  cost_cell_t* board;
} cost_t;

const char *get_option_string(const char *option_name, const char *default_value);
int get_option_int(const char *option_name, int default_value);
float get_option_float(const char *option_name, float default_value);
/* Our helper functions */


// overwrite the prev path with current path and randomly generate a new path
void new_rand_path(wire_t *wire);

#endif
