/**
 * Parallel VLSI Wire Routing via OpenMP
 * Jack Kasbeer (jkasbeer), Qifang Cai (qcai)
 */

#ifndef __WIREOPT_H__
#define __WIREOPT_H__

#include <omp.h>
#define WIRE_MAX 20
/* value_t struct is used to calculate the local minimum path
 */
typedef struct{
  int aggr_max;
  int m;
} value_t;

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
  int wire;
  int list[WIRE_MAX];
	omp_lock_t lock;
	int val;
} cost_cell_t;

/* cost_t *
 * the struct defines the board;
 * contains both the previous record and the current board
 */
typedef struct
{
  int dimX;
  int dimY;
  int prevMax;
  int prevAggrTotal;
  int currentMax;
  int currentAggrTotal;
  cost_cell_t* board;
} cost_t;

/* Command line helper functions */
const char *get_option_string(const char *option_name, const char *default_value);
int get_option_int(const char *option_name, int default_value);
float get_option_float(const char *option_name, float default_value);

/* Our helper functions */
void horizontalCost(cost_cell_t *C, int row, int startX, int endX, int dimY, int wire_n);
void verticalCost(cost_cell_t *C, int xCoord, int startY, int endY, int dimY, int wire_n);
void new_rand_path(wire_t *wire);
void incrCell(cost_cell_t *C, int x, int y, int dimY, int wire_n);
void updateBoard(cost_t* board);
inline int readBoard(cost_t* board, int x, int y, int wire_n);
value_t readVertical(cost_t* board, int x, int s_y, int e_y, int wire_n);
value_t readHorizontal(cost_t* board, int y, int s_x, int e_x, int wire_n);
value_t calculatePath(cost_t* board, int s_x, int s_y, int e_x, int e_y,
          int numBends, int b1_x, int b1_y, int b2_x, int b2_y, int wire_n);
value_t combineValue(value_t v1, value_t v2);
//void cleanUpWire( cost_t board, path_t * path);
//inline void decrValue(cost_t board, int x, int y);
//void copyBoard(cost_cell_t *dest, cost_cell_t *src, int dimX, int dimY);
#endif
