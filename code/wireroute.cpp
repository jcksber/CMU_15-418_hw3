/**
 * Parallel VLSI Wire Routing via OpenMP
 * Jack Kasbeer (jkasbeer), Qifang Cai (qcai)
 */

#include "wireroute.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "mic.h"

static int _argc;
static const char **_argv;

/////////////////////////////////////
// COMMAND LINE FUNCTIONS
////////////////////////////////////

const char *get_option_string(const char *option_name,
			      const char *default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return _argv[i + 1];
  return default_value;
}

int get_option_int(const char *option_name, int default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return atoi(_argv[i + 1]);
  return default_value;
}

float get_option_float(const char *option_name, float default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return (float)atof(_argv[i + 1]);
  return default_value;
}

static void show_help(const char *program_path)
{
    printf("Usage: %s OPTIONS\n", program_path);
    printf("\n");
    printf("OPTIONS:\n");
    printf("\t-f <input_filename> (required)\n");
    printf("\t-n <num_of_threads> (required)\n");
    printf("\t-p <SA_prob>\n");
    printf("\t-i <SA_iters>\n");
}
/////////////////////////////////////
/////////////////////////////////////


/////////////////////////////////////
// HELPER FUNCTIONS
/////////////////////////////////////

/* init_wires *
 * Clean up main routine, init wires here 
 */
void init_wires(FILE *input, wire_t *batch, int numWires)
{
  /* Read the grid dimension and wire information from file */
  int count = 0;
  while(count < numWires){
    int s_x, s_y, e_x, e_y;
    fscanf(input, "%d %d %d %d\n", &s_x, &s_y, &e_x, &e_y);
    batch[count].prevPath = (path_t*)calloc(1, sizeof(path_t));
    batch[count].currentPath = (path_t*)calloc(1, sizeof(path_t));
    batch[count].currentPath->numBends = 0;
    batch[count].currentPath->bounds[0] = s_x;
    batch[count].currentPath->bounds[1] = s_y;
    batch[count].currentPath->bounds[2] = e_x;
    batch[count].currentPath->bounds[3] = e_y;
    count++;
  }
}

/* init_cost_array *
 * Clean up main routine, init cost array here
 */
void init_cost_array(cost_t *arr, int numWires, int cols, int rows)
{
  costs->dimX = cols;
  costs->dimY = rows;
  costs->prevMax = num_of_wires;
  // no init for the prev_total of cost_t
  costs->board = (cost_cell_t *)calloc(cols * rows, sizeof(cost_cell_t));

  /* Initialize cost array */
  for( int y = 0; y < rows; y++){
    for( int x = 0; x < cols; x++){
      costs->board[y*rows + x].val = 0;
      omp_init_lock(costs->board[y*rows + x].lock);
    }
  }
}

/* new_rand_path *
 * Generate a random path in the space of delta_x + delta_y
 * 50% change pick x traversal  50% chance pick y traversal
 * random generate bends
 */
void new_rand_path(wire_t *wire)
{
  //overwrite previous path
  srand(time(NULL));
  int bend = 0;
  std::memcpy(wire->prevPath, wire->currentPath, sizeof(wire_t));
  int s_x, s_y, e_x, e_y, dy, yp, dx, xp;
  s_x = wire->currentPath->bounds[0];
  s_y = wire->currentPath->bounds[1];
  e_x = wire->currentPath->bounds[2];
  e_y = wire->currentPath->bounds[3];
  if (s_x == e_x || s_y == e_y){
    wire->currentPath->numBends = bend;
    return;
  }
  // not in the same line, need at least one bend
  bend +=1;
  dy = abs(e_y - s_y);
  dx = abs(e_x - s_x);
  // calculate random point on y axis
  if ((rand() % 10)> 5){
    // y first
    int ran_y = rand() % dy;
    if (ran_y == 0) ran_y = 1; // need to make progress
    if(s_y > e_y)  yp = s_y - ran_y;
    else yp = s_y + ran_y;
    // determind bends
    if(e_y != yp) bend +=1;
    // overwrite
    wire->currentPath->bends[0] = s_x;
    wire->currentPath->bends[1] = yp;
    wire->currentPath->bends[2] = e_x;
    wire->currentPath->bends[3] = yp;
    wire->currentPath->numBends = bend;
  }
  else{
    // x first traversal
    int ran_x = rand() % dx;
    if (ran_x == 0) ran_x = 1;
    if(s_x > e_x)  xp = s_x - ran_x;
    else xp = s_x + ran_x;
    // determind bends
    if(e_x != xp) bend +=1;
    // overwrite
    wire->currentPath->bends[0] = xp;
    wire->currentPath->bends[1] = s_y;
    wire->currentPath->bends[2] = xp;
    wire->currentPath->bends[3] = e_y;
    wire->currentPath->numBends = bend;
  }
}
/* incr_cell *
 * Increment value of cell (x,y)
 */
void incr_cell(cost_t *C, int x, int y)
{
  int dimY = C->dimY;
  /* Calculate the idx in board */
  cost_cell_t *c = C->board[y*dimY + x];
  /*### UPDATING CELL: CRITICAL REGION ###*/
  /* Lock -> increment -> unlock */
  omp_set_lock(&c->lock);
    c->val +=1;
  omp_unset_lock(&c->lock);
  /*######################################*/
}

/* path_cost *
 * Update cost array for given path determined by z (constant),
 * and the x OR y (one or the other) coordinates start & end
 */
void path_cost(cost_t *C, int z, int start, int end)
{
  int dir = start > end ? -1 : 1; //determine path direction
  int slip = start;

  /* Update cost array for given wire */
  while (slip != end)
  { 
    incr_cell(C, slip, z);
    slip += dir; 
  } 
  incr_cell(C, end, z); // update end point
}

/*DEAD CODE
  
  void horizontal_cost(cost_t *C, int y, int startX, int endX)
  {
    int dir = startX > endX ? -1 : 1; //determine path direction
    int x = startX;

     Update cost array for given wire 
    while (x != endX)
    { 
      incr_cell(C, x, y);
      x += dir; // add/subtract a column
    } 
    incr_cell(C, endX, y); // update end point
  }

  
  void vertical_cost(cost_t *C, int x, int startY, int endY)
  {
    cost_cell_t *c;
    // Determine path direction
    int dir = startY > endY ? -1 : 1;
    int y = startY;

    /* Update cost array for given wire 
    while (y != endY)
    { 
      incr_cell(C, x, y);
      y += dir;
    } 
    incr_cell(C, x, endY);
  }
*/

//////////////////////////////////////
//////////////////////////////////////


///////////////////////////////////////////////////////////
// MAIN ROUTINE
///////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
  using namespace std::chrono;
  typedef std::chrono::high_resolution_clock Clock;
  typedef std::chrono::duration<double> dsec;

  auto init_start = Clock::now();
  double init_time = 0;

  _argc = argc - 1;
  _argv = argv + 1;

  const char *input_filename = get_option_string("-f", NULL);
  int num_of_threads = get_option_int("-n", 1);
  double SA_prob = get_option_float("-p", 0.1f);
  int SA_iters = get_option_int("-i", 5);

  int error = 0;

  if (input_filename == NULL) {
    printf("Error: You need to specify -f.\n");
    error = 1;
  }

  if (error) {
    show_help(argv[0]);
    return 1;
  }

  printf("Number of threads: %d\n", num_of_threads);
  printf("Probability parameter for simulated annealing: %lf.\n", SA_prob);
  printf("Number of simulated anneling iterations: %d\n", SA_iters);
  printf("Input file: %s\n", input_filename);

  FILE *input = fopen(input_filename, "r");

  if (!input) {
    printf("Unable to open file: %s.\n", input_filename);
    return 1;
  }

  /* Parse for dimensions & num wires */
  int dim_x, dim_y;
  int num_of_wires;
  fscanf(input, "%d %d\n", &dim_x, &dim_y);
  fscanf(input, "%d\n", &num_of_wires);

  /* ALLOCATE for array of wires */
  wire_t *wires = (wire_t *)calloc(num_of_wires, sizeof(wire_t));
  /* INITIALIZE using input_file */
  init_wires(wires, num_of_wires, input);

  /* ALLOCATE for cost array struct */
  cost_t *costs = (cost_t *)calloc(1, sizeof(cost_t));
  /* INITIALIZE using dimensions of grid & num_of_wires */
  init_cost_array(costs, num_of_wires, dim_x, dim_y);
  printf("Complete initialize board\n");
  error = 0;

  init_time += duration_cast<dsec>(Clock::now() - init_start).count();
  printf("Initialization Time: %lf.\n", init_time);

  /**************************************
   **** START COMPUTATION 
   **************************************/
  auto compute_start = Clock::now();
  double compute_time = 0;
  omp_set_num_threads(num_of_threads);
#ifdef RUN_MIC /* Use RUN_MIC to distinguish between the target of compilation */

  /* This pragma means we want the code in the following block be executed in
   * Xeon Phi.
   */
   // ALGO
    // 1. With probability 1 - P, choose the current min path.  Otherwise, choose a
    //    a path uniformly at random from the space of delt_x + delt_y possible routes.
    // 2. Calculate cost of current path, if not known. This is the current min path.
    // 3. Consider all paths which first travel horizontally.
    //    If any costs less than the current min path, that is the new min path.
    // 4. Same as (2), using vertical paths.

    // Idea for later ?? Split up work of updating cost array by cells versus by wires
    //                   Structure to store "no touch points" (i.e. pt's with higher costs)??
    //                   Sort path points before updating cost array --> LOCALITY
#pragma offload target(mic) \
  inout(wires: length(num_of_wires) INOUT)    \
  inout(costs: length(dim_x*dim_y) INOUT)
#endif
  {
    // PRIVATE variables
    int i, j, x, y, w;
    path_t *mypath;
    int s_x, s_y, e_x, e_y;
    int b1_x, b1_y, b2_x, b2_y;

    // SHARED variables
    cost_t *C = costs

    /*@@@@@@@@@@@@@@ INIT LOOP @@@@@@@@@@@@@@*/
    /* Initialize all 'first' paths (create a start board) */
    #pragma omp parallel for default(shared)                       \
                         private(i) shared(wires) schedule(dynamic)
      for (i = 0; i < num_of_wires; i++) 
      {
        new_rand_path( &(wires[i]) );
      } /* implicit barrier */
    /* ############# END PRAGMA ############ */
    /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

    /*@@@@@@@@@@@@@@ MAIN LOOP @@@@@@@@@@@@@@*/
    for (i = 0; i < SA_iters; i++) // N iterations
    {
      /* Layout board */
      #pragma omp parallel for default(shared) \
          private(y, x) shared(B) schedule(dynamic)
        for( y = 0; y < dim_y; y++){
          for( x = 0; x < dim_x; x++){
            B[y*dim_y + x].val = 0;  // clean up board
          }
        }
      /* Update cost array */
      #pragma omp parallel for default(shared)                            \
              private(x,y) shared(B) schedule(dynamic)
        for (y = 0; y < dimY; j++) 
        {
          // Initialize wire private variables
          mypath    = wires[j].currentPath;
          num_bends = mypath->numBends;
          bends     = mypath->bends;
          bounds    = mypath->bounds;
          s_x = bounds[0];   // (start point)
          s_y = bounds[1];
          e_x = bounds[2];   // (end point)
          e_y = bounds[3];

          // Follow path & update cost array
          switch (num_bends) {
            case 0:
              if (s_y == e_y) // Horizontal path
                path_cost(costs, s_y, s_x, e_x);
              else            // Vertical path
                path_cost(costs, s_x, s_y, e_y);
            case 1: 
              b1_x = bends[0]; 
              b1_y = bends[1]; // Get bend coordinate
              if (s_y == b1_y) // Before bend is horizontal
              {
                path_cost(costs, s_y, s_x, b1_x);
                // horizontalCost(B, row, s_x, b1_x);
                // After bend must be vertical 
                path_cost(costs, s_x, s_y, b1_y);
                // verticalCost(B, row, b1_x, e_x, b1_y, e_y, dim_y);
              }
              else             // Before bend is vertical
              {
                path_cost(costs, b1_x, b1_y, e_y);
                // verticalCost(B, row, s_x, b1_x, s_y, b1_y, dim_y);
                // After bend must be horizontal 
                path_costs(costs, b1_y, b1_x, e_x);
                // horizontalCost(B, row, b1_x, e_x);
              }
            case 2: 
              b1_x = bends[0]; // Get both bend coordinates
              b1_y = bends[1];
              b2_x = bends[2];
              b2_y = bends[3]; 
              while (num_bends != 0)
              {
                if (s_y == b1_y) // Before bend is horizontal
                {
                  path_cost(costs, s_y, s_x, b1_x);
                  // horizontalCost(B, row, s_x, b1_x);
                  // After bend must be vertical 
                  path_cost(costs, s_x, s_y, b1_y);
                  // verticalCost(B, row, b1_x, e_x, b1_y, e_y, dim_y);
                }
                else             // Before bend is vertical
                {
                  path_cost(costs, b1_x, b1_y, e_y);
                  // verticalCost(B, row, s_x, b1_x, s_y, b1_y, dim_y);
                  // After bend must be horizontal 
                  path_costs(costs, b1_y, b1_x, e_x);
                  // horizontalCost(B, row, b1_x, e_x);
                }
                /* Update for next iteration */
                // start = bend1
                s_x = b1_x; s_y = b1_y;
                // bend1 = bend2
                b1_x = b2_x; b1_y = b2_y;
                // bend2 = end
                b2_x = e_x; b2_y = e_y;
                num_bends--;
              }
          }
        } /* implicit barrier */
      /* ############## END PRAGMA ############# */

      /* Parallel by wire, calculate cost of current path */

      /* Parallel by wire, determine NEW path */
      // With probability 1 - P, choose the current min path.

      // Otherwise, choose a path at random
    }
    /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  }
  /* #################### END PRAGMA ################### */

  compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
  printf("Computation Time: %lf.\n", compute_time);

  /* Write wires and costs to files */
  FILE *outputWire, *outputCost;
  char costFileName[128] = "costs_";
  char wireFileName[128] = "output_";
  strcat(costFileName, argv[0]);
  strcat(wireFileName, argv[0]);
  strcat(wireFileName, "_");
  strcat(costFileName, "_");
  char buf[5];
  snprintf(buf, sizeof(buf), "%d", num_of_threads);
  strcat(wireFileName, buf);
  strcat(costFileName, buf);
  strcat(wireFileName, ".txt");
  strcat(costFileName, ".txt");

  outputWire = fopen(wireFileName, "w");
  outputCost = fopen(costFileName, "w");
  if (outputCost == NULL || outputWire == NULL)
  {
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(outputWire, "%d %d\n", dim_x, dim_y);
  fprintf(outputCost, "%d %d\n", dim_x, dim_y);
  /*wrting to Cost */
  for(int row = 0 ; row < dim_y; row++){
    for(int col = 0; col < dim_x; col++){
      fprintf(outputCost, "%d ", costs->board[row*dim_y + col].val);
    }
    fprintf(outputCost, "\n");
  }
  /*wrting to wire */
  fprintf(outputWire, "%d\n", num_of_wires);
  for(int w_count = 0; w_count < num_of_wires; w_count++){
    fprintf(outputWire, "%d %d ", wires[w_count].currentPath->bounds[0],
                       wires[w_count].currentPath->bounds[1]);
    if(wires[w_count].currentPath->numBends == 2){
      fprintf(outputWire, "%d %d ", wires[w_count].currentPath->bends[0],
                       wires[w_count].currentPath->bends[1]);
      fprintf(outputWire, "%d %d ", wires[w_count].currentPath->bends[2],
                       wires[w_count].currentPath->bends[3]);
    }
    if(wires[w_count].currentPath->numBends == 1){
      fprintf(outputWire, "%d %d ", wires[w_count].currentPath->bends[0],
                       wires[w_count].currentPath->bends[1]);
    }
    fprintf(outputWire, "%d %d\n", wires[w_count].currentPath->bounds[2],
                       wires[w_count].currentPath->bounds[3]);
  }
  fclose(outputCost);
  fclose(outputWire);

  // free the allocated wire and DS
  for( int y = 0; y < dim_y; y++){
    for( int x = 0; x < dim_x; x++){
      omp_destroy_lock(costs->board[y*dim_y + x].lock);
    }
  }
  for(int i = 0; i < num_of_wires; i++){
    delete[] wires[i].currentPath;
    delete[] wires[i].prevPath;
  }
  delete[] wires;
  delete[] costs->board;
  delete[] costs;
  return 0;
}
