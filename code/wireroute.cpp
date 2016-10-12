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
// HELPER FUNCTION
////////////////////////////////////

void new_rand_path(wire_t *wire){
  //overwrite previous pathi
  std::memcpy(wire->prevPath, wire->currentPath, sizeof(wire_t));
  int s_x, s_y, e_x, e_y, dx, dy;
  s_x = wire->currentPath->bounds[0];
  s_y = wire->currentPath->bounds[1];
  e_x = wire->currentPath->bounds[2];
  e_y = wire->currentPath->bounds[3];
  dx
}

/////////////////////////////////////////

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

  // Parse for dimensions & num wires
  int dim_x, dim_y;
  int num_of_wires;
  fscanf(input, "%d %d\n", &dim_x, &dim_y);
  fscanf(input, "%d\n", &num_of_wires);
  // Allocate for array of wires
  wire_t *wires = (wire_t *)calloc(num_of_wires, sizeof(wire_t));

  /* Read the grid dimension and wire information from file */
  int count = 0;
  while(count < num_of_wires){
    int s_x, s_y, e_x, e_y;
    fscanf(input, "%d %d %d %d\n", &s_x, &s_y, &e_x, &e_y);
    wires[count].prevPath = (path_t*)calloc(1, sizeof(path_t));
    wires[count].currentPath = (path_t*)calloc(1, sizeof(path_t));
    wires[count].currentPath->numBends = 0;
    wires[count].currentPath->bounds[0] = s_x;
    wires[count].currentPath->bounds[1] = s_y;
    wires[count].currentPath->bounds[2] = e_x;
    wires[count].currentPath->bounds[3] = e_y;
    count++;
  }
  cost_t *costs = (cost_t *)calloc(1, sizeof(cost_t));
  costs->prev_max = num_of_wires;
  // no init for the prev_total of cost_t
  costs->board = (cost_cell_t *)calloc(dim_x * dim_y, sizeof(cost_cell_t));
  /* Initialize cost matrix */
  for( int y = 0; y < dim_y; y++){
    for( int x = 0; x < dim_x; x++){
      costs->board[y*dim_y + x].val = 0;
      omp_init_lock(&costs->board[y*dim_y + x].lock);
    }
  }

  /* Initailize additional data structures needed in the algorithm */
  // 1. Structure to store "no touch points" (i.e. pt's with higher costs)??
  error = 0;

  init_time += duration_cast<dsec>(Clock::now() - init_start).count();
  printf("Initialization Time: %lf.\n", init_time);

  // Init Param
  omp_set_num_threads(num_of_threads);

  auto compute_start = Clock::now();
  double compute_time = 0;
#ifdef RUN_MIC /* Use RUN_MIC to distinguish between the target of compilation */

  /* This pragma means we want the code in the following block be executed in
   * Xeon Phi.
   */
#pragma offload target(mic) \
  inout(wires: length(num_of_wires) INOUT)    \
  inout(costs: length(dim_x*dim_y) INOUT)
#endif
  {
    /* Implement the wire routing algorithm here
     * Feel free to structure the algorithm into different functions
     * Don't use global variables.
     * Use OpenMP to parallelize the algorithm. 
     */

    // PRIVATE variables
    int i;
    // SHARED variables
    cost_cell_t *B = costs->board;

    // ALGO
    // 1. With probability 1 - P, choose the current min path.  Otherwise, choose a
    //    a path uniformly at random from the space of delt_x + delt_y possible routes.
    // 2. Calculate cost of current path, if not known. This is the current min path.
    // 3. Consider all paths which first travel horizontally.
    //    If any costs less than the current min path, that is the new min path.
    // 4. Same as (2), using vertical paths.

    // Idea for later ?? Split up work of updating cost array by cells versus by wires

    /* INIT LOOP */
    /* Parallel by wire, initialize all wire 'first' paths (create a start board) */
    #pragma omp parallel for             \
                         default(shared) \
                         private(i)      \
                         shared(wires)   \
                         schedule(dynamic)
      for (i = 0; i < num_of_threads; i++) // 1 iteration
      {
        new_rand_path( &(wires[i]) );
      } /* implicit barrier */

    /* MAIN LOOP */
    for (i = 0; i < SA_iters; i++) // N iterations
    {
      /* Parallel by wire, update cost array */
      #pragma omp parallel for              \
                           default(shared)  \
                           private(j)       \
                           shared(wires, B) \
                           schedule(dynamic)
        for (j = 0; j < num_of_wires; j++) // 1 iteration
        {

        } /* implicit barrier */

      /* Parallel by wire, calculate cost of current path */

      /* Parallel by wire, determine NEW path */
      // With probability 1 - P, choose the current min path.

      // Otherwise, choose a path at random
    }
  }

  compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
  printf("Computation Time: %lf.\n", compute_time);

  /* Write wires and costs to files */
  // Output current wire set
  // Ouptut current cost set


  // free the alocated wire
  for(int i = 0; i < num_of_wires; i++){
    delete[] wires[i].currentPath;
    delete[] wires[i].prevPath;
  }
  delete[] wires;
  delete[] costs->board;
  delete[] costs;
  return 0;
}
