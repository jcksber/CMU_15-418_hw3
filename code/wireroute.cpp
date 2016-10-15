/**
 * Parallel VLSI Wire Routing via OpenMP
 * Jack Kasbeer (jkasbeer), Qifang Cai (qcai)
 */

#include "wireroute.h"
#include <chrono>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <assert.h>
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
// HELPER FUNCTIONS
/////////////////////////////////////

/* new_rand_path *
 * Generate a random path in the space of delta_x + delta_y
 * 50% change pick x traversal  50% chance pick y traversal
 * random generate bends
 */
void new_rand_path(wire_t *wire){
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

/* horizontalCost *
 * Update cost array for horizontal traversal
 * Input: ptr to board, y coord, starting x, ending x, dim_y
 */
void horizontalCost(cost_cell_t *C, int row, int startX, int endX, int dimY){
  int s_x = startX;
  // Determine path direction
  int dir = startX > endX ? -1 : 1;
  /* Update cost array for given wire */
  while (s_x != endX){
    /*### UPDATING CELL: CRITICAL REGION ###*/
      incrCell(C, s_x, row, dimY);
    /*######################################*/
    s_x += dir; // add/subtract a column
  }
}

/* verticalCost *
 * Update cost array for vertical traversal
 * Input: ptr to board, x coord, starting y, ending y, dim_y
 */
void verticalCost(cost_cell_t *C, int xCoord, int startY, int endY, int dimY){
  int s_y = startY;
  // Determine path direction
  int dir = startY > endY ? -1 : 1;
  /* Update cost array for given wire */
  while (s_y != endY){
    /*### UPDATING CELL: CRITICAL REGION ###*/
      incrCell(C, xCoord, s_y, dimY);
    /*######################################*/
    s_y += dir;
  }
}

// Use cell level lock to safely incre value by 1
// INPUT: ptr to board, x coord , y coord, dim_y
void incrCell(cost_cell_t *C, int x, int y, int dimY){
  cost_cell_t *c;
  c = &C[y*dimY + x]; // calculate the idx in board
  omp_set_lock(&c->lock);
    c->val +=1;
  omp_unset_lock(&c->lock);
}

/* use to run board statistic  */
void updateBoard(cost_t *board){
  // overwrite the previous data
  board->prevMax = board->currentMax;
  board->prevAggrTotal = board->currentAggrTotal;
  int Max = 0;
  int Total = 0;
  // traversal to count the board
  for (int row = 0; row < board->dimY; row++){
    for (int col = 0;  col < board->dimX; col++){
      int val = board->board[row* board->dimY + col].val;
      if(val > Max) Max = val;
      if(val > 1) Total += val;
    }
  }
  board->currentMax = Max;
  board->currentAggrTotal = Total;
}

// read a value in the board
inline int readBoard(cost_t *board, int x, int y){
  return board->board[y*board->dimY + x].val;
}

// get vertical cell values
value_t readVertical(cost_t* board, int x, int s_y, int e_y){
  value_t result;
  result.aggr_max = 0;
  result.m = 0;
  int dir = s_y > e_y ? -1:1;
  int c = s_y;
  while(c != e_y){
    int val = readBoard(board,x,c);
    if(result.m < val) result.m = val;
    if(val > 1) result.aggr_max += val;
    c += dir;
  }
  return result;
}

// get horizontal cell values
value_t readHorizontal(cost_t* board, int y, int s_x, int e_x){
  value_t result;
  result.aggr_max = 0;
  result.m = 0;
  int dir = s_x > e_x ? -1:1;
  int c = s_x;
  while(c != e_x){
    int val = readBoard(board,c,y);
    if(result.m < val) result.m = val;
    if(val > 1) result.aggr_max += val;
    c += dir;
  }
  return result;
}

// combine to value_t into one
value_t combineValue( value_t v1, value_t v2){
  value_t ret;
  ret.aggr_max = v1.aggr_max + v2.aggr_max;
  ret.m = (v1.m > v2.m) ? v1.m : v2.m;
  return ret;
}

/////// board cost calculation
value_t calculatePath(cost_t* board, int s_x, int s_y, int e_x, int e_y,
          int numBends, int b1_x, int b1_y, int b2_x, int b2_y){
  value_t result, temp, temp1, temp2;
  int tmp_val = readBoard(board, e_x, e_y);
  result.aggr_max = 0;
  result.m = 0;
  // Follow path & update cost array
  switch (numBends) {
    case 0:
      if (s_y == e_y){ // Horizontal path
        temp = readHorizontal(board, s_y, s_x, e_x);
        if (tmp_val > 1) result.aggr_max = temp.aggr_max + tmp_val;
        else result.aggr_max = temp.aggr_max;
        result.m = (temp.m > tmp_val) ? temp.m : tmp_val;
        break;
      }
      if (s_x == e_x){            // Vertical path
        temp = readVertical(board, s_x, s_y, e_y);
        if (tmp_val > 1) result.aggr_max = temp.aggr_max + tmp_val;
        else result.aggr_max = temp.aggr_max;
        result.m = (temp.m > tmp_val) ? temp.m : tmp_val;
        break;
      }
    case 1:
      if (s_y == b1_y) // Before bend is horizontal
      {
        temp1 = combineValue(readHorizontal(board, s_y, s_x, b1_x),
            // After bend must be vertical
            readVertical(board, e_x, b1_y, e_y));
        if (tmp_val > 1) result.aggr_max = temp1.aggr_max + tmp_val;
        else result.aggr_max = temp1.aggr_max;
        result.m = (temp1.m > tmp_val) ? temp1.m : tmp_val;
        break;
      }
      if (s_x == b1_x)           // Before bend is vertical
      {
        temp1 = combineValue(readVertical(board, s_x, s_y, b1_y),
            // After bend must be horizontal
              readHorizontal(board, e_y, b1_x, e_x));
        if (tmp_val > 1) result.aggr_max = temp1.aggr_max + tmp_val;
        else result.aggr_max = temp1.aggr_max;
        result.m = (temp1.m > tmp_val) ? temp1.m : tmp_val;
        break;
      }
    case 2:
      if (s_y == b1_y) // Before bend is horizontal
      {
        temp = combineValue(readHorizontal(board, s_y, s_x, b1_x),
                readVertical(board, b1_x, b1_y, b2_y));//after bend is vertical
        temp2 = combineValue(temp, readHorizontal(board, e_y, b2_x, e_x));
        if (tmp_val > 1) result.aggr_max = temp2.aggr_max + tmp_val;
        else result.aggr_max = temp2.aggr_max;
        result.m = (temp2.m > tmp_val) ? temp2.m : tmp_val;
        break;
      }
      if (s_x == b1_x) // Before bend is vertical
      {
        temp = combineValue(readVertical(board, s_x, s_y, b1_y),
            readHorizontal(board, b1_y, b1_x, b2_x));//after bend is horizontal
        temp2 = combineValue(temp, readVertical(board, b2_x, b2_y, e_y));
        if (tmp_val > 1) result.aggr_max = temp2.aggr_max + tmp_val;
        else result.aggr_max = temp2.aggr_max;
        result.m = (temp2.m > tmp_val) ? temp2.m : tmp_val;
        break;
      }
  }
  return result;
}

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

  /* Allocate for array of wires */
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
  printf("Complete read wires: %d\n", count);
  /* Allocate for cost array struct */
  cost_t *costs = (cost_t *)calloc(1, sizeof(cost_t));
  costs->dimX = dim_x;
  costs->dimY = dim_y;
  costs->currentMax = num_of_wires;
  costs->board = (cost_cell_t *)calloc(dim_x * dim_y, sizeof(cost_cell_t));
  printf("Complete allocate board\n");

  /* Initialize cell level locks */
  for( int y = 0; y < dim_y; y++){
    for( int x = 0; x < dim_x; x++){
      omp_init_lock(&(costs->board[y*dim_y + x].lock));
    }
  }
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
    int i, j, x, y, w, row, col;
    path_t *mypath;
    value_t localMax, tempMax;
    int s_x, s_y, e_x, e_y;
    int b1_x, b1_y, b2_x, b2_y;
    int n1_x, n1_y, n2_x, n2_y;
    int nBend, dir;
    // SHARED variables
    cost_cell_t *B = costs->board;
    /* ########## PARALLEL BY WIRE ##########*/
    /* Initialize all 'first' paths (create a start board) */
    #pragma omp parallel for default(shared)                       \
      private(w) shared(wires) schedule(dynamic)
    for (w = 0; w < num_of_wires; w++){
      new_rand_path( &(wires[w]) );
    } /* implicit barrier */

    /*@@@@@@@@@@@@@@ MAIN LOOP @@@@@@@@@@@@@@*/
    for (i = 0; i < SA_iters; i++){
      // Clean up the board
      #pragma omp parallel for default(shared) \
        private(y, x) shared(B) schedule(dynamic)
      for( y = 0; y < dim_y; y++){
        for( x = 0; x < dim_x; x++){
          B[y*dim_y + x].val = 0;  // clean up board
        }
      }
      /*  layout board */
      #pragma omp parallel for default(shared)                            \
        private(b1_x, b2_x, b1_y, b2_y,j,mypath,s_x,s_y,e_x,e_y)\
          shared(wires, B) schedule(dynamic)
      for (j = 0; j < num_of_wires; j++){
        // Initialize wire private variables
        mypath    = wires[j].currentPath;
        s_x = mypath->bounds[0];   // (start point)
        s_y = mypath->bounds[1];
        e_x = mypath->bounds[2];   // (end point)
        e_y = mypath->bounds[3];
        // Follow path & update cost array
        switch (mypath->numBends) {
          case 0:
            if (s_y == e_y){ // Horizontal path
              horizontalCost(B, s_y, s_x, e_x, dim_y);
              incrCell(B, e_x, e_y, dim_y);
              break;
            }
            if (s_x == e_x){            // Vertical path
              verticalCost(B, e_x, s_y, e_y, dim_y);
              incrCell(B, e_x, e_y, dim_y);
              break;
            }
          case 1:
            b1_x = mypath->bends[0];
            b1_y = mypath->bends[1]; // Get bend coordinate
            if (s_y == b1_y) // Before bend is horizontal
            {
              horizontalCost(B, s_y, s_x, b1_x, dim_y);
              // After bend must be vertical
              verticalCost(B, e_x, b1_y, e_y, dim_y);
              incrCell(B, e_x, e_y, dim_y);
              break;
            }
            if (s_x == b1_x)           // Before bend is vertical
            {
              verticalCost(B, s_x, s_y, b1_y, dim_y);
              // After bend must be horizontal
              horizontalCost(B, e_y, b1_x, e_x, dim_y);
              incrCell(B, e_x, e_y, dim_y);
              break;
            }
          case 2:
            b1_x = mypath->bends[0]; // Get both bend coordinates
            b1_y = mypath->bends[1];
            b2_x = mypath->bends[2];
            b2_y = mypath->bends[3];
            // Exam first bend
            if (s_y == b1_y) // Before bend is horizontal
            {
              horizontalCost(B, s_y, s_x, b1_x, dim_y);
              verticalCost(B, b1_x, b1_y, b2_y, dim_y);//after bend is vertical
              horizontalCost(B, e_y, b2_x, e_x,dim_y);
              incrCell(B, e_x, e_y, dim_y);
              break;
            }
            if (s_x == b1_x) // Before bend is vertical
            {
              verticalCost(B, s_x, s_y, b1_y, dim_y);
              horizontalCost(B, b1_y, b1_x, b2_x, dim_y);//after bend is horizontal
              verticalCost(B, b2_x, b2_y, e_y, dim_y);
              incrCell(B, e_x, e_y, dim_y);
              break;
            }
        }
      } /* implicit barrier */
      /* Parallel by wire, calculate cost of current path */
      /* Parallel by wire, determine NEW path */
      #pragma omp parallel for default(shared)       \
          private(w,row, col,  mypath, localMax, tempMax, s_x, s_y, e_x, e_y, b1_x, b2_x, \
                  b1_y, b2_y, n1_x, n1_y, n2_x, n2_y, nBend, dir) \
              shared(wires, costs) schedule(dynamic)
      for (w = 0; w < num_of_wires; w++){
        // With probability 1 - P, choose the current min path.
        srand(time(NULL));
        if((rand()%100) > int(SA_prob*100)){ // xx% chance pick the complicated  algo
          //printf("Running Complicated Path Generating Algo %d %d  \n", rand() %100, int(SA_prob*100));
          mypath = wires[w].currentPath;
          s_x = mypath->bounds[0];   // (start point)
          s_y = mypath->bounds[1];
          e_x = mypath->bounds[2];   // (end point)
          e_y = mypath->bounds[3];
          b1_x = mypath->bends[0];
          b1_y = mypath->bends[1];
          b2_x = mypath->bends[2];
          b2_y = mypath->bends[3];
          n1_x = b1_x;
          n1_y = b1_y;
          n2_x = b2_x;
          n2_y = b2_y;
          nBend = mypath->numBends;
          if ( s_x != e_x && s_y != e_y){
            localMax = calculatePath(costs, s_x, s_y, e_x, e_y, nBend, b1_x, b1_y, b2_x, b2_y);
            // case of one bend, at the end points
            // -> horizontal first:
            tempMax = calculatePath(costs, s_x, s_y, e_x, e_y, 1 ,e_x, s_y, 0, 0);
            if(tempMax.m < localMax.m && tempMax.aggr_max < localMax.aggr_max){
              localMax.m = tempMax.m;
              localMax.aggr_max = tempMax.aggr_max;
              nBend = 1;
              n1_x = e_x;
              n1_y = s_y;
            }
            // -> vertical one bend
            tempMax = calculatePath(costs, s_x, s_y, e_x, e_y, 1 ,s_x, e_y, 0, 0);
            if(tempMax.m < localMax.m && tempMax.aggr_max < localMax.aggr_max){
              localMax.m = tempMax.m;
              localMax.aggr_max = tempMax.aggr_max;
              nBend = 1;
              n1_x = s_x;
              n1_y = e_y;
            }
            // calculate horizontal path
            dir = (e_x > s_x) ? 1 : -1;
            for ( col = s_x + dir; col != e_x; col += dir){
              tempMax = calculatePath(costs, s_x, s_y, e_x, e_y,2, col, s_y, col , e_y);
              if(tempMax.m < localMax.m && tempMax.aggr_max < localMax.aggr_max){
                localMax.m = tempMax.m;
                localMax.aggr_max = tempMax.aggr_max;
                nBend = 2;
                n1_x = col;
                n1_y = s_y;
                n2_x = col;
                n2_y = e_y;
              }
            }
            // calculate vertical path
            dir = (e_y > s_y) ? 1 : -1;
            for(row = s_y + dir; row != e_y; row += dir){
              tempMax = calculatePath(costs, s_x, s_y, e_x, e_y, 2, s_x, row, e_x, row);
              if(tempMax.m < localMax.m && tempMax.aggr_max < localMax.aggr_max){
                localMax.m = tempMax.m;
                localMax.aggr_max = tempMax.aggr_max;
                nBend = 2;
                n1_x = s_x;
                n1_y = row;
                n2_x = e_x;
                n2_y = row;
              }
            }
          }
          // set new wire
          std::memcpy(wires[w].prevPath, mypath, sizeof(wire_t));
          mypath->numBends = nBend;
          mypath->bends[0] = n1_x;
          mypath->bends[1] = n1_y;
          mypath->bends[2] = n2_x;
          mypath->bends[3] = n2_y;
        }
        else{ // xx% chance take random path
          new_rand_path( &(wires[w]) );
        } /* implicit barrier */
      }
      // Finish picking the new path
    } /*  end iterations*/

    ////////////////////////////////////////////////////////////////////////////
     // clean up board
    #pragma omp parallel for default(shared) \
      private(y, x) shared(B) schedule(dynamic)
    for( y = 0; y < dim_y; y++){
      for( x = 0; x < dim_x; x++){
        B[y*dim_y + x].val = 0;
      }
    }
    /*  layout final result board  */
    #pragma omp parallel for default(shared)                            \
      private(j,mypath,s_x,s_y,e_x,e_y,b1_x, b2_x, b1_y, b2_y)\
        shared(wires, B) schedule(dynamic)
    for (j = 0; j < num_of_wires; j++){
      // Initialize wire private variables
      mypath    = wires[j].currentPath;
      s_x = mypath->bounds[0];   // (start point)
      s_y = mypath->bounds[1];
      e_x = mypath->bounds[2];   // (end point)
      e_y = mypath->bounds[3];
      // Follow path & update cost array
      switch (mypath->numBends) {
        case 0:
          if (s_y == e_y){ // Horizontal path
            horizontalCost(B, s_y, s_x, e_x, dim_y);
            incrCell(B, e_x, e_y, dim_y);
            break;
          }
          if (s_x == e_x){            // Vertical path
            verticalCost(B, e_x, s_y, e_y, dim_y);
            incrCell(B, e_x, e_y, dim_y);
            break;
          }
        case 1:
          b1_x = mypath->bends[0];
          b1_y = mypath->bends[1]; // Get bend coordinate
          if (s_y == b1_y) // Before bend is horizontal
          {
            horizontalCost(B, s_y, s_x, b1_x, dim_y);
            // After bend must be vertical
            verticalCost(B, e_x, b1_y, e_y, dim_y);
            incrCell(B, e_x, e_y, dim_y);
            break;
          }
          if (s_x == b1_x)           // Before bend is vertical
          {
            verticalCost(B, s_x, s_y, b1_y, dim_y);
            // After bend must be horizontal
            horizontalCost(B, e_y, b1_x, e_x, dim_y);
            incrCell(B, e_x, e_y, dim_y);
            break;
          }
        case 2:
          b1_x = mypath->bends[0]; // Get both bend coordinates
          b1_y = mypath->bends[1];
          b2_x = mypath->bends[2];
          b2_y = mypath->bends[3];
          // Exam first bend
          if (s_y == b1_y) // Before bend is horizontal
          {
            horizontalCost(B, s_y, s_x, b1_x, dim_y);
            verticalCost(B, b1_x, b1_y, b2_y, dim_y);//after bend is vertical
            horizontalCost(B, e_y, b2_x, e_x,dim_y);
            incrCell(B, e_x, e_y, dim_y);
            break;
          }
          if (s_x == b1_x) // Before bend is vertical
          {
            verticalCost(B, s_x, s_y, b1_y, dim_y);
            horizontalCost(B, b1_y, b1_x, b2_x, dim_y);//after bend is horizontal
            verticalCost(B, b2_x, b2_y, e_y, dim_y);
            incrCell(B, e_x, e_y, dim_y);
            break;
          }
      }
    } /* implicit barrier */
    ///////////////////////////////////////////////////////
  }
  /* #################### END PRAGMA ################### */

  compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
  printf("Computation Time: %lf.\n", compute_time);
  // update board statistic ////////////////
  updateBoard(costs);
  /////////////////////////////
  /* Write wires and costs to files */
  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != NULL)
    fprintf(stdout, "Current working dir: %s\n", cwd);
  else
    perror("getcwd() error");
  FILE *outputWire, *outputCost;
  char costFileName[1024];
  char wireFileName[1024];
  memset(costFileName, 0 , sizeof(costFileName));
  memset(wireFileName, 0 , sizeof(wireFileName));
  strcat(costFileName, cwd);
  strcat(wireFileName, cwd);
  strcat(costFileName, "/costs_");
  strcat(wireFileName, "/output_");
  // clean input_file_name
  while(strlen(input_filename) != 0 && strchr(input_filename, '/') != NULL){
    input_filename = strchr(input_filename, '/')+1;
  }
  memcpy(cwd, input_filename, sizeof(cwd));
  char* ptr = strchr(cwd,'.');
  *ptr = '\0';
  strcat(costFileName, cwd);
  strcat(wireFileName, cwd);
  strcat(wireFileName, "_");
  strcat(costFileName, "_");
  char buf[5];
  snprintf(buf, sizeof(buf), "%d", num_of_threads);
  strcat(wireFileName, buf);
  strcat(costFileName, buf);
  strcat(wireFileName, ".txt");
  strcat(costFileName, ".txt");
  // print stat
  printf("Input File: %s has total aggregated cost: [%d] and max layers: [%d]\n", cwd,
              costs->currentAggrTotal, costs->currentMax);
  outputWire = fopen(wireFileName, "w");
  outputCost = fopen(costFileName, "w");
  if (outputCost == NULL || outputWire == NULL)
  {
    perror("Error opening file!\n");
    printf("filename : %s\n", wireFileName);
    printf("filename : %s\n", costFileName);
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
      omp_destroy_lock(&(costs->board[y*dim_y + x].lock));
    }
  }
  for(int i = 0; i < num_of_wires; i++){
    free(wires[i].currentPath);
    free(wires[i].prevPath);
  }
  free(wires);
  free(costs->board);
  free(costs);
  return 0;
}
