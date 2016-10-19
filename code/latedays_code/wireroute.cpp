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

void new_rand_path(wire_t* wire, int *rand, int wire_n){
  //overwrite previous path
  int bend = 0;
  int s_x, s_y, e_x, e_y, dy, yp, dx, xp;
  s_x = wire->bound0;
  s_y = wire->bound1;
  e_x = wire->bound2;
  e_y = wire->bound3;
  if (s_x == e_x || s_y == e_y){
    return;
  }
  // not in` the same line, need at least one bend
  bend +=1;
  dy = abs(e_y - s_y);
  dx = abs(e_x - s_x);
  // calculate random point on y axis
  if ((rand[wire_n]% 10)> 5){
    // y first
    int ran_y = rand[wire_n] % dy;
    if (ran_y == 0) ran_y = 1; // need to make progress
    if(s_y > e_y)  yp = s_y - ran_y;
    else yp = s_y + ran_y;
    // determind bends
    if(e_y != yp) bend +=1;
    // overwrite
    wire->bends0 = s_x;
    wire->bends1 = yp;
    wire->bends2 = e_x;
    wire->bends3 = yp;
    wire->numBends = bend;
  }
  else{
    // x first traversal
    int ran_x = rand[wire_n] % dx;
    if (ran_x == 0) ran_x = 1;
    if(s_x > e_x)  xp = s_x - ran_x;
    else xp = s_x + ran_x;
    // determind bends
    if(e_x != xp) bend +=1;
    // overwrite
    wire->bends0 = xp;
    wire->bends1 = s_y;
    wire->bends2 = xp;
    wire->bends3 = e_y;
    wire->numBends = bend;
  }
}

/* horizontalCost *
 * Update cost array for horizontal traversal
 * Input: ptr to board, y coord, starting x, ending x, dim_y
 */
void horizontalCost(cost_cell_t *C, omp_lock_t *locks, int row, int startX, int endX, int dimY, int wire_n){
  int s_x = startX;
  // Determine path direction
  int dir = startX > endX ? -1 : 1;
  /* Update cost array for given wire */
  while (s_x != endX){
    /*### UPDATING CELL: CRITICAL REGION ###*/
      incrCell(C, locks, s_x, row, dimY, wire_n);
    /*######################################*/
    s_x += dir; // add/subtract a column
  }
}

/* verticalCost *
 * Update cost array for vertical traversal
 * Input: ptr to board, x coord, starting y, ending y, dim_y
 */
void verticalCost(cost_cell_t *C, omp_lock_t *locks, int xCoord, int startY, int endY, int dimY, int wire_n){
  int s_y = startY;
  // Determine path direction
  int dir = startY > endY ? -1 : 1;
  /* Update cost array for given wire */
  while (s_y != endY){
    /*### UPDATING CELL: CRITICAL REGION ###*/
      incrCell(C, locks, xCoord, s_y, dimY, wire_n);
    /*######################################*/
    s_y += dir;
  }
}

// Use cell level lock to safely incre value by 1
// INPUT: ptr to board, x coord , y coord, dim_y
inline void incrCell(cost_cell_t *C, omp_lock_t *lock, int x, int y, int dimY, int wire_n){
  cost_cell_t c = C[y*dimY + x]; // calculate the idx in board
  omp_set_lock(&lock[y*dimY + x]);
  //printf("here is the c val %d\n", C[y*dimY + x].val);
  /*
  switch(c.wire){
    case 0:
      c.list0 = wire_n;
      c.wire += 1;
      break;
    case 1:
      c.list1 = wire_n;
      c.wire += 1;
      break;
    case 2:
      c.list2 = wire_n;
      c.wire += 1;
      break;
    default:
      break;
  }
  */
  omp_unset_lock(&lock[y*dimY + x]);
}

/* use to run board statistic  */
void updateBoard(cost_t *board, cost_cell_t *B){
  // overwrite the previous data
  board->prevMax = board->currentMax;
  board->prevAggrTotal = board->currentAggrTotal;
  int Max = 0;
  int Total = 0;
  // traversal to count the board
  for (int row = 0; row < board->dimY; row++){
    for (int col = 0;  col < board->dimX; col++){
      int val = B[row* board->dimY + col].val;
      if(val > Max) Max = val;
      if(val > 1) Total += val;
    }
  }
  board->currentMax = Max;
  board->currentAggrTotal = Total;
}

// read a value in the board
int readBoard(cost_cell_t *board, int x, int y, int wire_n, int dimY){
  cost_cell_t *c =&board[y*dimY + x];
  if(wire_n == c->list0 || wire_n == c->list1 ||wire_n == c->list2 ){
      return c->val-1;
  }
  return c->val;
}

// get vertical cell values
value_t readVertical(cost_t* board, int x, int s_y, int e_y, int wire_n, cost_cell_t *B){
  value_t result;
  result.aggr_max = 0;
  result.m = 0;
  int dir = s_y > e_y ? -1:1;
  int c = s_y;
  while(c != e_y){
    int val = readBoard(B,x,c, wire_n, board->dimY);
    if(result.m < val) result.m = val;
    if(val > 1) result.aggr_max += val;
    c += dir;
  }
  return result;
}

// get horizontal cell values
value_t readHorizontal(cost_t* board, int y, int s_x, int e_x, int wire_n, cost_cell_t *B){
  value_t result;
  result.aggr_max = 0;
  result.m = 0;
  int dir = s_x > e_x ? -1:1;
  int c = s_x;
  while(c != e_x){
    int val = readBoard(B,c,y, wire_n, board->dimY);
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
value_t calculatePath(cost_t* board, cost_cell_t *B , int s_x, int s_y, int e_x, int e_y,
          int numBends, int b1_x, int b1_y, int b2_x, int b2_y, int wire_n){
  value_t result, temp, temp1, temp2;
  int tmp_val = readBoard(B, e_x, e_y, wire_n, board->dimY);
  result.aggr_max = 0;
  result.m = 0;
  // Follow path & update cost array
  switch (numBends) {
    case 0:
      if (s_y == e_y){ // Horizontal path
        temp = readHorizontal(board, s_y, s_x, e_x,wire_n, B);
        if (tmp_val > 1) result.aggr_max = temp.aggr_max + tmp_val;
        else result.aggr_max = temp.aggr_max;
        result.m = (temp.m > tmp_val) ? temp.m : tmp_val;
        break;
      }
      if (s_x == e_x){            // Vertical path
        temp = readVertical(board, s_x, s_y, e_y,wire_n, B);
        if (tmp_val > 1) result.aggr_max = temp.aggr_max + tmp_val;
        else result.aggr_max = temp.aggr_max;
        result.m = (temp.m > tmp_val) ? temp.m : tmp_val;
        break;
      }
    case 1:
      if (s_y == b1_y) // Before bend is horizontal
      {
        temp1 = combineValue(readHorizontal(board, s_y, s_x, b1_x,wire_n, B),
            // After bend must be vertical
            readVertical(board, e_x, b1_y, e_y,wire_n, B));
        if (tmp_val > 1) result.aggr_max = temp1.aggr_max + tmp_val;
        else result.aggr_max = temp1.aggr_max;
        result.m = (temp1.m > tmp_val) ? temp1.m : tmp_val;
        break;
      }
      if (s_x == b1_x)           // Before bend is vertical
      {
        temp1 = combineValue(readVertical(board, s_x, s_y, b1_y,wire_n, B),
            // After bend must be horizontal
              readHorizontal(board, e_y, b1_x, e_x,wire_n, B));
        if (tmp_val > 1) result.aggr_max = temp1.aggr_max + tmp_val;
        else result.aggr_max = temp1.aggr_max;
        result.m = (temp1.m > tmp_val) ? temp1.m : tmp_val;
        break;
      }
    case 2:
      if (s_y == b1_y) // Before bend is horizontal
      {
        temp = combineValue(readHorizontal(board, s_y, s_x, b1_x,wire_n, B),
                readVertical(board, b1_x, b1_y, b2_y,wire_n, B));//after bend is vertical
        temp2 = combineValue(temp, readHorizontal(board, e_y, b2_x, e_x,wire_n, B));
        if (tmp_val > 1) result.aggr_max = temp2.aggr_max + tmp_val;
        else result.aggr_max = temp2.aggr_max;
        result.m = (temp2.m > tmp_val) ? temp2.m : tmp_val;
        break;
      }
      if (s_x == b1_x) // Before bend is vertical
      {
        temp = combineValue(readVertical(board, s_x, s_y, b1_y,wire_n, B),
            readHorizontal(board, b1_y, b1_x, b2_x,wire_n, B));//after bend is horizontal
        temp2 = combineValue(temp, readVertical(board, b2_x, b2_y, e_y,wire_n, B));
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
  int SA_iters = get_option_int("-i", 5); // Need to pass to Phi

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
    wires[count].numBends = 0;
    wires[count].bends0 = -1;
    wires[count].bends1 = -1;
    wires[count].bends2 = -1;
    wires[count].bends3 = -1;
    wires[count].bound0 = s_x;
    wires[count].bound1 = s_y;
    wires[count].bound2 = e_x;
    wires[count].bound3 = e_y;
    count++;
  }
  printf("Complete read wires: %d\n", count);
  /* Allocate for cost array struct */
  cost_t *costs = (cost_t *)calloc(1, sizeof(cost_t));
  costs->dimX = dim_x;
  costs->dimY = dim_y;
  costs->currentMax = num_of_wires;
  cost_cell_t *B =  (cost_cell_t *)calloc(dim_x * dim_y, sizeof(cost_cell_t));
  printf("Complete allocate board\n");

  printf("Complete initialize board\n");
  error = 0;

  init_time += duration_cast<dsec>(Clock::now() - init_start).count();
  printf("Initialization Time: %lf.\n", init_time);

  // make random
  int *rand_iter = (int *)calloc(SA_iters, sizeof(int));
  for (int i = 0; i < SA_iters; i++){
    srand(time(0));
    rand_iter[i] = (rand()%100);
  }
  int *rand_wire = (int *)calloc(num_of_wires, sizeof(int));
  for (int i = 0; i < num_of_wires; i++){
    srand(time(0));
    rand_wire[i] = (rand());
  }
  /**************************************
   **** START COMPUTATION
   **************************************/
  auto compute_start = Clock::now();
  double compute_time = 0;
  	printf(">>> START computing \n");
  omp_set_num_threads(num_of_threads);
  	printf(">>> set number of threads to : %d \n", num_of_threads);
#ifdef RUN_MIC /* Use RUN_MIC to distinguish between the target of compilation */
	printf(">>>Runing on Xeon Phi ");
  /* This pragma means we want the code in the following block be executed in
   * Xeon Phi.
   */
#pragma offload target(mic) \
   	inout(wires: length(num_of_wires) INOUT)    \
    in(rand_wire:length(num_of_wires)) \
    in(rand_iter: length(SA_iters))   inout(B: length(dim_x*dim_y) INOUT) \
    inout(costs: length(1) INOUT)
   // in(SA_prob:length(1))in(SA_iters: length(1)) in(dim_x: length(1))  in(dim_y: length(1)) in(num_of_wires: length(1))
#else
	printf(">>> Running on Host ");
#endif
  {

  	printf(">>> Inside Phi ");
    // PRIVATE variables
    int i, j, x, y, w, row, col;
    wire_t mypath;
    value_t localMax, tempMax;
    int s_x, s_y, e_x, e_y;
    int b1_x, b1_y, b2_x, b2_y;
    int n1_x, n1_y, n2_x, n2_y;
    int nBend, dir;
    omp_lock_t *locks = (omp_lock_t *)calloc(dim_x*dim_y, sizeof(omp_lock_t));
    #pragma omp parallel for default(shared)\
      private (i, j) shared(locks) schedule(dynamic)
    for ( i = 0; i < dim_y; i++){
      for (j = 0; j < dim_x; j++){
        omp_init_lock(&(locks[i*dim_y+j]));
      }
    }
    // SHARED variables
    //Initialize all 'first' paths (create a start board)
    #pragma omp parallel for default(shared)                       \
      private(w) shared(wires) schedule(dynamic)
    for (w = 0; w < num_of_wires; w++){
      new_rand_path( &(wires[w]), rand_wire, w);
    } // implicit barrier
  	printf(">>> Get wire initialized \n");
    //@@@@@@@@@@@@@@ MAIN LOOP @@@@@@@@@@@@@@
    for (i = 0; i < SA_iters; i++){
      // Clean up the board
      #pragma omp parallel for default(shared) \
        private(y, x) shared(B) schedule(dynamic)
      for( y = 0; y < dim_y; y++){
        for( x = 0; x < dim_x; x++){
          B[y*dim_y + x].val = 0;  // clean up board
          B[y*dim_y + x].wire = 0;  // clean up board
          B[y*dim_y + x].list0 = 0;
          B[y*dim_y + x].list1 = 0;
          B[y*dim_y + x].list2 = 0;
        }
      }
  	  printf(">>> CleanUp board \n");
        //  layout board
      #pragma omp parallel for default(shared)                            \
        private(b1_x, b2_x, b1_y, b2_y,j,mypath,s_x,s_y,e_x,e_y)\
          shared(wires, B, locks) schedule(dynamic)
      for (j = 0; j < num_of_wires; j++){
        // Initialize wire private variables
        mypath = wires[j];
        s_x = mypath.bound0;   // (start point)
        s_y = mypath.bound1;
        e_x = mypath.bound2;   // (end point)
        e_y = mypath.bound3;
        // Follow path & update cost array
        switch (mypath.numBends) {
          /*
          case 0:
            if (s_y == e_y){ // Horizontal path
              horizontalCost(B, locks,s_y, s_x, e_x, dim_y, j);
              incrCell(B, locks, e_x, e_y, dim_y, j);
              break;
            }
            if (s_x == e_x){            // Vertical path
              verticalCost(B, locks, e_x, s_y, e_y, dim_y, j);
              incrCell(B, locks, e_x, e_y, dim_y, j);
              break;
            }
          case 1:
            b1_x = mypath.bends0;
            b1_y = mypath.bends1; // Get bend coordinate
            if (s_y == b1_y) // Before bend is horizontal
            {
              horizontalCost(B, locks, s_y, s_x, b1_x, dim_y,j);
              // After bend must be vertical
              verticalCost(B, locks, e_x, b1_y, e_y, dim_y , j);
              incrCell(B, locks, e_x, e_y, dim_y, j);
              break;
            }
            if (s_x == b1_x)           // Before bend is vertical
            {
              verticalCost(B, locks, s_x, s_y, b1_y, dim_y, j);
              // After bend must be horizontal
              horizontalCost(B, locks,  e_y, b1_x, e_x, dim_y,j);
              incrCell(B, locks, e_x, e_y, dim_y, j);
              break;
            }
            */

          case 2:
            b1_x = mypath.bends0;
            b1_y = mypath.bends1;
            b2_x = mypath.bends2;
            b2_y = mypath.bends3;
            // Exam first bend
            printf(">>>>>>>>>>>> %d %d %d %d\n", e_x, e_y, dim_y, j);
            if (s_y == b1_y) // Before bend is horizontal
            {
             //  horizontalCost(B, locks, s_y, s_x, b1_x, dim_y, j);
             // verticalCost(B, locks,  b1_x, b1_y, b2_y, dim_y, j);//after bend is vertical
             // horizontalCost(B, locks, e_y, b2_x, e_x,dim_y, j);
              incrCell(B, locks, e_x, e_y, dim_y, j);
              printf("TWO BEND\n");
              break;
            }
            if (s_x == b1_x) // Before bend is vertical
            {
             // verticalCost(B, locks, s_x, s_y, b1_y, dim_y, j);
             // horizontalCost(B, locks, b1_y, b1_x, b2_x, dim_y, j);//after bend is horizontal
             // verticalCost(B, locks, b2_x, b2_y, e_y, dim_y, j);
             // incrCell(B, locks, e_x, e_y, dim_y, j);
              break;
            }
            break;
        }
      } // implicit barrier
    }
      /*
  	  printf(">>> Layout new board \n");
      // Parallel by wire, determine NEW path
      #pragma omp parallel for default(shared)       \
          private(w,row, col,  mypath, localMax, tempMax, s_x, s_y, e_x, e_y, b1_x, b2_x, \
                  b1_y, b2_y, n1_x, n1_y, n2_x, n2_y, nBend, dir) \
              shared(wires, costs, rand_wire,  B) schedule(dynamic)
      for (w = 0; w < num_of_wires; w++){
        // With probability 1 - P, choose the current min path.
        if(rand_iter[i]> int(SA_prob*100)){
          // xx% chance pick the complicated  algo
          mypath = wires[w];
          s_x = mypath.bound0;   // (start point)
          s_y = mypath.bound1;
          e_x = mypath.bound2;   // (end point)
          e_y = mypath.bound3;
          b1_x = mypath.bends0;
          b1_y = mypath.bends1;
          b2_x = mypath.bends2;
          b2_y = mypath.bends3;
          n1_x = b1_x;
          n1_y = b1_y;
          n2_x = b2_x;
          n2_y = b2_y;
          nBend = mypath.numBends;
          if ( s_x != e_x && s_y != e_y){
            localMax = calculatePath(costs, B, s_x, s_y, e_x, e_y, nBend, b1_x, b1_y, b2_x, b2_y, -1);
            // case of one bend, at the end points
            // -> horizontal first:
            tempMax = calculatePath(costs, B, s_x, s_y, e_x, e_y, 1 ,e_x, s_y, 0, 0, w);
            if(tempMax.m < localMax.m && tempMax.aggr_max < localMax.aggr_max){
              localMax.m = tempMax.m;
              localMax.aggr_max = tempMax.aggr_max;
              nBend = 1;
              n1_x = e_x;
              n1_y = s_y;
            }
            // -> vertical one bend
            tempMax = calculatePath(costs, B, s_x, s_y, e_x, e_y, 1 ,s_x, e_y, 0, 0, w);
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
              tempMax = calculatePath(costs,B, s_x, s_y, e_x, e_y,2, col, s_y, col , e_y, w);
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
              tempMax = calculatePath(costs,B, s_x, s_y, e_x, e_y, 2, s_x, row, e_x, row, w);
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
          mypath.numBends = nBend;
          mypath.bends0 = n1_x;
          mypath.bends1 = n1_y;
          mypath.bends2 = n2_x;
          mypath.bends3 = n2_y;
        }
        else{ // xx% chance take random path
          new_rand_path( &(wires[w]), rand_wire, w );
        } // implicit barrier
      }
  	  printf(">>> Finish pick a wire set \n");
      // Finish picking the new path
    } //  end iterations

    ////////////////////////////////////////////////////////////////////////////
     // clean up board
    #pragma omp parallel for default(shared) \
      private(y, x) shared(B) schedule(dynamic)
    for( y = 0; y < dim_y; y++){
      for( x = 0; x < dim_x; x++){
        B[y*dim_y + x].val = 0;
        B[y*dim_y + x].wire = 0;  // clean up board
      }
    }
  	printf(">>> Final clean up \n");
    //  layout final result board
    #pragma omp parallel for default(shared)                            \
      private(j,mypath,s_x,s_y,e_x,e_y,b1_x, b2_x, b1_y, b2_y)\
        shared(wires, B, locks ) schedule(dynamic)
    for (j = 0; j < num_of_wires; j++){
      // Initialize wire private variables
      mypath = wires[j];
      s_x = mypath.bound0;   // (start point)
      s_y = mypath.bound1;
      e_x = mypath.bound2;   // (end point)
      e_y = mypath.bound3;
      // Follow path & update cost array
      switch (mypath.numBends) {
        case 0:
          if (s_y == e_y){ // Horizontal path
            horizontalCost(B, locks, s_y, s_x, e_x, dim_y, j);
            incrCell(B, locks, e_x, e_y, dim_y, j);
            break;
          }
          if (s_x == e_x){            // Vertical path
            verticalCost(B, locks, e_x, s_y, e_y, dim_y, j);
            incrCell(B,locks, e_x, e_y, dim_y, j);
            break;
          }
        case 1:
          b1_x = mypath.bends0;
          b1_y = mypath.bends1; // Get bend coordinate
          if (s_y == b1_y) // Before bend is horizontal
          {
            horizontalCost(B, locks, s_y, s_x, b1_x, dim_y, j);
            // After bend must be vertical
            verticalCost(B, locks, e_x, b1_y, e_y, dim_y, j);
            incrCell(B, locks, e_x, e_y, dim_y, j);
            break;
          }
          if (s_x == b1_x)           // Before bend is vertical
          {
            verticalCost(B, locks, s_x, s_y, b1_y, dim_y, j);
            // After bend must be horizontal
            horizontalCost(B, locks, e_y, b1_x, e_x, dim_y, j);
            incrCell(B, locks, e_x, e_y, dim_y, j);
            break;
          }
        case 2:
          b1_x = mypath.bends0; // Get both bend coordinates
          b1_y = mypath.bends1;
          b2_x = mypath.bends2;
          b2_y = mypath.bends3;
          // Exam first bend
          if (s_y == b1_y) // Before bend is horizontal
          {
            horizontalCost(B,locks, s_y, s_x, b1_x, dim_y, j);
            verticalCost(B,locks,  b1_x, b1_y, b2_y, dim_y, j);//after bend is vertical
            horizontalCost(B,locks, e_y, b2_x, e_x,dim_y, j);
            incrCell(B,locks, e_x, e_y, dim_y, j);
            break;
          }
          if (s_x == b1_x) // Before bend is vertical
          {
            verticalCost(B,locks, s_x, s_y, b1_y, dim_y, j);
            horizontalCost(B, b1_y, b1_x, b2_x, dim_y, j);//after bend is horizontal
            verticalCost(B, locks, b2_x, b2_y, e_y, dim_y, j);
            incrCell(B, locks, e_x, e_y, dim_y, j);
            break;
          }
      }
    } // implicit barrier //
  	printf(">>> Done PRAGMA \n");
    ///////////////////////////////////////////////////////
    */

    #pragma omp parallel for default(shared)\
      private (i, j) shared(locks) schedule(dynamic)
    // free the allocated wire and DS
    for(i = 0; i < dim_y; i++){
      for(j = 0; j < dim_x; j++){
        //omp_destroy_lock(&(locks[i*dim_y + j]));
      }
    }
    free(locks);
  }
  /* #################### END PRAGMA ################### */

  printf(" >>>>>>>>> DONE <<<<<<<<<<<<<<<<\n");
  compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
  printf("Computation Time: %lf.\n", compute_time);
  /*
  // update board statistic ////////////////
  updateBoard(costs);
  /////////////////////////////
  // Write wires and costs to files /
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
  // wrting to Cost
  for(int row = 0 ; row < dim_y; row++){
    for(int col = 0; col < dim_x; col++){
      fprintf(outputCost, "%d ", B[row*dim_y + col].val);
    }
    fprintf(outputCost, "\n");
  }
  // wrting to wire
  fprintf(outputWire, "%d\n", num_of_wires);
  for(int w_count = 0; w_count < num_of_wires; w_count++){
    fprintf(outputWire, "%d %d ", wires[w_count].bound0,
                       wires[w_count].bound1);
    if(wires[w_count].numBends == 2){
      fprintf(outputWire, "%d %d ", wires[w_count].bends0,
                       wires[w_count].bends1);
      fprintf(outputWire, "%d %d ", wires[w_count].bends2,
                       wires[w_count].bends3);
    }
    if(wires[w_count].numBends == 1){
      fprintf(outputWire, "%d %d ", wires[w_count].bends0,
                       wires[w_count].bends1);
    }
    fprintf(outputWire, "%d %d\n", wires[w_count].bound2,
                       wires[w_count].bound3);
  }
  fclose(outputCost);
  fclose(outputWire);
  free(wires);
  free(B);
  free(costs);
  */
  free(rand_wire);
  free(rand_iter);
  return 0;
}
