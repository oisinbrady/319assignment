#include <stdio.h>  /* needed for input, output, ... */
#include <stdlib.h> /* needed for EXIT_SUCCESS, ... */
#include <ctype.h>  /* needed for isdigit() */
#include <string.h> /* needed for memset() */
#include <glpk.h>   /* the linear programming toolkit */
#include <stddef.h>
#include <math.h>

/* global variables -- YOU SHOULD NOT CHANGE THIS! */
/* (You are allowed to add
 * your own if you want.) */
int numRows, numCols; /* number of rows and columns of the puzzle */
int *input;           /* array holding the input */
int *solution;        /* array holding the solution */
int debug;            /* flag for debug mode; 1 means debug mode, 0 means debug off */

/* prototypes of functions -- YOU SHOULD NOT CHANGE THIS! */
/* (Feel free to add your own as you like.) */
int readInput(char *filename); /* reads puzzle from file */
/* readInput creates and fills the global variables as needed */
/* it returns 0 if all is okay and 1 otherwise */
int computeSolution(void); /* computes a solution if possible */
/* the return value is 1 if there is a solution and 0 otherwise */

/* This is the function that actually solves the problem. */
/* It is currently basically empty and not functional. */
/* Your own implementation needs to go in here. */

bool different_color(int node_u, int node_w);

bool different_color(int node_u, int node_w){
  return (node_u != node_w && node_u != 0 && node_w != 0);
}

int count_edges(int t_nodes, int *adjaceny_matrix);

int count_edges(int t_nodes, int *adjaceny_matrix){
  // Traverse the adjaceny matrix and count the number of edges
  int edges = 0;
  for (int i = 0; i < t_nodes; i++)
  {
    for (int j = 0; j < t_nodes; j++)
    {
      printf("%d ", *((adjaceny_matrix+i*t_nodes) + j)); //DEBUG
      if (*((adjaceny_matrix+i*t_nodes) + j) == 1)
      {
        edges++;
      }
    }
    printf("\n");
  }
  return edges;
}

void print_edge_flows(int *solution, int *adjaceny_matrix, glp_prob *lp);

void print_edge_flows(int *solution, int *adjaceny_matrix, glp_prob *lp){
  // for debugging
  double flow;
  const int NODES = (numCols*numRows);
  for (int i = 0, edges = 0; i < NODES; i++)
  {
    for (int j = 0; j < NODES; j++)
    {
      if (*(adjaceny_matrix + i*NODES + j) > 0.0)
      {
        edges++;                            /* compute the number of the edge */
        flow = glp_get_col_prim(lp, edges); /* get the flow value */
        if (flow > 0.0)
        {
          fprintf(stdout, "flow %f on edge %d->%d\n", flow, i, j);
        }
      }
    }
  }
  printf("\n");
  glp_write_lp(lp, NULL, "output.txt"); // DEBUG: print the LP problem
}

int *determine_adjacent_nodes(int node, int *array);

int *determine_adjacent_nodes(int node, int *array)
{
  // find the nodes x and y coordinates
  // printf("\nNode offset:%i", node);
  int x = node % numCols;
  int y = node / numCols;
  // printf("(%i,%i)\n",x,y );

  if (y > 0)
  {
    int offset = int(y - 1) * int(numCols) + int(x);
    array[0] = offset; // top node offset for the 1d array
  }
  if (x > 0)
  {
    // left
    int offset = int(y) * int(numCols) + int(x - 1);
    array[1] = offset;
  }
  if (x < numCols - 1)
  {
    // right
    int offset = int(y) * int(numCols) + int(x + 1);
    array[2] = offset;
  }
  if (y < numRows - 1)
  {
    // bottom
    int offset = int(y + 1) * int(numCols) + int(x);
    array[3] = offset;
  }
  return array;
}

// used for the a linked list in printing the solution graph
struct Node {
  int index;
  int next_index;
};

struct Color {
  // the respective source and sink locations of each color
  int color;
  int source_location;
  int sink_location;
};

int computeSolution(void)
{
  const int NODES = numRows * numCols;
  glp_prob *lp;
  lp = glp_create_prob();
  glp_set_prob_name(lp, "MEDP"); // maximum edge disjoint paths - all or nothing


  // TODO convert to malloced 2d array so i can pass easily into functions
  int *adjaceny_matrix = (int *)malloc(NODES * NODES * sizeof(int));

  for (int i = 0; i < NODES; i++)
  {
    for (int j = 0; j < NODES; j++)
    {
      *(adjaceny_matrix + i*NODES + j) = 0;
    }
  }

  // map the input graph onto a 1d array for quick reference
  // !READ is this necessary?
  int input_1d[NODES];
  for (int i = 0, next_node = 0; i < numRows; i++)
  {
    for (int j = 0; j < numCols; j++)
    {
      input_1d[next_node] = input[i * numCols + j];
      next_node++;
    }
  }
  // count source & sink nodes
  int total_sources_and_sinks = 0;
  for (int i = 0; i < NODES; i++)
  {
    if (input_1d[i] != 0)
    {
      total_sources_and_sinks++;
    }
  }
  if (total_sources_and_sinks == 0)
  { // if there is no source or sink node
    fprintf(stdout, "No source/sink node(s) found");
    return 0;
  }

  // find source node(s)
  int source_sink_pair[total_sources_and_sinks/2][2];
  int pair_count = 0;

  const int PAIRS = total_sources_and_sinks/2;
  struct Color st_pairs[PAIRS];
  for (int i = 0, pair_counter = 0; i < NODES; i++)
  {
    for (int j = i+1; j < NODES; j++)
    {
      // if we have seen this color before
      if (input_1d[i] != 0 && input_1d[j] != 0)
      {
        if (input_1d[i] == input_1d[j] && input_1d[i] != 0){
        struct Color c;
        c.color = input_1d[i];
        c.source_location = i;
        c.sink_location = j;
        st_pairs[pair_counter] = c;
        pair_counter++;
        }
      }
    }
  }

  // make the adjaceny matrix
  // TODO rewrite
  for (int i = 0; i < NODES; i++)
  {
    if (input_1d[i] != 0)
    { // if we have a source or sink node
      // mark the node as a source if it has not been seen before
      // determine if node is source or sink
      bool is_source = false;
      struct Color matching_pair;
      for (int p = 0; p < PAIRS; p++)
      {
        if (st_pairs[p].color == input_1d[i])
        {
          matching_pair.color = st_pairs[p].color;
          matching_pair.source_location = st_pairs[p].source_location;
          matching_pair.sink_location = st_pairs[p].sink_location;
          if (st_pairs[p].source_location == i )
          {
            is_source = true;
          }
          // otherwise we are dealing with a sink node
        }
      }
      if (is_source)
      {
        // make the edges for the source/sink node
        // adjacent node array of the top, left, right, and bottom nodes
        // if the node does not have a certain directionally adjacent node
        // it will remain -1 - indicating that there is no edge in this case
        int source = i;
        int an[4] = {-1, -1, -1, -1};
        int *adjacent_nodes = determine_adjacent_nodes(source, an);
        // make all edges from source to its adjacent nodes
        const int MAX_ADJACENT = 4;
        for (int node = 0; node < MAX_ADJACENT; node++)
        {
          if (adjacent_nodes[node] != -1)
          {
            // if the edge does not go into another source/sink color
            int node_index = adjacent_nodes[node];
            int node_value = input_1d[adjacent_nodes[node]];
            if (node_index == matching_pair.sink_location || node_value == 0)
            {
              // make all edges for the source
              *(adjaceny_matrix + source*NODES + node_index) = 1;
            }
          }
        }
      }
      else
      { //otherwise construct edges going into the sink
        int sink = i;
        int an[4] = {-1, -1, -1, -1};
        int *adjacent_nodes = determine_adjacent_nodes(sink, an);
        // make all edges adjacent to the sink going into the sink
        const int MAX_ADJACENT = 4;
        for (int node = 0; node < MAX_ADJACENT; node++)
        {
          if (adjacent_nodes[node] != -1)
          {
            int node_index = adjacent_nodes[node];
            int node_value = input_1d[adjacent_nodes[node]];
            // if the edge does not go into a different colored source
            if (node_index == matching_pair.source_location || node_value == 0){
              // make all edges for the sink
              *(adjaceny_matrix + node_index*NODES + sink) = 1;
            }
          }
        }
      }
    }
    else
    { // for all other non-source/sink edges
      int an[4] = {-1, -1, -1, -1};
      int *adjacent_nodes = determine_adjacent_nodes(i, an);
      for (int ad_node = 0; ad_node < 4; ad_node++)
      {
        // if the node has an adjacent node in this direction
        if (adjacent_nodes[ad_node] != -1)
        {
          // if the adjacent node is not a source/sink
          if (input_1d[adjacent_nodes[ad_node]] == 0)
          {
            // make edges for all e in E not {s,t}
            *(adjaceny_matrix + i*NODES + adjacent_nodes[ad_node]) = 1;
            // a path in the converse direction is also possible
            *(adjaceny_matrix + adjacent_nodes[ad_node]*NODES + i) = 1;
          }
        }
      }
    }
  }

  int edges = count_edges(NODES, (int *)adjaceny_matrix);

  // set all node capacity constraints
  glp_add_cols(lp, edges);
  for (int i = 0, ad_edges = 0; i < NODES; i++)
  {
    for (int j = 0; j < NODES; j++)
    {
      if (*(adjaceny_matrix + i*NODES + j) > 0.0)
      {                                                                     // if we have an edge
        ad_edges++;                                                         /* compute the number of the edge */
        glp_set_col_bnds(lp, ad_edges, GLP_DB, 0.0, *(adjaceny_matrix + i*NODES + j)); /* 0<=e<=capacity*/
      }
    }
  }

  /* set objective function to max source edges */
  glp_set_obj_dir(lp, GLP_MAX); /* set maximisation as objective */
  for (int i = 0, ad_edges = 0; i < NODES; i++)
  {
    for (int j = 0; j < NODES; j++)
    {
      // if we come across any edge in the graph
      if (*(adjaceny_matrix + i*NODES + j) > 0.0 )
      {
        ad_edges++; /* compute the number of the edge */
        if (input_1d[i] > 0 || input_1d[j] > 0)  // if edge has source or sink
        {
          for (int p = 0; p < PAIRS; ++p)
          {
            if (i == st_pairs[p].source_location)
            {
              // add the edge to the max obj f
              printf("Adding e(%i,%i) to obj function\n", i, j);
              glp_set_obj_coef(lp, ad_edges, 1.0);
            }
          }
        }
        else
        {
          glp_set_obj_coef(lp, ad_edges, 0.0); /* other edges not involved */
        }
      }
    }
  }

  int index[NODES];    /* indices to define constraint coefficients */
  double value[NODES]; /* values to define constraint coefficients */
  int connected;             /* number of edges connected to specific node  */
  int i;
  int j;
  /* define constraints */
  glp_add_rows(lp, NODES); /* We have one constraint per node */
  for (int node = 0; node < NODES; node++)
  {
    printf("At node=%i\n", node);
    bool source = false;
    struct Color st_pair;
    /* constraint is sum of incoming - sum of outgoing = 0 */
    for (int p = 0; p < PAIRS; p++)
    {
      if (input_1d[node] == st_pairs[p].color)
      {
        st_pair.color = st_pairs[p].color;
        st_pair.source_location = st_pairs[p].source_location;
        st_pair.sink_location = st_pairs[p].sink_location;
        if (st_pairs[p].source_location == node )
          {
            source = true;
          }
      }
    }
    if (source)  // i.e., if we have the source node
    {
      glp_set_row_bnds(lp, 1 + node, GLP_UP, 0.0, 1.0); // set RHS
      for (i = 0, edges = 0, connected = 0; i < NODES; i++)
      {
        for (j = 0; j < NODES; j++)
        {
          if (*(adjaceny_matrix + i*NODES + j) > 0.0)
          {
            edges++; /* compute the number of the edge */
            if (i == node)
            {              /* edge is outgoing edge */
              connected++; /* count number of connected edges */
              index[connected] = edges;
              value[connected] = 1.0;
            }
            if (j == node)
            {              /* edge is incoming edge */
              connected++; /* count number of connected edges */
              index[connected] = edges;
              value[connected] = 1.0;
            }
          }
        }
      }
      glp_set_mat_row(lp, 1 + node, connected, index, value);
    }

    else if (node != st_pair.sink_location) // if the node is not a sink node
    {
      // if node is not adjacent to the sink node
      int an[4] = {-1, -1, -1, -1};
      int *adjacent_nodes = determine_adjacent_nodes(node, an);
      const int MAX_ADJACENT = 0;
      for (int adj_node = 0; adj_node < MAX_ADJACENT; adj_node++)
      {
        // if the node does not connect to a sink node
        if (input_1d[an[adj_node]] == st_pair.sink_location)
        {
            continue;
        }
      }

      glp_set_row_bnds(lp, 1 + node, GLP_FX, 0.0, 0.0); /* set RHS to 0 */
      for (i = 0, edges = 0, connected = 0; i < NODES; i++)
      {
        for (j = 0; j < NODES; j++)
        {
          if (*(adjaceny_matrix + i*NODES + j) > 0.0)
          {
            edges++; /* compute the number of the edge */
            if (i == node)
            {              /* edge is outgoing edge */
              connected++; /* count number of connected edges */
              index[connected] = edges;
              value[connected] = -1.0;
            }
            if (j == node)
            {              /* edge is incoming edge */
              connected++; /* count number of connected edges */
              index[connected] = edges;
              value[connected] = 1.0;
            }
          }
        }
      }
      /* now define the LHS of the constraint */
      glp_set_mat_row(lp, 1 + node, connected, index, value);
    }
  }


  glp_term_out(0);
  glp_simplex(lp, NULL); // solve LP via simplex algorithm

  printf("Maximal flow is %f\n", glp_get_obj_val(lp));

  printf("\n");
  double flow; /* flow on an edge */

  /* read and print actual flow */
  // TODO build the solution graph
  // IDEA store each flow in a linked list starting from the source node
  // i.e., A solvable graph with two st pairs would make an array of size 2
  // where array_i indicates the source node index and points to the linked list
  // of nodes to follow it

  print_edge_flows(solution, adjaceny_matrix, lp);
  /* house-keeping */
  glp_delete_prob(lp);
  // iff maximal flow == |source nodes| return 1
  return (glp_get_obj_val(lp) == PAIRS);
}

/* YOU SHOULD NOT CHANGE ANYTHING BELOW THIS LINE! */

/* printPuzzle(int *puzzle) prints either an input or a solution */
void printPuzzle(int *puzzle)
{
  int i, j; /* loop variables to go over rows and columns */
  for (i = 0; i < numRows; i++)
  { /* go over rows */
    for (j = 0; j < numCols; j++)
    {                                               /* go over columns */
      fputc('0' + puzzle[i * numCols + j], stdout); /* print the next char */
    }
    fprintf(stdout, "\n"); /* end the current line, start new one */
  }
}

int main(int argc, char **argv)
{
  int i; /* used to run over the command line parameters */

  if (argc < 2)
  { /* no command line parameter given */
    fprintf(stderr, "Usage: %s [file1] [file2] [file3] [...]\n"
                    "Where each [file] is the name of a file with a puzzle.\n",
            argv[0]);
    exit(EXIT_FAILURE);
  }

  if (argv[1][0] == '-' && argv[1][1] == 'd' && argv[1][2] == 0)
  {
    /* If the first parameter is -d we activate debug mode. */
    debug = 1;                                        /* switch debug mode on */
    fprintf(stdout, "DEBUG: Debug mode activated\n"); /* be explicit about it */
  }
  else
  {
    debug = 0; /* switch debug mode off */
  }

  for (i = 1 + debug; i < argc; i++)
  { /* go over remaining command line parameters */
    if (readInput(argv[i]))
    { /* try to read file */
      /* returned with error message */
      fprintf(stderr, "%s: Cannot read puzzle with filename %s. Skipping it.\n",
              argv[0], argv[i]);
    }
    else
    { /* input read successfully */
      fprintf(stdout, "%s: Looking at the following puzzle:\n", argv[i]);
      printPuzzle(input);
      if (computeSolution())
      { /* compute a solution if one exists */
        fprintf(stdout, "%s: Found the following solution:\n", argv[i]);
        printPuzzle(solution);
      }
      else
      {
        fprintf(stdout, "%s: Puzzle has no solution\n", argv[i]);
      }
      /* free memory for next input */
      free(input);
      free(solution);
    }
  }
  return EXIT_SUCCESS;
}

/* checkFile(FILE *fh) performs basic checks and sets numRows/numCols */
/* return value 1 indicates an error; otherwise 0 is returned */
int checkFile(FILE *fh)
{
  char c;
  int rows, cols; /* used to determine number of rows and columns */
  int read;       /* counts number of digits read in the current row */
  int firstRow;   /* indicates if we are reading the very first row */

  firstRow = 1; /* we start in the first row */
  rows = cols = read = 0;
  while (!feof(fh))
  {
    c = fgetc(fh); /* read the next char from the file */
    if (isdigit(c))
    {         /* normal character read */
      read++; /* count the digit we just read */
      if ((!firstRow) && (read > cols))
      {
        if (debug)
        {
          fprintf(stdout, "DEBUG: Row %d too long (%d, %d).\n", rows + 1, read, cols);
        }
        return 1; /* flag error because row is too long */
      }
    }
    else
    {
      if ((c == '\n') || (c == (char)-1))
      { /* end of line read */
        if (read > 0)
        {
          rows++; /* count the completed row if it was not empty */
        }
        if (firstRow)
        {               /* very first row read */
          cols = read;  /* accept number of characters as number of columns */
          firstRow = 0; /* not in the first row anymore after this */
          if (debug)
          {
            fprintf(stdout, "DEBUG: %d columns per row expected\n", cols);
          }
        }
        else
        {
          if ((read > 0) && (read != cols))
          { /* rows too short */
            if (debug)
            {
              fprintf(stdout, "DEBUG: Row %d too short.\n", rows + 1);
            }
            return 1; /* flag error because row is too short */
          }
        }
        read = 0; /* reset number of characters in current row */
      }
      else
      { /* illegal character found */
        if (debug)
        {
          fprintf(stdout, "DEBUG: Illegal character %c found.\n", c);
        }
        return 1; /* stop reading because of the error */
      }
    }
  }
  if (read > 0)
  {
    rows++; /* last row was not ended with newline */
  }
  /* use the determined size and prepare for reading */
  numRows = rows;
  numCols = cols;
  rewind(fh); /* reset to the beginning of the file to read the actual input */
  return 0;   /* signal all went well */
}

/* readInput(*char filename) reads the input and stores it */
/* return value 1 indicates an error; otherwise 0 is returned */
int readInput(char *filename)
{
  int i, j;        /* loop variables to go over the columns and rows of the input */
  int check[10];   /* array to check colours come in pairs */
  int keepReading; /* used to skip over newline */
  char c;          /* next char */
  FILE *fh;

  if ((fh = fopen(filename, "rt")) == NULL)
  {
    return 1;
  }

  /* perform basic checks and determine size of puzzle */
  if (checkFile(fh))
  { /* there was a problem */
    fclose(fh);
    return 1; /* signal error */
  }
  if ((input = (int *)malloc(sizeof(int) * numRows * numCols)) == NULL)
  {
    if (debug)
    {
      fprintf(stdout, "DEBUG: Unable to allocate %ld bytes of memory.\n",
              sizeof(int) * numRows * numCols);
    }
    fclose(fh);
    return 1;
  }
  if ((solution = (int *)malloc(sizeof(int) * numRows * numCols)) == NULL)
  {
    if (debug)
    {
      fprintf(stdout, "DEBUG: Unable to allocate %ld bytes of memory.\n",
              sizeof(int) * numRows * numCols);
    }
    free(input);
    fclose(fh);
    return 1;
  }
  memset(solution, 0, sizeof(int) * numRows * numCols); /*initialise solution empty*/
  if (debug)
  {
    fprintf(stdout, "DEBUG: Will read %dx%d sized puzzle\n", numRows, numCols);
  }
  /* prepare to count different digits */
  for (i = 0; i < 10; i++)
  {
    check[i] = 0;
  }
  /* Size is given in numRows, numCols; now we read */

  for (i = 0; i < numRows; i++)
  { /* go over rows */
    for (j = 0; j < numCols; j++)
    { /* go over columns */
      do
      {
        keepReading = 1; /* prepare to skip over newline */
        c = fgetc(fh);   /* get next digit */
        if (isdigit(c))
        {                                          /* store and count digit */
          input[i * numCols + j] = (int)(c - '0'); /* convert char to int */
          check[input[i * numCols + j]]++;
          keepReading = 0; /* mark digit as read */
        }
      } while (keepReading);
    }
  }
  for (i = 1; i < 10; i++)
  {
    if ((check[i] != 0) && (check[i] != 2))
    {
      if (debug)
      {
        fprintf(stdout, "DEBUG: Colour %d appears %d times.\n", i, check[i]);
        printPuzzle(input);
      }
      free(input);
      free(solution);
      fclose(fh);
      return 1;
    }
  }

  fclose(fh); /* close file after reading the input */
  return 0;   /* signal all went well */
}
