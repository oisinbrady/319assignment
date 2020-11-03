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

struct Color {
  // the respective source and sink locations of each color
  int color;
  int source_location;
  int sink_location;
};

int computeSolution(void)
{
  const int total_nodes = numRows * numCols;
  glp_prob *lp;
  lp = glp_create_prob();
  glp_set_prob_name(lp, "MEDP"); // maximum edge disjoint paths - all or nothing

  int adjaceny_matrix[total_nodes][total_nodes];
  for (int i = 0; i < total_nodes; i++)
  {
    for (int j = 0; j < total_nodes; j++)
    {
      adjaceny_matrix[i][j] = 0;
    }
  }

  // map the input graph onto a 1d array for quick reference
  // !READ is this necessary?
  int input_1d[total_nodes];
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
  for (int i = 0; i < total_nodes; i++)
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
  for (int i = 0; i < total_nodes; i++)
  {
    for (int j = i+1; j < total_nodes; j++)
    {
      // if we have seen this color before
      if (input_1d[i] == input_1d[j] && input_1d[i] != 0){
        struct Color c;
        c.color = input_1d[i];
        c.source_location = i;
        c.sink_location = j;
      }
    }
  }

  // make the adjaceny matrix
  // !TODO rewrite
  for (int i = 0; i < (numRows * numCols); i++)
  {
    if (input_1d[i] != 0)
    { // if we have a source or sink node
      // mark the node as a source if it has not been seen before
      // determine if node is source or sink
      bool is_source = false;
      for (int p = 0; p < PAIRS; p++)
      {
        if (st_pairs[p].color == input_1d[i])
        {
          if (st_pairs[p].source_location == i )
          {
            is_source = true;
          }
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
        for (int i = 0; i < MAX_ADJACENT; i++)
        {
          if (adjacent_nodes[i] != -1)
          {
            // make all edges for the source
            adjaceny_matrix[source][adjacent_nodes[i]] = 1;
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
        for (int i = 0; i < MAX_ADJACENT; i++)
        {
          if (adjacent_nodes[i] != -1)
          {
            // make all edges for the sink
            adjaceny_matrix[adjacent_nodes[i]][sink] = 1;
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
            adjaceny_matrix[i][adjacent_nodes[ad_node]] = 1;
            adjaceny_matrix[adjacent_nodes[ad_node]][i] = 1; // a path in the converse direction is also possible
          }
        }
      }
    }
  }

  int edges = count_edges(total_nodes, (int *)adjaceny_matrix);

  // set all node capacity constraints
  glp_add_cols(lp, edges);
  for (int i = 0, ad_edges = 0; i < total_nodes; i++)
  {
    for (int j = 0; j < total_nodes; j++)
    {
      if (adjaceny_matrix[i][j] > 0.0)
      {                                                                     // if we have an edge
        ad_edges++;                                                         /* compute the number of the edge */
        glp_set_col_bnds(lp, ad_edges, GLP_DB, 0.0, adjaceny_matrix[i][j]); /* 0<=e<=capacity*/
      }
    }
  }

  // !TODO Introduce a super node and super sink
  // super source = sum of all source involved edges
  // super sink = sum of all sink involved edges
  // obj func: max (super source - super sink)

  // loop through adjacency matrix:
  //    if (u or v in e(u,v) in {s,t}
  //      AND (u_color == v_color where u and v in {s,t}):
  //        add this edge to the obj function

  // N.B: I don't believe the sink currently has a constraint, therefore it can "catch flow"
  //      Does this have implications for super sink imlementation? dont think so

  // SHOULD I ADD NODES INBOUND TO SINK TO OBJ FUNCTION?!
  /* set objective function to maximise edges coming out of source */
  glp_set_obj_dir(lp, GLP_MAX); /* set maximisation as objective */
  for (int i = 0, ad_edges = 0; i < total_nodes; i++)
  {
    for (int j = 0; j < total_nodes; j++)
    {
      // if we come across any edge in the graph
      // !TODO also if edge e(s,t) where s_color != t_color

      if (adjaceny_matrix[i][j] > 0.0 )
      {
        ad_edges++; /* compute the number of the edge */
        if (input_1d[i] > 0 || input_1d[j] > 0)  // if edge has source or sink
        {
          // printf("here\n");
          // printf("i=%i\n", i);
          // printf("j=%i\n", j);
          // printf("input_1d[i]%i\n", input_1d[i]);
          // printf("input_1d[j]%i\n", input_1d[j]);
          // if (input_1d[i] == input_1d[j] || input_1d[i] == 0 || input_1d[j] ==0){
            printf("Adding edge(%i,%i) to obj function\n", i, j);
            glp_set_obj_coef(lp, ad_edges, 1.0);  // add the edge to the max obj f
          // }
        }
        else
        {
          glp_set_obj_coef(lp, ad_edges, 0.0); /* other edges not involved */
        }
      }
    }
  }

  int index[total_nodes];    /* indices to define constraint coefficients */
  double value[total_nodes]; /* values to define constraint coefficients */
  int connected;             /* number of edges connected to specific node  */
  int i;
  int j;
  /* define constraints */
  glp_add_rows(lp, total_nodes); /* We have one constraint per node */
  for (int node = 0; node <= total_nodes - 1; node++)
  {
    bool source = false;
    int source_index;
    /* constraint is sum of incoming - sum of outgoing = 0 */
    for (int i = 0; i < pair_count; i++)
    {
      if (node == source_sink_pair[i][0]){
        source_index = i;
        source = true;
      }
    }
    if (source)  // i.e., if we have the source node
    {
      glp_set_row_bnds(lp, 1 + node, GLP_UP, 0.0, 1.0); // set RHS
      for (i = 0, edges = 0, connected = 0; i < total_nodes; i++)
      {
        for (j = 0; j < total_nodes; j++)
        {
          if (adjaceny_matrix[i][j] > 0.0)
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
    else if (node != source_sink_pair[source_index][1])  // if node is not a sink node
    {
      glp_set_row_bnds(lp, 1 + node, GLP_FX, 0.0, 0.0); /* set RHS to 0 */
      for (i = 0, edges = 0, connected = 0; i < total_nodes; i++)
      {
        for (j = 0; j < total_nodes; j++)
        {
          if (adjaceny_matrix[i][j] > 0.0)
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
  for (int i = 0, edges = 0; i < total_nodes; i++)
  {
    for (int j = 0; j < total_nodes; j++)
    {
      if (adjaceny_matrix[i][j] > 0.0)
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
  /* house-keeping */
  glp_delete_prob(lp);

  // TODO if maximal flow == |source nodes| return 1
  return 0; /* this is not true for every puzzle, of course */
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
