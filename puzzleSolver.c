#include <stdio.h> /* needed for input, output, ... */
#include <stdlib.h> /* needed for EXIT_SUCCESS, ... */
#include <ctype.h> /* needed for isdigit() */
#include <string.h> /* needed for memset() */
#include <glpk.h> /* the linear programming toolkit */
#include <stddef.h>
#include <math.h>

/* global variables -- YOU SHOULD NOT CHANGE THIS! */
/* (You are allowed to add your own if you want.) */
int numRows, numCols; /* number of rows and columns of the puzzle */
int *input; /* array holding the input */
int *solution; /* array holding the solution */
int debug; /* flag for debug mode; 1 means debug mode, 0 means debug off */

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

int * determine_adjacent_nodes(int node, int *array);

struct Edge{  // edge for the directed graph
      int node_u;
      int node_v;
      int id[2]; // for debugging/printing, [node number][node number]
};
struct Edge construct_edge(int node_u, int node_v);
// construct edges for two given nodes
struct Edge construct_edge(int node_u, int node_v) {
   struct Edge e;
   e.node_u = node_u;
   e.node_v = node_v;
   e.id[0] = node_u;
   e.id[1] = node_v;
   return e;
}

int * determine_adjacent_nodes(int node, int *array){
    // find the nodes x and y coordinates
    // printf("\nNode offset:%i", node);
    int x = node % numCols;
    int y = node / numCols;
    // printf("(%i,%i)\n",x,y );

    if ( y > 0 ) {
      int offset = int(y-1)*int(numCols)+int(x);
      array[0] = offset;  // top node offset for the 1d array
    }
    if ( x > 0 ) {
        // left
        int offset = int(y)*int(numCols)+int(x-1);
        array[1] = offset;
    }
    if ( x < numCols - 1 ) {
        // right
        int offset = int(y)*int(numCols)+int(x+1);
        array[2] = offset;
    }
    if ( y < numRows - 1 ) {
        // bottom
        int offset = int(y+1)*int(numCols)+int(x);
        array[3] = offset;
    }
    return array;
}


int computeSolution(void) {
    glp_prob *lp;
    lp = glp_create_prob();
    int ia[1+1000], ja[1+1000];  // !READ
    double ar[1+1000];  // !READ
    glp_set_prob_name(lp, "MEDP");  // maximum edge disjoint paths - all or nothing
    glp_set_obj_dir(lp, GLP_MAX);  // this is a maximization problem

    /* !TODO
    construct a super-source and super-sink
    if more than one source and sink
    */
    const int total_nodes = numRows*numCols;

    // map the input graph onto a 1d array for quick reference
    // !READ is this necessary?
    int input_1d[total_nodes];
    int next_node = 0;
    for ( int i=0; i<numRows; i++ ) {
      for ( int j=0; j<numCols; j++ ) {
          input_1d[next_node] = input[i*numCols+j];
          next_node++;
      }
    }
    // find source & sink nodes
    int total_sources_and_sinks = 0;
    for ( int i=0; i<total_nodes; i++) {
        if (input_1d[i] != 0) {
            total_sources_and_sinks++;
        }
    }

    int const total_src_sink_pairs = total_sources_and_sinks / 2;
    if (total_src_sink_pairs == 0) {  // if there is no source or sink node
      fprintf(stdout, "No source/sink node(s) found");
      return 0;
    }
    int source_sink_nodes[total_src_sink_pairs][3];  // node numbers of source-sink pairs

    int pair_index = 0;
    for ( int i=0; i < total_nodes; i++ ) {
        if (input_1d[i] != 0) {
          int color = input_1d[i];
          for ( int j = i + 1; j < total_nodes; j++ ) {  // find a source/sink node
              if (input_1d[j] == color) {  // find its pair
                  source_sink_nodes[pair_index][0] = i;  // source node
                  source_sink_nodes[pair_index][1] = j;  // sink node
                  source_sink_nodes[pair_index][2] = color;  // the pair's color
                  pair_index++;
                  break;
              }
          }
        }
    }

    /*
      find all possible edges for this type of graph and
      make a new variable for each
    */
    int total_edges = (2*numRows + 2*(numCols-2)) + (3*((numRows*numCols)-(2*numRows + 2*(numCols-2))));
    Edge* edge_list = (Edge*)malloc(sizeof(Edge) * 40);
    printf("%d\n", total_edges );
    int edge_counter = 0;


    int ca_index = 0;
    int color_array[total_nodes];  // store the encountered colors
    for (int i = 0; i < total_nodes; i++) {
      color_array[i] = 0;
    }
    bool source = true;

    for ( int i=0; i<(numRows*numCols); i++ ) {
        if ( input_1d[i]!=0 ) {  // if we have a source or sink node
            // mark the node as a source if it has not been seen before
            for ( int c = 0; c < total_nodes; c++ ) {
              if ( color_array[c] == input_1d[i] ){
                source = false;
              }
            }
            if (source) {
                color_array[ca_index] = input_1d[i];
                ca_index++;
                // make the edges for the source/sink node
                int source = i;
                // adjacent node array of the top, left, right, and bottom nodes
                // if the node does not have a certain directionally adjacent node
                // it will remain -1 - indicating that there is no edge in this case
                int an[4] = {-1,-1,-1,-1};
                int *adjacent_nodes = determine_adjacent_nodes(source, an);

                // for ( int n=0; n<4; n++){
                //   printf("%i\n", adjacent_nodes[n]);
                // }
                // printf("\n");

                // make all edges from source to its adjacent nodes
                const int MAX_ADJACENT = 4;
                for ( int i=0; i<MAX_ADJACENT; i++) {
                  if (adjacent_nodes[i] != -1){
                    edge_list[edge_counter] = construct_edge(source, adjacent_nodes[i]);
                    edge_counter++;
                  }
                }
            }
            else {  //otherwise construct edges going into the sink
                int sink = i;
                int an[4] = {-1,-1,-1,-1};
                int *adjacent_nodes = determine_adjacent_nodes(sink, an);

                // for ( int n=0; n<4; n++){
                //   printf("%i\n", adjacent_nodes[n]);
                // }
                // printf("\n");

                // make all edges adjacent to the sink going into the sink
                const int MAX_ADJACENT = 4;
                for ( int i=0; i<MAX_ADJACENT; i++) {
                  if (adjacent_nodes[i] != -1){
                    edge_list[edge_counter] = construct_edge(adjacent_nodes[i], sink);
                    edge_counter++;
                  }
                }
            }
        }
        else {  // for all other non-source/sink edges

            // BUG: still producing edges into source node and not all adjacent nodes made
            // N.B: check if there are edges going from sink to adjacent node
            int an[4] = {-1,-1,-1,-1};
            int *adjacent_nodes = determine_adjacent_nodes(i, an);
            for (int an = 0; an < 4; an++) {
              for (int j = 0; j < total_src_sink_pairs; j++) {
                if (source_sink_nodes[j][3] != adjacent_nodes[an]) {
                    if (adjacent_nodes[an] != -1) {
                        edge_list[edge_counter] = construct_edge(i, adjacent_nodes[an]);
                        edge_counter++;
                    }
                }
              }
            }

            printf("Adjacent nodes to node offset, %i\n", i);
            for ( int n=0; n<4; n++){
                  printf("%i\n", adjacent_nodes[n]);
                }
            printf("\n");

            // const int MAX_ADJACENT = 4;

            // for ( int i=0; i<MAX_ADJACENT; i++) {
            //   if (adjacent_nodes[i] != -1) {
            //     edge_list[edge_counter] = construct_edge(i, adjacent_nodes[i]);
            //     edge_counter++;
            //   }
            // }
        }
    }


    for ( int i=0; i < total_edges; i++ ){
      printf("\n%i:", i);
      printf("\nid:[%i,%i]: %i,%i\n", edge_list[i].id[0], edge_list[i].id[1], edge_list[i].node_u, edge_list[i].node_v );
    }





    return 0;





    // set all node capacity constraints
    glp_add_cols(lp, 24);
    for ( int i = 1; i <= 24; i++){
      glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
    }

    // set all node flow constraints
    fprintf(stdout, "Total nodes = %d\n", total_nodes);
    glp_add_rows(lp, total_nodes);

    for ( int i = 1; i <= total_nodes; i++){
      if ( i == 4 ){
        printf("\nFound source node");

        glp_set_row_name(lp, i, "source");
        glp_set_row_bnds(lp, i, GLP_LO, 0.0, 1.0);

        // !TODO will need to contruct super-source/sink
        glp_set_col_name(lp, i, "source");  // objective function source node value
        glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
        glp_set_obj_coef(lp, i, 1.0);
      }
      else if( i == 13 ){
        printf("\nFound sink node");

        glp_set_row_name(lp, i, "sink");
        glp_set_row_bnds(lp, i, GLP_LO, 0.0, 1.0);

        glp_set_col_name(lp, i, "sink");  // objective function sink node value
        glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
        glp_set_obj_coef(lp, i, 1.0);
      }
      else{
        printf("\nelse");
        glp_set_row_bnds(lp, i, GLP_UP, 0.0, 1.0);  // set flow in = flow out (constraint)
      }
    }

    /*
    !TODO: find out how to change '+ source - sink <= 1' to '= 0'

    co-efficient matrix dummy data for hacking

    ia[i]: each constraint (flow law, in=out/source=sink),
    ja[i]: nodes belonging that constraint

    E.g:
    if node 4 = source node and 13 = sink node
    then:
    (constraint number)|(add edge to this constraint)   |c = the edge's coefficient
                       |(ja[x] where x is order of edge)|
    ---------------------------------------------------------------------------------
    ia[1] = 1,         | ja[1] = 4,                     | ar[1] = c
    */
    ia[1] = 1, ja[2] = 4, ar[1] =  -1.0; /* a[1,1] =  -1 */
    ia[2] = 1, ja[1] = 13, ar[2] =  1.0; /* a[1,2] =  1 */
    // ia[3] = 1, ja[3] = 3, ar[3] =  1.0; /* a[1,3] =  1 */
    // ia[4] = 2, ja[4] = 1, ar[4] =  10.0; /* a[2,1] = 10 */
    // ia[5] = 3, ja[5] = 1, ar[5] =  2.0; /* a[3,1] =  2 */
    // ia[6] = 2, ja[6] = 2, ar[6] =  1.0; /* a[2,2] =  4 */
    // ia[7] = 3, ja[7] = 2, ar[7] =  1.0; /* a[3,2] =  2 */
    // ia[8] = 2, ja[8] = 3, ar[8] =  1.0; /* a[2,3] =  5 */
    // ia[9] = 3, ja[9] = 3, ar[9] =  1.0; /* a[3,3] =  6 */

    glp_load_matrix(lp, 2, ia, ja, ar);  // coefficents

    glp_simplex(lp, NULL);  // solve LP via simplex algorithm
    /*
    z = objective function value
    I.e., The maximised value of source + sink
    Since f(source) = f(sink),
    we can 'reduce' the objective function to:
    Max(f(source))
    */
    double z = glp_get_obj_val(lp);
    for ( int i=1; i <= 24; i++){
      printf("\n%g" , glp_get_col_prim(lp, i));
    }
    printf("\n");

    printf("\n%g\n", z);

    glp_write_lp(lp, NULL, "/output_files/output_print");  // DEBUG: print the LP problem

    /* house-keeping */
    glp_delete_prob(lp);
    free(edge_list);
  	return 0; /* this is not true for every puzzle, of course */
}

/* YOU SHOULD NOT CHANGE ANYTHING BELOW THIS LINE! */

/* printPuzzle(int *puzzle) prints either an input or a solution */
void printPuzzle(int *puzzle) {
    int i, j; /* loop variables to go over rows and columns */
    for ( i=0; i<numRows; i++ ) { /* go over rows */
      for ( j=0; j<numCols; j++ ) { /* go over columns */
        fputc('0'+puzzle[i*numCols+j], stdout); /* print the next char */
      }
      fprintf(stdout, "\n"); /* end the current line, start new one */
    }
}

int main(int argc, char **argv) {
	int i; /* used to run over the command line parameters */

	if ( argc<2 ) { /* no command line parameter given */
		fprintf(stderr, "Usage: %s [file1] [file2] [file3] [...]\n"
      "Where each [file] is the name of a file with a puzzle.\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	if ( argv[1][0]=='-' && argv[1][1]=='d' && argv[1][2]==0 ) {
    /* If the first parameter is -d we activate debug mode. */
		debug=1; /* switch debug mode on */
		fprintf(stdout, "DEBUG: Debug mode activated\n"); /* be explicit about it */
	} else {
		debug=0; /* switch debug mode off */
	}

  for ( i=1+debug; i<argc; i++ ) { /* go over remaining command line parameters */
    if ( readInput(argv[i]) ) { /* try to read file */
      /* returned with error message */
      fprintf(stderr, "%s: Cannot read puzzle with filename %s. Skipping it.\n",
        argv[0], argv[i]);
    } else { /* input read successfully */
			fprintf(stdout, "%s: Looking at the following puzzle:\n", argv[i]);
      printPuzzle(input);
      if ( computeSolution() ) { /* compute a solution if one exists */
        fprintf(stdout, "%s: Found the following solution:\n", argv[i]);
        printPuzzle(solution);
      } else {
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
int checkFile(FILE *fh) {
  char c;
  int rows, cols; /* used to determine number of rows and columns */
  int read; /* counts number of digits read in the current row */
  int firstRow; /* indicates if we are reading the very first row */

  firstRow=1; /* we start in the first row */
  rows=cols=read=0;
  while ( !feof(fh) ) {
    c = fgetc(fh); /* read the next char from the file */
    if ( isdigit(c) ) { /* normal character read */
      read++; /* count the digit we just read */
      if ( ( !firstRow) && ( read>cols) ) {
        if ( debug ) {
          fprintf(stdout, "DEBUG: Row %d too long (%d, %d).\n", rows+1, read, cols);
        }
        return 1; /* flag error because row is too long */
      }
    } else {
      if ( ( c=='\n' ) || ( c==(char)-1 ) ) { /* end of line read */
        if ( read>0 ) {
          rows++; /* count the completed row if it was not empty */
        }
        if ( firstRow ) { /* very first row read */
          cols=read; /* accept number of characters as number of columns */
          firstRow=0; /* not in the first row anymore after this */
          if ( debug ) {
            fprintf(stdout, "DEBUG: %d columns per row expected\n", cols);
          }
        } else {
          if ( ( read>0 ) && ( read!=cols ) ) { /* rows too short */
            if ( debug ) {
              fprintf(stdout, "DEBUG: Row %d too short.\n", rows+1);
            }
            return 1; /* flag error because row is too short */
          }
        }
        read=0; /* reset number of characters in current row */
      } else { /* illegal character found */
        if ( debug ) {
          fprintf(stdout, "DEBUG: Illegal character %c found.\n", c);
        }
        return 1; /* stop reading because of the error */
      }
    }
  }
  if ( read>0 ) {
    rows++; /* last row was not ended with newline */
  }
  /* use the determined size and prepare for reading */
  numRows = rows;
  numCols = cols;
  rewind(fh); /* reset to the beginning of the file to read the actual input */
  return 0; /* signal all went well */
}

/* readInput(*char filename) reads the input and stores it */
/* return value 1 indicates an error; otherwise 0 is returned */
int readInput(char *filename) {
  int i, j; /* loop variables to go over the columns and rows of the input */
  int check[10]; /* array to check colours come in pairs */
  int keepReading; /* used to skip over newline */
  char c; /* next char */
  FILE *fh;

  if ( ( fh = fopen(filename, "rt") ) == NULL ) {
    return 1;
  }

  /* perform basic checks and determine size of puzzle */
  if ( checkFile(fh) ) { /* there was a problem */
    fclose(fh);
    return 1; /* signal error */
  }
  if ( ( input = (int *)malloc(sizeof(int)*numRows*numCols) ) == NULL ) {
    if ( debug ) {
      fprintf(stdout, "DEBUG: Unable to allocate %ld bytes of memory.\n",
        sizeof(int)*numRows*numCols);
    }
    fclose(fh);
    return 1;
  }
  if ( ( solution = (int *)malloc(sizeof(int)*numRows*numCols) ) == NULL ) {
    if ( debug ) {
      fprintf(stdout, "DEBUG: Unable to allocate %ld bytes of memory.\n",
        sizeof(int)*numRows*numCols);
    }
    free(input);
    fclose(fh);
    return 1;
  }
	memset(solution, 0, sizeof(int)*numRows*numCols);/*initialise solution empty*/
  if ( debug ) {
    fprintf(stdout, "DEBUG: Will read %dx%d sized puzzle\n", numRows, numCols);
  }
  /* prepare to count different digits */
  for ( i=0; i<10; i++ ) {
    check[i]=0;
  }
  /* Size is given in numRows, numCols; now we read */

  for ( i=0; i<numRows; i++ ) { /* go over rows */
    for ( j=0; j<numCols; j++ ) { /* go over columns */
      do {
        keepReading=1; /* prepare to skip over newline */
        c = fgetc(fh); /* get next digit */
        if ( isdigit(c) ) { /* store and count digit */
          input[i*numCols+j]=(int)(c-'0'); /* convert char to int */
          check[input[i*numCols+j]]++;
          keepReading=0; /* mark digit as read */
        }
      } while ( keepReading );
    }
  }
  for ( i=1; i<10; i++ ) {
    if ( ( check[i]!=0 ) && ( check[i]!=2 ) ) {
      if ( debug ) {
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
  return 0; /* signal all went well */
}
