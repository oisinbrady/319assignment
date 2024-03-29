#include <stdio.h>  /* needed for input, output, ... */
#include <stdlib.h> /* needed for EXIT_SUCCESS, ... */
#include <ctype.h>  /* needed for isdigit() */
#include <string.h> /* needed for memset() */
#include <glpk.h>   /* the linear programming toolkit */
#include <stdbool.h>

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

void print_edge_flows(glp_prob *lp, const int EDGES);

void print_edge_flows(glp_prob *lp, const int EDGES) {
    /* for debugging */
    double flow;
    for (int i = 1; i <= EDGES; i++) {
        const char *edge_name = glp_get_col_name(lp, i);
        flow = glp_get_col_prim(lp, i); /* get the flow value */
        if (flow > 0.0) {
            fprintf(stdout, "flow(%s)=%f\n", edge_name, flow);
        }
    }
    printf("\n");
}

int *determine_adjacent_nodes(int node, int *array);

int *determine_adjacent_nodes(int node, int *array) {
    /* Identify node's row & column index */
    int x = node % numCols;
    int y = node / numCols;
    /* Get all adjacent nodes */
    if (y > 0) {  /* top node */
        array[0] = (y - 1) * numCols + x;
    }
    if (x > 0) {  /* left node */
        array[1] = y * numCols + x - 1;
    }
    if (x < numCols - 1) {  /* right node */
        array[2] = y * numCols + x + 1;
    }
    if (y < numRows - 1) {  /* bottom node */
        array[3] = (y + 1) * numCols + x;
    }
    return array;
}

struct Color {
    /* Information for a Color - (i.e., a source/sink pair of nodes) */
    int color;
    int source_location;
    int sink_location;
    /* the node that the source/sink node connects to to form an edge */
    int source_edges[4];
    int sink_edges[4];
};

struct Node {
    int node_id;
    bool source;
    bool sink;
    int *incoming_edges;
    int *outgoing_edges;
};

void build_solution(glp_prob *lp, const int EDGES, int *solution);

void build_solution(glp_prob *lp, const int EDGES, int *solution) {
    for (int i = 1; i <= EDGES; i++) {
        const char *edge_name = glp_get_col_name(lp, i);
        if (glp_get_col_prim(lp, i) == 1.0) {
            int v1, v2, color;
            char regex[] = "%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d";
            sscanf(edge_name, regex, &v1, &v2, &color);
            solution[v1 - 1] = color;
            solution[v2 - 1] = color;
        }
    }
}
void find_st_pairs(const int NODES, const int PAIRS, struct Color *st_pairs, struct Node *node_info, int *input_1d);

void find_st_pairs(const int NODES, const int PAIRS, struct Color *st_pairs, struct Node *node_info, int *input_1d) {
    for (int i = 0, pair_counter = 0; i < NODES; i++) {
        for (int j = i + 1; j < NODES; j++) {
            // if we have seen this color before
            if (input_1d[i] != 0 && input_1d[j] != 0) {
                if (input_1d[i] == input_1d[j] && input_1d[i] != 0) {
                    struct Color c;
                    c.color = input_1d[i];
                    c.source_location = i;
                    c.sink_location = j;
                    st_pairs[pair_counter] = c;
                    pair_counter++;
                    node_info[i].source = true;
                    node_info[j].sink = true;
                }
            }
        }
    }
    // initialise all source's outbound edges
    for (int p = 0; p < PAIRS; p++) {
        for (int i = 0; i < 4; i++) {
            st_pairs[p].source_edges[i] = -1;
        }
    }
}

int get_color(int node_u, int node_v, int input[]);

int get_color(int node_u, int node_v, int input[]) {
    if (input[node_u] > 0) {
        return input[node_u];
    } else {
        return input[node_v];
    }
}

bool is_source(int node_u, struct Color *st_pairs, int PAIRS);

bool is_source(int node_u, struct Color *st_pairs, int PAIRS) {
    for (int p = 0; p < PAIRS; p++) {
        if (st_pairs[p].source_location == node_u) {
            return true;
        }
    }
    return false;
}

bool is_sink(int node_v, struct Color *st_pairs, int PAIRS);

bool is_sink(int node_v, struct Color *st_pairs, int PAIRS) {
    for (int p = 0; p < PAIRS; p++) {
        if (st_pairs[p].sink_location == node_v) {
            return true;
        }
    }
    return false;
}

void lp_make_col(glp_prob *lp, int col_id, int node_u, int node_v, int color);

void lp_make_col(glp_prob *lp, int col_id, int node_u, int node_v, int color) {
    char edge_name[256];
    sprintf(edge_name, "e(v%i,v%i,%i)", node_u + 1, node_v + 1, color);
    glp_set_col_bnds(lp, col_id, GLP_DB, 0.0, 1.0); /* 0<=e<=capacity*/
    glp_set_col_name(lp, col_id, edge_name);
}

int count_edges(const int NODES, const int PAIRS, int *adjacency_matrix , struct Node *node_info);

int count_edges(const int NODES, const int PAIRS, int *adjacency_matrix , struct Node *node_info) {
    /* Traverse the adjacency matrix and count the number of edges */
    int edges = 0;
    for (int i = 0; i < NODES; i++) {
        for (int j = 0; j < NODES; j++) {
            //printf("%d ", *((adjacency_matrix+i*t_nodes) + j)); //DEBUG print matrix
            if (*((adjacency_matrix + i * NODES) + j) == 1) {
                if (!node_info[i].sink && !node_info[i].source) {
                    for (int p = 0; p <= PAIRS; p++) {
                        edges++;
                    }
                } else {
                    edges++;
                }
            }
        }
        // printf("\n");
    }
    return edges;
}

void node_info_incoming_edge(struct Node *node_info, int node_id, int col_id, int NODES);

void node_info_incoming_edge(struct Node *node_info, int node_id, int col_id, int NODES) {
    for (int n = 0; n < NODES; n++) {
        if (node_info[node_id].incoming_edges[n] == -1) {
            node_info[node_id].incoming_edges[n] = col_id;
            break;
        }
    }
}

void node_info_outgoing_edge(struct Node *node_info, int node_id, int col_id, int NODES);

void node_info_outgoing_edge(struct Node *node_info, int node_id, int col_id, int NODES) {
    for (int n = 0; n < NODES; n++) {
        if (node_info[node_id].outgoing_edges[n] == -1) {
            node_info[node_id].outgoing_edges[n] = col_id;
            break;
        }
    }
}

void lp_make_row_flow_cons(glp_prob *lp, struct Node *node_info, int node, int NODES);

void lp_make_row_flow_cons(glp_prob *lp, struct Node *node_info, int node, int NODES) {
    /* Create the flow conservation constraint (3) */
    int index[NODES * NODES];
    double value[NODES * NODES];
    int row_id = 1 + node + NODES;
    int connected = 0;
    int o_edge_id = 0; // outgoing edge id
    while (node_info[node].outgoing_edges[o_edge_id] != -1) {  /* record edge id and set its value to +1.0 */
        // printf("%i\n",node_info[node].outgoing_edges[o_edge_id] );
        connected++;
        index[connected] = node_info[node].outgoing_edges[o_edge_id];
        value[connected] = -1.0;
        o_edge_id++;
    }

    int i_edge_id = 0;  // incoming edge id
    while (node_info[node].incoming_edges[i_edge_id] != -1) {  /* record edge id and set its value to +1.0 */
        connected++;
        index[connected] = node_info[node].incoming_edges[i_edge_id];
        value[connected] = 1.0;
        i_edge_id++;
    }

    glp_set_row_bnds(lp, row_id, GLP_FX, 0.0, 0.0); /* RHS = 0 */
    glp_set_row_name(lp, row_id, "F.C.L");  // "flow conservation constraint"
    /* LHS: sum(incoming edges) - sum(outgoing edges) */
    glp_set_mat_row(lp, row_id, connected, index, value);
}

int
lp_make_row_color_distinct(glp_prob *lp, struct Node *node_info, int node, int NODES, int PAIRS, struct Color *st_pairs,
                           int row);

int
lp_make_row_color_distinct(glp_prob *lp, struct Node *node_info, int node, int NODES, int PAIRS, struct Color *st_pairs,
                           int row) {
    // get all the edges for the node in each color
    // each time, add a constraint: incoming - outgoing = 0 (for only the current color edges)
    for (int p = 0; p < PAIRS; p++) {
        int index[NODES * NODES];
        double value[NODES * NODES];
        int connected = 0;
        int i_edge_id = 0;  // incoming edge id
        // get the incoming edges of the current color
        while (node_info[node].incoming_edges[i_edge_id] != -1) {
            const char *name = glp_get_col_name(lp, node_info[node].incoming_edges[i_edge_id]);
            int edge_color = (int) name[strlen(name) - 2] - 48;
            if (edge_color == st_pairs[p].color) {
                connected++;
                index[connected] = node_info[node].incoming_edges[i_edge_id];
                value[connected] = 1.0;
            }
            i_edge_id++;
        }


        int o_edge_id = 0; // outgoing edge id
        while (node_info[node].outgoing_edges[o_edge_id] != -1) {
            // printf("%i\n",node_info[node].outgoing_edges[o_edge_id] );
            const char *name = glp_get_col_name(lp, node_info[node].outgoing_edges[o_edge_id]);
            int edge_color = (int) name[strlen(name) - 2] - 48;
            if (edge_color == st_pairs[p].color) {
                connected++;
                index[connected] = node_info[node].outgoing_edges[o_edge_id];
                value[connected] = -1.0;
            }
            o_edge_id++;
        }
        // set RHS
        glp_set_row_bnds(lp, row, GLP_FX, 0.0, 0.0); /* RHS = 0 */
        glp_set_row_name(lp, row, "C.D.F");  // "Color distinct flow"
        /* LHS: sum(incoming edges) - sum(outgoing edges) = 0, where all edges of the same color */
        glp_set_mat_row(lp, row, connected, index, value);
        row++;
    }
    return row;
}

bool sink_shares_e_with_source(struct Node *node_info, int linked_source, int sink);

bool sink_shares_e_with_source(struct Node *node_info, int linked_source, int sink) {
    for (int source_e = 0; source_e < numCols * numRows; source_e++)
    {
        if (node_info[linked_source].outgoing_edges[source_e] != -1) {
            for (int sink_e = 0; sink_e < numCols * numRows; sink_e++) {
                if (node_info[sink].incoming_edges[sink_e] != -1) {
                    if (node_info[sink].incoming_edges[sink_e] == node_info[linked_source].outgoing_edges[source_e]) {
                        return true;
                    }
                } else {
                    break;
                }
            }
        } else {
            return false;
        }
    }
}


int count_colors(const int NODES, int input_1d[]);

int count_colors(const int NODES, int input_1d[]) {
    int color = 0;
    for (int i = 0; i < NODES; i++) {
        if (input_1d[i] != 0) {
            color++;
        }
    }
    return color;
}

void init_node_info(const int NODES, struct Node *node_info);

void init_node_info(const int NODES, struct Node *node_info) {
    /* Initialise Node information */
    for (int node = 0; node < NODES; node++) {
        struct Node n;
        n.incoming_edges = (int *) malloc((numCols * numRows) * sizeof(int));
        n.outgoing_edges = (int *) malloc((numCols * numRows) * sizeof(int));
        for (int i = 0; i < (numCols * numRows); i++) {
            n.incoming_edges[i] = -1;
            n.outgoing_edges[i] = -1;
        }
        n.sink = false;
        n.source = false;
        node_info[node] = n;
    }
}

void init_adjacency_matrix(int *adjacency_matrix, const int NODES);

void init_adjacency_matrix(int *adjacency_matrix, const int NODES) {
    for (int i = 0; i < NODES; i++) {
        for (int j = 0; j < NODES; j++) {
            *(adjacency_matrix + i * NODES + j) = 0;
        }
    }
}

void
build_adjacency_matrix(int input_1d[], const int NODES, const int PAIRS, struct Color *st_pairs, struct Node *node_info,
                       int *adjacency_matrix);

void
build_adjacency_matrix(int input_1d[], const int NODES, const int PAIRS, struct Color *st_pairs, struct Node *node_info,
                       int *adjacency_matrix) {
    const int MAX_ADJACENT = 4;
    for (int i = 0; i < NODES; i++) {
        if (input_1d[i] != 0) { // if we have a source or sink node
            /* determine if the node is a source node*/
            bool is_source = false;
            struct Color matching_pair;
            for (int p = 0; p < PAIRS; p++) {
                if (st_pairs[p].color == input_1d[i]) {
                    matching_pair.color = st_pairs[p].color;
                    matching_pair.source_location = st_pairs[p].source_location;
                    matching_pair.sink_location = st_pairs[p].sink_location;
                    if (st_pairs[p].source_location == i) {
                        is_source = true;  /* otherwise we are dealing with a sink node */
                        break;
                    }
                }
            }
            if (is_source) {
                // make the edges for the source/sink node
                // adjacent node array of the top, left, right, and bottom nodes
                // if the node does not have a certain directionally adjacent node
                // it will remain -1 - indicating that there is no edge in this case
                node_info[i].source = true;
                int source = i;
                int an[4] = {-1, -1, -1, -1};
                int *adjacent_nodes = determine_adjacent_nodes(source, an);
                // make all edges from source to its adjacent nodes
                for (int node = 0; node < MAX_ADJACENT; node++) {
                    if (adjacent_nodes[node] != -1) {

                        // if the edge does not go into another source/sink color
                        int node_index = adjacent_nodes[node];
                        int node_value = input_1d[adjacent_nodes[node]];
                        if (node_index == matching_pair.sink_location || node_value == 0) {
                            // make all edges for the source
                            *(adjacency_matrix + source * NODES + node_index) = 1;
                        }
                    }
                }
            } else { //otherwise construct edges going into the sink
                node_info[i].sink = true;

                int sink = i;
                int an[4] = {-1, -1, -1, -1};
                int *adjacent_nodes = determine_adjacent_nodes(sink, an);
                // make all edges adjacent to the sink going into the sink
                for (int node = 0; node < MAX_ADJACENT; node++) {
                    if (adjacent_nodes[node] != -1) {
                        int node_index = adjacent_nodes[node];
                        int node_value = input_1d[adjacent_nodes[node]];
                        // if the edge does not go into a different colored source
                        if (node_index == matching_pair.source_location || node_value == 0) {
                            // make all edges for the sink
                            *(adjacency_matrix + node_index * NODES + sink) = 1;
                        }
                    }
                }
            }
        } else { // for all other non-source/sink edges
            int an[4] = {-1, -1, -1, -1};
            int *adjacent_nodes = determine_adjacent_nodes(i, an);
            for (int ad_node = 0; ad_node < 4; ad_node++) {
                // if the node has an adjacent node in this direction
                if (adjacent_nodes[ad_node] != -1) {
                    // if the adjacent node is not a source/sink
                    if (input_1d[adjacent_nodes[ad_node]] == 0) {
                        // make edges for all e in E not {s,t}
                        *(adjacency_matrix + i * NODES + adjacent_nodes[ad_node]) = 1;
                        // a path in the converse direction is also possible
                        *(adjacency_matrix + adjacent_nodes[ad_node] * NODES + i) = 1;
                    }
                }
            }
        }
    }
}

void lp_make_bounds(glp_prob *lp, const int NODES, const int EDGES, const int PAIRS, struct Color *st_pairs, int *adjacency_matrix,
                     struct Node *node_info, int input_1d[]);

void lp_make_bounds(glp_prob *lp, const int NODES, const int EDGES, const int PAIRS, struct Color *st_pairs, int *adjacency_matrix,
                     struct Node *node_info, int input_1d[]) {
    /*
    Create a bound for all edges in the graph.
    Additionally, store the node's incoming/outgoing edges for later use in determining
    what constraints are applied for each node's edges.
    */
    glp_add_cols(lp, EDGES);
    int col_id = 0;
    for (int node_u = 0; node_u < NODES; node_u++) {
        for (int node_v = 0; node_v < NODES; node_v++) {
            if (*(adjacency_matrix + node_u * NODES + node_v) == 1) {
                if (is_source(node_u, st_pairs, PAIRS)) {
                    col_id++;
                    int color = get_color(node_u, node_v, input_1d);
                    lp_make_col(lp, col_id, node_u, node_v, color);
                    /* store node's edges */
                    node_info_outgoing_edge(node_info, node_u, col_id, NODES);
                    /* store incoming edge to node_v's (u -> v) info */
                    node_info_incoming_edge(node_info, node_v, col_id, NODES);
                } else if (is_sink(node_v, st_pairs, PAIRS)) {
                    col_id++;
                    int color = get_color(node_u, node_v, input_1d);
                    lp_make_col(lp, col_id, node_u, node_v, color);
                    node_info_incoming_edge(node_info, node_v, col_id, NODES);
                    /* store outgoing edge to node_u's info */
                    node_info_outgoing_edge(node_info, node_u, col_id, NODES);
                } else {  /* if the edge does NOT involve a source or sink node */
                    for (int p = 0; p < PAIRS; p++) {  /* create the edge for each color */
                        col_id++;
                        int color = st_pairs[p].color;
                        lp_make_col(lp, col_id, node_u, node_v, color);
                        node_info_incoming_edge(node_info, node_v, col_id, NODES);
                        node_info_outgoing_edge(node_info, node_u, col_id, NODES);
                    }
                }
            }
        }
    }
}

void init_input_1d(int input_1d[]);

void init_input_1d(int input_1d[]) {
    for (int i = 0, next_node = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            input_1d[next_node] = input[i * numCols + j];
            next_node++;
        }
    }
}
void lp_make_constraints(glp_prob *lp, const int NODES, const int PAIRS, struct Color *st_pairs,
                          struct Node *node_info, int input_1d[]);

void lp_make_constraints(glp_prob *lp, const int NODES, const int PAIRS, struct Color *st_pairs,
                          struct Node *node_info, int input_1d[]) {
    int constraint_3_row = 1 + NODES + NODES;  // row constraint ids for the third batch of constraints assigned to each non-s/t node
    const int CONSTRAINTS = (NODES * PAIRS) + (NODES * (PAIRS*2)); // The total number of constraints
    glp_add_rows(lp, CONSTRAINTS);
    for (int node = 0; node < NODES; node++) {
        int index[NODES];    /* indices to define constraint coefficients */
        double value[NODES]; /* values to define constraint coefficients */

        if (node_info[node].source) {
            int edge_id = 0;
            int connected = 0;
            /* Add the source node's outgoing edges to the objective function */
            while (node_info[node].outgoing_edges[edge_id] != -1) {
                /* set objective function to max flow of source edges */
                connected++;  /* count number of connected edges */
                glp_set_obj_coef(lp, node_info[node].outgoing_edges[edge_id], 1.0);
                // record edge id and set its value to +1.0
                index[connected] = node_info[node].outgoing_edges[edge_id];
                value[connected] = 1.0;
                edge_id++;
            }
            /* set constraint for source: sum(outgoing edges) <= 1 */
            /* Write-up LP constraint (4) */
            glp_set_row_name(lp, 1 + node, "SOURCE");
            glp_set_row_bnds(lp, 1 + node, GLP_UP, 0.0, 1.0);
            glp_set_mat_row(lp, 1 + node, connected, index, value);
        }
        else if (node_info[node].sink) {
            int connected = 0; /* number of edges in the constraint */
            int linked_source;
            for (int p = 0; p < PAIRS; p++) {  /* find the sink's source */
                if (st_pairs[p].color == input_1d[node]) {
                    linked_source = st_pairs[p].source_location;
                    break;
                }
            }

            if (!sink_shares_e_with_source(node_info, linked_source, node)) {
            /* for all non-shared s-t edges, add edge to sink constraint */
                int source_edge_id = 0;
                while (node_info[linked_source].outgoing_edges[source_edge_id] !=-1) {
                    /* record edge id and set its value to +1.0 */
                    connected++;
                    index[connected] = node_info[linked_source].outgoing_edges[source_edge_id];
                    value[connected] = 1.0;
                    source_edge_id++;
                }
                int sink_edge_id = 0;
                while (node_info[node].incoming_edges[sink_edge_id] != -1) {
                    connected++;
                    /* record edge id and set its value to -1.0 */
                    index[connected] = node_info[node].incoming_edges[sink_edge_id];
                    value[connected] = -1.0;
                    sink_edge_id++;
                }
                glp_set_row_name(lp, 1 + node, "SINK");
                glp_set_row_bnds(lp, 1 + node, GLP_FX, 0.0, 0.0); /* RHS: = 0 */
                /* LHS: sum(sink's incoming_edges) - sum(source's outgoing edges) 1 */
                glp_set_mat_row(lp, 1 + node, connected, index, value);
            }
        }
        else { /* For all v \ {s,t} */
            int connected = 0;
            int edge_id = 0;
            while (node_info[node].incoming_edges[edge_id] != -1) {
                /* record edge id and set its value to +1.0 */
                connected++;
                index[connected] = node_info[node].incoming_edges[edge_id];
                value[connected] = 1.0;
                edge_id++;
            }
            /* Write-up LP constraint (4) */
            glp_set_row_bnds(lp, 1 + node, GLP_UP, 0.0, 1.0); /* RHS <=1 */
            glp_set_mat_row(lp, 1 + node, connected, index, value);  /* LHS: sum(incoming edges) */

            /* Write-up LP constraint (3) */
            lp_make_row_flow_cons(lp, node_info, node, NODES);

            /* Write-up LP constraint (5)  */
            constraint_3_row = lp_make_row_color_distinct(lp, node_info, node, NODES, PAIRS, st_pairs,
                                                          constraint_3_row);
        }
    }
}

int computeSolution(void) {
    const int NODES = numRows * numCols; /* total nodes in puzzle */
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX); /* set maximisation as objective */
    //glp_set_prob_name(lp, "Max Flow"); /* allows for a title; useful when writing the LP to a file */

    /* map the input graph onto a 1d array - easier for me to understand when referencing nodes */
    int input_1d[NODES];
    init_input_1d(input_1d);

    /* count the total number of source & sink nodes */
    const int COLORS = count_colors(NODES, input_1d);

    /* An array of Node structs */
    struct Node *node_info = (struct Node *) malloc(NODES * sizeof(struct Node));
    init_node_info(NODES, node_info);

    /*
      Find source-sink pair nodes.
      Necessary for row/col creation, as these nodes will have different contraints to non-st nodes
    */
    const int PAIRS = COLORS / 2;  /* The number of source sink pairs (where equal colors) */
    struct Color *st_pairs = (struct Color *) malloc(PAIRS * sizeof(struct Color));
    find_st_pairs(NODES, PAIRS, st_pairs, node_info, input_1d);

    /* make the adjacency matrix */
    int *adjacency_matrix = (int *) malloc(NODES * NODES * sizeof(int));
    init_adjacency_matrix(adjacency_matrix, NODES);
    build_adjacency_matrix(input_1d, NODES, PAIRS, st_pairs, node_info, adjacency_matrix);

    /* set all edge capacities (rows/bounds) */
    const int EDGES = count_edges(NODES, PAIRS, adjacency_matrix, node_info);
    lp_make_bounds(lp, NODES, EDGES, PAIRS, st_pairs, adjacency_matrix, node_info, input_1d);

    /* set all constraints (columns)*/
    lp_make_constraints(lp, NODES, PAIRS, st_pairs, node_info, input_1d);

    glp_term_out(0); /* disable terminal output from glpk routines */
    glp_simplex(lp, NULL); /* solve LP via Simplex algorithm */

    /* DEBUG functions */
    // glp_write_lp(lp, NULL, "ignore_files/output.txt"); // write the LP problem to a file
    // print_edge_flows(lp, EDGES);  /* print edge flow values */
    // printf("\nMaximal flow is %f\n\n", glp_get_obj_val(lp));

    /* Let s = all source nodes in input, return 1 if max flow = |s|, otherwise 0 */
    if (glp_get_obj_val(lp) == PAIRS)
    {
      /* build the solution graph */
      build_solution(lp, EDGES, solution);
      glp_delete_prob(lp); /* house-keeping */
      return 1;
    }
    else
    {
      glp_delete_prob(lp);
      return 0;
    }
}

/* YOU SHOULD NOT CHANGE ANYTHING BELOW THIS LINE! */

/* printPuzzle(int *puzzle) prints either an input or a solution */
void printPuzzle(int *puzzle) {
    int i, j; /* loop variables to go over rows and columns */
    for (i = 0; i < numRows; i++) { /* go over rows */
        for (j = 0; j < numCols; j++) {                                               /* go over columns */
            fputc('0' + puzzle[i * numCols + j], stdout); /* print the next char */
        }
        fprintf(stdout, "\n"); /* end the current line, start new one */
    }
}

int main(int argc, char **argv) {
    int i; /* used to run over the command line parameters */

    if (argc < 2) { /* no command line parameter given */
        fprintf(stderr, "Usage: %s [file1] [file2] [file3] [...]\n"
                        "Where each [file] is the name of a file with a puzzle.\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }

    if (argv[1][0] == '-' && argv[1][1] == 'd' && argv[1][2] == 0) {
        /* If the first parameter is -d we activate debug mode. */
        debug = 1;                                        /* switch debug mode on */
        fprintf(stdout, "DEBUG: Debug mode activated\n"); /* be explicit about it */
    } else {
        debug = 0; /* switch debug mode off */
    }

    for (i = 1 + debug; i < argc; i++) { /* go over remaining command line parameters */
        if (readInput(argv[i])) { /* try to read file */
            /* returned with error message */
            fprintf(stderr, "%s: Cannot read puzzle with filename %s. Skipping it.\n",
                    argv[0], argv[i]);
        } else { /* input read successfully */
            fprintf(stdout, "%s: Looking at the following puzzle:\n", argv[i]);
            printPuzzle(input);
            if (computeSolution()) { /* compute a solution if one exists */
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
    int read;       /* counts number of digits read in the current row */
    int firstRow;   /* indicates if we are reading the very first row */

    firstRow = 1; /* we start in the first row */
    rows = cols = read = 0;
    while (!feof(fh)) {
        c = fgetc(fh); /* read the next char from the file */
        if (isdigit(c)) {         /* normal character read */
            read++; /* count the digit we just read */
            if ((!firstRow) && (read > cols)) {
                if (debug) {
                    fprintf(stdout, "DEBUG: Row %d too long (%d, %d).\n", rows + 1, read, cols);
                }
                return 1; /* flag error because row is too long */
            }
        } else {
            if ((c == '\n') || (c == (char) -1)) { /* end of line read */
                if (read > 0) {
                    rows++; /* count the completed row if it was not empty */
                }
                if (firstRow) {               /* very first row read */
                    cols = read;  /* accept number of characters as number of columns */
                    firstRow = 0; /* not in the first row anymore after this */
                    if (debug) {
                        fprintf(stdout, "DEBUG: %d columns per row expected\n", cols);
                    }
                } else {
                    if ((read > 0) && (read != cols)) { /* rows too short */
                        if (debug) {
                            fprintf(stdout, "DEBUG: Row %d too short.\n", rows + 1);
                        }
                        return 1; /* flag error because row is too short */
                    }
                }
                read = 0; /* reset number of characters in current row */
            } else { /* illegal character found */
                if (debug) {
                    fprintf(stdout, "DEBUG: Illegal character %c found.\n", c);
                }
                return 1; /* stop reading because of the error */
            }
        }
    }
    if (read > 0) {
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
int readInput(char *filename) {
    int i, j;        /* loop variables to go over the columns and rows of the input */
    int check[10];   /* array to check colours come in pairs */
    int keepReading; /* used to skip over newline */
    char c;          /* next char */
    FILE *fh;

    if ((fh = fopen(filename, "rt")) == NULL) {
        return 1;
    }

    /* perform basic checks and determine size of puzzle */
    if (checkFile(fh)) { /* there was a problem */
        fclose(fh);
        return 1; /* signal error */
    }
    if ((input = (int *) malloc(sizeof(int) * numRows * numCols)) == NULL) {
        if (debug) {
            fprintf(stdout, "DEBUG: Unable to allocate %ld bytes of memory.\n",
                    sizeof(int) * numRows * numCols);
        }
        fclose(fh);
        return 1;
    }
    if ((solution = (int *) malloc(sizeof(int) * numRows * numCols)) == NULL) {
        if (debug) {
            fprintf(stdout, "DEBUG: Unable to allocate %ld bytes of memory.\n",
                    sizeof(int) * numRows * numCols);
        }
        free(input);
        fclose(fh);
        return 1;
    }
    memset(solution, 0, sizeof(int) * numRows * numCols); /*initialise solution empty*/
    if (debug) {
        fprintf(stdout, "DEBUG: Will read %dx%d sized puzzle\n", numRows, numCols);
    }
    /* prepare to count different digits */
    for (i = 0; i < 10; i++) {
        check[i] = 0;
    }
    /* Size is given in numRows, numCols; now we read */

    for (i = 0; i < numRows; i++) { /* go over rows */
        for (j = 0; j < numCols; j++) { /* go over columns */
            do {
                keepReading = 1; /* prepare to skip over newline */
                c = fgetc(fh);   /* get next digit */
                if (isdigit(c)) {                                          /* store and count digit */
                    input[i * numCols + j] = (int) (c - '0'); /* convert char to int */
                    check[input[i * numCols + j]]++;
                    keepReading = 0; /* mark digit as read */
                }
            } while (keepReading);
        }
    }
    for (i = 1; i < 10; i++) {
        if ((check[i] != 0) && (check[i] != 2)) {
            if (debug) {
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
