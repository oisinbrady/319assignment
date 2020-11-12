<h1>CS31920 Advanced Algorithms Assignment</h1>

<h2>Problem Description</h2>

<p>"We are concerned with a small puzzle game for a single player, played on a rectangular grid. On the grid are pairs of squares marked with equal colours (each colour appearing exactly twice). You are asked to finda connection for all pairs simultaneously so that the path between each pair makes only use of horizontal and vertical steps (not diagonal steps) and each field of the grid is used by at most one connection or mark. [...] Note that possible solutions need not be unique [...] Not all inputs are necessarily solvable"[1]<p>

<h2>Implementation</h2>

This program reduces the puzzle into a maximum flow problem, before solving it via linear programming using the C GNU Linear Programming Kit (GLPK).

For graph, G(V,E,c):

MAXIMISE Σ f(e(s,.)), For s ∈ V [MAXIMISE SOURCE EDGE FLOW VALUES]

S.T: (GLPK: Columns [structural variables])

  [SOURCE-TO-SINK FLOW CONSERVATION]
  (1): Σ f(e(s,.)) - Σ f(e(.,t)) = 0, For v ∈ {s,t} where, color(s) = color(t)

  [FLOW CONSERVATION]
  (2): Σ f(e(v,.)) - Σ f(e(.,v)) = 0, For v ∈ V \{s,t}

  [NON-BRANCHING PATHS]
  (3): Σ f(e(.,v)) <= 1, For v ∈ V \{s}

  [EDGE DISJOINT PATHS B/W COLORS]
  (4): Σ f(e(.,v)) - Σ f(e(u,.)) = 0, For v ∈ V \{s,t} where, color(v) = color(u)

  Where,
    color(x): color of source/sink node.
    f(e): the flow value of an edge, e.

BOUNDS: (GLPK: Row [auxiliary variables])
  - 0 <= e <= 1 for e in E
 
[1]: Section 3.1 (Problem Description) of assignment brief

<h2>Compiling and Running</h2>

- Ensure GLPK is installed: https://www.gnu.org/software/glpk/

- ```gcc puzzleSolver.c -lglpk -o puzzleSolver_c```

- ```./puzzleSolver_c <input file> <-d (optional debug)>```
