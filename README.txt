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