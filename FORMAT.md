# Formats
We give descriptions of formats of configurations, rules, combined rules, auxiliary rules, cartwheels, and multi-boundary islands

## Configuration Format
(blank line) \
$N$ $R$ \
$R+1$ $d(v_{R+1})$ $a_{R+1}^1$ $\ldots$ $a_{R+1}^{d(v_{R+1})}$ \
$\vdots$ \
$N$ $d(v_{N})$ $a_{N}^1$ $\ldots$ $a_{N}^{d(v_{N})}$ \
$M$ \
$a_1$ $b_1$ \
$\vdots$ \
$a_M$ $b_M$

---

$N$ denotes the total number of vertices of a configuration.
$R$ denotes the ring size of a configuration.
The vertex set is $V=\{v_1, \ldots, v_{N}\}$, and vertex set in the ring is $V_R=\{v_1, \ldots, v_{R}\} \subseteq V$.
For each internal vertex $v_i$ ($i=R+1, \ldots, N$), the file lists the vertex index $i$, the number of adjacent vertices $d(v_i)$, and the clockwise rotation of indices of its neighbors $a_i^1, \ldots, a_i^{d(v_i)}$.
If $v_i$ is not incident to a digon, $\delta_K(v_i)=d(v_i)$, but otherwise, $\delta_K(v_i)=d(v_i) + m_i$, where $m_i$ is the number of incident digons (but $m_i=1$ for our configurations.)
Specifically, the vertices in the ring $v_1, \ldots, v_R$ are clockwisely ordered.
$M$ denotes the number of digons in the configuration.
$a_i$ $b_i$ ($i=1,\ldots,M$) are endpoints of $i$-th digon.
If $M=0$, we can omit $M$.

## Rule Format
(blank line) \
$N$ $s$ $t$ $r$ \
$1$ $\delta^-(v_1)$ $\delta^+(v_1)$ $a_1^1$ $\ldots$ $a_1^{d(v_1)}$ \
$\vdots$ \
$N$ $\delta^-(v_N)$ $\delta^+(v_N)$ $a_N^1$ $\ldots$ $a_N^{d(v_N)}$ \
$M$ \
$a_1$ $b_1$ \
$\vdots$ \
$a_M$ $b_M$

---

$N$ denotes the total number of vertices of a rule.
The vertex set is $V=\{v_1, \ldots, v_N\}$.
$s, t$ denotes the index of the vertex that charge sends/receives respecitvely.
$r$ denotes the amount of charge.
For each vertex $v_i$ ($i=1,\ldots,N$), the file lists the vertex index $i$, its degree-range $\delta^-(v_i), \delta^+(v_i)$, and the clockwise rotation of indices of its neighbors $a_i^1, \ldots, a_i^{d(v_i)}$.
The degree $\infty$ is represented by $0$.
When $a_i^j=-1$, it representes the boundary.
$M$ denotes the number of diongs in the configuration.
$a_i$ $b_i$ ($i=1,\ldots,M$) are endpoints of $i$-th digon.
If $M=0$, we can omit $M$.

## Combined Rule Format
(blank line) \
$N$ $s$ $t$ $r$ \
$1$ $\delta^-(v_1)$ $\delta^+(v_1)$ $a_1^1$ $\ldots$ $a_1^{d(v_1)}$ \
$\vdots$ \
$N$ $\delta^-(v_N)$ $\delta^+(v_N)$ $a_N^1$ $\ldots$ $a_N^{d(v_N)}$ \
$M$ \
$a_1$ $b_1$ \
$\vdots$ \
$a_M$ $b_M$ \
$F$

---

Except for the last line, the format is identical to the Rule format, but we must write $M$ even if $M=0$.
The last line contains a $01$ sequence $F$ representing which rules are combined.
When combining rules, we order rules by their filenames, resulting $R_0, \ldots, R_{k-1}$.
The length of $F$ is $k$.
If the $i$-th character $F[i]$ is $1$, it means the rule $R_i$ is included in the combination.

## Auxiliary Rule Format
$R$ \
(blank line) \
$k$ \
$R_1$ \
$\vdots$ \
$R_k$ \

---
$R,R_1,\ldots,R_k$ are rules represented by the above rule format.
$k$ is the number of rules that consists of homomorphic cover of the rule $R$.
We must write the number of digons $M$ even if $M=0$ for each rule $R, R_1, \ldots, R_k$.
It is necessary to parse the next rule correctly.

## CartWheel Format
(blank line) \
$N$ $c$ \
$1$ $\delta^-(v_1)$ $\delta^+(v_1)$ $a^1_1$ $\ldots$ $a^{d(v_1)}_1$ \
$\vdots$ \
$N$ $\delta^-(v_N)$ $\delta^+(v_N)$ $a^1_N$ $\ldots$ $a^{d(v_1)}_N$ \
$M$ \
$a_1$ $b_1$ \
$\vdots$ \
$a_M$ $b_M$

---

$N$ denotes the total number of vertices of a cartwheel.
The vertex set is $V=\{v_1, \ldots, v_N\}$.
$c$ denotes the index of the center of a cartwheel.
The other format is identical to the Rule format.


## Multi Boundary Island Format
$N$ \
$R_1$ ... $R_k$ \
$M$ \
$a_1$ $b_1$ $c_1$ \
...               \
$a_N$ $b_N$ $c_N$

---
$N$ denotes the number of vertices of degree 2 or 3.
$k$ is the number of rings.
$R_i$ ($1 \leq i \leq k$) denotes the number of edges that are incident to the degree-one vertices in the $i$-th ring.
The index of these edges are $\{a_i, a_i + 1, \ldots, a_i + R_i - 1\}$ where $a_i = \sum_{1 \leq j < i} R_j$.
These indices are ordered in clockwise order of rotatinos.
$M$ denotes the number of vertices of degree 2.
We add a dummy edge to each vertex of degree 2 to make its degree 3.
The index of these edges are $\{a_{k+1}, a_{k+1} + 1, \ldots, a_{k+1} + M - 1\}$ where $a_{k+1} = \sum_{1 \leq j < k + 1} R_j$.
For the $i$-th vertex, $a_i,b_i,c_i$ is the indices of the three incident edges.
If the degree of this vertex is 3, they are listed in the cyclic order around the vertex in the planar embedding.
