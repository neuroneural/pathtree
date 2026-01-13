% reverse.lp
% Expand PathTree/PathForest facts (root/loop) into a unit-lag directed graph
% using the minimal-node construction from Section 3.2.
% Also expand bidirected edges into a hidden parent with unit chains whose lengths
% encode the required lag differences.

% Domain / nodes
node(U) :- observed(U).
node(U) :- root(U,_,_).
node(V) :- root(_,V,_).
node(U) :- loop(U,_,_).
node(V) :- loop(_,V,_).

% Latent symbols (helpers)
% Chain nodes for the root path of length A>1: p(U,V,A,1..A-1)
latent(p(U,V,A,K)) :- root(U,V,A), A > 1, K = 1..A-1.

% --- Back-compat adapter for unscoped loops loop(U,V,D) -> loop(U,V,A,D) ---

% Is there any root with A > 1?
has_root_gt1(U,V) :- root(U,V,A), A > 1.

% Pick the smallest A > 1 (so periodic trees like (2,<2>) bind to A=2)
amin_gt1(U,V,Am) :- has_root_gt1(U,V), Am = #min { A : root(U,V,A), A > 1 }.

% Map legacy loops to that A; otherwise fall back to A=1.
loop(U,V,Am,D) :- loop(U,V,D), amin_gt1(U,V,Am).
loop(U,V,1,D)  :- loop(U,V,D), not has_root_gt1(U,V).

% ---------- Bundle selection (per-root version) ----------
% If A>1, bundle is first chain node p(U,V,A,1)
bundle(U,V,A,p(U,V,A,1)) :- root(U,V,A), A > 1.

% Only create a bundle if THIS root itself has a loop
has_loop_pair(U,V,A) :- root(U,V,A), loop(U,V,A).
bundle(U,V,A,b(U,V)) :-
    root(U,V,A), has_loop_pair(U,V,A), A = 1,
    not bi_edge(U,V).     % <-- new guard: skip if bidirected


% If A=1 and there are NO loops
bundle(U,V,A,V) :- root(U,V,A), A = 1,
    not has_loop_pair(U,V,A), not bi_edge(U,V).
bundle(U,V,A,b(U,V)) :- root(U,V,A), A = 1,
    not has_loop_pair(U,V,A), bi_edge(U,V).

% Declare the latent bundle node only when used
latent(b(U,V)) :- bundle(U,V,_,b(U,V)).

% ---------- Build the root path U -> ... -> V ----------
% helper: mark (U,V) as synthetic-only if there is a bidirected edge
% but no reverse root(V,U,_)
synthetic_bidir_only(U,V) :-
    bi_edge(U,V),
    not root(U,V,_),
    not root(V,U,_).

% A = 1: keep the direct edge unless it's only a synthetic bidirected mirror
dir_unique(U,V,1) :-
    root(U,V,1),
    not has_loop_pair(U,V,1),
    not synthetic_bidir_only(U,V).
    
% Allow A=1 looped roots to keep direct edge if an A>1 sibling exists
dir_unique(U,V,1) :-
    root(U,V,1), has_loop_pair(U,V,1),
    has_root_gt1(U,V).

dir_unique(U,b(U,V),1)  :- root(U,V,1), has_loop_pair(U,V,1).
dir_unique(b(U,V),V,1)  :- root(U,V,1), has_loop_pair(U,V,1).

% A > 1: chain U -> p(...,1) -> ... -> p(...,A-1) -> V
dir_unique(U,p(U,V,A,1),1)            :- root(U,V,A), A > 1.
dir_unique(p(U,V,A,K),p(U,V,A,K+1),1) :- root(U,V,A), A > 2, K = 1..A-2.
dir_unique(p(U,V,A,A-1),V,1)          :- root(U,V,A), A > 1.

% For the cycle backbone we use positions 1..M-1 as new latent nodes;
% position 0 is the bundle X itself.
% M is the largest loop length among siblings for (U,V).
% Compute M per specific root(A)
mmax(U,V,A,M) :- root(U,V,A), M = #max { D : loop(U,V,A,D) }.

% Backbone latent cycle nodes for each A
cycpos(U,V,A,K) :- mmax(U,V,A,M), M > 1, K = 1..M-1.
latent(cyc(U,V,A,K)) :- cycpos(U,V,A,K).

% Close the length-M cycle for that A
dir_unique(X,cyc(U,V,A,1),1)                :- mmax(U,V,A,M), M > 1, bundle(U,V,A,X).
dir_unique(cyc(U,V,A,K),cyc(U,V,A,K+1),1)   :- mmax(U,V,A,M), M > 2, K = 1..M-2.
dir_unique(cyc(U,V,A,M-1),X,1)              :- mmax(U,V,A,M), M > 1, bundle(U,V,A,X).

% D=1 means a true 1-cycle on the bundle for this root(A)
% This is needed so a PathTree like (A,<1>) actually becomes a loop source.
dir_unique(X,X,1) :-
    loop(U,V,A,1),
    bundle(U,V,A,X).

dir_unique(cyc(U,V,A,D-1),X,1) :-
    loop(U,V,A,D), D > 1,
    bundle(U,V,A,X).

% ================================================================
% Expand bidirected lag sets using baseline distance (from hide_nodes)
% ================================================================
% Extend bi_diff(U,V,D0) using base_min(U,V,Dmin)
% D0 = |Lx - Ly| ; Dmin = distance from hidden parent to nearest observed child
% This yields all valid D = D0 + K - Dmin  up to (maxlag - Dmin)

lag(0..maxlag).

% Zero-diff pair: shift by -Dmin to preserve equality semantics
lag_expanded(U,V,D) :-
    bi_zero(U,V),
    base_min_final(U,V,Dmin),
    lag(Draw),
    D = Draw - Dmin,
    D >= 0,
    D <= maxlag - Dmin.

% Nonzero-diff pair: offset by both D0 and -Dmin
lag_expanded(U,V,D) :-
    bi_diff(U,V,D0),
    base_min_final(U,V,Dmin),
    lag(K), K >= 0,
    D = D0 + K - Dmin,
    D >= 0,
    D <= maxlag - Dmin.

% Fallback (no base_min available)
lag_expanded(U,V,D) :-
    bi_diff(U,V,D),
    not base_min_final(U,V,_).

% Bidirected edges -> introduce a hidden parent H with unit chains

% Create a zero-diff parent for every bidirected pair, independent of order
latent(h0(U,V)) :- bi_edge(U,V), bi_zero(U,V).

% Each child sits one unit from parent; Dmin affects lag-expansion only.
dir_unique(h0(U,V),U,1) :-
    bi_edge(U,V), bi_zero(U,V),
    base_min_final(U,V,_).

dir_unique(h0(U,V),V,1) :-
    bi_edge(U,V), bi_zero(U,V),
    base_min_final(U,V,_).

dir_unique(h0(U,V),U,1) :-
    bi_edge(U,V), bi_zero(U,V),
    not base_min_final(U,V,_).

dir_unique(h0(U,V),V,1) :-
    bi_edge(U,V), bi_zero(U,V),
    not base_min_final(U,V,_).

% --- attach any loops from (U,V) onto its zero-diff parent h0(U,V)
%   D = 1  → self-loop on the parent
%   D > 1  → chord from cyc(D-1) back to the parent through shared backbone
% --- Adapter: expose per-root loops to bidirected builders ---
loop(U,V,D) :- loop(U,V,A,D).
% Largest loop length among these siblings
mmax_bi(U,V,M) :- bi_edge(U,V), U < V, bi_zero(U,V), M = #max { D : loop(U,V,D) }.

% backbone cycle for h0(U,V)
mmax_bi(U,V,M) :- bi_edge(U,V), bi_zero(U,V), M = #max { D : loop(U,V,D) }.

bcycpos(U,V,K) :- mmax_bi(U,V,M), M > 1, K = 1..M-1.
latent(bcyc(U,V,K)) :- bcycpos(U,V,K).

dir_unique(h0(U,V), bcyc(U,V,1), 1)            :- bi_edge(U,V), bi_zero(U,V), mmax_bi(U,V,M), M > 1.
dir_unique(bcyc(U,V,K), bcyc(U,V,K+1), 1)      :- bi_edge(U,V), bi_zero(U,V), mmax_bi(U,V,M), M > 2, K = 1..M-2.
dir_unique(bcyc(U,V,M-1), h0(U,V), 1)          :- bi_edge(U,V), bi_zero(U,V), mmax_bi(U,V,M), M > 1.

dir_unique(h0(U,V), h0(U,V), 1)                :- bi_edge(U,V), bi_zero(U,V), loop(U,V,1).
dir_unique(bcyc(U,V,D-1), h0(U,V), 1)          :- bi_edge(U,V), bi_zero(U,V), loop(U,V,D), D > 1.

% Nonzero diff: parent->U length 1, parent->V length (1 + D) via unit chain
latent(h(U,V,D))          :- bi_edge(U,V), lag_expanded(U,V,D), U != V.
vlen(U,V,L)               :- bi_edge(U,V), lag_expanded(U,V,D), L = D + 1, U != V.
dir_unique(h(U,V,D),U,1)  :- bi_edge(U,V), lag_expanded(U,V,D), U != V.

latent(pu(h(U,V,D),V,L,K)) :-
    bi_edge(U,V), lag_expanded(U,V,D), U != V,
    vlen(U,V,L), L > 1, K = 1..L-1.

dir_unique(h(U,V,D), pu(h(U,V,D),V,L,1), 1) :-
    bi_edge(U,V), lag_expanded(U,V,D), U != V,
    vlen(U,V,L), L > 1.

dir_unique(pu(h(U,V,D),V,L,K), pu(h(U,V,D),V,L,K+1), 1) :-
    bi_edge(U,V), lag_expanded(U,V,D), U != V,
    vlen(U,V,L), L > 2, K = 1..L-2.

dir_unique(pu(h(U,V,D),V,L,L-1), V, 1) :-
    bi_edge(U,V), lag_expanded(U,V,D), U != V,
    vlen(U,V,L), L > 1.

% ================================================================
% Bidirected loops: attach loopset only to the canonical zero-diff parent
% ================================================================
% Only h0(U,V) participates in the backbone cycle
% Nonzero-diff parents h(U,V,D) do NOT connect to the ring or to 8 directly

% Canonical parent for the bidirected pair
parent_zero(U,V,h0(U,V)) :- bi_edge(U,V), bi_zero(U,V).

% Determine maximum loop length
mmax_bi(U,V,M) :- bi_edge(U,V), M = #max { R : loop(U,V,R) }.

% Create the backbone ring of M−1 latent nodes
bcycpos(U,V,K) :- mmax_bi(U,V,M), M > 1, K = 1..M-1.
latent(bcyc(U,V,K)) :- bcycpos(U,V,K).

% Connect only h0 into the ring
dir_unique(h0(U,V), bcyc(U,V,1), 1)            :- parent_zero(U,V,_), mmax_bi(U,V,M), M > 1.
dir_unique(bcyc(U,V,K), bcyc(U,V,K+1), 1)      :- mmax_bi(U,V,M), M > 2, K = 1..M-2.
dir_unique(bcyc(U,V,M-1), h0(U,V), 1)          :- parent_zero(U,V,_), mmax_bi(U,V,M), M > 1.

% True self-loop on h0 (D=1) or closure for longer loopsets (D>1)
dir_unique(h0(U,V), h0(U,V), 1)                :- bi_edge(U,V), bi_zero(U,V), loop(U,V,1).
dir_unique(bcyc(U,V,D-1), h0(U,V), 1)          :- bi_edge(U,V), bi_zero(U,V), loop(U,V,D), D > 1.

% All other parents h(U,V,D) remain simple (no backbone connections)


% Output node set (for visibility)
out_node(U) :- observed(U).
out_node(X) :- latent(X).
out_node(U) :- node(U).

% Show only the outputs we care about

#show dir_unique/3.
#show latent/1.
#show out_node/1.
#show node/1.
#show observed/1.
#show bundle/3.
#show loop/3.