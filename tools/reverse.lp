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

% ---------- Bundle selection ----------
% If A>1, bundle is first chain node p(U,V,A,1) (unchanged).
bundle(U,V,p(U,V,A,1)) :- root(U,V,A), A > 1.

% If A=1 and there is ANY loop on (U,V), create a dedicated bundle b(U,V)
% so loops don't land on V itself.
% Only treat loops as belonging to the *directed* edge when the pair is not bidirected.
has_loop_pair(U,V) :- loop(U,V,_), not bi_edge(U,V).
bundle(U,V,b(U,V)) :- root(U,V,1), has_loop_pair(U,V).

% If A=1 and there are NO loops:
%   - for non-bidirected pairs, we can use V directly (original behavior)
%   - for bidirected pairs, we never bind to the observed node V (avoid bundle(8,9,9))
bundle(U,V,V) :- root(U,V,1), not has_loop_pair(U,V), not bi_edge(U,V).
bundle(U,V,b(U,V)) :- root(U,V,1), not has_loop_pair(U,V), bi_edge(U,V).

% Declare the latent bundle node only when used.
latent(b(U,V)) :- bundle(U,V,b(U,V)).


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
    not has_loop_pair(U,V),
    not synthetic_bidir_only(U,V).

% A = 1: if there ARE loops, route via the bundle b(U,V)
dir_unique(U,b(U,V),1)  :- root(U,V,1), has_loop_pair(U,V).
dir_unique(b(U,V),V,1)  :- root(U,V,1), has_loop_pair(U,V).

% A > 1: chain U -> p(...,1) -> ... -> p(...,A-1) -> V
dir_unique(U,p(U,V,A,1),1)            :- root(U,V,A), A > 1.
dir_unique(p(U,V,A,K),p(U,V,A,K+1),1) :- root(U,V,A), A > 2, K = 1..A-2.
dir_unique(p(U,V,A,A-1),V,1)          :- root(U,V,A), A > 1.

% For the cycle backbone we use positions 1..M-1 as new latent nodes;
% position 0 is the bundle X itself.
% M is the largest loop length among siblings for (U,V).
mmax(U,V,M) :- root(U,V,_), M = #max { D : loop(U,V,D) }.

% Backbone latent cycle nodes for M>1: cyc(U,V,1..M-1)
cycpos(U,V,K) :- mmax(U,V,M), M > 1, K = 1..M-1.
latent(cyc(U,V,K)) :- cycpos(U,V,K).

% ================================================================
% Build the bundled cycles for siblings (leaf-level minimality)
% ================================================================
% Let X be the bundle node for (U,V).
% Build the longest cycle of length M using only M-1 new nodes cyc(...,1..M-1)
% and re-use that backbone for all shorter loop lengths.

% Close the length-M cycle:
% X -> cyc1 -> cyc2 -> ... -> cyc(M-1) -> X
dir_unique(X,cyc(U,V,1),1)                :- mmax(U,V,M), M > 1, bundle(U,V,X).
dir_unique(cyc(U,V,K),cyc(U,V,K+1),1)     :- mmax(U,V,M), M > 2, K = 1..M-2.
dir_unique(cyc(U,V,M-1),X,1)              :- mmax(U,V,M), M > 1, bundle(U,V,X).

% For each loop length D:
% - if D = 1  → add a self-loop on X
% - if D > 1  → add a "chord" from cyc(D-1) back to X
dir_unique(X,X,1)                         :- loop(U,V,1), bundle(U,V,X).
dir_unique(cyc(U,V,D-1),X,1)              :- loop(U,V,D), D > 1, bundle(U,V,X).

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
% ================================================================

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
latent(h(U,V,D))          :- bi_edge(U,V), lag_expanded(U,V,D).
vlen(U,V,L)               :- bi_edge(U,V), lag_expanded(U,V,D), L = D + 1.
dir_unique(h(U,V,D),U,1)  :- bi_edge(U,V), lag_expanded(U,V,D).

latent(pu(h(U,V,D),V,L,K)) :-
    bi_edge(U,V), lag_expanded(U,V,D),
    vlen(U,V,L), L > 1, K = 1..L-1.

dir_unique(h(U,V,D), pu(h(U,V,D),V,L,1), 1) :-
    bi_edge(U,V), lag_expanded(U,V,D),
    vlen(U,V,L), L > 1.

dir_unique(pu(h(U,V,D),V,L,K), pu(h(U,V,D),V,L,K+1), 1) :-
    bi_edge(U,V), lag_expanded(U,V,D),
    vlen(U,V,L), L > 2, K = 1..L-2.

dir_unique(pu(h(U,V,D),V,L,L-1), V, 1) :-
    bi_edge(U,V), lag_expanded(U,V,D),
    vlen(U,V,L), L > 1.

% ================================================================
% Bidirected loops: attach loopset to the hidden parent(s)
% ================================================================
% Collect all parents introduced for the bidirected pair:
parent_bi(U,V,h0(U,V))      :- bi_edge(U,V), bi_zero(U,V).
parent_bi(U,V, h(U,V,D))    :- bi_edge(U,V), bi_diff(U,V,D).

mmax_bi(U,V,M) :- bi_edge(U,V), M = #max { R : loop(U,V,R) }.

bcycpos(U,V,K) :- mmax_bi(U,V,M), M > 1, K = 1..M-1.
latent(bcyc(U,V,K)) :- bcycpos(U,V,K).

dir_unique(P, bcyc(U,V,1), 1)             :- parent_bi(U,V,P), mmax_bi(U,V,M), M > 1.
dir_unique(bcyc(U,V,K), bcyc(U,V,K+1), 1) :- mmax_bi(U,V,M), M > 2, K = 1..M-2.
dir_unique(bcyc(U,V,M-1), P, 1)           :- parent_bi(U,V,P), mmax_bi(U,V,M), M > 1.

dir_unique(P, P, 1)                       :- parent_bi(U,V,P), loop(U,V,1).
dir_unique(bcyc(U,V,D-1), P, 1)           :- parent_bi(U,V,P), loop(U,V,D), D > 1.

% Output node set (for visibility)
% ================================================================
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