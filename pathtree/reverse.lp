% ============================================================
% reverse.lp  (abduction over hide_nodes.lp)
% Given targets (target_dir/3, target_bi/3) over observed nodes,
% synthesize an input graph (edges + hidden nodes) such that
% running hide_nodes.lp reproduces exactly those targets.
% ============================================================

#const maxlag = 1000.
#const k_aux  = 6.
lag(0..maxlag).

aux(1..k_aux).
cand_node(X) :- observed(X).
cand_node(aux(I)) :- aux(I).

hidden(aux(I)) :- aux(I).

% ------------------------------------------------------------
% Targets provided as input facts:
%   observed(X).
%   target_dir(X,Y,L).
%   target_bi(X,Y,D).
% ------------------------------------------------------------

% ------------------------------------------------------------
% Mandatory edges: if target_dir(X,Y,1) then we must have edge(X,Y)
% because dir_unique(X,Y,1) in hide_nodes comes from edge(X,Y).
% ------------------------------------------------------------
edge(X,Y) :- target_dir(X,Y,1).

% ------------------------------------------------------------
% reverse_z3 idea: preserve observed self-loops if demanded
% ------------------------------------------------------------
edge(X,X) :- observed(X), target_dir(X,X,1).

% ------------------------------------------------------------
% Detect "consecutive run" directed lagsets for a pair (U,V).
% This is the common pattern like {2,3,4,...,17}.
% We use it to force a tiny latent motif that tends to generate
% all intermediate lags.
% ------------------------------------------------------------

% U,V has at least 3 lags, and they include a and a+1
has_pair(U,V) :- target_dir(U,V,_).
minlag(U,V,A) :- A = #min { L : target_dir(U,V,L) } , has_pair(U,V).
maxlagp(U,V,B) :- B = #max { L : target_dir(U,V,L) } , has_pair(U,V).

% Count how many lags are in the target for (U,V)
ctlag(U,V,N) :- N = #count { L : target_dir(U,V,L) }, has_pair(U,V).

% "looks_consecutive" is a cheap necessary condition:
% at least 3 lags, and includes min and min+1, and includes max
looks_consecutive(U,V) :-
    observed(U), observed(V), U != V,
    ctlag(U,V,N), N >= 3,
    minlag(U,V,A),
    target_dir(U,V,A),
    A1 = A + 1,
    target_dir(U,V,A1),
    maxlagp(U,V,B),
    target_dir(U,V,B).

% ------------------------------------------------------------
% Motif forcing: two hidden nodes forming a 2-cycle, plus U->both, H0->V.
% This mirrors the successful reverse_z3 construction for consecutive lags.
% We use aux(1) and aux(2) as the motif latents.
%
% If you do not want this heuristic, comment out this block.
% ------------------------------------------------------------
edge(U,aux(1)) :- looks_consecutive(U,V).
edge(U,aux(2)) :- looks_consecutive(U,V).
edge(aux(1),aux(2)) :- looks_consecutive(U,V).
edge(aux(2),aux(1)) :- looks_consecutive(U,V).
edge(aux(1),V) :- looks_consecutive(U,V).

% ------------------------------------------------------------
% Candidate edge generation space.
% We keep it structured:
% - observed -> observed edges are only allowed if target_dir(...,1)
%   or if we really need them (you can relax this later).
% - observed -> hidden and hidden -> observed are allowed.
% - hidden -> hidden are allowed.
% ------------------------------------------------------------

% Allow observed->observed only if it is explicitly required at lag 1
{ edge(U,V) } :-
    observed(U), observed(V), U != V,
    target_dir(U,V,1).

% Allow observed->hidden and hidden->observed freely (search space)
{ edge(U,H) } :- observed(U), hidden(H).
{ edge(H,V) } :- hidden(H), observed(V).

% Allow hidden->hidden edges (including cycles)
{ edge(H1,H2) } :- hidden(H1), hidden(H2), H1 != H2.
{ edge(H,H) } :- hidden(H).   % allow hidden self-loops, useful for periodicity

% Optional: forbid edges into observed from observed unless target_dir(...,1)
% Keep this ON initially to avoid huge search.
:- edge(U,V), observed(U), observed(V), U != V, not target_dir(U,V,1).

% ------------------------------------------------------------
% Include forward semantics (hide_nodes)
% ------------------------------------------------------------
#include "hide_nodes.lp".

% ------------------------------------------------------------
% Exact matching constraints
% ------------------------------------------------------------

% All targets must be produced
:- target_dir(X,Y,L), not dir_unique(X,Y,L).
:- target_bi(X,Y,D),  not bi_diff(X,Y,D).

% No extras allowed
:- dir_unique(X,Y,L), not target_dir(X,Y,L).
:- bi_diff(X,Y,D),    not target_bi(X,Y,D).

% ------------------------------------------------------------
% Used auxiliaries and edges
% ------------------------------------------------------------
used_hidden(H) :- hidden(H), edge(H,_).
used_hidden(H) :- hidden(H), edge(_,H).
used_edge(U,V) :- edge(U,V).

% ------------------------------------------------------------
% Optimization
% ------------------------------------------------------------
#minimize { 1,H : used_hidden(H) }.
#minimize { 1,U,V : used_edge(U,V) }.

% ------------------------------------------------------------
% Output
% ------------------------------------------------------------
#show edge/2.
#show used_hidden/1.
#show observed/1.
#show hidden/1.