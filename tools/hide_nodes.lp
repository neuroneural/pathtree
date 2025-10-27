% Graph marginalization with lags

#const maxlag=1000.
lag(0..maxlag).
% ----- BASIC SANITY (optional) ---
:- observed(X), hidden(X).
:- node(X), not observed(X), not hidden(X).

% Paths whose intermediate nodes are hidden
path_h(Z,Y,1) :-
  hidden(Z), observed(Y), edge(Z,Y).
path_h(Z,Y,L) :-
  lag(L), L>1,
  hidden(Z), edge(Z,Z2), hidden(Z2),
  path_h(Z2,Y,L-1).

pathvv(X,Y,1) :-
  observed(X), observed(Y), edge(X,Y), not hidden(X), not hidden(Y).
pathvv(X,Y,L) :-
  lag(L), L>1,
  observed(X), edge(X,Z), hidden(Z),
  path_h(Z,Y,L-1).

% Directed edges in G' and lag set
dir_edge(X,Y) :-
  pathvv(X,Y,_).
dir_lag(X,Y,L) :-
  pathvv(X,Y,L).

% Base_h: loop-free hidden→observed distance counting each edge exactly once
base_h(H,Y,1) :-
    hidden(H), observed(Y), edge(H,Y).

base_h(H,Y,L) :-
    hidden(H), edge(H,Z), hidden(Z), H != Z,
    base_h(Z,Y,L2),
    L = L2 + 1,
    lag(L).

% Bridge base_h results into path_h for hidden ancestors with loops
path_h(H,Y,L) :-
    base_h(H,Y,L),
    has_loop(H).

% Also inherit one level up: if H1->H2 and H2 has loop, propagate distance
path_h(H1,Y,L) :-
    hidden(H1), edge(H1,H2), hidden(H2),
    path_h(H2,Y,L2),
    L = L2 + 1,
    lag(L).

% Find the first looped hidden ancestor (source of periodicity)
% ---------------------------------
% Case 1: the hidden node itself has a loop
looped_ancestor(H0,Y,L) :-
  hidden(H0), has_loop(H0),
  base_h(H0,Y,L).
% Periodic expansion from looped ancestor (adds multiples of loop distance)
looped_ancestor(H0,Y,L) :-
  hidden(H1), edge(H1,H2), hidden(H2),
  base_h(H2,Y,L2),
  not has_loop(H1),
  looped_ancestor(H0,H1,L3),
  L = L2 + L3.

% true self-loop (edge to itself)
path_hidden(H,H,1) :-
    hidden(H), edge(H,H).

% direct hidden→hidden link (no self-loop)
path_hidden(H1,H2,1) :-
    hidden(H1), hidden(H2),
    edge(H1,H2), H1 != H2.

% recursive hidden→hidden path (forward only, no observed reachability)
path_hidden(H1,H3,L) :-
    lag(L), L>1, L<maxlag,
    hidden(H1), edge(H1,H2), hidden(H2),
    hidden(H3),
    path_hidden(H2,H3,L2),
    L = L2 + 1.

% loop length: only if hidden node truly closes within hidden subgraph
loop_len(H,L) :-
    hidden(H),
    path_hidden(H,H,L).

has_loop(H) :- loop_len(H,_).

% Bidirected edges (common hidden ancestor)
% ---------------------------------
bi_lag(X,Y,Lx,Ly) :-
  observed(X), observed(Y), X != Y,
  hidden(H),
  path_h(H,X,Lx),
  path_h(H,Y,Ly),
  Lx <= Ly.

bi_edge(X,Y) :- bi_lag(X,Y,_,_).

% absolute difference form
bi_diff(X,Y,D) :-
  observed(X), observed(Y), X != Y,
  hidden(H),
  base_h(H,X,Lx0), base_h(H,Y,Ly0),
  Ly0 >= Lx0, D = Ly0 - Lx0.

bi_diff(X,Y,D2) :-
  observed(X), observed(Y), X != Y,
  looped_ancestor(H0,X,Lx0),
  looped_ancestor(H0,Y,Ly0),
  has_loop(H0), loop_len(H0,Lloop), lag(Lloop),
  lag(K), K>0,
  Lx0 <= Ly0,
  D0 = Ly0 - Lx0,
  D2 = D0 + K * Lloop,
  D2 <= maxlag, lag(D2).

bi_zero(X,Y)   :- bi_diff(X,Y,0).

dir_unique(X,Y,L) :- dir_lag(X,Y,L).
bi_unique(X,Y,Lx,Ly) :- bi_lag(X,Y,Lx,Ly).

% ================================
% Bi-directed self-loops (A <-> A)
% ================================

% --- Helpers: observed -> hidden distances (mirror of base_h/path_h) ---
base_vh(X,H,1) :-
  observed(X), hidden(H), edge(X,H).

% one observed->hidden step then purely hidden steps to H
base_vh(X,H,L) :-
  lag(L), L>1,
  observed(X), edge(X,Z), hidden(Z),
  path_hidden(Z,H,L2),
  L = L2 + 1.

% observed->hidden paths that can profit from loops upstream
path_vh(X,H,L) :- base_vh(X,H,L).
path_vh(X,H,L) :-
  base_vh(X,H0,L0),
  H = H0,
  loop_len(H0,Lloop),
  lag(K), K>0,
  L = L0 + K * Lloop, lag(L).

% -------------------------------
% CASE 1: R -> A on both ends
% -------------------------------
% differences of two (possibly loop-augmented) R->A paths

% base differences without loops
bi_self_base(A,D0) :-
  observed(A),
  hidden(R),
  base_h(R,A,Lx), base_h(R,A,Ly),
  Ly >= Lx,
  D0 = Ly - Lx.

% close under multiples of R's loop length
bi_self_diff(A,D) :-
  bi_self_base(A,D0),
  hidden(R), has_loop(R), loop_len(R,Lloop),
  lag(N), N >= 0,
  D = D0 + N * Lloop, lag(D).

% if R has no loop, we still keep the D0
bi_self_diff(A,D0) :-
  bi_self_base(A,D0),
  hidden(R), no_loop(R).

% mark the presence of a bidirected self-loop if any D exists
bi_self_edge(A) :- bi_self_diff(A,_).

% -------------------------------
% CASE 2 (optional): A <-> R on one end, R -> A on the other
% -------------------------------
% Provide this fact in input if the compressed pre-marginal graph had it:
%    bi_obs_hidden(A,R).
%
% Then the induced self-bidir lags include differences between (A<->R) and (R->A),
% also closed under R's loop multiples, with the usual no-collider restriction
% (we never allow A<->R on both sides).

bi_self_mixed(A,D) :-
  observed(A), hidden(R),
  bi_obs_hidden(A,R),           % input fact if applicable
  base_h(R,A,Lra),
  D = abs(Lra - 0),             % treat A<->R as contributing '0' base distance to A
  lag(D).

bi_self_mixed(A,D2) :-
  observed(A), hidden(R),
  bi_obs_hidden(A,R),
  base_h(R,A,Lra),
  has_loop(R), loop_len(R,Lloop),
  lag(K), K>0,
  D2 = abs(Lra - K * Lloop),
  lag(D2).

% union of both sources
bi_self_all(A,D) :- bi_self_diff(A,D).
bi_self_all(A,D) :- bi_self_mixed(A,D).

% expose as "bi_unique" using X=Y arity-3 form (store as (A,A,D))
% emit only if there exists at least one nonzero lag
bi_unique(A,A,D,D) :-
    bi_self_all(A,D),
    D != 0.

bi_diff(A,A,D) :-
    bi_unique(A,A,D,D).

% convenient boolean: only true if a nonzero lag exists
bi_edge(A,A) :-
    bi_unique(A,A,_,_).


% forbid spurious unit-lag between observed nodes without a direct edge
:- dir_unique(X,Y,1), observed(X), observed(Y), not edge(X,Y).

% classify hidden nodes by loop presence
no_loop(H) :- hidden(H), not has_loop(H).

% ensure total knowledge: each hidden node is either has_loop or no_loop
:- hidden(H), not has_loop(H), not no_loop(H).
% --- D_min candidates = min(Lx, Ly) for hidden ancestors that have loops
base_min_cand(X,Y,Lx) :-
  observed(X), observed(Y), X < Y,
  hidden(H), has_loop(H),          % <--- already correct but now guaranteed true/false
  base_h(H,X,Lx), base_h(H,Y,Ly),
  Lx <= Ly, Lx > 0, Ly > 0.

base_min_cand(X,Y,Ly) :-
  observed(X), observed(Y), X < Y,
  hidden(H), has_loop(H),
  base_h(H,X,Lx), base_h(H,Y,Ly),
  Ly < Lx, Lx > 0, Ly > 0.

% Pick the minimal D among those candidates
base_min_final(X,Y,D) :-
    D = #min { D1 : base_min_cand(X,Y,D1) },
    base_min_cand(X,Y,D).

#show base_min_cand/3.
#show base_min_final/3.
#show dir_edge/2.
#show dir_unique/3.
#show bi_edge/2.
#show bi_unique/4.
#show bi_diff/3.
#show bi_zero/2.
#show base_h/3.
#show has_loop/1.