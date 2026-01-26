#const maxlag = 17.
lag(0..maxlag).

% Inputs:
% observed(X).
% hidden(X).           (only needed for forward mode facts; reverse will define its own hidden ids)
% edge(U,V,1).         (edges in the candidate or given graph)

% PATH: reachability with total lag
path(U,V,1) :- edge(U,V,1), observed(U), observed(V).
ohpath(U,H,1) :- observed(U), hidden(H), edge(U,H,1).
ohpath(U,H,L2) :- ohpath(U,H1,L1), edge(H1,H,1), hidden(H1), hidden(H), L2 = L1 + 1, lag(L2).

% Close hidden chain to observed endpoint
path(U,V,L2) :- observed(U), observed(V), ohpath(U,H,L1), edge(H,V,1), hidden(H), L2 = L1 + 1, lag(L2).
% Hidden to hidden chain
hhpath(H1,H2,1) :- hidden(H1), hidden(H2), edge(H1,H2,1).
hhpath(H1,H2,L2) :- hhpath(H1,Hm,L1), edge(Hm,H2,1), hidden(Hm), hidden(H2), L2 = L1 + 1, lag(L2).

% Hidden to observed, with hidden-only intermediates
hpath(H,V,1) :- hidden(H), observed(V), edge(H,V,1).
hpath(H,V,L2) :- hidden(H), observed(V), hhpath(H,Hm,L1), edge(Hm,V,1), hidden(Hm), L2 = L1 + 1, lag(L2).
% Prune forks that have observed interior nodes on an arm
blocked(H,X) :- hidden(H), observed(Z), Z != X, hpath(H,Z,L1), path(Z,X,L2), L = L1 + L2, lag(L), hpath(H,X,L).

% FORK: oriented bidirected lag, hidden ancestor H, no observed interior nodes on either arm
bidirected_lag(X,Y,D) :- hidden(H), observed(X), observed(Y), hpath(H,X,Lx), hpath(H,Y,Ly),
    lag(Lx), lag(Ly), Lx < Ly, D = Ly - Lx, lag(D), not blocked(H,X), not blocked(H,Y).
bidirected_lag(Y,X,D) :- hidden(H), observed(X), observed(Y), hpath(H,X,Lx), hpath(H,Y,Ly),
    lag(Lx), lag(Ly), Ly < Lx, D = Lx - Ly, lag(D), not blocked(H,X), not blocked(H,Y).
bidirected_lag(X,Y,0) :- hidden(H), observed(X), observed(Y), X < Y, hpath(H,X,L), hpath(H,Y,L), lag(L),
    not blocked(H,X), not blocked(H,Y).
