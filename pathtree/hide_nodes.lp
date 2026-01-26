#const maxlag = 17.
#const hhmax  = 10.
#const inverse = 0.  % 0 = forward, 1 = inverse
lag(0..maxlag).
% Inputs (provided externally) are:
% observed(X).
% In forward mode: edge(U,V,1).
% In inverse mode: target_dir(X,Y,L).  target_bi(X,Y,D).
% Node domain
hidden(1000+I) :- inverse=1, I=1..hhmax.
node(X) :- observed(X).
node(H) :- hidden(H).
% Candidate edge guessing (inverse mode)
% Unit-lag, directed, self-loops allowed
{ edge(U,V,1) } :- inverse=1, node(U), node(V).

% PATH: reachability with total lag
path(U,V,1) :- edge(U,V,1).
path(U,V,L2) :- path(U,Z,L1), edge(Z,V,1), L2 = L1 + 1, lag(L2).
% Only keep observed endpoints
cand_dir(X,Y,L) :- observed(X), observed(Y), path(X,Y,L), lag(L).

% FORK: oriented bidirected lag
% (X is closer to the hidden ancestor than Y)
cand_bi(X,Y,D) :- hidden(H), observed(X), observed(Y), path(H,X,Lx), path(H,Y,Ly),
    lag(Lx), lag(Ly), Lx < Ly, D = Ly - Lx, lag(D), not blocked(H,X), not blocked(H,Y).
cand_bi(Y,X,D) :- hidden(H), observed(X), observed(Y), path(H,X,Lx), path(H,Y,Ly),
    lag(Lx), lag(Ly), Ly < Lx, D = Lx - Ly, lag(D), not blocked(H,X), not blocked(H,Y).

% Equal distance case (D = 0), canonical orientation
cand_bi(X,Y,0) :- hidden(H), observed(X), observed(Y), X < Y, path(H,X,L), path(H,Y,L), lag(L),
    not blocked(H,X), not blocked(H,Y).
blocked(H,X) :- observed(Z), Z != X, path(H,Z,L1), path(Z,X,L2), lag(L1), lag(L2).
% FORWARD OUTPUT
dir_lag(X,Y,L) :- inverse=0, cand_dir(X,Y,L).
bi_lag(X,Y,D) :- inverse=0, cand_bi(X,Y,D).

% INVERSE CONSTRAINTS
% Exact lagset equality with gset
:- inverse=1, out_dir(X,Y,L), not cand_dir(X,Y,L).
:- inverse=1, cand_dir(X,Y,L), not out_dir(X,Y,L).

:- inverse=1, out_bi(X,Y,D), not cand_bi(X,Y,D).
:- inverse=1, cand_bi(X,Y,D), not out_bi(X,Y,D).

% Node minimality. Minimize number of hidden nodes that participate in edges
active_hidden(H) :- inverse=1, hidden(H), edge(H,_,1).
active_hidden(H) :- inverse=1, hidden(H), edge(_,H,1).

#minimize { 1,H : active_hidden(H) }.

out_dir(X,Y,L) :- inverse=0, dir_lag(X,Y,L).
out_bi(X,Y,D)  :- inverse=0, bi_lag(X,Y,D).
out_edge(U,V,L) :- inverse=1, edge(U,V,L).
out_active(H)   :- inverse=1, active_hidden(H).

#show out_dir/3.
#show out_bi/3.
#show out_edge/3.
#show out_active/1.