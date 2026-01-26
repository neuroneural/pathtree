#const hhmax = 10.
hidden(1000+I) :- I=1..hhmax.

node(X) :- observed(X).
node(H) :- hidden(H).

% Guess additional unit-lag directed edges among all nodes (self loops allowed)
{ edge(U,V,1) } :- node(U), node(V).

% Exact lagset equality with provided facts out_dir/out_bi
:- out_dir(X,Y,L), not path(X,Y,L).
:- path(X,Y,L), not out_dir(X,Y,L).

:- out_bi(X,Y,D), not bidirected_lag(X,Y,D).
:- bidirected_lag(X,Y,D), not out_bi(X,Y,D).

% Node minimality: minimize number of hidden nodes that participate in any edge
active_hidden(H) :- hidden(H), edge(H,_,1).
active_hidden(H) :- hidden(H), edge(_,H,1).
#minimize { 1,H : active_hidden(H) }.

% What to output from reverse
out_edge(U,V,1) :- edge(U,V,1).
out_active(H)   :- active_hidden(H).
