# from copy import deepcopy
# from bclique import bcliques
from latents import hide_nodes
from latentsb import bpts, convert_bclique
from pathtree import PathTree, assign_alpha_labels
from pathtreetools import update_edge_lags
from pprint import pprint

def test_print_normalized_path_forests():
    print("=== Test Case: Graph -> PathForest for Induced Edges ===\n")
    
    # Original graph: each edge is stored as: {edge_type: {edge-lag set}}
    graph = {
        1: {2: {1: {1}}},
        2: {2: {1: {1}}, 3: {1: {1}}},
        3: {}
    }
    print("Original Graph:")
    pprint(graph)
    
    # Hide vertex 2
    hidden_graph = hide_nodes(graph, [2])
    print("\nGraph after hiding vertex 2:")
    pprint(hidden_graph)
    print("\n--------------------------------------------------\n")
    
    # For each edge in the hidden graph, show its normalized PathForest.
    print("The path forests of the new edges:")
    for u in sorted(hidden_graph.keys()):
        for v in sorted(hidden_graph[u].keys()):
            print(f"\nEdge ({u} -> {v}):")
            # Form a b-clique for the edge as a set with one tuple.
            bc = {(u, v)}
            # Convert the b-clique to the expected tuple format.
            bc_converted = convert_bclique(bc)
            # Build the PathForest for the edge.
            pts = bpts(hidden_graph, bc)
            if not pts:
                print("  No PathTree found.")
                continue
            # For each PathTree, update edge-lags (i.e. normalize the tree).
            for pt in pts:
                normalized_graph = update_edge_lags(hidden_graph, bc_converted, pt)
                # Extract the normalized PathTree from the updated graph.
                v1 = next(iter(bc_converted[0]))
                v2 = next(iter(bc_converted[1]))
                normalized_pt = normalized_graph[v1][v2][1]
                assign_alpha_labels(normalized_pt)
                print("  Normalized PathTree:", repr(normalized_pt))
                print("  Number of children of normalized_pt:", len(normalized_pt.children))
                try:
                    expr = normalized_pt.get_overall_delay_expr()
                    print("  Overall delay expression:", expr)
                except Exception as e:
                    print("  Error computing overall delay expression:", e)

if __name__ == "__main__":
    test_print_normalized_path_forests()

# graph = {
#     1: {2: {1: {1}}},
#     2: {3: {1: {2}}},
#     3: {}
# }

# # The node we want to hide (make latent)
# latent_node = 2
# print("Original Graph:", graph)

# # res_hide = hide_node(deepcopy(graph), latent_node)
# # print("\nResult from hide_node(...):")
# # print(res_hide)
# res_sfi = sequential_forward_inference(deepcopy(graph), {latent_node})
# print("\nResult from sequential_forward_inference(...):")
# print(res_sfi)


# graph_multiple_parents = {
#     1: {4: {1: {1,2}}},
#     2: {4: {1: {3}}},
#     3: {},
#     4: {3: {1: {5}}, 4: {1: {10}}}
# }

# print("Original Graph (Multiple Parents for 4):")
# print(graph_multiple_parents)

# # Hide node 4
# res_multiple_parents = hide_node(deepcopy(graph_multiple_parents), 4)

# print("\nAfter hide_node(...), removing node 4:")
# print(res_multiple_parents)



# graph_multichild_loop = {
#     1: {2: {1: {1}}},
#     2: {
#         2: {1: {10}},  # self-loop
#         3: {1: {5}},
#         4: {1: {8}}
#     },
#     3: {},
#     4: {}
# }

# print("Original Graph (Hidden Node=2, multiple children + self-loop):")
# print(graph_multichild_loop)

# res_multichild_loop = hide_node(deepcopy(graph_multichild_loop), 2)

# print("\nAfter hide_node(...), removing node 2:")
# print(res_multichild_loop)


# graph_two_hidden = {
#     1: {2: {1: {1}}},
#     2: {
#         2: {1: {10}},  # Self-loop on 2
#         3: {1: {2}},
#         5: {1: {3}}
#     },
#     3: {
#         3: {1: {10}},  # Self-loop on 3
#         4: {1: {4}},
#         6: {1: {5}}
#     },
#     4: {},
#     5: {},
#     6: {}
# }

# print("Original Graph (Two Hidden Nodes in Sequence, 2 and 3):")
# print(graph_two_hidden)

# # First, remove node 2
# graph_after_2 = hide_node(deepcopy(graph_two_hidden), 2)
# print("\nAfter hide_node(...), removing node 2:")
# print(graph_after_2)

# # Then, remove node 3
# graph_final = hide_node(graph_after_2, 3)
# print("\nAfter hide_node(...), removing node 3:")
# print(graph_final)

# ggg = hide_nodes(graph_two_hidden, [2, 3])
# print("ggg: ", ggg)


# graph_multiple_paths = {
#     1: {2: {1: {1, 2}}},
#     4: {2: {1: {3}}},
#     6: {2: {1: {4}}},
#     2: {
#         2: {1: {10}},  # Self-loop on 2
#         3: {1: {5}},
#         5: {1: {6}},
#         7: {1: {7}}
#     },
#     3: {},
#     5: {},
#     7: {}
# }

# print("Original Graph (Multiple Paths Through Node 2):")
# print(graph_multiple_paths)

# # Remove node 2
# graph_after_removal = hide_node(deepcopy(graph_multiple_paths), 2)

# print("\nAfter hide_node(...), removing node 2:")
# print(graph_after_removal)


# graph = {
#     1: {3: {1: {5}}, 2: {1: {4}}},
#     2: {2: {1: {6}}, 3: {1: {2}}, 4: {1: {3}}},
#     3: {5: {1: {1}}},
#     4: {5: {1: {2}}},
#     5: {}
# }

# print(hide_node(graph, 2))



# # Find CRD Paths
# print("CRD Paths (1 → 6):", find_crd_paths(graph_with_cycle, 1, 5))

# # Find H-Treks
# print("H-Treks (1 → 6):", find_h_treks(graph, 4, 6))

###
# graph_directed = {
#     1: {2: {1: {1}}},       # 1 → 2
#     2: {3: {1: {1}}, 4: {1: {2}}},  # 2 → 3, 2 → 4
#     3: {},                  # 3 is a sink
#     4: {}                   # 4 is a sink
# }

# graph_bidirected = {
#     1: {2: {2: {1}}},       # 1 ↔ 2
#     2: {1: {2: {1}}, 3: {2: {1}}},  # 2 ↔ 3
#     3: {2: {2: {1}}}        # 3 ↔ 2
# }

# # Test Case 1: Directed Edges Update
# print("Original Directed Graph:")
# print(graph_directed)
# updated_directed = update_directed_edges(graph_directed, 2)
# print("\nUpdated Directed Graph After Removing Node 2:")
# print(updated_directed)

# # Test Case 2: Bi-Directed Edges Update
# print("\nOriginal Bi-Directed Graph:")
# print(graph_bidirected)
# updated_bidirected = update_bidirected_edges(graph_bidirected, 2)
# print("\nUpdated Bi-Directed Graph After Removing Node 2:")
# print(updated_bidirected)

###
# graph_directed = {
#     1: {2: {1: {1, 2, 3}}},       # 1 → 2 with lags {1, 2, 3}
#     2: {3: {1: {4, 5}}, 4: {1: {6}}},  # 2 → 3 with lags {4, 5}, 2 → 4 with lag {6}
#     3: {},                        # 3 is a sink
#     4: {}                         # 4 is a sink
# }

# # Test update_directed_edges for node 2
# print("Original Directed Graph:")
# print(graph_directed)
# updated_directed = update_directed_edges(graph_directed, 2)
# print("\nUpdated Directed Graph After Removing Node 2:")
# print(updated_directed)

###
# graph_bidirected = {
#     1: {2: {2: {1, 2, 3}}},       # 1 ↔ 2 with lags {1, 2, 3}
#     2: {1: {2: {1, 2, 3}}, 3: {2: {4, 5}}},  # 2 ↔ 3 with lags {4, 5}
#     3: {2: {2: {4, 5}}}           # 3 ↔ 2 with lags {4, 5}
# }

# # Test update_bidirected_edges for node 2
# print("\nOriginal Bi-Directed Graph:")
# print(graph_bidirected)
# updated_bidirected = update_bidirected_edges(graph_bidirected, 2)
# print("\nUpdated Bi-Directed Graph After Removing Node 2:")
# print(updated_bidirected)\

###
# from itertools import product
# def test_single_edge_pathtree():
#     # Construct a simple graph with one edge from '1' to '2'.
#     # For this edge, we assume the observed edge-lag set is {3, 5}.
#     # That is, the full graph G would have g['1']['2'] = {1: {3, 5}}.
#     graph = {
#        '1': {'2': {1: {3, 5}}},
#        '2': {}
#     }
    
#     # The corresponding b-clique is a pair of sets: V1 = {'1'}, V2 = {'2'}.
#     V1 = {'1'}
#     V2 = {'2'}
#     bc = list(product(V1, V2))
    
#     # Now, call our bpts function (which builds PathTrees for the b-clique)
#     pts = bpts(graph, bc)
    
#     print("PathTree(s) for edge (1,2):")
#     for pt in pts:
#         assign_alpha_labels(pt)
#         print(pt)
#         print("Overall delay expression:")
#         print(pt.get_overall_delay_expr())
    

# if __name__ == "__main__":
#     test_single_edge_pathtree()

###
# def test_complex_pathtree():
#     # Create the root node with preset = 10.
#     root = PathTree(preset=10)
    
#     # Create Child A: a leaf node with preset = 3.
#     A = PathTree(preset=3)
    
#     # Create Child B: preset = 5, with two children:
#     #   B1: preset = 2 (a leaf)
#     #   B2: preset = 4 (a leaf)
#     B = PathTree(preset=5)
#     B1 = PathTree(preset=2)
#     B2 = PathTree(preset=4)
#     B.add_child(B1)
#     B.add_child(B2)
    
#     # Create Child C: preset = 7, with two children:
#     #   C1: preset = 2 (a leaf)
#     #   C2: preset = 3 (a leaf)
#     C = PathTree(preset=7)
#     C1 = PathTree(preset=2)
#     C2 = PathTree(preset=3)
#     C.add_child(C1)
#     C.add_child(C2)
    
#     # Create Child D: preset = 8, with one child D1;
#     # D1 has one child D1a with preset = 2.
#     D = PathTree(preset=8)
#     D1 = PathTree(preset=1)
#     D1a = PathTree(preset=2)
#     D1.add_child(D1a)
#     D.add_child(D1)
    
#     # Create Child E: preset = 9, with two children:
#     #   E1: preset = 2 (a leaf)
#     #   E2: preset = 1 (a leaf)
#     E = PathTree(preset=9)
#     E1 = PathTree(preset=2)
#     E2 = PathTree(preset=1)
#     E.add_child(E1)
#     E.add_child(E2)
    
#     # Attach children A, B, C, D, E to the root.
#     root.add_child(A)
#     root.add_child(B)
#     root.add_child(C)
#     root.add_child(D)
#     root.add_child(E)
    
#     # Now assign alpha labels recursively.
#     assign_alpha_labels(root)
    
#     # Print the PathTree structure.
#     print("Complex PathTree structure:")
#     print(root)
    
#     # Print the overall delay expression.
#     print("\nOverall delay expression:")
#     expr = root.get_overall_delay_expr()
#     pprint(expr)
    
#     # Print the assigned alpha labels for the root and its children.
#     print("\nAlpha labels:")
#     print("Root label:", root.label)
#     for idx, child in enumerate(root.children, start=1):
#         print(f"Child {idx} label:", child.label)
#         if child.children:
#             for j, gc in enumerate(child.children, start=1):
#                 print(f"  Child {idx}.{j} label:", gc.label)

# if __name__ == "__main__":
#     test_complex_pathtree()

#!/usr/bin/env python
# from sympy import pprint
# from pathtree import PathTree, assign_alpha_labels
# from latentsb import bpts, learn_path_tree
# from pathtreetools import ptelements, update_edge_lags

# def test_graph_to_pathtree():
#     print("=== Test Case 1: Graph -> PathTree ===")
#     # Create a simple graph with a single edge (from '1' to '2')
#     # and an observed edge-lag set of {3, 5}.
#     graph = {
#         '1': {'2': {1: {3, 5}}},
#         '2': {}
#     }
#     # For a single edge, the b-clique is represented as two sets:
#     # V1 = {'1'} and V2 = {'2'}.
#     bc = {("1", "2")}
    
#     # Use our bpts function to build the PathTree(s) (a PathForest)
#     pts = bpts(graph, bc)
    
#     print("PathForest for edge (1,2):")
#     for pt in pts:
#         # Now assign alpha labels to every node in the PathTree.
#         assign_alpha_labels(pt)
#         print(pt)
#         print("Overall delay expression:")
#         pprint(pt.get_overall_delay_expr())
#         print("-" * 40)

# def test_pathtree_to_graph():
#     print("=== Test Case 2: PathTree -> Graph ===")
#     # Manually construct a PathTree for an edge.
#     # Suppose the observed edge-lag set is {3, 5}. One way to represent this is:
#     # a flat tree with preset=3 (covering delay 3) extended by a child with preset=2
#     # so that overall delay = 3 + 2*1 = 5.
#     pt = PathTree(preset=3)
#     child = PathTree(preset=2)
#     pt.add_child(child)
#     assign_alpha_labels(pt)
    
#     print("Constructed PathTree for edge (1,2):")
#     print(pt)
#     print("Overall delay expression:")
#     pprint(pt.get_overall_delay_expr())
    
#     # Now update a graph from this PathTree.
#     # We simulate the conversion by taking a graph with an edge (1,2) and
#     # replacing its edge-lag set with the children from the PathTree.
#     graph = {
#         '1': {'2': {1: set()}},
#         '2': {}
#     }
#     # For a single edge, we can represent the b-clique as V1={'1'}, V2={'2'}.
#     bc = ({'1'}, {'2'})
#     update_edge_lags(graph, bc, pt)
    
#     print("Graph after updating edge-lag set from the PathTree:")
#     print(graph)

# def test_edge_lag_to_pathtree():
#     print("=== Test Case 3: Edge-lag list -> PathTree ===")
#     # Given an observed edge-lag set, e.g., {3, 5}, build a PathTree.
#     edge_lag_set = {3, 5}
#     # Start by initializing a flat PathTree with the smallest element.
#     pt = learn_path_tree(PathTree(preset=min(edge_lag_set)), edge_lag_set)
#     assign_alpha_labels(pt)
    
#     print("PathTree built from observed edge-lag set {3, 5}:")
#     print(pt)
#     print("Overall delay expression:")
#     pprint(pt.get_overall_delay_expr())

# def test_pathtree_to_edge_lag():
#     print("=== Test Case 4: PathTree -> Edge-lag list ===")
#     # Build a PathTree from an observed edge-lag set {3, 5} as in Test Case 3.
#     edge_lag_set = {3, 5}
#     pt = learn_path_tree(PathTree(preset=min(edge_lag_set)), edge_lag_set)
#     assign_alpha_labels(pt)
    
#     # Now, using ptelements, extract a list of delay values that this PathTree can generate.
#     delays = ptelements(pt, seqlen=5, verbose=True)
    
#     print("Edge-lag values computed from the PathTree:")
#     print(delays)

# if __name__ == "__main__":
#     test_graph_to_pathtree()
#     print("\n" + "=" * 50 + "\n")
#     test_pathtree_to_graph()
#     print("\n" + "=" * 50 + "\n")
#     test_edge_lag_to_pathtree()
#     print("\n" + "=" * 50 + "\n")
#     test_pathtree_to_edge_lag()

###
# def test_nested_pathtree():
#     print("Testing Nested PathTree Construction\n" + "="*50)
    
#     # Example 1: A simple (flat) node with no children.
#     pt1 = PathTree(preset=5, label=0)
#     print("Example 1: One Node (no children)")
#     print("PathTree:", pt1)
#     print("Overall delay expression:")
#     pprint(pt1.get_overall_delay_expr())
#     print("\n" + "-"*50 + "\n")
    
#     # Example 2: A node with one child.
#     # This corresponds to a two-node structure: the root and one cycle (child node).
#     root2 = PathTree(preset=3, label=0)
#     child2 = PathTree(preset=4, label=0)  # A cycle with delay 4.
#     root2.add_child(child2)
#     print("Example 2: Two Nodes (one child attached to the root)")
#     print("PathTree:", root2)
#     print("Overall delay expression:")
#     pprint(root2.get_overall_delay_expr())
#     print("\n" + "-"*50 + "\n")
    
#     # Example 3: A nested structure (three nodes in a chain).
#     # The root has a child, and that child has its own child.
#     root3 = PathTree(preset=2, label=0)
#     child3 = PathTree(preset=3, label=0)
#     nested_child3 = PathTree(preset=4, label=0)
#     child3.add_child(nested_child3)
#     root3.add_child(child3)
#     print("Example 3: Three Nodes (nested structure)")
#     print("PathTree:", root3)
#     print("Overall delay expression:")
#     pprint(root3.get_overall_delay_expr())
#     print("\n" + "-"*50 + "\n")

# if __name__ == "__main__":
#     test_nested_pathtree()

### Two children, one has its own child
# def test_complex_pathtree():
#     # Construct the nested PathTree:
#     # Root: preset = 10
#     # Child A: preset = 4, with one child A1: preset = 2
#     # Child B: preset = 6
#     root = PathTree(preset=10, label=0)
    
#     child_A = PathTree(preset=4, label=0)
#     child_A1 = PathTree(preset=2, label=0)
#     child_A.add_child(child_A1)
    
#     child_B = PathTree(preset=6, label=0)
    
#     root.add_child(child_A)
#     root.add_child(child_B)
    
#     # Compute the overall delay symbolically.
#     overall_expr = root.get_overall_delay_expr()
    
#     print("Complex PathTree:")
#     print(root)
#     print("\nOverall delay expression:")
#     pprint(overall_expr)

# if __name__ == "__main__":
#     test_complex_pathtree()

### ----- Test Cases for Nested PathTree with Alpha Labels -----

# def test_nested_pathtree_with_alpha():
#     print("Testing Nested PathTree with Alpha Labels\n" + "="*50)

#     # Example 1: A simple node with no children.
#     pt1 = PathTree(preset=5)
#     assign_alpha_labels(pt1)
#     print("Example 1: One Node (no children)")
#     print("PathTree:", pt1)
#     print("Overall delay expression:")
#     pprint(pt1.get_overall_delay_expr())
#     print("Alpha label:", pt1.label)
#     print("\n" + "-"*50 + "\n")
    
#     # Example 2: Two nodes: a root with one child.
#     root2 = PathTree(preset=3)
#     child2 = PathTree(preset=4)
#     root2.add_child(child2)
#     assign_alpha_labels(root2)
#     print("Example 2: Two Nodes (root with one child)")
#     print("PathTree:", root2)
#     print("Overall delay expression:")
#     pprint(root2.get_overall_delay_expr())
#     print("Alpha labels: Root =", root2.label, ", Child =", child2.label)
#     print("\n" + "-"*50 + "\n")
    
#     # Example 3: Three nodes (nested structure): root -> child -> nested child.
#     root3 = PathTree(preset=2)
#     child3 = PathTree(preset=3)
#     nested_child3 = PathTree(preset=4)
#     child3.add_child(nested_child3)
#     root3.add_child(child3)
#     assign_alpha_labels(root3)
#     print("Example 3: Three Nodes (nested structure)")
#     print("PathTree:", root3)
#     print("Overall delay expression:")
#     pprint(root3.get_overall_delay_expr())
#     print("Alpha labels: Root =", root3.label,
#           ", Child =", child3.label,
#           ", Nested Child =", nested_child3.label)
#     print("\n" + "-"*50 + "\n")
    
#     # Example 4: More complex structure with two children at the root.
#     root4 = PathTree(preset=10)
#     child_A = PathTree(preset=4)
#     child_A1 = PathTree(preset=2)
#     child_A.add_child(child_A1)
#     child_B = PathTree(preset=6)
#     root4.add_child(child_A)
#     root4.add_child(child_B)
#     assign_alpha_labels(root4)
#     print("Example 4: Complex structure with two children at the root")
#     print("PathTree:", root4)
#     print("Overall delay expression:")
#     pprint(root4.get_overall_delay_expr())
#     print("Alpha labels: Root =", root4.label,
#           ", Child A =", child_A.label,
#           ", Child A1 =", child_A1.label,
#           ", Child B =", child_B.label)
#     print("\n" + "-"*50 + "\n")

# if __name__ == "__main__":
#     test_nested_pathtree_with_alpha()

###
# def test_nodes_same_and_distinct_labels():
#     # Construct a tree with 5 nodes:
#     # Root: preset = 10
#     #  - Child1: preset = 5, no children
#     #  - Child2: preset = 5, no children (identical to Child1)
#     #  - Child3: preset = 8, with one child:
#     #         Child3a: preset = 3, no children
#     root = PathTree(preset=10)
#     child1 = PathTree(preset=5)
#     child2 = PathTree(preset=5)
#     child3 = PathTree(preset=8)
#     child3a = PathTree(preset=3)
    
#     # Build the tree structure:
#     root.add_child(child1)
#     root.add_child(child2)
#     root.add_child(child3)
#     child3.add_child(child3a)
    
#     # Recursively assign alpha labels based on local cycle structure.
#     assign_alpha_labels(root)
    
#     # Print out the tree structure and the labels of each node.
#     print("Tree structure:")
#     print(root)
#     print("\nAssigned alpha labels:")
#     print("Root label:", root.label)
#     print("Child1 label:", child1.label)
#     print("Child2 label:", child2.label)
#     print("Child3 label:", child3.label)
#     print("Child3a label:", child3a.label)

# if __name__ == "__main__":
#     test_nodes_same_and_distinct_labels()

###
# def test_hide_nodes_and_print_pathtrees():
#     # Example graph: keys are strings for nodes.
#     # Each edge is stored as: g[u][v] = {1: {lag_set}}
#     graph = {
#         '1': {'2': {1: {4}}, '3': {1: {5}}},
#         '2': {'3': {1: {2}}, '4': {1: {3}}},
#         '3': {'5': {1: {1}}},
#         '4': {'5': {1: {2}}},
#         '5': {}
#     }
#     graph2={'1': {'2': {1: {1}}},
#         '2': {'2': {1: {1}}, '3': {1: {1}}},
#         '3': {}

#     }
#     # Hide multiple nodes (e.g., hide node '2')
#     # You can provide a list of nodes to hide.
#     hidden_graph = hide_nodes(graph, ['2'])
    
#     print("Graph after hiding node '2':")
#     pprint(hidden_graph)
#     print("\n--------------------------------------------------\n")
    
#     # Extract b-cliques from the hidden graph
#     bcliques_list = bcliques(hidden_graph, verbose=True)
#     print("B-cliques in the hidden graph:")
#     for bc in bcliques_list:
#         print(bc)
#     print("\n--------------------------------------------------\n")
    
#     # For each b-clique, compute the associated PathTrees.
#     # The function bpts (in latentsb.py) returns a set of PathTree objects for the b-clique.
#     for bc in bcliques_list:
#         pts = bpts(hidden_graph, bc)
#         print("PathTree(s) for b-clique", bc, ":")
#         for pt in pts:
#             # Print the top-level representation of the PathTree.
#             # If the PathTree is nested, its __repr__ will show the children as part of it.
#             print(pt)
#         print("\n--------------------------------------------------\n")

# if __name__ == "__main__":
#     test_hide_nodes_and_print_pathtrees()


# def test_multiple_pathtrees_single_edge_manual():
#     # We want to represent the edge from '1' to '3' with an observed edge-lag set {5, 6}.
#     #
#     # Option 1: A simple (flat) PathTree that directly represents delay 5.
#     T1 = PathTree(preset=5)
    
#     # Option 2: A nested PathTree that represents delay 6.
#     # For instance, let the root have preset=2, and it has one child with preset=4.
#     # Then the overall delay is 2 + 4 = 6.
#     T2 = PathTree(preset=2)
#     child = PathTree(preset=4)
#     T2.add_child(child)
    
#     # Combine these into a path-forest for the edge (1,3).
#     path_forest = {T1, T2}
    
#     print("PathForest for edge (1,3):")
#     for pt in path_forest:
#         print(pt)
#         print("Overall delay expression:")
#         # Assuming get_overall_delay_expr() returns a symbolic expression,
#         # we call it and print its result.
#         expr = pt.get_overall_delay_expr()
#         pprint(sympify(expr))
#         print("-" * 40)

# if __name__ == "__main__":
#     test_multiple_pathtrees_single_edge_manual()