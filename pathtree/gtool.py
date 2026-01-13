import graph_tool as gt
from graph_tool import draw as gtd
import numpy as np

def lg2gt_forest(g):
    """
    Build a graph_tool Graph from our graph structure (g).
    Collect all vertices (both keys and those appearing as neighbors)
    and add one edge per occurrence.
    If an edge’s delay is stored as a list (i.e. a path-forest), add one edge per preset.
    """
    gr = gt.Graph()
    vlabel = gr.new_vertex_property("string")
    verts = {}
    # Gather all vertices from keys and neighbors.
    all_vertices = set(g.keys())
    for v in g:
        for w in g[v]:
            all_vertices.add(w)
    for v in sorted(all_vertices):
        verts[v] = gr.add_vertex()
        vlabel[verts[v]] = str(v)
    gr.vertex_properties["label"] = vlabel

    # Add edges: if the delay is a list (multiple presets), add one edge for each.
    for v in g:
        for w in g[v]:
            edge_data = g[v][w]
            if 1 in edge_data:
                data = edge_data[1]
                if isinstance(data, list):
                    for _ in data:
                        gr.add_edge(verts[v], verts[w])
                else:
                    gr.add_edge(verts[v], verts[w])
            elif 2 in edge_data:
                data = edge_data[2]
                if isinstance(data, list):
                    for _ in data:
                        gr.add_edge(verts[v], verts[w])
                else:
                    gr.add_edge(verts[v], verts[w])
    return gr, verts

def plotg_forest(g, layout='sfdp'):
    """
    Plot the graph constructed from g using graph_tool.
    This version simply draws the base graph, so if an edge's delay is a list,
    you will see multiple edges (arrows) drawn between the same vertices.
    """
    gg, verts = lg2gt_forest(g)
    pos_prop = gtd.sfdp_layout(gg)
    gtd.graph_draw(gg, pos=pos_prop,
                   vertex_text=gg.vertex_properties['label'],
                   vertex_font_size=32,
                   edge_pen_width=1,
                   edge_marker_size=15,
                   vertex_pen_width=1,
                   vertex_fill_color=[0.62109375,
                                      0.875,
                                      0.23828125,
                                      1],
                   output_size=(600,600))
