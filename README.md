# pathtree

Tools to explore dynamic causal graphs in the case of undersampled data, helping to unfold the apparent structure into the underlying truth.  

**Based on gunfolds** — a framework for latent variable dynamic causality.

---

## Dependencies

1. **igraph**  
2. **LaTeX** with TikZ (and `pdflatex`)  
3. **scipy**  
4. **networkx**  

---

## Gunfolds: Latent Variable Dynamic Causality via PathTrees

`gunfolds` (the underlying approach) is a Python project that implements a framework for **dynamic causal inference** in graphs with latent variables. The approach transforms a given graph into a set of **PathForests** (collections of **PathTrees**), which represent edge‑lag sets and the cycle structure underlying observed delays. This enables marginalizing over latent nodes and obtaining a **normalized, symbolic representation** of causal dynamics.

### Overview

- **PathTrees & PathForests**  
  Each edge delay is modeled as a **PathTree**. A PathTree node is defined as a triple `(l, α, T)` where:
  - `l` *(preset)*: Represents the constant, cycle‑free delay (base delay).
  - `α` *(label)*: An alpha label assigned based on the node’s local cycle structure.
  - `T` *(children)*: Represents the cycle contributions (if any) as child nodes.  
  A collection of PathTrees for a given edge forms a **PathForest**.

- **Latent Marginalization**  
  The project implements procedures to “hide” (marginalize) latent vertices in the original graph. The resulting graph has new edges whose delay information is represented by refined PathForests, separating cycle‑free (base) delay from cyclic (child) delay contributions.

- **Symbolic Delay Expressions**  
  Using **Sympy** and **OR‑Tools**, we compute overall delay expressions symbolically. This helps compare delays across edges and supports dynamic causality analysis.

---

## Features

- **Graph to PathForest Transformation**  
  Convert an input graph (with edge‑lag sets) into a normalized representation by marginalizing latent nodes.

- **PathTree Construction and Refinement**  
  Build and refine PathTrees using a Batch PTS algorithm; only extend a tree when an observed candidate delay is missing.

- **Alpha Label Assignment**  
  Assign alpha labels to nodes in a PathTree based on local cycle structure, facilitating equivalence checking.

- **Symbolic Delay Computation**  
  Compute overall delay expressions symbolically, separating constant (cycle‑free) delay from cyclic contributions.

---

## Installation

1. **Clone** the repository:
   ```bash
   git clone https://github.com/neuroneural/pathtree.git
   cd pathtree