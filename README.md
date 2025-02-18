# pathtree

Tools to explore dynamic causal graphs in the case of undersampled data, helping to unfold the apparent structure into the underlying truth.  

**Based on gunfolds** — a framework for latent variable dynamic causality.

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
2. **Create and activate** a Python virtual environment (optional but recommended):
   ```bash

    # Create a new virtual environment in a folder named "venv"
    python -m venv venv

    # Activate it on Linux/macOS
    source venv/bin/activate

    # On Windows (Command Prompt):
    venv\Scripts\activate

3. **Install** the required packages:

   ```bash

    pip install -r requirements.txt
Requirements typically include OR‑Tools, Sympy, Numpy, and SortedContainers.

Usage
Run the test cases to see how the framework converts graphs into PathForests, refines PathTrees, and computes symbolic delay expressions:

    ```bash

    python tests_selim.py
A typical test run will show:

The original graph.
The graph after marginalizing latent vertices.
The induced PathForest for a specific edge, with its normalized PathTree, overall delay expression, and alpha labels.
Example
Consider a simple graph:

python

graph = {
    1: {2: {1: {1}}},
    2: {2: {1: {1}}, 3: {1: {1}}},
    3: {}
}
After hiding vertex 2, the edge from 1 to 3 might be represented by a PathTree with:

A root preset of 2 (cycle‑free/base delay).
A child node representing a cycle contribution of 1*a (symbolic).
The overall delay expression then becomes a + 2.