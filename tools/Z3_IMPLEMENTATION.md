# Z3-Based PathTree Implementation

This document describes the pure Python implementation of the PathTree forward/reverse passes, replacing the clingo/ASP constraint solver.

## Overview

The PathTree system performs causal inference with latent variables by:
1. **Forward pass (hide_nodes)**: Marginalize hidden nodes from a graph, computing lag sets for paths between observed nodes
2. **Backward pass (reverse)**: Given a marginalized graph, find a latent structure that produces it

## Files

| File | Description |
|------|-------------|
| `hide_nodes_z3.py` | Forward pass - marginalizes hidden nodes |
| `reverse_z3.py` | Backward pass - finds latent structure via constraint satisfaction |
| `test_z3_vs_clingo.py` | Unit tests for forward pass |
| `test_roundtrip.py` | Round-trip tests (forward → reverse → forward) |

## Forward Pass: `hide_nodes_z3`

### Function Signature

```python
def hide_nodes_z3(
    graph: Dict[Any, Dict[Any, Dict[int, Set[int]]]],
    hidden: Set[Any],
    maxlag: int = 17,
    verbose: bool = False
) -> Dict[Any, Dict[Any, Dict[int, Set[int]]]]
```

### Input Format

Graph format: `{u: {v: {etype: {lags}}}}`
- `etype=1`: directed edge
- `etype=2`: bidirected edge
- `lags`: set of integer lag values

Example:
```python
graph = {
    5: {6: {1: {1}}, 7: {1: {1}}, 5: {1: {1}}},  # 5→6, 5→7, 5→5 (self-loop)
    6: {7: {1: {1}}},                              # 6→7
    7: {8: {1: {1}}, 6: {1: {1}}},                 # 7→8, 7→6 (cycle)
    8: {}
}
hidden = {6, 7}
```

### Output

Marginalized graph with only observed nodes:
```python
{
    5: {5: {1: {1}}, 8: {1: {2, 3, 4, ..., 17}}},
    8: {8: {2: {2, 4, 6, ..., 16}}}  # bidirected self-loop
}
```

### Algorithm

The implementation uses fixed-point iteration (not Z3 constraint solving) to compute:

1. **`path_h(Z, Y, L)`**: Paths from hidden node Z to observed node Y with length L
   - Base case: direct edges from hidden to observed
   - Recursive: extend through hidden-to-hidden edges
   - Self-loops: expand with multiples of loop length

2. **`pathvv(X, Y, L)`**: Paths from observed X to observed Y
   - Direct observed-to-observed edges
   - Paths through hidden subgraph: X → (hidden) → Y

3. **`base_h(H, Y, L)`**: Loop-free distances (for bidirected edges)
   - BFS avoiding node revisits

4. **Bidirected edges**: From common hidden ancestors
   - If hidden H reaches both X and Y, creates X ↔ Y
   - Difference set: `{|Ly - Lx| : H reaches X at Lx, Y at Ly}`

## Backward Pass: `reverse_z3`

### Function Signature

```python
def reverse_z3(
    observed_graph: Dict[Any, Dict[Any, Dict[int, Set[int]]]],
    maxlag: int = 17,
    max_latents: int = 5,
    verbose: bool = False
) -> Optional[Dict[Any, Dict[Any, Dict[int, Set[int]]]]]
```

### Constraint Satisfaction Approach

Given marginalized graph Y, find X such that `forward(X) = Y`.

The implementation uses **generate-and-test**:
1. Generate candidate latent structures based on target patterns
2. Test each candidate: does `hide_nodes_z3(candidate, latents) == target`?

### Heuristics for Candidate Generation

**Consecutive lag progressions** (e.g., {2,3,4,5,...}):
- Two latents H0, H1 forming a cycle
- Source connects to BOTH latents (for odd + even path lengths)
- One latent connects to target

```
     ┌──→ H0 ──┐
     │    ↑↓   │
u ───┤         ├──→ v
     │    H1   │
     └──→ ↑↓ ──┘
```

**Bidirected self-loops** (e.g., X ↔ X with diffs {2,4,6,...}):
- Hidden ancestor reaches X via multiple path lengths
- Cycle in hidden subgraph creates periodic differences

## Round-Trip Property

The key invariant:
```
forward(reverse(forward(G))) == forward(G)
```

**Note**: `reverse(forward(G)) ≠ G` in general, because marginalization is not injective—multiple latent structures can produce the same observed graph.

## Usage Example

```python
from hide_nodes_z3 import hide_nodes_z3
from reverse_z3 import reverse_z3, check_roundtrip

# Original graph with hidden nodes
graph = {
    5: {6: {1: {1}}, 7: {1: {1}}, 5: {1: {1}}},
    6: {7: {1: {1}}},
    7: {8: {1: {1}}, 6: {1: {1}}},
    8: {}
}
hidden = {6, 7}

# Forward pass
G_forward = hide_nodes_z3(graph, hidden, maxlag=17)
# Result: {5: {5: {1: {1}}, 8: {1: {2,3,4,...,17}}}, 8: {8: {2: {2,4,6,...,16}}}}

# Reverse pass (find equivalent latent structure)
G_reverse = reverse_z3(G_forward, maxlag=17)
# Result: {5: {H0: ..., H1: ...}, H0: {H1: ..., 8: ...}, H1: {H0: ...}}

# Verify round-trip
assert check_roundtrip(graph, hidden, maxlag=17)  # True
```

## Installation

```bash
# Clone the repository
git clone https://github.com/spikedoanz/pathtree.git
cd pathtree

# Install with uv (recommended)
uv pip install -e ".[dev]"

# Or with pip
pip install -e ".[dev]"
```

### Optional Dependencies

```bash
# Install clingo support (optional, for comparison with ASP implementation)
uv pip install -e ".[clingo]"

# Install everything
uv pip install -e ".[all]"
```

## Running Tests

```bash
# Run all tests with pytest
pytest

# Run with verbose output
pytest -v

# Run specific test file
pytest tools/tests/test_hide_nodes.py

# Run with coverage
pytest --cov=tools
```

## Dependencies

- Python 3.10+
- numpy>=1.24.0
- pytest>=7.0.0 (dev)

## Comparison with Clingo Implementation

| Aspect | Clingo (ASP) | Python (this implementation) |
|--------|--------------|------------------------------|
| Forward pass | Declarative rules in `hide_nodes.lp` | Fixed-point iteration |
| Reverse pass | `reverse.lp` (broken for some cases) | Constraint satisfaction via generate-and-test |
| Dependencies | Requires clingo binary | Pure Python |
| Performance | Good for complex constraints | Good for reachability problems |
