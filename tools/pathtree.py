import sys
from infiniteset import InfiniteSet
from sympy import symbols, Add, sympify
sys.path.append('./tools/')


class InfiniteExpression:
    """Represents a sum of independent infinite terms like 20*k + 33*m."""
    def __init__(self, base=0, terms=None):
        self.base = base  # Constant term (e.g., 23 in 23 + 20*k + 33*m)
        self.terms = terms if terms is not None else {}  # Dictionary: {variable: coefficient}

    def __add__(self, other):
        if isinstance(other, int):
            return InfiniteExpression(self.base + other, self.terms)
        if isinstance(other, InfiniteExpression):
            new_base = self.base + other.base
            new_terms = self.terms.copy()
            for var, coef in other.terms.items():
                new_terms[var] = new_terms.get(var, 0) + coef
            return InfiniteExpression(new_base, new_terms)
        raise TypeError(f"Unsupported operand type for +: '{type(other)}'")

    def __repr__(self):
        # Always include the multiplication operator.
        term_parts = [f"{coef}*{var}" for var, coef in sorted(self.terms.items())]
        if not term_parts:
            return f"{self.base}"
        if self.base:
            return f"{self.base} + " + " + ".join(term_parts)
        else:
            return " + ".join(term_parts)

    def __eq__(self, other):
        return (isinstance(other, InfiniteExpression) and
                self.base == other.base and
                self.terms == other.terms)

    def __hash__(self):
        return hash((self.base, frozenset(self.terms.items())))

    def add_term(self, coefficient, variable):
        if variable in self.terms:
            self.terms[variable] += coefficient
        else:
            self.terms[variable] = coefficient
    def to_sympy(self):
        expr = sympify(self.base)
        for var, coef in self.terms.items():
            expr += coef * symbols(var)
        return expr

def osumnum(s, num):
    if isinstance(s, InfiniteSet):
        return InfiniteSet(s.base + num, s.step)  # Shift the InfiniteSet
    return set(num + x for x in s)

def osumset(*sets):
    """
    Compute the sum of multiple sets or InfiniteExpressions.
    Ensures finite numbers are summed first before merging with multiple independent InfiniteExpressions.
    """
    finite_sum = 0
    infinite_expression = InfiniteExpression(0)  # Initialize for infinite terms

    for s in sets:
        print(f"Processing set: {s}")  # Debugging step
        if isinstance(s, InfiniteExpression):
            infinite_expression += s  # Correctly sum infinite expressions
        elif isinstance(s, set):
            for x in s:
                if isinstance(x, InfiniteExpression):
                    infinite_expression += x  # Add infinite terms properly
                elif isinstance(x, int):  # Sum finite terms normally
                    finite_sum += x
                elif isinstance(x, str):  # ❌ Catch early string conversion errors
                    raise TypeError(f"osumset() received an unexpected string: {x}")
                else:
                    raise TypeError(f"osumset() received an unexpected type: {type(x)} ({x})")

    # Return the final expression as an InfiniteExpression object
    final_expr = InfiniteExpression(finite_sum) + infinite_expression
    return {final_expr}  # Return as an object, not a string

def osumset_full(*sets):
    """
    Perform a full cartesian sum on integer sets, and merge infinite expressions symbolically.
    
    If passing multiple sets:
      - Each set can contain integers or InfiniteExpressions (or be empty).
      - We do a cartesian product among them.
      - Summation of an int + int is a plain int,
        an int + InfiniteExpression is a new InfiniteExpression, 
        or an InfiniteExpression + InfiniteExpression merges them symbolically.
        
    Returns:
      A set of integers and/or InfiniteExpression objects representing all possible sums.
    """

    # Start with a result that has one element: integer 0
    # so we can sum the next set with it cartesian-style.
    result = {0}
    
    # Helper to sum two “elements,” each either an int or an InfiniteExpression.
    def merge_two(e1, e2):
        if isinstance(e1, int) and isinstance(e2, int):
            return e1 + e2
        elif isinstance(e1, int) and isinstance(e2, InfiniteExpression):
            return e2 + e1  # same as e2 + int
        elif isinstance(e1, InfiniteExpression) and isinstance(e2, int):
            return e1 + e2
        elif isinstance(e1, InfiniteExpression) and isinstance(e2, InfiniteExpression):
            return e1 + e2
        else:
            raise TypeError(f"osumset_full can only handle ints or InfiniteExpressions, got {e1} + {e2}")
    
    for s in sets:
        # If s is not a set, convert it to one so we can iterate
        # (some code might pass a single int or single InfiniteExpression)
        if isinstance(s, set) and all(isinstance(x, int) for x in s):
            result = result.union(s)
            continue  # Skip the Cartesian product below
        if isinstance(s, (int, InfiniteExpression)):
            s = {s}
        elif not isinstance(s, set):
            raise TypeError(f"osumset_full expects sets, ints, or InfiniteExpressions, got {type(s)}: {s}")

        new_result = set()
        for r_elem in result:
            for s_elem in s:
                merged_val = merge_two(r_elem, s_elem)
                new_result.add(merged_val)
        result = new_result  # move on to the next set

    return result
from sympy import sympify

class PathTree:
    def __init__(self, preset=0, label=0, children=None):
        """
        Initialize a PathTree node.
        
        :param preset: Base delay (l) for this node.
        :param label: The alpha label (default 0 means unlabeled).
        :param children: A list of child PathTree objects (for recursive nesting).
        """
        self.preset = preset          # Base delay (ideally the constant part only)
        self.label = label            # Alpha label for equivalence
        self.children = children if children is not None else []  # List of child nodes
        self.signature = None         # Will be computed from preset and children

    def add_child(self, child):
        """Add a child node (a nested PathTree)."""
        if not isinstance(child, PathTree):
            raise TypeError("Child must be a PathTree instance.")
        self.children.append(child)

    def get_overall_delay_expr(self, counter=None):
        if counter is None:
            counter = [1]  # mutable counter to ensure unique weight symbols

        if isinstance(self.preset, InfiniteExpression):
            base_expr = self.preset.to_sympy()
        else:
            base_expr = sympify(self.preset)

        terms = []
        for child in self.children:
            w = symbols(f'w{counter[0]}')
            counter[0] += 1
            child_expr = child.get_overall_delay_expr(counter)
            terms.append(w * child_expr)
        if terms:
            return base_expr + Add(*terms)
        else:
            return base_expr

    def __repr__(self):
        """
        Custom representation:
        - The root (base) shows only the constant (cycle‑free) part.
        - The children list shows any nested PathTrees (which capture cycle contributions).
        """
        from sympy import sympify
        try:
            # Convert preset to a sympy expression and split it into constant and nonconstant parts.
            expr = sympify(self.preset)
            const_part, _ = expr.as_coeff_Add()
        except Exception:
            const_part = self.preset
        return f"PathTree(base={const_part}, label={self.label}, children={self.children})"


    
def assign_alpha_labels(pt, label_map=None, next_label=None):
    """
    Recursively assign an alpha label to each node in the path tree based on its local cycle structure.
    
    The signature for a node is defined as a tuple:
         (preset, (sorted(child_signature_1, child_signature_2, ...)))
    Nodes with the same signature receive the same label.
    
    :param pt: The PathTree node (root) to process.
    :param label_map: A dictionary mapping signatures to assigned labels.
    :param next_label: A mutable counter (e.g., a list with one integer) for the next new label.
    :return: The signature for pt.
    """
    if label_map is None:
        label_map = {}
    if next_label is None:
        next_label = [1]  # mutable counter for labels

    # Recursively assign labels for children and compute their signatures.
    child_signatures = []
    for child in pt.children:
        sig = assign_alpha_labels(child, label_map, next_label)
        child_signatures.append(sig)
    
    # Sort the child signatures to ensure a canonical order.
    child_signatures = tuple(sorted(child_signatures))
    # The signature for this node is its preset plus the ordered child signatures.
    signature = (pt.preset, child_signatures)
    pt.signature = signature

    # Assign the label based on the signature.
    if signature in label_map:
        pt.label = label_map[signature]
    else:
        pt.label = next_label[0]
        label_map[signature] = next_label[0]
        next_label[0] += 1

    return signature