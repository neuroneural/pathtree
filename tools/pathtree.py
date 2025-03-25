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

    def __gt__(self, other):
        if isinstance(other, InfiniteExpression):
            return self.base > other.base
        elif isinstance(other, (int, float)):
            return self.base > other
        else:
            raise TypeError(f"Cannot compare InfiniteExpression with {type(other)}")
    
    def __lt__(self, other):
        if isinstance(other, InfiniteExpression):
            return self.base < other.base
        elif isinstance(other, (int, float)):
            return self.base < other
        else:
            raise TypeError(f"Cannot compare InfiniteExpression with {type(other)}")
    
    def __eq__(self, other):
        if isinstance(other, InfiniteExpression):
            return self.base == other.base and self.terms == other.terms
        elif isinstance(other, (int, float)):
            return self.base == other and not self.terms
        else:
            return False

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
                elif isinstance(x, str):
                    raise TypeError(f"osumset() received an unexpected string: {x}")
                else:
                    raise TypeError(f"osumset() received an unexpected type: {type(x)} ({x})")

    # Return the final expression as an InfiniteExpression object
    final_expr = InfiniteExpression(finite_sum) + infinite_expression
    return {final_expr}  # Return as an object, not a string

def osumset_full(*sets):
    # Start with a result that contains only 0.
    result = {0}
    
    # Helper function to merge two elements.
    def merge_two(e1, e2):
        if isinstance(e1, int) and isinstance(e2, int):
            return e1 + e2
        # If you have InfiniteExpression types, handle them accordingly.
        elif isinstance(e1, int) and isinstance(e2, InfiniteExpression):
            return e2 + e1
        elif isinstance(e1, InfiniteExpression) and isinstance(e2, int):
            return e1 + e2
        elif isinstance(e1, InfiniteExpression) and isinstance(e2, InfiniteExpression):
            return e1 + e2
        else:
            raise TypeError(f"osumset_full can only handle ints or InfiniteExpressions, got {e1} + {e2}")
    
    for s in sets:
        # Ensure s is a set.
        if isinstance(s, (int, InfiniteExpression)):
            s = {s}
        elif not isinstance(s, set):
            raise TypeError(f"osumset_full expects sets, ints, or InfiniteExpressions, got {type(s)}: {s}")
        
        # Always do the Cartesian sum.
        new_result = set()
        for r_elem in result:
            for s_elem in s:
                merged_val = merge_two(r_elem, s_elem)
                new_result.add(merged_val)
        result = new_result

    return result
class PathTree:
    """
    A simpler node structure where:
      - preset: the base delay (int or set of ints).
      - loopset: a set containing either:
          (a) integers (the minimal single-pass delay of a loop), or
          (b) other PathTree objects (nested loops).
    """
    def __init__(self, preset=0, loopset=None):
        self.preset = preset
        # loopset can contain ints or PathTree objects
        self.loopset = loopset if loopset is not None else set()
    
    def add_loop(self, loop_delay):
        """
        Add a single-pass delay integer to loopset.
        """
        self.loopset.add(loop_delay)

    def add_nested_loop(self, nested_loop):
        """
        Add a nested PathTree to loopset.
        """
        if not isinstance(nested_loop, PathTree):
            raise TypeError("Nested loop must be a PathTree.")
        self.loopset.add(nested_loop)

    def __repr__(self):
        # If preset is a set with exactly one element, use that element.
        if isinstance(self.preset, set) and len(self.preset) == 1:
            preset_str = str(next(iter(self.preset)))
        else:
            preset_str = str(self.preset)
        if self.loopset:
            loopset_str = "<" + ', '.join([str(x) for x in self.loopset]) + ">"
            return f"({preset_str}, {loopset_str})"
        else:
            return f"({preset_str})"
    
    # def __repr__(self):
    #     return "(" + str(self.preset) + ", <" + ', '.join([str(x) for x in self.loopset])+">)" if self.loopset else "("+str(self.preset) + ")"
        #return f"PathTree(preset={self.preset}, loopset={self.loopset})"
        #return f"({self.preset}, <{self.loopset}>)" if self.loopset else f"({self.preset})"
def assign_alpha_labels(pt, label_map=None, next_label=None):
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
    # If the preset is a set, convert it to a frozenset so that it's hashable.
    base = frozenset(pt.preset) if isinstance(pt.preset, set) else pt.preset
    signature = (base, child_signatures)
    pt.signature = signature

    if signature in label_map:
        pt.label = label_map[signature]
    else:
        pt.label = next_label[0]
        label_map[signature] = next_label[0]
        next_label[0] += 1

    return signature