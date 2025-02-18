class InfiniteSet:
    def __init__(self, base, step):
        self.base = base  # Smallest value in the series
        self.step = step  # Incremental step (for cycles)

    def __contains__(self, value):
        """Check if a value is in the infinite set."""
        if value < self.base:
            return False
        return (value - self.base) % self.step == 0

    def __repr__(self):
        return f"{{ {self.base} + {self.step}k | k ≥ 0 }}"

    def union(self, other):
        """Union with another InfiniteSet or finite set."""
        if isinstance(other, InfiniteSet):
            return self._union_with_inf(other)
        else:
            return {self.base + k * self.step for k in range(10)}.union(other)

    def _union_with_inf(self, other):
        """Handles the union of two InfiniteSets (approximate)."""
        if self.step % other.step == 0:
            return self
        if other.step % self.step == 0:
            return other

        return {self.base + k * self.step for k in range(10)}.union(
            {other.base + k * other.step for k in range(10)}
        )

    def intersect(self, other):
        """Intersection with another InfiniteSet (if step sizes align)."""
        if self.step % other.step == 0 or other.step % self.step == 0:
            return InfiniteSet(max(self.base, other.base), max(self.step, other.step))
        return set()  # Otherwise, intersection is difficult