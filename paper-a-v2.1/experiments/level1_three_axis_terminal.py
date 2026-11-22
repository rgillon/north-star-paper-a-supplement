"""Exact terminal exclusion on small three-chain shapes, over all prime values."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import prod

import sympy as sp

from level1_four_prime_symbolic_exclusion import VARIABLES
from level1_four_prime_terminal import necessary

X = VARIABLES[:3]


def below(left, right):
    return all(x <= y for x, y in zip(left, right))


def row(prime, kernel, order):
    if order == 0:
        return 1
    if order <= kernel:
        return prime**(order-1) * (prime-1)
    return -prime**kernel if order == kernel+1 else 0


def certify(shape):
    cells = list(product(*(range(a+1) for a in shape)))
    size = len(cells)
    indices = range(size)
    upper = [sum(1 << j for j, y in enumerate(cells) if below(x, y)) for x in cells]
    ideals = [0]
    for j, y in enumerate(cells):
        predecessors = sum(1 << i for i, x in enumerate(cells) if x != y and below(x, y))
        ideals += [h | (1 << j) for h in ideals if h & predecessors == predecessors]
    masks = {tuple(sorted((a & ~b, b & ~a))) for a, b in combinations(ideals, 2)
             if a & ~b and b & ~a}
    primes = [x+1 for x in X]
    rows = [[prod(row(p, h, v) for p, h, v in zip(primes, c, s)) for s in cells]
            for c in cells]

    @lru_cache(None)
    def dependencies(up, zeros):
        physical = [cells[i] for i in indices if up >> i & 1]
        matrix = sp.Matrix([[prod(1 if h == v else -1 if h+1 == v else 0
                                  for h, v in zip(c, u)) for u in physical]
                            for l, c in enumerate(cells) if zeros >> l & 1])
        return matrix.nullspace()

    counts = Counter()
    unresolved = []
    for positive, negative in sorted(masks):
        support = positive | negative
        delta = [((positive >> i) & 1)-((negative >> i) & 1) for i in indices]
        polys = [sp.expand(sum(y*v for y, v in zip(r, delta) if v)) for r in rows]
        if necessary(polys[-1]) == 1:
            continue
        counts["root"] += 1
        zeros = sum(1 << i for i, f in enumerate(polys)
                    if necessary(f) != 1 and int(f.subs(dict.fromkeys(X, 0))) % 2 == 0)
        cubes = []
        up = 0
        for i, s in enumerate(cells):
            if not delta[i]:
                continue
            up |= upper[i]
            if upper[i] & support != 1 << i:
                continue
            choices = [[(v-1, v)] if v else [(h,) for h in range(axis+1)]
                       for v, axis in zip(s, shape)]
            for box in product(*choices):
                cubes.append(sum(1 << cells.index(c) for c in product(*box)))
        if not all(c & zeros for c in cubes):
            continue
        counts["cube"] += 1
        values = []
        for i, s in enumerate(cells):
            if up >> i & 1:
                contrast = sum((-1)**sum(y-x for x, y in zip(s, v)) * abs(delta[j])
                               for j, v in enumerate(cells)
                               if all(0 <= y-x <= 1 for x, y in zip(s, v)))
                values.append(contrast * prod(p**v for p, v in zip(primes, s)))

        def exclude(possible, required):
            if not all(c & possible for c in cubes):
                return True
            equations = [necessary(f) for i, f in enumerate(polys) if required >> i & 1]
            equations.extend(necessary(sum(x*y for x, y in zip(vector, values)))
                             for vector in dependencies(up, possible))
            if 1 in equations:
                return True
            basis = sp.groebner([f for f in equations if f != 0], *X,
                                order="grevlex", domain=sp.QQ)
            if any(necessary(f.as_expr()) == 1 for f in basis.polys):
                return True
            undecided = possible & ~required
            if not undecided:
                return False
            bit = undecided & -undecided
            return exclude(possible ^ bit, required) and exclude(possible, required | bit)

        if not exclude(zeros, 1 << (size-1)):
            unresolved.append((positive, negative))
    print(shape, len(ideals), len(masks), dict(counts), unresolved, flush=True)
    return unresolved


if __name__ == "__main__":
    from itertools import permutations

    for ordered_shape in sorted({p for s in [(2, 1, 1), (3, 1, 1), (2, 2, 1)]
                                 for p in permutations(s)}):
        if certify(ordered_shape):
            raise RuntimeError(f"Unresolved differences for {ordered_shape}")

