"""Exact four-prime terminal Fourier test. Run from the repository root."""
from collections import Counter
from functools import lru_cache
from math import prod

import sympy as sp

from level1_four_prime_symbolic_exclusion import (
    CELLS, VARIABLES, arm_candidates, factors, necessary_factor,
)


def prime_obstruction(polynomial):
    for i, variable in enumerate(VARIABLES):
        residue = polynomial.subs(variable, -1)
        if residue == 0:
            continue
        content, terms = sp.factor_list(residue)
        units = {x+1 for x in VARIABLES if x != variable} | set(VARIABLES[:i])
        if all(f in units for f, _ in terms) and all(p <= 2 for p in sp.factorint(abs(content))):
            return True
    return False


@lru_cache(None)
def necessary(poly):
    reduced = necessary_factor(poly)
    if reduced == 0:
        return 0
    return prod(f for f, _ in sp.factor_list(reduced)[1] if not prime_obstruction(f))


def inconsistent(equations):
    retained = [necessary(f) for f in equations]
    if 1 in retained:
        return True
    basis = sp.groebner([f for f in retained if f != 0], *VARIABLES,
                        domain=sp.QQ, order="grevlex")
    return any(necessary(row.as_expr()) == 1 for row in basis.polys)


@lru_cache(None)
def dependencies(missing, zero_mask):
    rows = [[int(c & u == u) for u in missing] for c in CELLS if zero_mask >> c & 1]
    return sp.Matrix(rows).nullspace()


def exclude(positive, negative, zeros, cubes, polynomial):
    support = positive | negative
    missing = tuple(u for u in CELLS if any(support >> s & 1 and not u & s for s in CELLS))
    contrasts = [sum((-1)**j.bit_count() * (support >> (15 ^ j) & 1)
                     for j in CELLS if j & u == j) for u in missing]
    weights = [prod(1+x for i, x in enumerate(VARIABLES) if not u >> i & 1)
               for u in missing]
    rows = [polynomial.subs({x: -1 for i, x in enumerate(VARIABLES) if c >> i & 1})
            for c in CELLS]

    def descend(possible, required):
        if not all(possible & cube for cube in cubes):
            return True
        equations = [polynomial] + [rows[i] for i in CELLS if required >> i & 1]
        for vector in dependencies(missing, possible):
            equations.append(sum(a*b*w for a, b, w in zip(vector, contrasts, weights)))
        if inconsistent(equations):
            return True
        undecided = possible & ~required
        if not undecided:
            return False
        bit = undecided & -undecided
        return descend(possible ^ bit, required) and descend(possible, required | bit)

    return descend(zeros, 1)


if __name__ == "__main__":
    _, arms = arm_candidates()
    counts = Counter()
    for positive, negative, zeros, cubes, _, polynomial in arms:
        if len(polynomial.free_symbols) < 4 or len(factors(positive, negative)) != 1:
            counts["smaller vanishing block"] += 1
        elif exclude(positive, negative, zeros, cubes, polynomial):
            counts["four-variable exclusion"] += 1
        else:
            raise AssertionError(("unresolved", positive, negative))
    print(dict(counts))
    assert counts == {"smaller vanishing block": 44, "four-variable exclusion": 46}

