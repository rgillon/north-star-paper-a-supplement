"""Exclude exclusive-arm pairs for all strictly ordered odd prime quadruples."""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import prod

import numpy as np
import sympy as sp

from level1_four_prime_arm_sieve import factors, four_axis_ideals

N = 4
CELLS = range(1 << N)
VARIABLES = sp.symbols("x:4")
MONOMIALS = [prod(VARIABLES[i] for i in range(N) if m >> i & 1) for m in CELLS]


@lru_cache(None)
def expand_monomial(exponents):
    """Substitute x_i=2(i+1)+g_0+...+g_i, with nonnegative gaps."""
    terms = Counter({(0, 0, 0, 0): 1})
    for i, degree in enumerate(exponents):
        for _ in range(degree):
            following = Counter()
            for powers, value in terms.items():
                following[powers] += 2*(i+1)*value
                for j in range(i+1):
                    new = list(powers)
                    new[j] += 1
                    following[tuple(new)] += value
            terms = following
    return terms


def shell_rows():
    polys = [expand_monomial(tuple(m >> i & 1 for i in range(N))) for m in CELLS]
    basis = sorted(set().union(*(set(p) for p in polys)))
    assert basis[0] == (0, 0, 0, 0)
    return np.array([[[(-1)**((c & m).bit_count()) * polys[m & ~c][k]
                       for k in basis] for m in CELLS] for c in CELLS], dtype=np.int64)


def fixed_sign(row):
    return (row[0] > 0 and min(row) >= 0) or (row[0] < 0 and max(row) <= 0)


def excluded_joint_zero(a, b):
    if fixed_sign(b):
        return True
    return any(fixed_sign(a*y-b*x) for x, y in zip(a, b) if x and y)


def arm_candidates():
    ideals = four_axis_ideals()
    masks = sorted({tuple(sorted((a & ~b, b & ~a))) for a, b in combinations(ideals, 2)
                    if a & ~b and b & ~a})
    rows = shell_rows().transpose(0, 2, 1)
    arms = []
    for pos, neg in masks:
        delta = [((pos >> i) & 1) - ((neg >> i) & 1) for i in CELLS]
        support = [i for i in CELLS if delta[i]]
        maximal = [s for s in support if not any(t != s and t & s == s for t in support)]
        if any(s.bit_count() == 1 for s in maximal):
            continue
        transformed = rows @ delta
        root = transformed[0]
        if fixed_sign(root):
            continue
        zeros = sum(1 << c for c in CELLS
                    if sum(delta[m] for m in support if m & c == m) % 2 == 0
                    and not excluded_joint_zero(root, transformed[c]))
        cubes = [sum(1 << c for c in CELLS if c & ~s == outside)
                 for s in maximal for outside in CELLS if not outside & s]
        if all(zeros & cube for cube in cubes):
            polynomial = sum(v*m for v, m in zip(delta, MONOMIALS))
            arms.append((pos, neg, zeros, cubes, root, polynomial))
    return masks, arms


@lru_cache(None)
def positive_gap_certificate(polynomial):
    expanded = Counter()
    for exponents, coefficient in sp.Poly(polynomial, *VARIABLES).terms():
        assert coefficient.is_Integer
        for powers, value in expand_monomial(exponents).items():
            expanded[powers] += int(coefficient)*value
    return (expanded[(0, 0, 0, 0)] != 0
            and (min(expanded.values()) >= 0 or max(expanded.values()) <= 0))


@lru_cache(None)
def necessary_factor(polynomial):
    polynomial = sp.expand(polynomial)
    if not polynomial:
        return 0
    # Clear rational denominators exactly; use squarefree factors for zero tests.
    integral = sp.Poly(polynomial, *VARIABLES, domain=sp.QQ).clear_denoms(convert=True)[1]
    return prod(f.as_expr() for f, _ in integral.factor_list()[1]
                if not positive_gap_certificate(f.as_expr()))


def main():
    masks, arms = arm_candidates()
    counts = Counter()
    for left, right in combinations(arms, 2):
        a, b, z, cubes, root, f = left
        c, d, zz, other_cubes, other_root, h = right
        if not all(z & zz & cube for cube in cubes + other_cubes):
            continue
        if excluded_joint_zero(root, other_root) or factors(a, b) & factors(c, d):
            continue
        counts["joint candidates"] += 1
        equations = [necessary_factor(f), necessary_factor(h)]
        for i, j in combinations(VARIABLES, 2):
            minor = sp.diff(f, i)*sp.diff(h, j) - sp.diff(f, j)*sp.diff(h, i)
            equations.append(necessary_factor(minor))
            if equations[-1] == 1:
                break
        if 1 in equations:
            counts["factor sign"] += 1
            continue
        basis = sp.groebner([g for g in equations if g != 0],
                            *VARIABLES, order="grevlex", domain=sp.QQ)
        if any(necessary_factor(g.as_expr()) == 1 for g in basis.polys):
            counts["basis sign"] += 1
            continue
        raise AssertionError(("unresolved", a, b, c, d))
    print("signed differences", len(masks), "individual candidates", len(arms), dict(counts))
    assert (len(masks), len(arms), counts["joint candidates"]) == (3801, 90, 799)


if __name__ == "__main__":
    main()

