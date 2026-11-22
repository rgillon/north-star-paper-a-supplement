"""Check all bounded fourth primes for selected three-prime cores."""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from math import gcd, prod

from sympy import factor_list, primerange, symbols

AXES = symbols("x:4")
MONOMIALS = [prod(AXES[j] for j in range(4) if mask >> j & 1)
             for mask in range(16)]


def four_axis_ideals():
    ideals = [0]
    for cell in range(16):
        required = sum(1 << (cell ^ (1 << j)) for j in range(4) if cell >> j & 1)
        ideals += [ideal | (1 << cell) for ideal in ideals
                   if ideal & required == required]
    assert len(ideals) == 168
    return ideals


@lru_cache(None)
def factors(positive, negative):
    poly = sum(term * (((positive >> i) & 1) - ((negative >> i) & 1))
               for i, term in enumerate(MONOMIALS))
    return {factor for factor, _ in factor_list(poly)[1]}


def surviving_ratios(primes, ideals):
    weights = [prod(primes[j]-1 for j in range(4) if cell >> j & 1)
               for cell in range(16)]
    populations = defaultdict(list)
    for ideal in ideals:
        populations[sum(w for i, w in enumerate(weights) if ideal >> i & 1)].append(ideal)
    directions = defaultdict(set)
    for group in populations.values():
        for first, second in combinations(group, 2):
            delta = [((first >> i) & 1) - ((second >> i) & 1) for i in range(16)]
            moments = [sum(w*d for i, (w, d) in enumerate(zip(weights, delta)) if i >> j & 1)
                       for j in range(4)]
            content = gcd(*moments)
            assert content  # FM1 on four axes
            if next(value for value in moments if value) < 0:
                moments = [-value for value in moments]
                first, second = second, first
            directions[tuple(value//content for value in moments)].add(
                (content, first & ~second, second & ~first))
    ratios = set()
    for group in directions.values():
        for (a, pos, neg), (b, pos2, neg2) in combinations(group, 2):
            if a == b or factors(pos, neg) & factors(pos2, neg2):
                continue
            g = gcd(a, b)
            ratios.add(tuple(sorted((a//g, b//g))))
    return ratios


if __name__ == "__main__":
    ideals = four_axis_ideals()
    for triple in combinations((3, 5, 7, 11, 13), 3):
        candidates = [p for p in primerange(3, prod(triple)+2) if p not in triple]
        surviving = [(p, r) for p in candidates
                     if (r := surviving_ratios((*triple, p), ideals))]
        print(triple, len(candidates), surviving)

