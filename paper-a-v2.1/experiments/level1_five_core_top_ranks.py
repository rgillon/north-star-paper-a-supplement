"""Exact sign, projection and prime-divisor tests for missing-pair graphs."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import lcm

import numpy as np
from sympy import Matrix


def ordered_rows():
    monomials = [m for size in range(5)
                 for m in combinations_with_replacement(range(5), size)]
    index = {m: j for j, m in enumerate(monomials)}
    rows = np.zeros((31, 32, len(index)), dtype=np.int64)
    for shell in range(31):
        for conductor in range(32):
            coeff = {(): 1}
            for i in range(5):
                if not shell & (1 << i):
                    continue
                if conductor & (1 << i):
                    coeff = {m: -v for m, v in coeff.items()}
                    continue
                expanded = {m: (2*i+2)*v for m, v in coeff.items()}
                for m, v in coeff.items():
                    for k in range(i+1):
                        term = tuple(sorted((*m, k)))
                        expanded[term] = expanded.get(term, 0) + v
                coeff = expanded
            for term, value in coeff.items():
                rows[shell, conductor, index[term]] = value
    return rows


def sign_patterns(edges):
    for rest in product((1, -1), repeat=len(edges)-1):
        signs = (1, *rest)
        choices = []
        for vertex in range(5):
            incident = {sign for edge, sign in zip(edges, signs) if vertex in edge}
            assert incident
            choices.append((0, *incident) if len(incident) == 1 else (0,))
        for coats in product(*choices):
            if -1 in (*signs, *coats):
                yield (*signs, *coats)


@lru_cache(maxsize=2048)
def dependencies(faces, zeros):
    matrix = Matrix([[int(c & s == s) for s in faces] for c in zeros])
    result = []
    for vector in matrix.nullspace():
        scale = lcm(*(int(n.q) for n in vector))
        result.append(tuple(int(n*scale) for n in vector))
    return result


def excluded_by_prime(faces, zeros, contrasts):
    for relation in dependencies(faces, zeros):
        terms = [(s, w*b) for s, w, b in zip(faces, relation, contrasts) if w*b]
        for i in range(5):
            residues = [v for s, v in terms if s & (1 << i)]
            if len(residues) == 1 and abs(residues[0]) < 2*i+3:
                return True
    return False


def certify(edges, rows):
    counts = Counter()
    patterns = np.array(list(sign_patterns(edges)), dtype=np.int64)
    for order in permutations(range(5)):
        pairs = [(order[i], order[j]) for i, j in edges]
        masks = [(1 << i) | (1 << j) for i, j in pairs]
        cells = [31 ^ s for s in masks] + [31 ^ (1 << i) for i in order]
        faces = (0, 1, 2, 4, 8, 16, *masks)
        leaves = {1 << i for i in range(5) if sum(i in e for e in pairs) == 1}
        polynomial = (patterns @ rows[cells].reshape(len(cells), -1)).reshape(
            len(patterns), 32, -1)
        constant = polynomial[:, :, 0]
        fixed_sign = ((polynomial.min(axis=2) >= 0) & (constant > 0)) | (
            (polynomial.max(axis=2) <= 0) & (constant < 0))
        possible = ~fixed_sign & (constant % 2 == 0)
        for pattern, allowed in zip(patterns, possible):
            if not allowed[0]:
                continue
            counts["root"] += 1
            zeros = tuple(int(c) for c in np.flatnonzero(allowed) if c not in leaves)
            h = [abs(int(pattern[len(edges)+order.index(i)])) for i in range(5)]
            for (i, j), mask in zip(pairs, masks):
                required = {0}
                if h[i] or not h[j]:
                    required.add(1 << i)
                if h[j] or not h[i]:
                    required.add(1 << j)
                if h[i] == h[j]:
                    required.add(mask)
                if not required <= {c & mask for c in zeros}:
                    break
            else:
                counts["projection"] += 1
                contrasts = [0, *(-v for v in h), *(1-h[i]-h[j] for i, j in pairs)]
                if not excluded_by_prime(faces, zeros, contrasts):
                    counts["unresolved"] += 1
    return counts


def graph_representatives(size):
    seen = set()
    for edges in combinations(combinations(range(5), 2), size):
        if edges in seen or len({v for pair in edges for v in pair}) != 5:
            continue
        for order in permutations(range(5)):
            seen.add(tuple(sorted(tuple(sorted((order[i], order[j]))) for i, j in edges)))
        yield edges


if __name__ == "__main__":
    rows = ordered_rows()
    for size in range(3, 11):
        total = Counter()
        for edges in graph_representatives(size):
            total["graphs"] += 1
            total.update(certify(edges, rows))
        print(size, dict(total), flush=True)
        if total["unresolved"]:
            raise AssertionError(total)

