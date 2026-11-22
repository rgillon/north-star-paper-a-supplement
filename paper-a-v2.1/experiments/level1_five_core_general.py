"""Exact terminal screens on the signed differences of five-axis downsets."""
from collections import Counter
from functools import lru_cache
from math import prod

import numpy as np
import sympy as sp

from level1_five_core_top_ranks import excluded_by_prime, ordered_rows


def signed_differences():
    ideals = [0]
    for shell in range(32):
        lower = sum(1 << t for t in range(shell) if t & shell == t)
        ideals.extend([m | (1 << shell) for m in ideals if m & lower == lower])
    ideals = np.array(ideals, dtype=np.uint64)
    chunks = []
    for j, first in enumerate(ideals):
        second = ideals[:j]
        pos, neg = first & ~second, second & ~first
        mask = (pos != 0) & (neg != 0)
        low, high = np.minimum(pos[mask], neg[mask]), np.maximum(pos[mask], neg[mask])
        chunks.append((low << 32) | high)
    return np.unique(np.concatenate(chunks))


def possible_zeros(poly):
    constant = poly[..., 0]
    fixed = ((poly.min(axis=-1) >= 0) & (constant > 0)) | (
        (poly.max(axis=-1) <= 0) & (constant < 0))
    return ~fixed & (constant % 2 == 0)


def byte_tables(rows):
    tables = np.zeros((4, 256, *rows.shape[1:]), dtype=np.int32)
    for j in range(4):
        for byte in range(1, 256):
            bit = byte & -byte
            shell = 8*j + bit.bit_length()-1
            tables[j, byte] = tables[j, byte-bit]
            if shell < 31:
                tables[j, byte] += rows[shell].astype(np.int32)
    return tables


def evaluate(tables, pos, neg):
    return sum(tables[j, (pos >> (8*j)) & 255] - tables[j, (neg >> (8*j)) & 255]
               for j in range(len(tables)))


X = sp.symbols("a:e")


@lru_cache(None)
def nonvanishing_factor(factor):
    if factor in X or factor in tuple(x+1 for x in X):
        return True
    for i, x in enumerate(X):
        residue = factor.subs(x, -1)
        if residue == 0:
            continue
        constant, terms = sp.factor_list(residue)
        units = set(X[:i]) | {y+1 for y in X if y != x}
        if all(q == 2 for q in sp.factorint(abs(constant))) and all(
                term in units for term, _ in terms):
            return True
    return False


def nonvanishing_row(poly):
    if poly == 0:
        return False
    return all(nonvanishing_factor(f) for f, _ in sp.factor_list(poly)[1])


def certify():
    rows = ordered_rows()
    tables = byte_tables(rows)
    patterns = signed_differences()
    upper = [sum(1 << t for t in range(32) if t & s == s) for s in range(31)]

    @lru_cache(None)
    def needed_projections(shell, prescribed):
        local = [t ^ shell for t in range(31) if prescribed >> t & 1]
        required = ~possible_zeros(rows[local].sum(axis=0))
        return tuple((c, 31 ^ shell) for c in range(32) if not c & shell and required[c])

    counts = Counter(total=len(patterns))
    residual = []
    for start in range(0, len(patterns), 256):
        batch = patterns[start:start+256]
        pos, neg = batch >> 32, batch & 0xffffffff
        root_ok = possible_zeros(evaluate(tables[:, :, 0], pos, neg))
        pos, neg = pos[root_ok], neg[root_ok]
        if not len(pos):
            continue
        possible = possible_zeros(evaluate(tables, pos, neg))
        for p, n, allowed in zip(pos, neg, possible):
            p, n = int(p), int(n)
            zeros = tuple(map(int, np.flatnonzero(allowed)))
            support = p | n
            cells = [s for s in range(31) if support >> s & 1]
            counts["root"] += 1
            if any(not any(z & mask == c for z in zeros) for s in cells
                   for c, mask in needed_projections(s, support & upper[s])):
                continue
            counts["projection"] += 1
            faces = tuple(u for u in range(32) if any(not (u & s) for s in cells))
            contrasts = [sum((-1)**(31 ^ s).bit_count() for s in cells
                             if (31 ^ s) & u == (31 ^ s)) for u in faces]
            if excluded_by_prime(faces, zeros, contrasts):
                continue
            counts["factor_test"] += 1
            poly = sum((((p >> s) & 1) - ((n >> s) & 1)) *
                       prod(X[i] for i in range(5) if s & (1 << i)) for s in cells)
            if nonvanishing_row(poly):
                counts["root_factor"] += 1
                continue
            zeros = tuple(c for c in zeros if not nonvanishing_row(
                poly.subs({x: -1 for i, x in enumerate(X) if c & (1 << i)})))
            if any(not any(z & mask == c for z in zeros) for s in cells
                   for c, mask in needed_projections(s, support & upper[s])):
                counts["projection_factor"] += 1
                continue
            residual.append((p, n, zeros))
    counts["unresolved"] = len(residual)
    return counts, residual


if __name__ == "__main__":
    counts, residual = certify()
    print(dict(counts), flush=True)
    if residual:
        raise AssertionError(residual)

