"""Exact ordered-prime terminal tests on grids with at most 36 cells."""
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from math import lcm, prod

import numpy as np
import sympy as sp

from level1_five_core_general import possible_zeros, nonvanishing_row, evaluate, X
from level1_three_axis_terminal import below, row


def differences(cells):
    ideals = [0]
    for i, s in enumerate(cells):
        lower = sum(1 << j for j, v in enumerate(cells[:i]) if below(v, s))
        ideals += [m | (1 << i) for m in ideals if m & lower == lower]
    ideals = np.array(ideals, dtype=np.uint64)
    # Distinct signed differences have distinct balanced ternary encodings.
    # This packs two 36-bit masks into a single sortable integer.
    wide = len(cells) > 32
    if wide:
        signed_keys = np.array([sum(3**i for i in range(len(cells)) if int(m) >> i & 1)
                                for m in ideals], dtype=np.int64)
    chunks = []
    for j, first in enumerate(ideals):
        second = ideals[:j]
        pos, neg = first & ~second, second & ~first
        keep = (pos != 0) & (neg != 0)
        if wide:
            chunks.append(np.abs(signed_keys[j] - signed_keys[:j][keep]))
        else:
            chunks.append((np.minimum(pos[keep], neg[keep]) << 32) |
                          np.maximum(pos[keep], neg[keep]))
    return np.unique(np.concatenate(chunks))


def largest_prime_excludes(poly, axes):
    # At the largest prime, even a smaller prime plus one is a unit.
    residue = poly.subs(X[axes-1], -1)
    if residue == 0:
        return False
    constant, factors = sp.factor_list(residue)
    units = {x+i for x in X[:axes-1] for i in range(3)}
    return all(q == 2 for q in sp.factorint(abs(constant))) and all(
        f in units for f, _ in factors)


def certify(shape):
    cells = list(product(*(range(a+1) for a in shape)))
    size, axes = len(cells), len(shape)
    assert size <= 36 and sum(shape) <= 6 and axes <= 4
    gaps = sp.symbols(f"g:{axes}")
    primes = [3+2*i+sum(gaps[:i+1]) for i in range(axes)]
    mons = sorted(m for m in product(range(sum(shape)+1), repeat=axes)
                  if sum(m) <= sum(shape))
    rows = np.zeros((size, size, len(mons)), dtype=np.int64)
    raw = []
    for i, s in enumerate(cells):
        raw.append([])
        for j, h in enumerate(cells):
            expression = prod(row(pr, k, v) for pr, k, v in zip(primes, h, s))
            coefficients = sp.Poly(expression, *gaps).as_dict()
            rows[i, j] = [int(coefficients.get(m, 0)) for m in mons]
            raw[i].append(prod(row(x+1, k, v) for x, k, v in zip(X, h, s)))
    upper = [sum(1 << k for k, v in enumerate(cells) if below(s, v)) for s in cells]
    adjacent = [[(j, (-1)**sum(vv-uu for vv, uu in zip(v, u)))
                 for j,v in enumerate(cells) if all(0 <= vv-uu <= 1 for vv,uu in zip(v,u))]
                for u in cells]
    incidence = [[prod(int(k == v)-int(k+1 == v) for k,v in zip(h,u))
                  for u in cells] for h in cells]

    @lru_cache(None)
    def cubes(i, supported):
        s = cells[i]
        boxes = []
        for h in cells:
            if any(k for k,v in zip(h,s) if v):
                continue
            coeff = np.zeros(len(mons), dtype=np.int64)
            for t,v in enumerate(cells):
                if not supported >> t & 1:
                    continue
                if any(b-a > 1 for a,b in zip(s,v) if a):
                    continue
                projected = tuple(0 if a else b for a,b in zip(s,v))
                sign = (-1)**sum(b-a for a,b in zip(s,v) if a)
                coeff += sign*rows[cells.index(projected),cells.index(h)]
            if not possible_zeros(coeff):
                choices = [(a-1,a) if a else (k,) for a,k in zip(s,h)]
                boxes.append(sum(1 << cells.index(v) for v in product(*choices)))
        return boxes

    @lru_cache(maxsize=4096)
    def relations(indices, zero_indices):
        matrix = sp.Matrix([[incidence[j][i] for i in indices] for j in zero_indices])
        out = []
        for v in matrix.nullspace():
            scale = lcm(*(int(x.q) for x in v))
            out.append(tuple(int(x*scale) for x in v))
        return out

    def excluded(support, zeros, boxes):
        if any(not (box & zeros) for box in boxes):
            return "cube"
        indices = tuple(i for i in range(size) if any(support >> j & 1 and
                        upper[j] >> i & 1 for j in range(size)))
        zero_indices = tuple(i for i in range(size) if zeros >> i & 1)
        b = [sum(sign for j,sign in adjacent[i] if support >> j & 1) for i in indices]
        for v in relations(indices, zero_indices):
            terms = [(cells[i], t*u) for i,t,u in zip(indices,v,b) if t*u]
            if not terms:
                continue
            for axis in range(axes):
                power = min(s[axis] for s,_ in terms)
                digits = [w for s,w in terms if s[axis] == power]
                if len(digits) == 1 and abs(digits[0]) < 2*axis+3:
                    return "incidence"
        return None

    # Lookup signed shell sums in blocks of eight. A byte contributes at
    # most eight times 13**6 in absolute coefficient sum for these shapes.
    nbytes = (size+7)//8
    tables = np.zeros((nbytes, 256, size, len(mons)), dtype=np.int32)
    for j in range(nbytes):
        for byte in range(1, 256):
            bit = byte & -byte
            shell = 8*j + bit.bit_length() - 1
            tables[j, byte] = tables[j, byte-bit]
            if shell < size:
                tables[j, byte] += rows[shell].astype(np.int32)

    patterns = differences(cells)
    counts = Counter(total=len(patterns))
    residual = []
    for start in range(0, len(patterns), 256):
        batch = patterns[start:start+256]
        if size <= 32:
            pos, neg = batch >> 32, batch & 0xffffffff
        else:
            digits = batch + (3**size-1)//2
            pos, neg = (np.zeros(len(batch), dtype=np.uint64) for _ in range(2))
            for i in range(size):
                digit = digits % 3
                pos |= (digit == 2).astype(np.uint64) << i
                neg |= (digit == 0).astype(np.uint64) << i
                digits = digits // 3
        shifts = np.arange(size, dtype=np.uint64)
        values = (((pos[:,None] >> shifts) & 1).astype(np.int64) -
                  ((neg[:,None] >> shifts) & 1).astype(np.int64))
        keep = possible_zeros(evaluate(tables[:, :, -1], pos, neg))
        if not keep.any():
            continue
        polys = evaluate(tables, pos[keep], neg[keep])
        for p,n,vals,ok in zip(pos[keep],neg[keep],values[keep],possible_zeros(polys)):
            counts["root"] += 1
            support = int(p | n)
            zeros = sum(1 << int(j) for j in np.flatnonzero(ok))
            boxes = [box for i in range(size) if support >> i & 1
                     for box in cubes(i,support & upper[i])]
            why = excluded(support,zeros,boxes)
            if why:
                counts[why] += 1
                continue
            def equation(j):
                return sp.expand(sum(int(vals[i])*raw[i][j] for i in range(size) if vals[i]))
            if nonvanishing_row(equation(size-1)):
                counts["root_factor"] += 1
                continue
            zeros = sum(1 << j for j in range(size) if zeros >> j & 1 and
                        not nonvanishing_row(equation(j)))
            if excluded(support,zeros,boxes):
                counts["row_factor"] += 1
            elif largest_prime_excludes(equation(size-1), axes):
                counts["largest_prime"] += 1
            else:
                residual.append((int(p),int(n),zeros))
    return counts, residual


if __name__ == "__main__":
    shapes = sorted({p for shape in [(2,2,1,1), (3,1,1,1), (2,1,1,1), (4,1,1), (3,2,1), (2,2,2)]
                     for p in permutations(shape)})
    for shape in shapes:
        counts, residual = certify(shape)
        print(shape, dict(counts), residual, flush=True)
        if residual:
            raise AssertionError(residual)

