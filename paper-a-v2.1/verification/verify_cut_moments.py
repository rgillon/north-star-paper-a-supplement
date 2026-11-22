#!/usr/bin/env python3
"""Exact controls for the source-independent cut coefficient criterion.

Standard library only. Run: python3 verify_cut_moments.py
This does not verify the unrecovered LEVEL1_FACTOR_CUT_MOMENTS.md.
The positive control deliberately enlarges the class to radial binary sets.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from math import gcd
import json


@lru_cache(None)
def divisors(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


@lru_cache(None)
def factor(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            a = 0
            while n % p == 0:
                n //= p
                a += 1
            out.append((p, a))
        p += 1
    if n > 1:
        out.append((n, 1))
    return tuple(out)


def mu(n):
    f = factor(n)
    return 0 if any(a > 1 for _, a in f) else (-1) ** len(f)


@lru_cache(None)
def ramanujan(d, k):
    return sum(h * mu(d // h) for h in divisors(gcd(d, k)))


@lru_cache(None)
def ideals(n):
    ds = divisors(n)
    predecessors = [sum(1 << j for j, e in enumerate(ds[:i]) if d % e == 0)
                    for i, d in enumerate(ds)]
    result = []

    def visit(i, mask):
        if i == len(ds):
            result.append(mask)
            return
        visit(i + 1, mask)
        if mask & predecessors[i] == predecessors[i]:
            visit(i + 1, mask | (1 << i))

    visit(0, 0)
    return tuple(result)


def spectrum(n, mask):
    ds = divisors(n)
    return tuple(sum(ramanujan(d, n // c) for i, d in enumerate(ds)
                     if mask >> i & 1) for c in ds)


def literal_correlation(n, mask):
    ds = divisors(n)
    order_bit = {d: i for i, d in enumerate(ds)}
    xs = [x for x in range(n) if mask >> order_bit[n // gcd(x, n)] & 1]
    values = [0] * n
    for x in xs:
        for y in xs:
            values[(y - x) % n] += 1
    return tuple(values)


def square_signature(v):
    return tuple(x * x for x in v)


def polynomial(vectors):
    return Counter(square_signature(v) for v in vectors)


def cut_moment(plus, minus):
    inner = sum(count * minus.get(sig, 0) for sig, count in plus.items())
    support = set(plus) | set(minus)
    norm_difference = sum((plus.get(s, 0) + minus.get(s, 0)) ** 2
                          - (plus.get(s, 0) - minus.get(s, 0)) ** 2
                          for s in support)
    assert norm_difference == 4 * inner
    assert inner >= 0
    return inner


def transfer(n, p, a, cut, bit):
    """Recurrence Q, retaining an actual ideal-mask witness for every monomial."""
    m = n // p ** a
    ds, ms, cs = divisors(n), divisors(m), divisors(m)
    states = ideals(m)
    index_n = {d: i for i, d in enumerate(ds)}
    cut_shell = divisors(n)[cut]
    cut_t = 0
    cut_m = cut_shell
    while cut_m % p == 0:
        cut_t += 1
        cut_m //= p
    cut_m_bit = ms.index(cut_m)
    # Use global conductor ordering in exponent vectors.
    coords = [(next(j for j in range(a + 1) if c % p ** (j + 1) != 0),
               c // p ** next(j for j in range(a + 1) if c % p ** (j + 1) != 0))
              for c in ds]
    fj = {j: {c: spectrum(m, j)[ms.index(c)] for c in cs} for j in states}
    q = {}
    for t in range(a + 1):
        next_q = {}
        for state in states:
            if t == cut_t and (state >> cut_m_bit & 1) != bit:
                continue
            contribution = tuple(ramanujan(p ** t, p ** (a - j)) * fj[state][c]
                                 for j, c in coords)
            shell_mask = sum(1 << index_n[p ** t * d] for k, d in enumerate(ms)
                             if state >> k & 1)
            entries = {}
            if t == 0:
                entries[contribution] = shell_mask
            else:
                for upper, monomials in q.items():
                    if state & ~upper:
                        continue
                    for vector, witness in monomials.items():
                        vector = tuple(x + y for x, y in zip(vector, contribution))
                        new_witness = witness | shell_mask
                        assert vector not in entries or entries[vector] == new_witness
                        entries[vector] = new_witness
            if entries:
                next_q[state] = entries
        q = next_q
    result = {}
    for entries in q.values():
        for vector, witness in entries.items():
            assert vector not in result or result[vector] == witness
            result[vector] = witness
    return result


def verify_modulus(n):
    ds = divisors(n)
    ideal_list = ideals(n)
    full = {mask: spectrum(n, mask) for mask in ideal_list}
    literal = {mask: literal_correlation(n, mask) for mask in ideal_list}
    # Compare equivalence relations by both directions, not merely class counts.
    sigma_to_c = {}
    c_to_sigma = {}
    for mask in ideal_list:
        sig, corr = square_signature(full[mask]), literal[mask]
        transformed_corr = tuple(sum(ramanujan(d, n // c) * corr[(n // d) % n]
                                     for d in ds) for c in ds)
        assert sig == transformed_corr
        reconstructed = tuple(sum(full[mask][j] * ramanujan(c, n // d)
                                  for j, c in enumerate(ds)) for d in ds)
        assert reconstructed == tuple(n * (mask >> i & 1) for i in range(len(ds)))
        assert sig not in sigma_to_c or sigma_to_c[sig] == corr
        assert corr not in c_to_sigma or c_to_sigma[corr] == sig
        sigma_to_c[sig] = corr
        c_to_sigma[corr] = sig
    transfers = 0
    moments = []
    for cut in range(len(ds)):
        partitions = [{v: mask for mask, v in full.items() if (mask >> cut & 1) == bit}
                      for bit in (0, 1)]
        for p, a in factor(n):
            for bit in (0, 1):
                assert transfer(n, p, a, cut, bit) == partitions[bit]
                transfers += 1
        pm = [polynomial(partitions[bit]) for bit in (0, 1)]
        k = cut_moment(pm[1], pm[0])
        lc = [Counter(literal[mask] for mask in ideal_list
                      if (mask >> cut & 1) == bit) for bit in (0, 1)]
        assert k == sum(count * lc[0].get(corr, 0) for corr, count in lc[1].items())
        moments.append(k)
    return dict(N=n, ideals=len(ideal_list), cuts=len(ds), transfer_comparisons=transfers,
                distinct_signatures=len(sigma_to_c), moments=moments)


def radial_positive_control():
    n = 12
    ds = divisors(n)
    masks = range(1 << len(ds))
    spec = {mask: spectrum(n, mask) for mask in masks}
    literal = {mask: literal_correlation(n, mask) for mask in masks}
    for mask in masks:
        assert square_signature(spec[mask]) == tuple(
            sum(ramanujan(d, n // c) * literal[mask][(n // d) % n] for d in ds)
            for c in ds)
    results = []
    witness = None
    for cut, d in enumerate(ds):
        pm = [polynomial(v for mask, v in spec.items() if (mask >> cut & 1) == bit)
              for bit in (0, 1)]
        k = cut_moment(pm[1], pm[0])
        lc = [Counter(literal[mask] for mask in masks if (mask >> cut & 1) == bit)
              for bit in (0, 1)]
        assert k == sum(count * lc[0].get(corr, 0) for corr, count in lc[1].items())
        results.append(k)
        if k and witness is None:
            for left in masks:
                for right in masks:
                    if (left >> cut & 1) and not (right >> cut & 1):
                        if square_signature(spec[left]) == square_signature(spec[right]):
                            assert literal[left] == literal[right]
                            assert left not in ideals(n) or right not in ideals(n)
                            witness = dict(cut=d, left_orders=[e for i, e in enumerate(ds) if left >> i & 1],
                                           right_orders=[e for i, e in enumerate(ds) if right >> i & 1])
                            break
                if witness:
                    break
    assert witness is not None
    return dict(N=n, class_name="all radial binary sets, not subgroup unions",
                cuts=list(ds), moments=results, witness=witness)


if __name__ == "__main__":
    moduli = list(range(1, 73)) + [84, 90, 96, 108, 120, 126, 132, 140, 144, 150, 168, 180, 210]
    rows = [verify_modulus(n) for n in moduli]
    print(json.dumps(dict(status="PASS", exact_arithmetic=True, moduli=len(rows),
                          total_ideals=sum(x["ideals"] for x in rows),
                          total_cut_moments=sum(x["cuts"] for x in rows),
                          total_transfer_comparisons=sum(x["transfer_comparisons"] for x in rows),
                          divisor_union_cases=rows,
                          radial_positive_control=radial_positive_control()), indent=2))
