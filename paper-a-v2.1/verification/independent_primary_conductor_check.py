#!/usr/bin/env python3
"""Exact independent checks of the two reported three-prime cancellations.

For all eight conductor classes, compare the shell-product expression with
literal residue enumeration reduced modulo the cyclotomic polynomial. No
North Star verifier or third-party package is imported. These finite examples
illustrate a separately proved character identity, not a proof by enumeration.
"""
from functools import lru_cache
from fractions import Fraction
from math import gcd, prod
import json

def trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a

def remainder(a, b):
    a = trim(list(a))
    assert b[-1] == 1
    while len(a) >= len(b) and any(a):
        shift = len(a) - len(b)
        lead = a[-1]
        for j, coeff in enumerate(b):
            a[shift+j] -= lead*coeff
        trim(a)
    return a

def monic_quotient(a, b):
    a = trim(list(a))
    assert b[-1] == 1
    q = [0]*max(1, len(a)-len(b)+1)
    while len(a) >= len(b) and any(a):
        shift = len(a)-len(b)
        lead = a[-1]
        q[shift] = lead
        for j, coeff in enumerate(b):
            a[shift+j] -= lead*coeff
        trim(a)
    assert a == [0]
    return tuple(trim(q))

@lru_cache(None)
def cyclotomic(n):
    f = [-1]+[0]*(n-1)+[1]
    for d in range(1, n):
        if n % d == 0:
            f = monic_quotient(f, cyclotomic(d))
    return tuple(f)

def run(primes, rank, missing_signs):
    n = prod(primes)
    if rank == 1:
        profile = {p:s for p,s in zip(primes, missing_signs)}
    else:
        assert rank == 2
        profile = {n//p:s for p,s in zip(primes, missing_signs)}
    rows = []
    for mask in range(8):
        nontrivial_axes = [i for i in range(3) if (mask>>i)&1]
        conductor = prod(primes[i] for i in nontrivial_axes)
        xs = [-1 if i in nontrivial_axes else p-1
              for i,p in enumerate(primes)]
        expression = (sum(s*x for s,x in zip(missing_signs, xs))
            if rank == 1 else
            sum(missing_signs[i]*prod(xs[j] for j in range(3) if j != i)
                for i in range(3)))
        # Literal physical points; their character is exp(2*pi*i*x/c).
        poly = [0]*conductor
        for x in range(n):
            order = n//gcd(n, x)
            poly[x % conductor] += profile.get(order, 0)
        residue = remainder(poly, cyclotomic(conductor))
        assert residue == [expression], (primes, conductor, residue, expression)
        assert (conductor == 1) == (xs == [p-1 for p in primes])
        assert (conductor > 1) == any(x == -1 for x in xs)
        rows.append({
            'conductor': conductor,
            'nontrivial_axis_primes': [primes[i] for i in nontrivial_axes],
            'kernel_exponents': [0 if i in nontrivial_axes else 1 for i in range(3)],
            'singleton_shell_sums': xs,
            'Fourier_value': expression,
            'literal_cyclotomic_remainder': residue,
        })
    zero_conductors = [r['conductor'] for r in rows if r['Fourier_value'] == 0]
    assert zero_conductors == [1]
    return {
        'primes': primes, 'modulus': n, 'rank': rank,
        'signs_by_singleton_or_missing_axis': missing_signs,
        'physical_shell_profile': profile,
        'zero_conductors': zero_conductors,
        'rows': sorted(rows, key=lambda r:r['conductor'])
    }

if __name__ == '__main__':
    rank1 = run([3,5,7], 1, [-1,-1,1])
    rank2 = run([5,7,13], 2, [1,-1,-1])
    assert Fraction(1,4)-Fraction(1,6)-Fraction(1,12) == 0
    print(json.dumps({
        'status':'PASS',
        'scope':'Two explicit dangerous triples, all eight conductor classes each; literal cyclic residue polynomials checked exactly by cyclotomic reduction.',
        'third_party_imports':False,
        'examples':[rank1,rank2],
        'conclusion':'Both displayed cancellations occur only at conductor 1, the trivial character. Every nontrivial conductor has a singleton-shell sum -1 and a nonzero displayed polynomial.'
    }, indent=2))
