#!/usr/bin/env python3
"""Exact check of the two referee cancellations, at all eight conductors.

No construction or existing verification module is imported.
The shell-factor calculation is checked against a direct residue polynomial
reduced modulo the appropriate cyclotomic polynomial, using SymPy exactly.
This finite check is a regression, not the proof of the general implication.
"""
import itertools
import json
import math
import platform
from pathlib import Path
import sympy as sp

X = sp.Symbol("X")
cases = [
    ((3, 5, 7), {1: -1, 2: -1, 4: 1}, "rank_one"),
    ((5, 7, 13), {6: 1, 5: -1, 3: -1}, "rank_two"),
]
records = []
for primes, shell_values, label in cases:
    K = math.prod(primes)
    profile = []
    for y in range(K):
        order = K // math.gcd(K, y)
        mask = sum(1 << i for i, p in enumerate(primes) if order % p == 0)
        profile.append(shell_values.get(mask, 0))
    spectra = []
    for mask in range(8):
        conductor = math.prod(p for i, p in enumerate(primes) if mask >> i & 1)
        kernel = tuple(0 if mask >> i & 1 else 1 for i in range(3))
        factors = tuple(-1 if mask >> i & 1 else p - 1 for i, p in enumerate(primes))
        value = sum(
            sign * math.prod(factors[i] for i in range(3) if physical_mask >> i & 1)
            for physical_mask, sign in shell_values.items()
        )
        coefficients = [0] * conductor
        for y, f in enumerate(profile):
            coefficients[y % conductor] += f
        residue_polynomial = sp.Poly.from_list(list(reversed(coefficients)), X)
        cyclotomic = sp.Poly(sp.cyclotomic_poly(conductor, X), X)
        remainder = residue_polynomial.rem(cyclotomic)
        assert remainder == sp.Poly(value, X)
        assert (conductor != 1) == any(v == 0 for v in kernel)
        assert (conductor != 1) == any(v == -1 for v in factors)
        spectra.append(
            dict(conductor=conductor, kernel_exponents=kernel,
                 local_shell_factors=factors, exact_fourier_value=value,
                 direct_cyclotomic_remainder=str(remainder.as_expr()))
        )
    zeros = [v["conductor"] for v in spectra if v["exact_fourier_value"] == 0]
    assert zeros == [1]
    records.append(dict(case=label, primes=primes, modulus=K,
                        signed_shell_masks=shell_values, zero_conductors=zeros,
                        conductors=spectra))

out = dict(status="PASS", python=platform.python_version(), sympy=sp.__version__,
           scope="Exact eight-conductor controls for the two referee cancellations; not an exhaustive rigidity search.",
           cases=records)
print(json.dumps(out, indent=2))
